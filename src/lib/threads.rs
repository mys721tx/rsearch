//! Running a tool's work on several threads without changing its output.
//!
//! The tools that search or align give every record to one of the workers
//! of a pool and write the results in input order. A worker keeps its own
//! state (a searcher, a classifier) from one record to the next. With one
//! thread no pool is built and the work runs inline, as a plain loop.

use std::sync::{Mutex, MutexGuard};

use clap::{Arg, ArgMatches, Command, value_parser};
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};

/// adds `-t/--threads`, the threads to run on, 0 (the default) for one per
/// logical core.
pub fn threads_args(command: Command) -> Command {
    command.arg(
        Arg::new("threads")
            .short('t')
            .long("threads")
            .value_name("COUNT")
            .value_parser(value_parser!(usize))
            .default_value("0")
            .help("Run on this many threads, 0 for one per core"),
    )
}

/// the threads `-t/--threads` asks for, at least one.
pub fn threads(matches: &ArgMatches) -> usize {
    match *matches.get_one::<usize>("threads").unwrap() {
        0 => std::thread::available_parallelism().map_or(1, |n| n.get()),
        n => n,
    }
}

/// A pool of workers, each with its own state `S`.
///
/// # Example
/// ```rust
/// # use rsearch::threads::Workers;
/// let workers = Workers::new(4, || 0u64);
/// let squares = workers.map((1..=5u64).collect(), |calls, x| {
///     *calls += 1;
///     x * x
/// });
/// assert_eq!(squares, [1, 4, 9, 16, 25]);
///
/// let mut items = (0..1000u64).into_iter();
/// let mut written = Vec::new();
/// workers.stream(64, || items.next(), |_, x| x + 1, |y| written.push(y));
/// assert_eq!(written, (1..=1000).collect::<Vec<_>>());
/// ```
pub struct Workers<S> {
    pool: Option<ThreadPool>,
    states: Vec<Mutex<S>>,
}

impl<S: Send> Workers<S> {
    /// `threads` workers, each with a state made by `init`. With one, the
    /// work runs on the calling thread.
    pub fn new(threads: usize, mut init: impl FnMut() -> S) -> Workers<S> {
        let threads = threads.max(1);
        let pool = (threads > 1).then(|| {
            ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("a thread pool")
        });
        Workers {
            pool,
            states: (0..threads).map(|_| Mutex::new(init())).collect(),
        }
    }

    /// the number of workers.
    pub fn threads(&self) -> usize {
        self.states.len()
    }

    /// the state of the first worker, for serial work between batches.
    pub fn state(&self) -> MutexGuard<'_, S> {
        self.states[0].lock().unwrap()
    }

    /// `work` on each item, by the workers, the results in the items' order.
    pub fn map<T, R>(&self, items: Vec<T>, work: impl Fn(&mut S, T) -> R + Sync) -> Vec<R>
    where
        T: Send,
        R: Send,
    {
        match &self.pool {
            None => {
                let mut state = self.state();
                items.into_iter().map(|t| work(&mut state, t)).collect()
            }
            Some(pool) => pool.install(|| {
                items
                    .into_par_iter()
                    .map(|t| {
                        let index = rayon::current_thread_index().unwrap_or(0);
                        work(&mut self.states[index].lock().unwrap(), t)
                    })
                    .collect()
            }),
        }
    }

    /// reads items with `next` until it returns `None`, gives each to `work`,
    /// and hands the results to `emit` in the items' order. Items go to the
    /// workers `batch` at a time. While the workers run on one batch, the
    /// calling thread emits the results of the one before and reads the one
    /// after. `next` and `emit` run only on the calling thread. `next` is not
    /// called again after it returns `None`.
    ///
    /// An item that stands for an error (a record that could not be read)
    /// should make `next` return `None` after it, and `emit` fail on its
    /// result. Every result before it is then emitted first, as in a serial
    /// loop.
    pub fn stream<T, R>(
        &self,
        batch: usize,
        mut next: impl FnMut() -> Option<T>,
        work: impl Fn(&mut S, T) -> R + Sync,
        mut emit: impl FnMut(R),
    ) where
        T: Send,
        R: Send,
    {
        let Some(pool) = &self.pool else {
            let mut state = self.state();
            while let Some(t) = next() {
                emit(work(&mut state, t));
            }
            return;
        };
        let batch = batch.max(1);
        let mut exhausted = false;
        let mut read = || {
            let mut items = Vec::with_capacity(batch);
            while !exhausted && items.len() < batch {
                match next() {
                    Some(t) => items.push(t),
                    None => exhausted = true,
                }
            }
            items
        };
        let mut current = read();
        let mut done: Vec<R> = Vec::new();
        while !current.is_empty() {
            let mut results: Vec<R> = Vec::new();
            let mut following = Vec::new();
            let work = &work;
            pool.in_place_scope(|scope| {
                let items = std::mem::take(&mut current);
                let results = &mut results;
                scope.spawn(move |_| *results = self.map(items, work));
                for r in done.drain(..) {
                    emit(r);
                }
                following = read();
            });
            done = results;
            current = following;
        }
        for r in done {
            emit(r);
        }
    }
}
