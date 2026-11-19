//! Scrambling the bases of a read as vsearch's `--scramble` does, drawn
//! step for step as vsearch draws it. It is a shuffle of the bases or, for
//! words of k bases, a uniformly random sequence with the same count of
//! every word of up to k bases (uShuffle, Jiang et al. 2008).

use std::collections::HashMap;

use crate::random::{Generator, bounded, shuffle};

/// The de Bruijn multigraph of a read: one vertex per distinct word of
/// k - 1 bases, numbered by first appearance, and one edge per word of k
/// bases, grouped by source in the order of the read.
#[derive(Default)]
struct Graph {
    /// Where each vertex's edges start in `targets`, and one past the last.
    offsets: Vec<usize>,
    /// The target vertex of each edge.
    targets: Vec<u32>,
    /// The last base of each vertex's word.
    last_bases: Vec<u8>,
    /// The vertex of the word at each position.
    positions: Vec<u32>,
    ids: HashMap<u64, u32>,
}

impl Graph {
    fn build(&mut self, seq: &[u8], width: usize) {
        self.ids.clear();
        self.last_bases.clear();
        self.positions.clear();
        // words packed into a rolling key, a byte per base
        let mask = if width == 8 {
            u64::MAX
        } else {
            (1u64 << (8 * width)) - 1
        };
        let mut key = 0u64;
        for (i, &base) in seq.iter().enumerate() {
            key = ((key << 8) | u64::from(base)) & mask;
            if i + 1 < width {
                continue;
            }
            let next = self.last_bases.len() as u32;
            let vertex = *self.ids.entry(key).or_insert(next);
            if vertex == next {
                self.last_bases.push(base);
            }
            self.positions.push(vertex);
        }
        let vertices = self.last_bases.len();
        self.offsets.clear();
        self.offsets.resize(vertices + 1, 0);
        for &source in &self.positions[..self.positions.len() - 1] {
            self.offsets[source as usize + 1] += 1;
        }
        for v in 0..vertices {
            self.offsets[v + 1] += self.offsets[v];
        }
        self.targets.clear();
        self.targets.resize(self.positions.len() - 1, 0);
        let mut cursors = self.offsets[..vertices].to_vec();
        for pair in self.positions.windows(2) {
            let source = pair[0] as usize;
            self.targets[cursors[source]] = pair[1];
            cursors[source] += 1;
        }
    }

    fn edges(&self, vertex: u32) -> std::ops::Range<usize> {
        self.offsets[vertex as usize]..self.offsets[vertex as usize + 1]
    }
}

/// Scrambles reads, reusing its buffers from one read to the next.
#[derive(Default)]
pub struct Scrambler {
    graph: Graph,
    last_exits: Vec<u32>,
    in_tree: Vec<bool>,
    cursors: Vec<usize>,
}

impl Scrambler {
    /// scrambles `seq` in place keeping the count of every word of up to
    /// `k` bases (from 1 to 9); a read of k bases or fewer is left as it
    /// is when k is 2 or more.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::random::SplitMix64;
    /// # use rsearch::scramble::Scrambler;
    /// let mut seq = b"ACGTTGCAAC".to_vec();
    /// let mut scrambler = Scrambler::default();
    /// scrambler.scramble(&mut seq, 2, &mut SplitMix64::new(42));
    /// let mut pairs: Vec<_> = seq.windows(2).collect();
    /// pairs.sort();
    /// let mut before: Vec<_> = b"ACGTTGCAAC".windows(2).collect();
    /// before.sort();
    /// assert_eq!(pairs, before);
    /// ```
    pub fn scramble(&mut self, seq: &mut [u8], k: usize, generator: &mut impl Generator) {
        assert!((1..=9).contains(&k), "k must be from 1 to 9");
        if k == 1 {
            shuffle(seq, generator);
            return;
        }
        if seq.len() <= k {
            return;
        }
        self.graph.build(seq, k - 1);
        self.select_last_exits(generator);
        self.shuffle_edges(generator);
        self.walk(seq, k - 1);
    }

    /// Wilson's algorithm: from each vertex, a loop-erased random walk to
    /// the tree rooted at the read's last vertex, which picks each vertex's
    /// last exit.
    fn select_last_exits(&mut self, generator: &mut impl Generator) {
        let graph = &self.graph;
        let vertices = graph.last_bases.len();
        let end = *graph.positions.last().unwrap();
        self.last_exits.clear();
        self.last_exits.resize(vertices, 0);
        self.in_tree.clear();
        self.in_tree.resize(vertices, false);
        self.in_tree[end as usize] = true;
        for vertex in 0..vertices as u32 {
            let mut walker = vertex;
            while !self.in_tree[walker as usize] {
                let edges = graph.edges(walker);
                let pick = bounded(generator, edges.len() as u64) as usize;
                let target = graph.targets[edges.start + pick];
                self.last_exits[walker as usize] = target;
                walker = target;
            }
            walker = vertex;
            while !self.in_tree[walker as usize] {
                self.in_tree[walker as usize] = true;
                walker = self.last_exits[walker as usize];
            }
        }
    }

    /// keeps one edge to each vertex's last exit at the end of its edges,
    /// and shuffles the others; the last vertex keeps none back.
    fn shuffle_edges(&mut self, generator: &mut impl Generator) {
        let graph = &mut self.graph;
        let end = *graph.positions.last().unwrap();
        for vertex in 0..graph.last_bases.len() as u32 {
            let range = graph.edges(vertex);
            let edges = &mut graph.targets[range];
            if vertex == end {
                shuffle(edges, generator);
            } else {
                let exit = self.last_exits[vertex as usize];
                let kept = edges.iter().position(|&t| t == exit).unwrap();
                let last = edges.len() - 1;
                edges.swap(kept, last);
                shuffle(&mut edges[..last], generator);
            }
        }
    }

    /// walks from the first vertex, taking each vertex's edges in order,
    /// and writes the last base of each vertex reached.
    fn walk(&mut self, seq: &mut [u8], width: usize) {
        let graph = &self.graph;
        self.cursors.clear();
        self.cursors.resize(graph.last_bases.len(), 0);
        let mut current = graph.positions[0];
        for base in &mut seq[width..] {
            let edges = graph.edges(current);
            let target = graph.targets[edges.start + self.cursors[current as usize]];
            self.cursors[current as usize] += 1;
            *base = graph.last_bases[target as usize];
            current = target;
        }
    }
}
