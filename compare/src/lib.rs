//! Compares rsearch with vsearch.
//!
//! The crate translates a vsearch command line into the rsearch pipeline
//! that gives the same results ([`translate`]), runs both and compares
//! their outputs, exit codes, errors and warnings through named
//! comparators ([`execute`], [`comparators`], [`cases`]), makes its cases
//! from seeded generators of random and malformed input ([`generators`],
//! [`sources`]), and runs them in parallel, recording each as a line of
//! JSON ([`runner`]).

pub mod options;
pub mod translate;

pub mod generators;
pub mod rng;

pub mod cases;
pub mod comparators;
pub mod execute;
pub mod runner;
pub mod sources;
