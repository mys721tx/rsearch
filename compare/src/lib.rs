//! Compares rsearch with vsearch.
//!
//! The crate translates a vsearch command line into the rsearch pipeline
//! that gives the same results ([`translate`]). It runs both and compares
//! their outputs, exit codes, errors and warnings through named comparators
//! ([`execute`], [`comparators`], [`cases`]). It makes its cases from
//! seeded generators of random and malformed input ([`generators`],
//! [`sources`]). It runs them in parallel and records each as a line of
//! JSON ([`runner`]). It also fetches real datasets ([`fetch`]), builds
//! vsearch ([`vsearch_build`]) and times two builds ([`timing`]). The
//! `rscompare` binary is its command line.

pub mod options;
pub mod translate;

pub mod fetch;
pub mod generators;
pub mod rng;
pub mod timing;
pub mod vsearch_build;

pub mod cases;
pub mod comparators;
pub mod execute;
pub mod runner;
pub mod sources;
