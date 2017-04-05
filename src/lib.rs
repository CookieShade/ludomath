//! A crate containing various math & number functions,
//! intended for light 2D graphics & game programming.
//!
//! Includes utilities for cheap random number generation,
//! 2D vector math & transformation, as well as a few
//! basic numeric functions.

#![warn(missing_docs)]

extern crate rand;

pub mod vec2d;
pub mod num;
pub mod rng;
pub mod consts;

#[cfg(test)]
mod tests;
