//! # Symmetria
//! 
//! Symmetria is a small toy implementation of useful numerical
//! algorithms, implemented without using complicated data types
//! that add to your workload. Instead, everything (excepting a few
//! cases) is simply an `f64`. 
//! 
//! Currently this crate provides:
//! - *n*-point Gauss-Legendre quadratures
//! - Adaptive 4/7/13-point Gauss-Lobatto-Kronrod quadrature,
//!   computing to machine precision
//! - Lagrange interpolation
//! - Bisection method root finding
//! - Cubic formula returning complex-valued roots
//! - Quadratic formula doing the same
//! 
//! ---
//! Eventually, I hope this crate will provide high coverage of the
//! routines implemented in the NSWC Library of Mathematics Subroutines,
//! (1990 ver.) published by the United States' Naval Warfare Center. It's
//! a well-established FORTRAN library, and that would be a fairly complete
//! selection. However, that's a lot of work.

pub mod functions;
pub mod interpolation;
pub mod quadrature;
