# Symmetria

I intend for this crate to become a useful source of ready-to-use numerical algorithms.
Many crates implement their own confusing data types, which take time and mental effort
to learn. Instead of abstracting, this crate will use simple `f64` variables, and will
duplicate work if/when `f32` versions are needed.

Currently, this crate provides the following procedures:
- Gaussian quadrature
- Gaussian-Lobatto-Kronrod adaptive quadrature
- Batch evaluation of functions

I hope to eventually provide all the routines found in the NSWC Library of Mathematics Subroutines,
published by the United States Naval Surface Warfare Center (1990 version), which is an established
FORTRAN library. That would be a complete selection. However, that's an enormous amount of work.

Don't hold your breath. :-)
