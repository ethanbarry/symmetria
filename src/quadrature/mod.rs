mod gaussian;
mod interpolation;

// Expose methods to callers outside this module.
pub use gaussian::gaussian_kronrod_quad;
pub use gaussian::gaussian_quad;
