/// Evaluate a one-dimensional function at a slice of points.
/// ## Inputs
/// - f: A function to be evaluated.
/// - batch: A slice of `f64` values at which to evaluate f.
/// ## Outputs
/// A `Vec<f64>` containing the evaluated results.
/// ```rust
/// use symmetria::functions::batch_eval;
/// use std::f64::consts;
/// let f = |x: f64| x + 1.0;
/// let x_coords = vec![0.0, 0.5, 1.0];
/// let y_coords = batch_eval(&f, &x_coords);
/// assert_eq!(vec![1.0, 1.5, 2.0], y_coords);
/// ```
pub fn batch_eval<F>(f: &F, batch: &[f64]) -> Vec<f64>
where
    F: Fn(f64) -> f64,
{
    batch.into_iter().map(|val| f(*val)).collect()
}
