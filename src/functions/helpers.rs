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
    batch.iter().map(|val| f(*val)).collect()
}

/// A String-wrapping error type for iterative root-finding
/// methods, such as the bisection method.
#[non_exhaustive]
#[derive(Debug)]
pub enum RootFindError {
    TooManyIterations(String),
    NonBoundedInputs(String),
}

/// Return a root of a function given known positive and
/// negative bounding values, and the desired precision.
/// ## Inputs
/// - f: The function f(x) we wish to solve.
/// - pos & neg: The x-values of the positive and negative
///   bounds on the root.
/// - precision: The precision to which the root will be evaluated.
///   Usually something like 1.0e-9.
/// - max_iter: The maximum number of iterations to carry out.
fn bisection<F>(
    f: F,
    pos: f64,
    neg: f64,
    precision: f64,
    max_iter: u16,
) -> Result<f64, RootFindError>
where
    F: Fn(f64) -> f64,
{
    let mut fpos = f(pos);
    let fneg = f(neg);
    let mut dx = 0.0;
    let mut xmid = 0.0;

    // Handle invalid input.
    if fpos.is_sign_negative() && fneg.is_sign_negative()
        || fpos.is_sign_positive() && fneg.is_sign_positive()
    {
        return Err(RootFindError::NonBoundedInputs(
            "Inputs `pos` and `neg` must have different signs to bracket a root.".to_string(),
        ));
    }

    let mut rtb = if fneg < 0.0 {
        dx = pos - neg;
        neg
    } else {
        dx = neg - pos;
        pos
    };

    for i in 0..max_iter {
        dx *= 0.5;
        xmid = rtb + dx;
        fpos = f(xmid);
        if fpos <= 0.0 {
            rtb = xmid;
        }
        if dx.abs() < precision || fpos == 0.0 {
            return Ok(rtb);
        }
    }

    Err(RootFindError::TooManyIterations(format!(
        "Did not converge in {max_iter} iterations."
    )))
}

#[cfg(test)]
mod test {
    use crate::functions::helpers::bisection;
    #[test]
    fn bisection_success() {
        let f = |x: f64| x * x + x - 2.0;
        let rt1 = bisection(f, 1.3, 0.8, 1.0e-9, 50).unwrap();
        let rt2 = bisection(f, -3.0, -1.0, 1.0e-9, 50).unwrap();
        assert!((rt1 - 1.0).abs() < 1.0e-9);
        assert!((rt2 - -2.0).abs() < 1.0e-9);
    }
}
