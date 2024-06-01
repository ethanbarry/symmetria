// Basic Lagrange implementation to interpolate points.

/// A hilariously-bad version of the Lagrange interpolation process.
/// We can improve this by precomputing and reusing the table of points,
/// instead of calling it every time, and using Neville's algorithm
/// for the interpolation process.
///
/// ## Inputs
/// - `f` is an integrand of type `Fn(f64) -> f64`
/// - `x` is the value to interpolate
/// - `a` is the lower bound
/// - `b` is the upper bound
/// - `steps` is the number of points to use in interpolation.
///
/// ## Conditions
/// Inputs must satisfy the following:
/// - a < b
pub fn lagrange_interpolate<F>(f: F, x: f64, a: f64, b: f64, steps: u32) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    if a >= b {
        return Err(String::from("Input bounds invalid - a greater than b."));
    }
    let table = fn_eval(f, a, b, steps);

    let n = table.len();
    let mut sum = 0.0;

    for i in 0..n {
        if let Some((x_i, y_i)) = table.get(i) {
            // Evaluate the ith Lagrange poly'n. at xval and add it to sum.
            let mut product = 1.0;

            // Iterator grabs all but i.
            for j in (0..n)
                .enumerate()
                .filter(|&(pos, _)| (pos != i))
                .map(|(_, e)| e)
            {
                let x_j = table.get(j).unwrap().0;
                product *= (x - x_j) / (x_i - x_j);
            } // product is now the value of L_i(xval).

            sum += y_i * product;
        }
    }
    Ok(sum)
}

/// Generates a table of points to use in the interpolation.
fn fn_eval<F>(f: F, a: f64, b: f64, steps: u32) -> Vec<(f64, f64)>
where
    F: Fn(f64) -> f64,
{
    let mut res: Vec<(f64, f64)> = vec![];
    let delta_x = (b - a) / steps as f64;

    for i in 0..steps {
        let x_val = a + delta_x * i as f64;
        let y_val = f(x_val);
        res.push((x_val, y_val));
    }

    res
}

#[cfg(test)]
mod test {
    use super::fn_eval;
    use super::lagrange_interpolate;

    #[test]
    fn check_fn_eval() {
        let f = |x: f64| x + 2_f64;
        let a = 0_f64;
        let b = 3_f64;
        let expected = vec![(0.0, 2.0), (1.0, 3.0), (2.0, 4.0)];
        assert_eq!(expected, fn_eval(f, a, b, 3));
    }

    #[test]
    fn check_lagrange_interpolate() {
        let f = |x: f64| x.sin();
        let a = 0.0;
        let b = 0.5;
        let x = 0.2;
        let steps = 4;
        let res = lagrange_interpolate(f, x, a, b, steps).unwrap();

        assert!((res - x.sin()).abs() <= 1e-3);
    }
}
