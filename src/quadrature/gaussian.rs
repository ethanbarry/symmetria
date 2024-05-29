use std::{f64::consts, vec};

/// An *n*-point Gauss-Lagrange integration method,
/// using Newton-Raphson root-finding for the Lagrange polynomials.
/// Accuracy rises with the square of n. Choose a low value.
pub fn gaussian_quad<F>(f: F, a: f64, b: f64, n: u32) -> f64
where
    F: Fn(f64) -> f64,
{
    const EPS: f64 = 1.0e-10; // TODO: Adjust this as needed.
    let mut res = 0.0;
    // Scale the bounds of integration a and b to [-1, 1].
    // This converts
    // \int_a^b f(t)dt to \int_{-1}^1 f(x)dx
    // using the change of variables
    // x = mt + c for t ∈ [-1, 1]
    // where
    // m = (b - a) / 2
    // and
    // c = (b + a) / 2.
    let m = (b - a) / 2.0;
    let c = (b + a) / 2.0;
    // We have converted to the form
    // \int_{-1}^1 f(mx + c)mdx.
    // Now all we need to do is use this approximation
    // \int_{-1}^1 f(x)dx ≈ \sum_{i=1}^n c_i * f(x_i)
    // to calculate the integral. We can represent f(x) by using
    // Lagrange polynomials like so
    // y = \sum_{i=1}^n f(x_i) * L_i(x).
    //
    // Now, we can see that the coefficient c_i in the approximation
    // is actually \int_{-1}^1 L_i(x) dx.
    // Next, we'll find the roots of the Legendre polynomials.
    let m = (n + 1) / 2;
    let mut roots = vec![0.0; n as usize];
    let mut weights = vec![0.0; n as usize];
    for i in 0..m {
        // Approximate the ith root of this Legendre polynomial.
        let mut z = f64::cos(consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5));
        let mut z1 = z - 1.0; // Ensure difference.
        let mut p_prime = 0.0;
        while (z - z1).abs() > EPS {
            let mut p1 = 1.0;
            let mut p2 = 0.0;
            for j in 0..n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j as f64 + 1.0) * z * p2 - j as f64 * p3) / (j + 1) as f64;
            }
            // p1 is now the desired Legendre polynomial.
            // let p_prime be its derivative, computed by a known relation.
            p_prime = n as f64 * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / p_prime; // Newton-Raphson.
        }
        dbg!(&z);
        // z is now the ith root of the polynomial.
        roots[i as usize] = z;
        roots[(n - 1 - i) as usize] = -z;
        weights[i as usize] = 2.0 / ((1.0 - z * z) * p_prime * p_prime);
        weights[(n - 1 - i) as usize] = weights[i as usize];
        dbg!(&weights[i as usize]);
    }

    // \int_{-1}^1 f(x)dx = \sum_0^n w_i * f(x_i) where w_i is the weight and x_i
    // is the ith root of the nth Legendre polynomial.
    for i in 0..n {
        res += weights[i as usize] * f(roots[i as usize]);
    }

    dbg!(&res);
    res
}

#[cfg(test)]
mod test {
    use super::gaussian_quad;

    #[test]
    fn check_gaussian_quad() {
        let f = |x: f64| x.cos();
        let v = gaussian_quad(f, 0.0, 2.0, 3);
        assert!(v > 0.0);
    }
}
