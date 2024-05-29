use std::{f64::consts, vec};

/// An *n*-point Gauss-Legendre integration method,
/// using Newton-Raphson root-finding on Legendre polynomials.
/// Accuracy rises with the square of *n*. Don't go crazy.
///
/// ## Inputs
/// - f: The integrand.
/// - a: Lower bound of integration.
/// - b: Upper bound of integration.
/// - n: Number of points to interpolate.
/// ## Examples
/// ```rust
/// use symmetria::quadrature::gaussian_quad;
/// use std::f64::consts;
/// let f = |x: f64| x.sin();
/// let val = gaussian_quad(f, 0.0, consts::PI, 7);
/// assert!((2.0 - val).abs() < 1.0e-9);
/// ```
pub fn gaussian_quad<F>(f: F, a: f64, b: f64, n: u32) -> f64
where
    F: Fn(f64) -> f64,
{
    const NEWTON_PRECISION: f64 = 1.0e-10;
    let mut res = 0.0;
    // Scale the bounds of integration a and b to [-1, 1].
    // This converts
    // \int_a^b f(x)dx to \int_{-1}^1 f(sx + c)sdx
    // using the change of variables
    // x = st + c for t ∈ [-1, 1]
    // where
    // s = (b - a) / 2
    // and
    // c = (b + a) / 2.
    let s = (b - a) / 2.0;
    let c = (b + a) / 2.0;
    let m = (n + 1) / 2;
    let mut roots = vec![0.0; n as usize];
    let mut weights = vec![0.0; n as usize];
    for i in 0..m {
        // Approximate the ith root of this Legendre polynomial.
        let mut z = f64::cos(consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5));
        let mut z1 = z - 1.0; // Ensure an initial difference.
        let mut p_prime = 0.0;
        while (z - z1).abs() > NEWTON_PRECISION {
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
        // z is now the ith root of the polynomial.
        roots[i as usize] = z;
        roots[(n - 1 - i) as usize] = -z;
        weights[i as usize] = 2.0 / ((1.0 - z * z) * p_prime * p_prime);
        weights[(n - 1 - i) as usize] = weights[i as usize];
    }

    // \int_{-1}^1 f(sx + c)sdx = \sum_0^n w_i * f(x_i) where w_i is the weight and x_i
    // is the ith root of the nth Legendre polynomial.
    for i in 0..n {
        res += weights[i as usize] * f(s * roots[i as usize] + c) * s;
    }

    res
}

#[cfg(test)]
mod test {
    use super::gaussian_quad;

    #[test]
    fn check_gaussian_quad() {
        // Test case from here: <https://phys.libretexts.org/@go/page/8094?pdf>
        // \int_0^1 \frac{ x^4 }{ \sqrt{ 2(1+x) }} dx ≈ 0.108 709 465.
        let f = |x: f64| (x * x * x * x) / (2.0 * (1.0 + x * x)).sqrt();
        let v = gaussian_quad(f, 0.0, 1.0, 7);
        assert!((v - 0.108_709_465).abs() < 1.0e-9);
    }
}
