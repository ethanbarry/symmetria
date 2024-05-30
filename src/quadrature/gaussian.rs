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

/// An *n*-point adaptive Gaussian quadrature routine which
/// returns an estimate of the error along with the result.
/// Kronrod's Rule is used to reduce the computations needed,
/// while using a 2n+1-point Kronrod quadrature.
///
/// ## Inputs
/// - f: The integrand (must be integrable on \[a,b\])
/// - a: The lower bound of integration.
/// - b: The upper bound of integration.
/// - n: The order of the Gaussian quadrature.
/// ## Outputs
/// A result of `(f64, f64)` or an error message inside a `String`.
/// The tuple is (answer, error) where error is the absolute difference
/// between the Gaussian *n*-point and Kronrod 2n+1-point results.
/// ## Errors
/// - When Newton-Raphson iteration cannot reliably find the Gaussian or Kronrod abscissae.
/// ## Examples
/// ```rust
/// use symmetria::quadrature::gaussian_kronrod_quad;
/// use std::f64::consts;
/// let f = |x: f64| x.cos();
/// let val = gaussian_kronrod_quad(f, 0.0, (consts::PI / 2.0), 9).unwrap().0; // Answer is first elem.
/// assert!((0.0 - val).abs() < 1.0e-9);
/// ```
pub fn gaussian_kronrod_quad<F>(f: F, a: f64, b: f64, n: u32) -> Result<(f64, f64), String>
where
    F: Fn(f64) -> f64,
{
    let m = (n + 1) / 2;
    let even = if m % 2 == 0 { true } else { false };
    // The Gauss-Kronrod rule will include 2n+1 points.  However, by symmetry,
    // only n + 1 of them need to be listed.
    //
    // The arrays abscissa, kweights and gweights contain the nonnegative
    // abscissas in decreasing order, and the weights of each abscissa in
    //  the Gauss-Kronrod and Gauss rules respectively.  This means that
    // about half the entries in W2 are zero.

    // First task, compute the abscissae.
    let mut cheby = vec![0.0_f64; (m + 1) as usize]; // `b` in jburkardt's Python.
    let mut tau = vec![0.0; m as usize];
    let mut kweights = vec![0.0; (n + 1) as usize];
    let mut gweights = vec![0.0; (n + 1) as usize];
    let mut abscissa = vec![0.0; (n + 1) as usize];

    let mut d = 2.0;
    let mut a_n = 0.0;
    for _ in 1..n + 1 {
        a_n = a_n + 1.0;
        d = d * a_n / (a_n + 0.5);
    }
    // Calculating Chebyshev coefficients of orthogonal poly'n.
    tau[0] = (a_n + 2.0) / (a_n + a_n + 3.0);
    cheby[(m - 1) as usize] = tau[0] - 1.0;
    let mut a_k = a_n;

    for l in 1..m {
        a_k += 2.0;
        tau[l as usize] =
            ((a_k - 1.0) * a_k - a_n * (a_n + 1.0)) * (a_k + 2.0) * tau[(l - 1) as usize]
                / (a_k * ((a_k + 3.0) * (a_k + 2.0) - a_n * (a_n + 1.0)));
        cheby[(m - l - 1) as usize] = tau[(l + 1 - 1) as usize];

        for ll in 1..l + 1 {
            cheby[(m - l - 1) as usize] = cheby[(m - l - 1) as usize]
                + tau[(ll - 1) as usize] * cheby[(m - l + ll - 1) as usize];
        }
    }
    cheby[m as usize] = 1.0;

    // Now calculate approximate values for the abscissae.
    let mut bb = (0.5 * consts::PI / (a_n + a_n + 1.0)).sin();
    let mut x1 = (1.0 - bb * bb).sqrt();
    let s = 2.0 * bb * x1;
    let c = (1.0 - s * s).sqrt();
    let coef = 1.0 - (1.0 - 1.0 / a_n) / (8.0 * a_n * a_n);
    let mut xx = coef * x1;
    // Coefficient needed for weights.
    let mut coef2 = 2.0 / (2.0 * n as f64 + 1.0);
    for i in 1..n + 1 {
        coef2 = coef2 * 4.0 * i as f64 / (n + i) as f64;
    }

    // Calculate the kth Kronrod abscissa and corresponding weight.
    for k in (1..n + 1).step_by(2) {
        (xx, kweights[(k - 1) as usize]) = abwe1(n, m, 1.0e-9, coef2, even, &cheby, xx)?;
        gweights[(k - 1) as usize] = 0.0;

        abscissa[(k - 1) as usize] = xx;
        let mut y = x1;
        x1 = y * c - bb * s;
        bb = y * s + bb * c;

        xx = if k == n { 0.0 } else { coef * x1 };

        // Calculate the k+1th abscissa and its weight.
        (xx, kweights[k as usize], gweights[k as usize]) =
            abwe2(n, m, 1.0e-9, coef2, even, &cheby, xx)?;
        abscissa[k as usize] = xx;
        y = x1;
        x1 = y * c - bb * s;
        bb = y * s + bb * c;
        xx = coef * x1;
    }

    // If n was even, we have one more Kronrod abscissa to compute -- the origin.
    if even {
        xx = 0.0;
        (xx, kweights[n as usize]) = abwe1(n, m, 1.0e-9, coef2, even, &cheby, xx)?;
        gweights[n as usize] = 0.0;
        abscissa[n as usize] = xx;
    }

    // The vectors abscissa, kweights, and gweights now contain the correct values for
    // the interval [-1, 1]. Now we need to scale & sum.
    // Remember that the Gaussian quadrature formula is
    // \int_a^b f(x)dx = \sum_1^n w_i * f(x_i)
    for i in 0..=n {
        abscissa[i as usize] =
            ((1.0 - abscissa[i as usize]) * a + (1.0 + abscissa[i as usize]) * b) / 2.0;
        kweights[i as usize] = ((b - a) / 2.0) * kweights[i as usize];
        gweights[i as usize] = ((b - a) / 2.0) * gweights[i as usize];
    }

    // Now we'll compute the estimates.
    let mut kronrod_integral = kweights[n as usize] * f(abscissa[n as usize]);
    let mut gauss_integral = gweights[n as usize] * f(abscissa[n as usize]);

    for i in 0..n {
        gauss_integral +=
            gweights[i as usize] * (f(-abscissa[i as usize]) + f(abscissa[i as usize]));
        kronrod_integral +=
            kweights[i as usize] * (f(-abscissa[i as usize]) + f(abscissa[i as usize]));
    }

    // The last steps here are summing scaled weights times f(scaled roots),
    // and subtracting the Gaussian result from the Gauss-Kronrod result.
    let err = (gauss_integral - kronrod_integral).abs();
    if err > 1.0e-9 {
        dbg!(&gauss_integral, &kronrod_integral, &err);

        Err(format!(
            "Bad estimate with {n} points -- try a larger number.",
        ))
    } else {
        Ok((kronrod_integral, err))
    }
}

/// Helper function for Kronrod Rule quadrature.
/// Calculates a Kronrod abscissa and its weight.
/// ## Input
/// - n: Order of the Gaussian quadrature.
/// - m: Order of the Kronrod quadrature (2n + 1).
/// - tol: Desired absolute accuracy of the abscissae.
/// - coef2: Value needed to compute weights.
/// - even: `true` if n is even.
/// - cheby: A `vec![f64; m + 1]` of Chebyshev coefficients.
/// - xx: An estimate of the abscissa.
/// ## Outputs
/// Result of `(f64, f64)` with form (abscissa, weight).
/// ## Errors
/// - When Newton-Raphson iteration cannot reliably find the Gaussian or Kronrod abscissae.
fn abwe1(
    n: u32,
    m: u32,
    tol: f64,
    coef2: f64,
    even: bool,
    b: &[f64],
    xx: f64,
) -> Result<(f64, f64), String> {
    let mut ka = if xx == 0.0 { 1 } else { 0 };
    let mut delta = 0.0;
    let mut x = xx;
    let mut fd = 0.0;

    let mut b1 = 0.0;
    let mut b0 = 0.0;
    let mut b2 = b[m as usize];
    let mut yy = 4.0 * x * x - 2.0;
    let mut d1 = 0.0;

    let mut ai = 0.0;
    let mut d2 = 0.0;
    let mut dif = 0.0;
    let mut i = 0;

    // Iterative process to find the abscissa.
    for iter in 1..51 {
        b1 = 0.0;
        b2 = b[m as usize];
        yy = 4.0 * x * x - 2.0;
        d1 = 0.0;

        if even {
            ai = (m + m + 1) as f64;
            d2 = ai as f64 * b[m as usize];
            dif = 2.0;
        } else {
            ai = (m + 1) as f64;
            d2 = 0.0;
            dif = 1.0;
        }

        for k in 1..m + 1 {
            ai = ai - dif;
            i = m - k + 1;
            b0 = b1;
            b1 = b2;
            let d0 = d1;
            d1 = d2;
            b2 = yy * b1 - b0 + b[(i - 1) as usize];

            if !even {
                i += 1;
            }
            d2 = yy * d1 - d0 + ai * b[(i - 1) as usize];
        }

        let mut f = 0.0;
        if even {
            f = x * (b1 - b2);
            fd = d2 + d1;
        } else {
            f = 0.5 * (b2 - b0);
            fd = 4.0 * x * d2;
        }
        // Newton correction...
        delta = f / fd;
        x = x - delta;
        if ka == 1 {
            break;
        }
        if delta.abs() <= tol {
            ka = 1;
        }
    }
    // Handle non-convergence.
    if ka != 1 {
        return Err(format!(
            "Function abwe1() - Fatal error.\nIteration limit reached.\nLast delta was {}",
            delta
        ));
    }

    // We're OK; compute the weight.
    let mut d0 = 1.0;
    d1 = x;
    ai = 0.0;
    for k in 2..n + 1 {
        ai = ai + 1.0;
        d2 = ((ai + ai + 1.0) * x * d1 - ai * d0) / (ai + 1.0);
        d0 = d1;
        d1 = d2;
    }

    let w = coef2 / (fd * d2);

    Ok((x, w))
}

/// Helper function for Kronrod Rule quadrature.
/// Calculates a Gaussian abscissa and its weight.
/// ## Input
/// - n: Order of the Gaussian quadrature.
/// - m: Order of the Kronrod quadrature (2n + 1).
/// - tol: Desired absolute accuracy of the abscissae.
/// - coef2: Value needed to compute weights.
/// - even: `true` if n is even.
/// - cheby: A `vec![f64; m + 1]` of Chebyshev coefficients.
/// - xx: An estimate of the abscissa.
/// ## Outputs
/// Result of `(f64, f64, f64)` with form (abscissa, Gauss-Kronrod weight, Gauss weight).
/// ## Errors
/// - When Newton-Raphson iteration cannot reliably find the Gaussian or Kronrod abscissae.
fn abwe2(
    n: u32,
    m: u32,
    tol: f64,
    coef2: f64,
    even: bool,
    b: &[f64],
    mut x: f64,
) -> Result<(f64, f64, f64), String> {
    let mut ka = if x == 0.0 { 1 } else { 0 };
    let mut delta = 0.0;
    let mut p0 = 1.0;
    let mut p1 = x;
    let mut pd0 = 0.0;
    let mut pd1 = 1.0;
    let mut p2 = 0.0;
    let mut pd2 = 0.0;
    // Iterative process for computing a Gaussian abscissa.
    for _ in 1..51 {
        p0 = 1.0;
        p1 = x;
        pd0 = 0.0;
        pd1 = 1.0;

        // If n is one, initialize p2 and pd2 to avoid problems.
        if n <= 1 {
            if x != 0.0 {
                p2 = (3.0 * x * x - 1.0) / 2.0;
                pd2 = 3.0 * x;
            } else {
                p2 = 3.0 * x;
                pd2 = 3.0;
            }
        }

        let mut ai = 0.0;
        for k in 2..n + 1 {
            ai += 1.0;
            p2 = ((ai + ai + 1.0) * x * p1 - ai * p0) / (ai + 1.0);
            pd2 = ((ai + ai + 1.0) * (p1 + x * pd1) - ai * pd0) / (ai + 1.0);
            p0 = p1;
            p1 = p2;
            pd0 = pd1;
            pd1 = pd2;
        }

        // Newton correction.
        delta = p2 / pd2;
        x = x - delta;

        if ka == 1 {
            break;
        }
        if delta.abs() <= tol {
            ka = 1;
        }
    }

    // Handle non-convergence.
    if ka != 1 {
        return Err(format!(
            "Function abwe2() - Fatal error.\nIteration limit reached.\nLast delta was {}",
            delta
        ));
    }

    let an = n;

    let w2 = 2.0 / (an as f64 * pd2 * p0);

    p1 = 0.0;
    p2 = b[m as usize];
    let yy = 4.0 * x * x - 2.0;
    for k in 1..m + 1 {
        let i = m - k + 1;
        p0 = p1;
        p1 = p2;
        p2 = yy * p1 - p0 + b[(i - 1) as usize];
    }

    let w1 = if even {
        w2 + coef2 / (pd2 * x * (p2 - p1))
    } else {
        w2 + 2.0 * coef2 / (pd2 * (p2 - p0))
    };

    Ok((x, w1, w2))
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
