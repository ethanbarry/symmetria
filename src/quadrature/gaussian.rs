use std::{f64::consts, vec};

use crate::functions::batch_eval;

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
    for i in 0..n as usize {
        res += weights[i] * f(s * roots[i] + c) * s;
    }

    res
}

/// An adaptive 4-point Gauss-Lobatto quadrature routine.
/// Kronrod's Rule is used to reduce the computations needed,
/// while using a 13-point Kronrod quadrature for error checking.
///
/// ## Inputs
/// - f: The integrand (must be integrable on \[a,b\])
/// - a: The lower bound of integration.
/// - b: The upper bound of integration.
///
/// NOTE: The algorithm does not assume a < b.
/// ## Outputs
/// A result of `f64` or an error message inside a `String`.
/// The value should be an approximation of the integral within the given
/// tolerance.
/// ## Examples
/// ```rust
/// use symmetria::quadrature::gaussian_kronrod_quad;
/// use std::f64::consts;
/// let f = |x: f64| x.cos();
/// let val = gaussian_kronrod_quad(f, 0.0, (consts::PI / 2.0), 7); // Answer is first elem.
/// assert!((1.0 - val).abs() < 1.0e-9);
/// ```
pub fn gaussian_kronrod_quad<F>(func: F, a: f64, b: f64, n: u32) -> f64
where
    F: Fn(f64) -> f64,
{
    // This algorithm was adopted from the paper
    // "Adaptive Quadrature --- Revisited" by Walter Gander and Walter Gautschi.
    // published in BIT Numerical Mathematics, vol. 40, No. 1, pp. 84-101.

    // Scale the bounds of integration a and b to [-1, 1].
    // This converts
    // \int_a^b f(x)dx to \int_{-1}^1 f(sx + c)sdx
    // using the change of variables
    // x = scale * t + coscale for t ∈ [-1, 1]
    // where
    // scale = (b - a) / 2
    // and
    // coscale = (b + a) / 2.
    let scale = (b - a) / 2.0;
    let coscale = (b + a) / 2.0;
    // Now, f must always be evaluated at scale * x + coscale, not simply x, and the answer multiplied by scale.

    let gl4_weights = vec![1.0 / 6.0, 5.0 / 6.0];
    let gl4_fnevals = vec![
        func(scale * -1.0 + coscale),
        func(scale * 1.0 + coscale),
        func(scale * -1.0 / 5.0_f64.sqrt() + coscale),
        func(scale * 1.0 / 5.0_f64.sqrt() + coscale),
    ];
    // Next, use the 4-point Gauss-Lobatto quadrature given by
    // \int_{-1}^1 f(x)dx ≈ \frac{1}{6}[f(-1) + f(1)]
    //                      + \frac{5}{6}[f(-\frac{1}{\sqrt{5}}) + f(\frac{1}{\sqrt{5}})]
    let gl4 = gl4_weights[0] * ((gl4_fnevals[0] + gl4_fnevals[1]) * scale)
        + gl4_weights[1] * ((gl4_fnevals[2] + gl4_fnevals[3]) * scale);

    let kron7_fnevals = vec![
        func(scale * -((2. / 3.0_f64).sqrt()) + coscale),
        func(scale * (2. / 3.0_f64).sqrt() + coscale),
        func(coscale), // scale * 0.0 is 0.
    ];
    // Then we compare this to its 7-point Kronrod extension, given by
    // \int_{-1}^1 f(x)dx ≈ \frac{11}{210}[f(-1) + f(1)]
    //                      + \frac{72}{245}[f(-\sqrt{\frac{2}{3}}) + f(\sqrt{\frac{2}{3}})]
    //                      + \frac{125}{294}[f(-\frac{1}{\sqrt{5}}) + f(\frac{1}{\sqrt{5}})]
    //                      + \frac{16}{35}[f(0)]
    let kron7 = (11. / 210.) * ((gl4_fnevals[0] + gl4_fnevals[1]) * scale)
        + (72. / 245.) * ((kron7_fnevals[0] + kron7_fnevals[1]) * scale)
        + (125. / 294.) * ((gl4_fnevals[2] + gl4_fnevals[3]) * scale)
        + (16. / 35.) * (kron7_fnevals[2] * scale);

    // Finally, we compare all of them to the 13-point Kronrod extension, given by
    // \int_{-1}^1 f(x)dx ≈ A * [f(-1) + f(1)]
    //                      + B * [f(-x_1) + f(x_1)]
    //                      + C * [f(-\sqrt{\frac{2}{3}}) + f(\sqrt{\frac{2}{3}})]
    //                      + D * [f(-x_2) + f(x_2)]
    //                      + E * [f(-\frac{1}{\sqrt{5}}) + f(\frac{1}{\sqrt{5}})]
    //                      + F * [f(-x_3) + f(x_3)]
    //                      + G * [f(0)]
    // where
    // A = .015827191973480183087169986733305510591,
    // B = .094273840218850045531282505077108171960,
    // C = .15507198733658539625363597980210298680,
    // D = .18882157396018245442000533937297167125,
    // E = .19977340522685852679206802206648840246,
    // F = .22492646533333952701601768799639508076,
    // G = .24261107190140773379964095790325635233,
    // and
    // x_1 = .94288241569547971905635175843185720232,
    // x_2 = .64185334234578130578123554132903188354,
    // x_3 = .23638319966214988028222377349205292599.
    let a_weight = 0.015827191973480183087169986733305510591;
    let b_weight = 0.094273840218850045531282505077108171960;
    let c = 0.15507198733658539625363597980210298680;
    let d = 0.18882157396018245442000533937297167125;
    let e = 0.19977340522685852679206802206648840246;
    let f = 0.22492646533333952701601768799639508076;
    let g = 0.24261107190140773379964095790325635233;

    let x_1 = 0.94288241569547971905635175843185720232;
    let x_2 = 0.64185334234578130578123554132903188354;
    let x_3 = 0.23638319966214988028222377349205292599;

    let kron13 = a_weight * (gl4_fnevals[0] + gl4_fnevals[1]) * scale
        + b_weight * (func(scale * -x_1 + coscale) + func(scale * x_1 + coscale)) * scale
        + c * (kron7_fnevals[0] + kron7_fnevals[1]) * scale
        + d * (func(scale * -x_2 + coscale) + func(scale * x_2 + coscale)) * scale
        + e * (gl4_fnevals[2] + gl4_fnevals[3]) * scale
        + f * (func(scale * -x_3 + coscale) + func(scale * x_3 + coscale)) * scale
        + g * kron7_fnevals[2] * scale;

    // Now we need to:
    // - a: implement the helper function that does the recursion,
    // - b: connect that function to this one,
    // - c: and finally return the correct answer.

    if a < b {
        return kron13;
    } else {
        return -kron13;
    }
}

fn gkquad_recursive<F>(f: &F, a: f64, b: f64, fa: f64, fb: f64, iest: f64) -> f64
where
    F: Fn(f64) -> f64,
{
    let h = (b - a) / 2.0;
    let m = (b + a) / 2.0;
    let alpha = (2.0 / 3.0_f64).sqrt();
    let beta = 1.0 / 5.0_f64.sqrt();
    // Our new evaluation points...
    let mll = m - alpha * h;
    let ml = m - beta * h;
    let mr = m + beta * h;
    let mrr = m + alpha * h;
    let x = vec![mll, ml, m, mr, mrr];
    let y = batch_eval(&f, &x);

    let fmll = y[0];
    let fml = y[1];
    let fm = y[2];
    let fmr = y[3];
    let fmrr = y[4];

    let i2 = (h / 6.0) * (fa + fb + 5.0 * (fml + fmr));
    let i1 =
        (h / 1470.) * (77. * (fa + fb) + 432. * (fmll + fmrr) + 625. * (fml + fmr) + 672. * fm);

    let q = if iest + (i1 - i2) == iest || mll <= a || b <= mrr {
        if m <= a || b <= m {
            // At this point, we have exhausted the machine's precision.
            // We should handle this somehow...
        }
        i1
    } else {
        // Sum of the integrals on the ranges [a, mll], [mll, ml], [ml, m], [m, mr], [mr, mrr], [mrr, b];
        gkquad_recursive(f, a, mll, fa, fmll, iest)
            + gkquad_recursive(f, mll, ml, fmll, fml, iest)
            + gkquad_recursive(f, ml, m, fml, fm, iest)
            + gkquad_recursive(f, m, mr, fm, fmr, iest)
            + gkquad_recursive(f, mr, mrr, fmr, fmrr, iest)
            + gkquad_recursive(f, mrr, b, fmrr, fb, iest)
    };
    println!("{}", q);
    q
}

#[cfg(test)]
mod test {
    use super::{gaussian_kronrod_quad, gaussian_quad};

    #[test]
    fn check_gaussian_quad() {
        // Test case from here: <https://phys.libretexts.org/@go/page/8094?pdf>
        // \int_0^1 \frac{ x^4 }{ \sqrt{ 2(1+x) }} dx ≈ 0.108 709 465.
        let f = |x: f64| (x * x * x * x) / (2.0 * (1.0 + x * x)).sqrt();
        let v = gaussian_quad(f, 0.0, 1.0, 7);
        assert!((v - 0.108_709_465).abs() < 1.0e-9);
    }

    #[test]
    fn check_gaussian_kronrod_quad() {
        // Test case from here: <https://phys.libretexts.org/@go/page/8094?pdf>
        // \int_0^1 \frac{ x^4 }{ \sqrt{ 2(1+x) }} dx ≈ 0.108 709 465.
        let f = |x: f64| (x * x * x * x) / (2.0 * (1.0 + x * x)).sqrt();
        let v = gaussian_kronrod_quad(f, 1.0, 0.0, 7);
        assert!((v - 0.108_709_465).abs() < 1.0e-9);
    }
}
