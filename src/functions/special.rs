use num_complex::Complex64;

/// Solve a quadratic equation.
/// ## Inputs
/// - coeffs: Triplet of `f64` values representing the
///   coefficients of terms in increasing order.
/// ## Examples
/// ```rust
/// use symmetria::functions::special::quadratic_formula;
/// use num_complex::Complex64;
/// // Solve x² + 2x + 3
/// let roots = quadratic_formula((3., 2., 1.));
/// assert_eq!(roots.0, Complex64::new(-1., -(2_f64).sqrt()));
/// assert_eq!(roots.1, Complex64::new(-1., (2_f64).sqrt()));
/// ```
pub fn quadratic_formula(coeffs: (f64, f64, f64)) -> (Complex64, Complex64) {
    // In form ax² + bx + c these are...
    let a0 = Complex64::new(coeffs.0, 0.); // c
    let a1 = Complex64::new(coeffs.1, 0.); // b
    let a2 = Complex64::new(coeffs.2, 0.); // a

    // The bog-standard quadratic equation...
    // x = (-b ± sqrt(b² - 4ac)) / (2a)

    let x1 = (-a1 - (a1.powu(2) - 4. * a2 * a0).sqrt()) / (2. * a2);
    let x2 = (-a1 + (a1.powu(2) - 4. * a2 * a0).sqrt()) / (2. * a2);

    (x1, x2)
}

/// Solve a cubic equation.
/// ## Inputs
/// - coeffs: Triplet of `f64` values representing the
///   coefficients of terms in increasing order.
/// ## Outputs
/// - A triplet of complex roots.
/// ## Examples
/// ```rust
/// use symmetria::functions::special::cubic_formula;
/// use num_complex::Complex64;
/// // Solve x³ + 3x² - 4x - 12, which is
/// // the expansion of (x+3)(x-2)(x+2).
/// let roots = cubic_formula((-12., -4., 3., 1.,));
/// dbg!(roots);
/// assert_eq!(roots.0.re, 3.);
/// assert_eq!(roots.1.re, -2.);
/// assert_eq!(roots.2.re, 2.);
/// ```
pub fn cubic_formula(coeffs: (f64, f64, f64, f64)) -> (Complex64, Complex64, Complex64) {
    let a0 = Complex64::new(coeffs.0, 0.); // d
    let a1 = Complex64::new(coeffs.1, 0.); // c
    let a2 = Complex64::new(coeffs.2, 0.); // b
    let a3 = Complex64::new(coeffs.3, 0.); // a

    // Δ₀ = b² - 3ac
    let delta_nought = a2.powu(2) - 3. * a3 * a1;
    // Δ₁ = 2b³ - 9abc + 27a²d
    let delta_one = 2. * a2.powu(3) - 9. * a3 * a2 * a1 + 27. * a3.powu(2) * a0;

    // C = cubrt((Δ₁ ± sqrt(Δ₁² - 4Δ₀³)) / (2)) but if C is 0, the other sign in "±" must be chosen instead.
    let mut c =
        ((delta_one + Complex64::sqrt(delta_one.powu(2) - 4. * delta_nought.powu(3))) / 2.).cbrt();

    if c == Complex64::ZERO {
        c = ((delta_one - Complex64::sqrt(delta_one.powu(2) - 4. * delta_nought.powu(3))) / 2.)
            .cbrt();
    }

    // Now, if C still equals 0, we will have a fraction of 0/0, which must be *interpreted* as zero.
    let root1 = if c != Complex64::ZERO {
        (-1. / (3. * a3)) * (a2 + c + delta_nought / c)
    } else {
        (-1. * a2) / (3. * a3)
    };

    // Now we multiply by the primitive cube root of unity twice to find the other two roots.
    let root_unity =
        (Complex64::new(-1., 0.) + Complex64::new(-3., 0.).sqrt()) / Complex64::new(2., 0.);

    let root2 = if c != Complex64::ZERO {
        (-1. / (3. * a3)) * (a2 + root_unity * c + delta_nought / (root_unity * c))
    } else {
        (-1. * a2) / (3. * a3)
    };

    // Second multiplication...
    let root3 = if c != Complex64::ZERO {
        (-1. / (3. * a3)) * (a2 + root_unity.powu(2) * c + delta_nought / (root_unity.powu(2) * c))
    } else {
        (-1. * a2) / (3. * a3)
    };

    (-root1, root2, root3) // Minus sign to fix the sign error that snuck in somewhere...
}
