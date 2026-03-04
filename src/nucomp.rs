//! NUCOMP and NUDUPL form composition.
//!
//! Port of chiavdf/src/nucomp.h (William Hart's algorithm).

use crate::form::Form;
use crate::integer::{divexact, fast_extended_gcd, fast_gcd_coeff_b, fdiv_q, fdiv_r, tdiv_r};
use crate::xgcd_partial::xgcd_partial;
use num_bigint::BigInt;
use num_traits::{Signed, Zero};

/// Compose two forms: result = f * g.
/// This is qfb_nucomp.
pub fn nucomp(f: &Form, g: &Form, d: &BigInt, l: &BigInt) -> Form {
    // Ensure a1 <= a2 (swap if needed)
    if f.a > g.a {
        return nucomp(g, f, d, l);
    }

    let a1 = f.a.clone();
    let a2 = g.a.clone();
    let c2 = g.c.clone();

    // ss = (f.b + g.b) / 2
    let ss = (&f.b + &g.b) >> 1usize;
    // m = (f.b - g.b) / 2
    let m = (&f.b - &g.b) >> 1usize;

    // Compute sp = gcd(a2 mod a1, a1) and v1
    let t = fdiv_r(&a2, &a1);
    let (sp, v1) = if t.is_zero() {
        (a1.clone(), BigInt::from(0u32))
    } else {
        let (gcd, x, _) = fast_extended_gcd(&t, &a1);
        (gcd, x)
    };

    // k = m * v1 mod a1
    let mut k = fdiv_r(&(&m * &v1), &a1);

    let (a1_new, a2_new, c2_new);

    if sp != BigInt::from(1u32) {
        let (s, v2, u2) = gcd_ext3(&ss, &sp);

        // k = k * u2 - v2 * c2
        k = &k * &u2 - &v2 * &c2;

        if s != BigInt::from(1u32) {
            a1_new = divexact(&a1, &s);
            a2_new = divexact(&a2, &s);
            c2_new = &c2 * &s;
        } else {
            a1_new = a1.clone();
            a2_new = a2.clone();
            c2_new = c2.clone();
        }

        k = fdiv_r(&k, &a1_new);
    } else {
        a1_new = a1.clone();
        a2_new = a2.clone();
        c2_new = c2.clone();
    }

    if a1_new < *l {
        // Short path: no partial GCD needed
        let t = &a2_new * &k;
        let ca = &a2_new * &a1_new;
        let cb = (&t << 1usize) + &g.b;
        let cc_num = (&g.b + &t) * &k + &c2_new;
        let cc = divexact(&cc_num, &a1_new);

        Form::new(ca, cb, cc)
    } else {
        // Long path: use xgcd_partial
        let mut r2 = a1_new.clone();
        let mut r1 = k;
        let mut co2 = BigInt::zero();
        let mut co1 = BigInt::zero();

        xgcd_partial(&mut co2, &mut co1, &mut r2, &mut r1, l);

        // m1 = (m * co1 + a2_new * r1) / a1_new
        let m1 = divexact(&(&m * &co1 + &a2_new * &r1), &a1_new);

        // m2 = (ss * r1 - c2_new * co1) / a1_new
        let m2 = divexact(&(&ss * &r1 - &c2_new * &co1), &a1_new);

        // ca = -sgn(co1) * (r1*m1 - co1*m2)
        // Note: ca may be negative here; cb must be computed before ca is made positive.
        let ca_unsigned = &r1 * &m1 - &co1 * &m2;
        let mut ca = if co1.is_negative() {
            ca_unsigned // -sgn(co1) = +1
        } else {
            -ca_unsigned // -sgn(co1) = -1
        };

        // t = a2_new * r1 (preserved for cb computation, matching C's local `t`)
        let t_val = &a2_new * &r1;

        // cb = (2*(t - ca*co2)/co1 - g.b) mod (2*ca)  [ca may be negative]
        let cb_inner = &t_val - &ca * &co2;
        let cb_scaled = &cb_inner << 1usize;
        let cb_divided = divexact(&cb_scaled, &co1);
        let cb_shifted = cb_divided - &g.b;
        let ca2 = &ca << 1usize; // may be negative; fdiv_r honors sign of modulus
        let cb = fdiv_r(&cb_shifted, &ca2);

        // cc = (cb^2 - D) / (4*ca)
        // Use 4*ca as the divisor (exact since b^2 ≡ D mod 4a for valid class-group forms).
        // ca may be negative here; fdiv_r already used the signed ca above.
        let cc_num = &cb * &cb - d;
        let cc_denom = &ca << 2usize; // 4*ca, may be negative
        let mut cc = divexact(&cc_num, &cc_denom);

        // Make ca positive; negate cc to keep the form equation b^2 - 4ac = D valid.
        if ca.is_negative() {
            ca = -ca;
            cc = -cc;
        }

        Form::new(ca, cb, cc)
    }
}

/// Extended GCD returning (gcd, coeff_a, coeff_b) where gcd = coeff_a * a + coeff_b * b.
fn gcd_ext3(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    fast_extended_gcd(a, b)
}

/// Duplicate a form: result = f^2.
/// This is qfb_nudupl.
pub fn nudupl(f: &Form, d: &BigInt, l: &BigInt) -> Form {
    let a1 = f.a.clone();
    let c1 = f.c.clone();

    // s = gcd(|b|, a), v2 = coefficient for b in gcd = v2*|b| + ?*a
    // Use half-GCD since we only need the coefficient of |b|, not both Bezout coefficients.
    let b_abs = f.b.abs();
    let (s, v2) = {
        let (gcd, coeff_b) = fast_gcd_coeff_b(&a1, &b_abs);
        // fast_gcd_coeff_b(a1, b_abs) returns coeff of b_abs in gcd = ?*a1 + coeff_b*b_abs
        let v2 = if f.b.is_negative() { -coeff_b } else { coeff_b };
        (gcd, v2)
    };

    // k = -(c * inv(b)) mod a = -c * v2 mod a
    let k_raw = -&c1 * &v2;
    // Use truncating remainder and fix up negative results
    let mut k = tdiv_r(&k_raw, &a1);
    if k.is_negative() {
        k += &a1;
    }

    let a1_new;
    let c1_new;

    let s_is_1 = s == BigInt::from(1u32);
    if !s_is_1 {
        a1_new = fdiv_q(&a1, &s);
        c1_new = &c1 * &s;
    } else {
        a1_new = a1.clone();
        c1_new = c1.clone();
    }

    if a1_new < *l {
        // Short path
        let t = &a1_new * &k;
        let new_a = &a1_new * &a1_new;
        let cb = (&t << 1usize) + &f.b;
        let cc_num = (&f.b + &t) * &k + &c1_new;
        let cc = fdiv_q(&cc_num, &a1_new);

        Form::new(new_a, cb, cc)
    } else {
        // Long path: xgcd_partial
        let mut r2 = a1_new.clone();
        let mut r1 = k;
        let mut co2 = BigInt::zero();
        let mut co1 = BigInt::zero();

        xgcd_partial(&mut co2, &mut co1, &mut r2, &mut r1, l);

        // m2 = (b * r1 - c1_new * co1) / a1_new
        let m2_num = &f.b * &r1 - &c1_new * &co1;
        let m2 = divexact(&m2_num, &a1_new);

        // new_a = r1^2 - co1*m2, then negate if co1 >= 0
        // (matches C: mpz_submul(r->a, co1, m2); if sgn(co1)>=0: neg(r->a))
        let mut new_a = &r1 * &r1 - &co1 * &m2;
        if !co1.is_negative() {
            new_a = -new_a;
        }
        // new_a may be negative here — keep it signed for cb and cc computation

        // cb = 2*(a1*r1 - new_a*co2)/co1 - f.b  (mod 2*new_a)
        // Mirrors C: cb=new_a*co2; submul(a1*r1); neg; *=2; divexact(co1); sub(b); fdiv_r(2*new_a)
        let cb_tmp = &new_a * &co2 - &a1_new * &r1; // = new_a*co2 - a1*r1
        let cb_neg = -cb_tmp; // = a1*r1 - new_a*co2
        let cb_doubled = &cb_neg << 1usize;
        let cb_div = divexact(&cb_doubled, &co1);
        let cb_pre = cb_div - &f.b;
        // fdiv_r mod 2*new_a (new_a may be negative — GMP fdiv_r honors sign of modulus)
        let two_new_a = &new_a << 1usize;
        let cb = fdiv_r(&cb_pre, &two_new_a);

        // cc = (cb^2 - D) / new_a  (divexact, new_a may be negative)
        //   then >> 2 (truncating, per mpz_tdiv_q_2exp)
        let cc_num = &cb * &cb - d;
        let cc_pre = divexact(&cc_num, &new_a);
        let cc = &cc_pre >> 2usize;

        // Fix signs: if new_a < 0, negate both a and c
        let (final_a, final_c) = if new_a.is_negative() {
            (-new_a, -cc)
        } else {
            (new_a, cc)
        };

        Form::new(final_a, cb, final_c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::form::Form;
    use num_bigint::BigInt;

    fn discriminant_ok(f: &Form, d: &BigInt) -> bool {
        let disc = &f.b * &f.b - BigInt::from(4) * &f.a * &f.c;
        &disc == d
    }

    #[test]
    fn test_nucomp_preserves_discriminant() {
        let d = BigInt::from(-47i64);
        let l = Form::compute_l(&d);
        let f = Form::new(BigInt::from(2), BigInt::from(1), BigInt::from(6));
        let g = Form::new(BigInt::from(3), BigInt::from(1), BigInt::from(4));
        assert!(discriminant_ok(&f, &d));
        assert!(discriminant_ok(&g, &d));
        let result = nucomp(&f, &g, &d, &l);
        assert!(
            discriminant_ok(&result, &d),
            "nucomp result has wrong discriminant: a={}, b={}, c={}",
            result.a,
            result.b,
            result.c
        );
    }

    #[test]
    fn test_nudupl_preserves_discriminant() {
        let d = BigInt::from(-47i64);
        let l = Form::compute_l(&d);
        let f = Form::new(BigInt::from(2), BigInt::from(1), BigInt::from(6));
        assert!(discriminant_ok(&f, &d));
        let result = nudupl(&f, &d, &l);
        assert!(
            discriminant_ok(&result, &d),
            "nudupl result has wrong discriminant: a={}, b={}, c={}",
            result.a,
            result.b,
            result.c
        );
    }
}
