use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::math::core::{LinearSpace, ScalarSpace};

pub trait ComplexSpace: ScalarSpace {
    type Real: ScalarSpace;

    fn i() -> Self;
    fn from_real(re: Self::Real) -> Self;
    fn from_polar(r: Self::Real, phase: Self::Real) -> Self;
    fn cis(phase: Self::Real) -> Self;
    fn re(&self) -> Self::Real;
    fn im(&self) -> Self::Real;
    fn phase(&self) -> Self::Real;
    fn arg(&self) -> Self::Real {
        self.phase()
    }
}

#[derive(Clone, Copy, PartialEq)]
#[repr(C)]
pub struct Complex<S: ScalarSpace> {
    re: S,
    im: S,
}

impl<S: ScalarSpace> Complex<S> {
    pub fn new(re: S, im: S) -> Self {
        Self { re: re, im: im }
    }
}

impl<S: ScalarSpace> LinearSpace for Complex<S> {
    type Data = [S; 2];

    fn new(data: Self::Data) -> Self {
        Self {
            re: data[0],
            im: data[1],
        }
    }
    fn zero() -> Self {
        Self {
            re: S::zero(),
            im: S::zero(),
        }
    }
    fn size(&self) -> usize {
        1
    }
    fn get_data(&self) -> Self::Data {
        [self.re, self.im]
    }
}

impl ScalarSpace for Complex<f64> {
    fn one() -> Self {
        Self { re: 1.0, im: 0.0 }
    }
    fn abs_square(&self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    fn abs(&self) -> f64 {
        self.re.hypot(self.im)
    }
    fn conj(&self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }
    fn powi(&self, n: i32) -> Self {
        Self::from_polar(self.abs().powi(n), self.phase() * (n as f64))
    }
    fn powf(&self, n: f64) -> Self {
        Self::from_polar(self.abs().powf(n), self.phase() * (n as f64))
    }
}

impl ComplexSpace for Complex<f64> {
    type Real = f64;

    fn i() -> Self {
        Self { re: 0.0, im: 1.0 }
    }
    fn cis(phase: Self::Real) -> Self {
        Self {
            re: phase.cos(),
            im: phase.sin(),
        }
    }
    fn from_real(re: Self::Real) -> Self {
        Self { re: re, im: 0.0 }
    }
    fn from_polar(r: Self::Real, phase: Self::Real) -> Self {
        Self {
            re: r * phase.cos(),
            im: r * phase.sin(),
        }
    }
    fn re(&self) -> Self::Real {
        self.re
    }
    fn im(&self) -> Self::Real {
        self.im
    }
    fn phase(&self) -> Self::Real {
        self.im.atan2(self.re)
    }
}

impl<S: ScalarSpace> fmt::Display for Complex<S> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Complex(re: {:?}, im: {:?})", self.re, self.im)
    }
}

impl<S: ScalarSpace> fmt::Debug for Complex<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Complex(re: {:?}, im: {:?})", self.re, self.im)
    }
}

impl<S: ScalarSpace> Neg for Complex<S> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            re: -self.re,
            im: -self.im,
        }
    }
}

impl<S: ScalarSpace> Add for Complex<S> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl<S: ScalarSpace> Sub for Complex<S> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl<S: ScalarSpace> Mul<f64> for Complex<S> {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            re: self.re * rhs,
            im: self.im * rhs,
        }
    }
}

impl<S: ScalarSpace> Div<f64> for Complex<S> {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            re: self.re / rhs,
            im: self.im / rhs,
        }
    }
}

impl<S: ScalarSpace> Mul for Complex<S> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl<S: ScalarSpace> Div for Complex<S> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self {
            re: (self.re * rhs.re + self.im * rhs.im) / (rhs.re * rhs.re + rhs.im * rhs.im),
            im: (self.im * rhs.re - self.re * rhs.im) / (rhs.re * rhs.re + rhs.im * rhs.im),
        }
    }
}

impl<S: ScalarSpace> Mul<Complex<S>> for f64 {
    type Output = Complex<S>;

    fn mul(self, rhs: Complex<S>) -> Self::Output {
        Self::Output {
            re: rhs.re * self,
            im: rhs.im * self,
        }
    }
}

impl<S: ScalarSpace> Div<Complex<S>> for f64 {
    type Output = Complex<S>;

    fn div(self, rhs: Complex<S>) -> Self::Output {
        Self::Output {
            re: (rhs.re * self) / (rhs.re * rhs.re + rhs.im * rhs.im),
            im: (-rhs.im * self) / (rhs.re * rhs.re + rhs.im * rhs.im),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    const EPS: f64 = 1e-10;

    #[test]
    fn test_complex_op() {
        let z1 = Complex::new(1.0, 2.0);
        let z2 = Complex::new(3.0, 4.0);
        let r = 3.0;

        assert_eq!(z1 + z2, Complex::new(4.0, 6.0), "Complex + Complex");
        assert_eq!(z1 - z2, Complex::new(-2.0, -2.0), "Complex - Complex");

        assert_eq!(-z1, Complex::new(-1.0, -2.0), "-Complex");

        assert_eq!(z1 * r, Complex::new(3.0, 6.0), "Complex * f64");
        assert_eq!(r * z1, Complex::new(3.0, 6.0), "f64 * Complex");

        assert_eq!(z1 / r, Complex::new(1.0 / 3.0, 2.0 / 3.0), "Complex / f64");
        assert_eq!(r / z1, Complex::new(0.6, -1.2), "f64 / Complex");

        assert_eq!(z1 * z2, Complex::new(-5.0, 10.0), "Complex * Complex");

        let i = Complex::i();
        assert_eq!(i * i, Complex::new(-1.0, 0.0), "i * i");

        assert_eq!(
            z1 / z2,
            Complex::new(11.0 / 25.0, 2.0 / 25.0),
            "Complex / Complex"
        );
        assert_eq!(z1 / z1, Complex::new(1.0, 0.0), "Complex / Complex");
    }

    #[test]
    fn test_complex_methods() {
        let z = Complex::new(3.0, 4.0);

        assert_eq!(z.abs(), 5.0, "Complex Absolute");
        assert_eq!(z.conj(), Complex::new(3.0, -4.0), "Complex Conjugate");
        assert_eq!(
            Complex::<f64>::one(),
            Complex::new(1.0, 0.0),
            "Complex::one()"
        );
        assert_eq!(z.size(), 1, "Complex.size() Dimension (1)");
        assert_eq!(z.get_data(), [3.0, 4.0]);
    }

    #[test]
    fn test_polar_conversion() {
        let z1 = Complex::from_polar(1.0, 0.0);
        assert!(
            (z1 - Complex::new(1.0, 0.0)).abs() < EPS,
            "(r=1, theta=0) -> 1 + 0i"
        );

        let z2 = Complex::from_polar(1.0, FRAC_PI_2);
        assert!(
            (z2 - Complex::new(0.0, 1.0)).abs() < EPS,
            "(r=1, theta=pi/2) -> 0 + 1i"
        );

        let sqrt_2 = 2.0_f64.sqrt();
        let z3 = Complex::from_polar(sqrt_2, FRAC_PI_4);
        assert!(
            (z3 - Complex::new(1.0, 1.0)).abs() < EPS,
            "(r=sqrt(2), theta=pi/4) -> 1 + 1i"
        );

        let z4 = Complex::cis(PI);
        assert!(
            (z4 - Complex::new(-1.0, 0.0)).abs() < EPS,
            "cis(PI) = e^(i*pi) = -1"
        );
    }

    #[test]
    fn test_phase_calculation() {
        let z1 = Complex::new(1.0, 1.0);
        assert!(
            (z1.phase() - FRAC_PI_4).abs() < EPS,
            "1 + i = sqrt(2) * e^(pi/4)"
        );
        assert!(
            (z1.arg() - FRAC_PI_4).abs() < EPS,
            "1 + i = sqrt(2) * e^(pi/4)"
        ); // alias 확인

        let z2 = Complex::new(-1.0, 0.0);
        assert!((z2.phase() - PI).abs() < EPS, "-1 = e^pi");

        let z3 = Complex::new(0.0, -1.0);
        assert!((z3.phase() - -FRAC_PI_2).abs() < EPS, "-i = e^(-pi)");
    }

    #[test]
    fn test_polar_round_trip() {
        let original = Complex::new(3.0, -4.0);

        let r = original.abs(); // 5.0
        let phi = original.phase(); // atan2(-4, 3)

        let restored = Complex::from_polar(r, phi);

        assert!(
            (original - restored).abs() < EPS,
            "Complex Ortho -> Polar -> Ortho Error is too Large!"
        );
    }
}
