use std::f64::consts::PI;

use crate::math::{
    complex::{Complex, ComplexSpace},
    core::{LinearSpace, Vector},
};

pub fn dft1d<const N: usize>(x_n: Vector<Complex, N>) -> Vector<Complex, N> {
    Vector::new(std::array::from_fn(|k| {
        (0..N).fold(Complex::zero(), |acc, n| {
            acc + Complex::cis(-2.0 * PI * (k as f64) * (n as f64) / (N as f64)) * x_n[n]
        })
    }))
}

pub fn idft1d<const N: usize>(x_k: Vector<Complex, N>) -> Vector<Complex, N> {
    Vector::new(std::array::from_fn(|n| {
        (0..N).fold(Complex::zero(), |acc, k| {
            acc + Complex::cis(2.0 * PI * (k as f64) * (n as f64) / (N as f64)) * x_k[k]
        }) / (N as f64)
    }))
}

pub fn fft1d<const N: usize>(x: Vector<Complex, N>) -> Vector<Complex, N> {
    Vector::zero() // TODO: Implement
}

pub fn ifft1d<const N: usize>(x: Vector<Complex, N>) -> Vector<Complex, N> {
    Vector::zero() // TODO: Implement
}

#[cfg(test)]
mod tests {
    use crate::math::core::ScalarSpace;

    use super::*;

    const EPS: f64 = 1e-14;

    fn assert_complex_vector_eq<const N: usize>(
        a: Vector<Complex, N>,
        b: Vector<Complex, N>,
        msg: &str,
    ) {
        let result = a
            .get_data()
            .iter()
            .zip(b.get_data().iter())
            .map(|(a, b)| (a.re() - b.re()).abs() < EPS && (a.im() - b.im()).abs() < EPS)
            .fold(true, |acc, val| acc && val);
        assert!(result, "{msg}\n left: {:?}\nright: {:?}", a, b)
    }

    #[test]
    fn test_dft1d() {
        let x = Vector::new([
            Complex::zero(),
            Complex::one(),
            Complex::zero(),
            -Complex::one(),
        ]);
        let frequancy = Vector::new([
            Complex::zero(),
            Complex::new(0.0, -2.0),
            Complex::zero(),
            Complex::new(0.0, 2.0),
        ]);
        let output = dft1d(x);
        assert_complex_vector_eq(frequancy, output, "1D DFT")
    }

    #[test]
    fn test_idft1d() {
        let frequancy = Vector::new([
            Complex::zero(),
            Complex::new(0.0, -2.0),
            Complex::zero(),
            Complex::new(0.0, 2.0),
        ]);
        let x = Vector::new([
            Complex::zero(),
            Complex::one(),
            Complex::zero(),
            -Complex::one(),
        ]);
        let output = idft1d(frequancy);
        assert_complex_vector_eq(x, output, "1D IDFT");
    }

    #[test]
    fn test_fft1d() {
        let x = Vector::new([
            Complex::zero(),
            Complex::one(),
            Complex::zero(),
            -Complex::one(),
        ]);
        let frequancy = Vector::new([
            Complex::zero(),
            Complex::new(0.0, -2.0),
            Complex::zero(),
            Complex::new(0.0, 2.0),
        ]);
        let output = fft1d(x);
        assert_complex_vector_eq(frequancy, output, "1D FFT")
    }

    #[test]
    fn test_ifft1d() {
        let frequancy = Vector::new([
            Complex::zero(),
            Complex::new(0.0, -2.0),
            Complex::zero(),
            Complex::new(0.0, 2.0),
        ]);
        let x = Vector::new([
            Complex::zero(),
            Complex::one(),
            Complex::zero(),
            -Complex::one(),
        ]);
        let output = ifft1d(frequancy);
        assert_complex_vector_eq(x, output, "1D IFFT");
    }
}
