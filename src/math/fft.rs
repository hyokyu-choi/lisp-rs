use std::f64::consts::PI;

use crate::math::{
    complex::{Complex, ComplexSpace},
    core::{LinearSpace, ScalarSpace, Vector},
};

pub fn dft1d<const N: usize>(x_n: Vector<Complex, N>) -> Vector<Complex, N> {
    if N == 0 {
        return Vector::zero();
    }
    Vector::new(std::array::from_fn(|k| {
        (0..N).fold(Complex::zero(), |acc, n| {
            acc + Complex::cis(-2.0 * PI * (k as f64) * (n as f64) / (N as f64)) * x_n[n]
        })
    }))
}

pub fn idft1d<const N: usize>(x_k: Vector<Complex, N>) -> Vector<Complex, N> {
    if N == 0 {
        return Vector::zero();
    }
    Vector::new(std::array::from_fn(|n| {
        (0..N).fold(Complex::zero(), |acc, k| {
            acc + Complex::cis(2.0 * PI * (k as f64) * (n as f64) / (N as f64)) * x_k[k]
        }) / (N as f64)
    }))
}

/// FFT with Cooley-Tukey algorithm on a slice
pub fn fft1d_slice(x: &mut [Complex]) {
    let n = x.len();
    if n == 0 || (n & (n - 1)) != 0 {
        panic!("FFT length N must be a power of 2");
    }
    
    // Bit-reversal permutation
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;

        if i < j {
            x.swap(i, j);
        }
    }

    // Butterfly diagram
    let mut len = 2;
    while len <= n {
        let w_step = Complex::cis(-2.0 * PI / (len as f64));
        for i in (0..n).step_by(len) {
            let mut w = Complex::one();
            for j in 0..(len / 2) {
                let u = x[i + j];
                let v = x[i + j + len / 2] * w;
                x[i + j] = u + v;
                x[i + j + len / 2] = u - v;
                w = w * w_step;
            }
        }
        len <<= 1;
    }
}

/// IFFT with Cooley-Tukey algorithm on a slice
pub fn ifft1d_slice(x: &mut [Complex]) {
    let n = x.len();
    if n == 0 || (n & (n - 1)) != 0 {
        panic!("IFFT length N must be a power of 2");
    }

    // Bit-reversal permutation
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;

        if i < j {
            x.swap(i, j);
        }
    }

    // Butterfly diagram
    let mut len = 2;
    while len <= n {
        let w_step = Complex::cis(2.0 * PI / (len as f64));
        for i in (0..n).step_by(len) {
            let mut w = Complex::one();
            for j in 0..(len / 2) {
                let u = x[i + j];
                let v = x[i + j + len / 2] * w;
                x[i + j] = u + v;
                x[i + j + len / 2] = u - v;
                w = w * w_step;
            }
        }
        len <<= 1;
    }

    // Normalization
    let inv_n = 1.0 / (n as f64);
    for val in x.iter_mut() {
        *val = *val * inv_n;
    }
}

pub fn fft1d<const N: usize>(x_n: &Vector<Complex, N>) -> Vec<Complex> {
    let n_pow2 = N.next_power_of_two();
    let mut data = vec![Complex::zero(); n_pow2];
    data[..N].copy_from_slice(x_n.as_slice());
    fft1d_slice(&mut data);
    data
}

pub fn ifft1d<const N: usize>(x_k: &Vector<Complex, N>) -> Vec<Complex> {
    let n_pow2 = N.next_power_of_two();
    let mut data = vec![Complex::zero(); n_pow2];
    data[..N].copy_from_slice(x_k.as_slice());
    ifft1d_slice(&mut data);
    data
}

#[cfg(test)]
mod tests {
    use crate::math::core::ScalarSpace;

    use super::*;

    const EPS: f64 = 1e-14;

    fn assert_complex_slice_eq(
        a: &[Complex],
        b: &[Complex],
        msg: &str,
    ) {
        assert_eq!(a.len(), b.len(), "{} (length mismatch)", msg);
        let result = a
            .iter()
            .zip(b.iter())
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
        let frequancy = [
            Complex::zero(),
            Complex::new(0.0, -2.0),
            Complex::zero(),
            Complex::new(0.0, 2.0),
        ];
        let output = dft1d(x);
        assert_complex_slice_eq(&frequancy, output.as_slice(), "1D DFT")
    }

    #[test]
    fn test_idft1d() {
        let frequancy = Vector::new([
            Complex::zero(),
            Complex::new(0.0, -2.0),
            Complex::zero(),
            Complex::new(0.0, 2.0),
        ]);
        let x = [
            Complex::zero(),
            Complex::one(),
            Complex::zero(),
            -Complex::one(),
        ];
        let output = idft1d(frequancy);
        assert_complex_slice_eq(&x, output.as_slice(), "1D IDFT");
    }

    #[test]
    fn test_fft1d() {
        let x = Vector::new([
            Complex::zero(),
            Complex::one(),
            Complex::zero(),
            -Complex::one(),
        ]);
        let frequancy = [
            Complex::zero(),
            Complex::new(0.0, -2.0),
            Complex::zero(),
            Complex::new(0.0, 2.0),
        ];
        let output = fft1d(&x);
        assert_complex_slice_eq(&frequancy, &output, "1D FFT")
    }

    #[test]
    fn test_ifft1d() {
        let frequancy = Vector::new([
            Complex::zero(),
            Complex::new(0.0, -2.0),
            Complex::zero(),
            Complex::new(0.0, 2.0),
        ]);
        let x = [
            Complex::zero(),
            Complex::one(),
            Complex::zero(),
            -Complex::one(),
        ];
        let output = ifft1d(&frequancy);
        assert_complex_slice_eq(&x, &output, "1D IFFT");
    }

    #[test]
    fn test_fft1d_non_pow2() {
        let x = Vector::new([
            Complex::one(),
            Complex::one(),
            Complex::one(),
        ]);
        let output = fft1d(&x);
        assert_eq!(output.len(), 4); // Next power of 2
        
        // Manual calculation for padded FFT of [1, 1, 1, 0]
        // X[0] = 3, X[1] = -i, X[2] = 1, X[3] = i
        let expected = [
            Complex::new(3.0, 0.0),
            Complex::new(0.0, -1.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 1.0),
        ];
        assert_complex_slice_eq(&expected, &output, "1D FFT non-power-of-2 (padded)");
    }

    #[test]
    fn test_dft1d_non_pow2() {
        let x = Vector::new([
            Complex::one(),
            Complex::one(),
            Complex::one(),
        ]);
        let output = dft1d(x);
        // DFT of [1, 1, 1] is [3, 0, 0]
        let expected = [
            Complex::new(3.0, 0.0),
            Complex::zero(),
            Complex::zero(),
        ];
        assert_complex_slice_eq(&expected, output.as_slice(), "1D DFT non-power-of-2");
    }

    #[test]
    fn test_ifft1d_non_pow2() {
        // Input from the non-power-of-2 FFT test: FFT of [1, 1, 1] padded to 4
        let frequancy = Vector::new([
            Complex::new(3.0, 0.0),
            Complex::new(0.0, -1.0),
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 1.0), // This was the missing part
        ]);
        let output = ifft1d(&frequancy);
        assert_eq!(output.len(), 4);

        // Expected result is the original signal [1, 1, 1] zero-padded to 4
        let expected = [
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0),
            Complex::zero(),
        ];
        
        assert_complex_slice_eq(&expected, &output, "1D IFFT non-power-of-2 (padded)");
    }
}
