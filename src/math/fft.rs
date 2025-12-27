use std::f64::consts::PI;

use crate::math::{
    complex::{Complex, ComplexSpace},
    core::{LinearSpace, ScalarSpace, Vector},
};

pub fn dft1d<const N: usize>(x_n: Vector<Complex, N>) -> Vector<Complex, N> {
    if N == 0 || (N & (N - 1)) != 0 {
        panic!("DFT length N must be a power of 2"); // TODO: Implement zero psort
    }
    Vector::new(std::array::from_fn(|k| {
        (0..N).fold(Complex::zero(), |acc, n| {
            acc + Complex::cis(-2.0 * PI * (k as f64) * (n as f64) / (N as f64)) * x_n[n]
        })
    }))
}

pub fn idft1d<const N: usize>(x_k: Vector<Complex, N>) -> Vector<Complex, N> {
    if N == 0 || (N & (N - 1)) != 0 {
        panic!("IDFT length N must be a power of 2"); // TODO: Implement zero psort
    }
    Vector::new(std::array::from_fn(|n| {
        (0..N).fold(Complex::zero(), |acc, k| {
            acc + Complex::cis(2.0 * PI * (k as f64) * (n as f64) / (N as f64)) * x_k[k]
        }) / (N as f64)
    }))
}

/// FFT with Cooley-Tukey algorithm
///
/// using bit reverse for Radix-2 DIT divides
pub fn fft1d<const N: usize>(mut x_n: Vector<Complex, N>) -> Vector<Complex, N> {
    if N == 0 || (N & (N - 1)) != 0 {
        panic!("FFT length N must be a power of 2"); // TODO: Implement zero psort
    }
    // Radix-2 DIT divides using bit reverse sort
    // 반대방향 bit 덧셈
    let mut j = 0; // j: 000 부터 시작
    for i in 1..N {
        let mut bit = N >> 1; // bit: 100 부터 시작
        while j & bit != 0 {
            // j의 bit 자리가 1 이라면
            j ^= bit; // j의 bit 자리를 1에서 0으로 변경
            bit >>= 1; // bit 자리가 오른쪽으로 한칸 이동
        }
        j ^= bit; // j의 bit 자리를 0에서 1로 변경

        if i < j {
            // Index swap x_n.swap(i, j)
            let temp = x_n[i];
            x_n[i] = x_n[j];
            x_n[j] = temp;
        }
    }
    // Butterfly diagram
    let mut len = 2;
    while len <= N {
        let w_step = Complex::cis(-2.0 * PI / (len as f64)); // twiddle factor
        for i in (0..N).step_by(len) {
            let mut w = Complex::one(); // twiddle factor
            for j in 0..(len / 2) {
                let u = x_n[i + j]; // even
                let v = x_n[i + j + len / 2] * w; // odd

                // 1D DFT
                x_n[i + j] = u + v;
                x_n[i + j + len / 2] = u - v;

                w = w * w_step;
            }
        }
        len <<= 1;
    }

    x_n
}

/// IFFT with Cooley-Tukey algorithm
///
/// using bit reverse for Radix-2 DIT divides
pub fn ifft1d<const N: usize>(mut x_k: Vector<Complex, N>) -> Vector<Complex, N> {
    if N == 0 || (N & (N - 1)) != 0 {
        panic!("IFFT length N must be a power of 2"); // TODO: Implement zero psort
    }
    // Divide and Conquer using bit reverse sort
    // 반대방향 bit 덧셈
    let mut j = 0; // j: 000 부터 시작
    for i in 1..N {
        let mut bit = N >> 1; // bit: 100 부터 시작
        while j & bit != 0 {
            // j의 bit 자리가 1 이라면
            j ^= bit; // j의 bit 자리를 1에서 0으로 변경
            bit >>= 1; // bit 자리가 오른쪽으로 한칸 이동
        }
        j ^= bit; // j의 bit 자리를 0에서 1로 변경

        if i < j {
            // Index swap x_k.swap(i, j)
            let temp = x_k[i];
            x_k[i] = x_k[j];
            x_k[j] = temp;
        }
    }
    // Butterfly diagram
    let mut len = 2;
    while len <= N {
        let w_step = Complex::cis(2.0 * PI / (len as f64)); // twiddle factor
        for i in (0..N).step_by(len) {
            let mut w = Complex::one(); // twiddle factor
            for j in 0..(len / 2) {
                let u = x_k[i + j]; // even
                let v = x_k[i + j + len / 2] * w; // odd

                // 1D DFT
                x_k[i + j] = u + v;
                x_k[i + j + len / 2] = u - v;

                w = w * w_step;
            }
        }
        len <<= 1;
    }

    x_k / (N as f64)
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
