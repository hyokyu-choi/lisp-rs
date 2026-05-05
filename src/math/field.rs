use crate::math::{
    complex::Complex,
    fft::{fft1d_slice, ifft1d_slice},
};
use std::ops::{Index, IndexMut};

pub struct Field3D {
    pub n: usize,
    pub data: Vec<Complex>,
}

impl Field3D {
    pub fn new(n: usize) -> Self {
        let n_pow2 = n.next_power_of_two();
        Self {
            n: n_pow2,
            data: vec![Complex::new(0.0, 0.0); n_pow2 * n_pow2 * n_pow2],
        }
    }

    pub fn zero(n: usize) -> Self {
        Self::new(n)
    }

    #[inline]
    pub fn get_idx(&self, x: usize, y: usize, z: usize) -> usize {
        (x * self.n + y) * self.n + z
    }
}

impl Index<[usize; 3]> for Field3D {
    type Output = Complex;
    fn index(&self, idx: [usize; 3]) -> &Self::Output {
        &self.data[self.get_idx(idx[0], idx[1], idx[2])]
    }
}

impl IndexMut<[usize; 3]> for Field3D {
    fn index_mut(&mut self, idx: [usize; 3]) -> &mut Self::Output {
        let i = self.get_idx(idx[0], idx[1], idx[2]);
        &mut self.data[i]
    }
}

pub trait Field {
    fn fft(&mut self) -> &mut Self;
    fn ifft(&mut self) -> &mut Self;
}

impl Field for Field3D {
    fn fft(&mut self) -> &mut Self {
        let n = self.n;

        // 1. FFT along Z axis (contiguous)
        for x in 0..n {
            for y in 0..n {
                let start = (x * n + y) * n;
                fft1d_slice(&mut self.data[start..start + n]);
            }
        }

        // 2. FFT along Y axis (stride n)
        let mut temp = vec![Complex::new(0.0, 0.0); n];
        for x in 0..n {
            for z in 0..n {
                for y in 0..n {
                    temp[y] = self[[x, y, z]];
                }
                fft1d_slice(&mut temp);
                for y in 0..n {
                    self[[x, y, z]] = temp[y];
                }
            }
        }

        // 3. FFT along X axis (stride n^2)
        for y in 0..n {
            for z in 0..n {
                for x in 0..n {
                    temp[x] = self[[x, y, z]];
                }
                fft1d_slice(&mut temp);
                for x in 0..n {
                    self[[x, y, z]] = temp[x];
                }
            }
        }
        self
    }

    fn ifft(&mut self) -> &mut Self {
        let n = self.n;

        // 1. IFFT along Z axis (contiguous)
        for x in 0..n {
            for y in 0..n {
                let start = (x * n + y) * n;
                ifft1d_slice(&mut self.data[start..start + n]);
            }
        }

        // 2. IFFT along Y axis (stride n)
        let mut temp = vec![Complex::new(0.0, 0.0); n];
        for x in 0..n {
            for z in 0..n {
                for y in 0..n {
                    temp[y] = self[[x, y, z]];
                }
                ifft1d_slice(&mut temp);
                for y in 0..n {
                    self[[x, y, z]] = temp[y];
                }
            }
        }

        // 3. IFFT along X axis (stride n^2)
        for y in 0..n {
            for z in 0..n {
                for x in 0..n {
                    temp[x] = self[[x, y, z]];
                }
                ifft1d_slice(&mut temp);
                for x in 0..n {
                    self[[x, y, z]] = temp[x];
                }
            }
        }
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::complex::{Complex, ComplexSpace};

    #[test]
    fn test_field3d_reversibility() {
        let n = 4;
        let mut field = Field3D::new(n);

        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    field[[x, y, z]] = Complex::new((x + y + z) as f64, 1.0);
                }
            }
        }

        let original_val = field[[1, 2, 1]];
        field.fft().ifft();

        let diff = (field[[1, 2, 1]].re() - original_val.re()).abs();
        assert!(diff < 1e-10, "Reversibility failed: diff {}", diff);
    }

    #[test]
    fn test_field3d_auto_padding() {
        let field = Field3D::new(10); // 10 is not power of 2
        assert_eq!(field.n, 16); // Should pad to 16
        assert_eq!(field.data.len(), 16 * 16 * 16);
    }
}
