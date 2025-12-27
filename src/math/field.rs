use crate::math::{
    complex::Complex,
    core::Vector,
    fft::{fft1d, ifft1d},
};

pub type Field1D<const N: usize> = Vector<Complex, N>;
pub type Field2D<const N: usize> = Vector<Vector<Complex, N>, N>;
pub type Field3D<const N: usize> = Vector<Vector<Vector<Complex, N>, N>, N>;

pub trait Field {
    fn fft(&mut self) -> &mut Self;
    fn ifft(&mut self) -> &mut Self;
}

impl<const N: usize> Field for Field1D<N> {
    fn fft(&mut self) -> &mut Self {
        fft1d(self);
        self
    }
    fn ifft(&mut self) -> &mut Self {
        ifft1d(self);
        self
    }
}

impl<const N: usize> Field for Field2D<N> {
    fn fft(&mut self) -> &mut Self {
        for x in 0..N {
            fft1d(&mut self[x]);
        }

        // Transpose
        for x in 0..N {
            for y in (x + 1)..N {
                let temp = self[x][y];
                self[x][y] = self[y][x];
                self[y][x] = temp;
            }
        }

        for x in 0..N {
            fft1d(&mut self[x]);
        }

        // Transpose
        for x in 0..N {
            for y in (x + 1)..N {
                let temp = self[x][y];
                self[x][y] = self[y][x];
                self[y][x] = temp;
            }
        }
        self
    }
    fn ifft(&mut self) -> &mut Self {
        for x in 0..N {
            ifft1d(&mut self[x]);
        }

        // Transpose
        for x in 0..N {
            for y in (x + 1)..N {
                let temp = self[x][y];
                self[x][y] = self[y][x];
                self[y][x] = temp;
            }
        }

        for x in 0..N {
            ifft1d(&mut self[x]);
        }

        // Transpose
        for x in 0..N {
            for y in (x + 1)..N {
                let temp = self[x][y];
                self[x][y] = self[y][x];
                self[y][x] = temp;
            }
        }
        self
    }
}

impl<const N: usize> Field for Field3D<N> {
    fn fft(&mut self) -> &mut Self {
        // 2D FFT in YZ surface
        for x in 0..N {
            self[x].fft();
        }

        // Transpose xz
        for x in 0..N {
            for y in 0..N {
                for z in (x + 1)..N {
                    let temp = self[x][y][z];
                    self[x][y][z] = self[z][y][x];
                    self[z][y][x] = temp;
                }
            }
        }

        // 1D FFT in X axis
        for z in 0..N {
            for y in 0..N {
                self[z][y].fft();
            }
        }

        // Transpose XZ
        for x in 0..N {
            for y in 0..N {
                for z in (x + 1)..N {
                    let temp = self[x][y][z];
                    self[x][y][z] = self[z][y][x];
                    self[z][y][x] = temp;
                }
            }
        }
        self
    }
    fn ifft(&mut self) -> &mut Self {
        // 2D IFFT in YZ surface
        for x in 0..N {
            self[x].ifft();
        }

        // Transpose xz
        for x in 0..N {
            for y in 0..N {
                for z in (x + 1)..N {
                    let temp = self[x][y][z];
                    self[x][y][z] = self[z][y][x];
                    self[z][y][x] = temp;
                }
            }
        }

        // 1D IFFT in X axis
        for z in 0..N {
            for y in 0..N {
                self[z][y].ifft();
            }
        }

        // Transpose XZ
        for x in 0..N {
            for y in 0..N {
                for z in (x + 1)..N {
                    let temp = self[x][y][z];
                    self[x][y][z] = self[z][y][x];
                    self[z][y][x] = temp;
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
    use crate::math::core::{LinearSpace, ScalarSpace};

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
    fn test_field1d_reversibility() {
        const N: usize = 8;
        let mut data = [Complex::zero(); N];
        for i in 0..N {
            data[i] = Complex::new(i as f64, (i as f64) * 0.5);
        }

        let mut field = Field1D::<N>::new(data);
        let original = field.clone();

        field.fft().ifft();

        assert_complex_vector_eq(field, original, "Field1D.fft() reversibility failed.");
    }

    #[test]
    fn test_field2d_reversibility() {
        const N: usize = 4;
        let mut field = Field2D::<N>::zero();
        let mut original = Field2D::<N>::zero();

        for x in 0..N {
            for y in 0..N {
                let val = Complex::new((x + y) as f64, (x * y) as f64);
                field[x][y] = val;
                original[x][y] = val;
            }
        }

        field.fft().ifft();

        for x in 0..N {
            assert_complex_vector_eq(
                field[x],
                original[x],
                "Field2D.fft() reversibility failed. Check transpose logic",
            );
        }
    }

    #[test]
    fn test_field3d_reversibility() {
        const N: usize = 4;
        let mut field = Field3D::<N>::zero();
        let mut original = Field3D::<N>::zero();

        for x in 0..N {
            for y in 0..N {
                for z in 0..N {
                    let val = Complex::new((x * 10 + y * 5 + z) as f64, 1.0);
                    field[x][y][z] = val;
                    original[x][y][z] = val;
                }
            }
        }

        field.fft().ifft();

        for x in 0..N {
            for y in 0..N {
                assert_complex_vector_eq(
                    field[x][y],
                    original[x][y],
                    "Field3D.fft() reversibility failed. Check transpose logic",
                );
            }
        }
    }

    #[test]
    fn test_impulse_response_check() {
        // FFT(Delta function) = Const
        // x[0] = 1, otherwise 0 -> FFT -> X[k] = 1
        const N: usize = 4;

        let mut field = Field1D::<N>::zero();
        field[0] = Complex::one(); // Delta function

        field.fft();

        let freq = Vector::new(std::array::from_fn(|_| Complex::one()));

        assert_complex_vector_eq(field, freq, "Field1D.fft() failed. Delta fuction check.");
    }

    #[test]
    fn test_field3d_axis_ordering() {
        const N: usize = 4;
        let mut field = Field3D::<N>::zero();

        field[1][2][3] = Complex::one();

        field.fft();

        field.ifft();

        assert_eq!(
            field[1][2][3],
            Complex::one(),
            "Field3D transpose test are failed."
        );
        assert_eq!(
            field[0][0][0],
            Complex::zero(),
            "Field3D transpose test are failed."
        );
        assert_eq!(
            field[3][2][1],
            Complex::zero(),
            "Field3D transpose test are failed."
        );
    }
}
