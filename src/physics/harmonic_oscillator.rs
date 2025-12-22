use crate::math::core::{Scalar, ScalarSpace};
use crate::math::integrate::{EulerMethod, RK4Method, System};

pub struct SimpleHarmonicOscillator {
    pub omega_square: Scalar,
}

pub struct DampedHarmonicOscillator {
    pub k: Scalar,
    pub b: Scalar,
}

pub struct DrivenHarmonicOscillator {
    pub k: Scalar,
    pub b: Scalar,
    pub f0: Scalar,
    pub omega: Scalar,
}

impl System for SimpleHarmonicOscillator {
    type Vector = Scalar;
    fn derivative(&self, _t: Scalar, y: Self::Vector, _y_prime: Self::Vector) -> Self::Vector {
        -self.omega_square * y
    }
}

impl System for DampedHarmonicOscillator {
    type Vector = Scalar;
    fn derivative(&self, _t: Scalar, y: Self::Vector, y_prime: Self::Vector) -> Self::Vector {
        -self.k * y - self.b * y_prime
    }
}

impl System for DrivenHarmonicOscillator {
    type Vector = Scalar;
    fn derivative(&self, t: Scalar, y: Self::Vector, y_prime: Self::Vector) -> Self::Vector {
        self.f0 * (self.omega * t).cos() - self.k * y - self.b * y_prime
    }
}

#[cfg(test)]
mod tests {
    use crate::math::integrate::Solver;

    use super::*;

    #[test]
    pub fn sho_test() {
        type Method = RK4Method;

        let sho_ode = SimpleHarmonicOscillator {
            omega_square: Scalar(1.0),
        };

        let y0 = Scalar(0.0);
        let y0_prime = Scalar(1.0);
        let h = Scalar(0.1);
        let steps = 32;

        let mut sho_solver = Solver::new(RK4Method, sho_ode, y0, y0_prime);
        sho_solver.run(h, steps);
        let (ts, ys, ys_prime) = sho_solver.get_results_f64();

        for ((t, y), y_prime) in ts.iter().zip(ys).zip(ys_prime) {
            println!("({:.2}, {:.6}, {:.6})", t, y, y_prime);
        }

        assert!(false);
    }

    #[test]
    pub fn dho_test() {
        type Method = RK4Method;

        let dho_ode = DampedHarmonicOscillator {
            k: Scalar(1.0),
            b: Scalar(10.0),
        };

        let y0 = Scalar(0.0);
        let y0_prime = Scalar(1.0);
        let h = Scalar(0.1);
        let steps = 32;

        let mut dho_solver = Solver::new(RK4Method, dho_ode, y0, y0_prime);
        dho_solver.run(h, steps);
        let (ts, ys, ys_prime) = dho_solver.get_results_f64();

        for ((t, y), y_prime) in ts.iter().zip(ys).zip(ys_prime) {
            println!("({:.2}, {:.6}, {:.6})", t, y, y_prime);
        }

        assert!(false);
    }
}
