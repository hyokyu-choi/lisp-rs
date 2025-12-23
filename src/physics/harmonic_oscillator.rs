use crate::math::core::LinearSpace;
use crate::math::integrate::{EulerMethod, RK4Method, System};

pub struct SimpleHarmonicOscillator {
    pub omega_square: f64,
}

pub struct DampedHarmonicOscillator {
    pub k: f64,
    pub b: f64,
}

pub struct DrivenHarmonicOscillator {
    pub k: f64,
    pub b: f64,
    pub f0: f64,
    pub omega: f64,
}

impl System for SimpleHarmonicOscillator {
    type Vector = f64;

    fn derivative(&self, _t: f64, y: Self::Vector, _y_prime: Self::Vector) -> Self::Vector {
        -self.omega_square * y
    }
}

impl System for DampedHarmonicOscillator {
    type Vector = f64;

    fn derivative(&self, _t: f64, y: Self::Vector, y_prime: Self::Vector) -> Self::Vector {
        -self.k * y - self.b * y_prime
    }
}

impl System for DrivenHarmonicOscillator {
    type Vector = f64;

    fn derivative(&self, t: f64, y: Self::Vector, y_prime: Self::Vector) -> Self::Vector {
        Self::Vector::new(self.f0 * (self.omega * t).cos()) - self.k * y - self.b * y_prime
    }
}

#[cfg(test)]
mod tests {
    use crate::math::integrate::Solver;

    use super::*;

    #[test]
    pub fn sho_test() {
        let method = RK4Method;

        let sho_ode = SimpleHarmonicOscillator { omega_square: 1.0 };

        let y0 = 0.0;
        let y0_prime = 1.0;
        let h = 0.1;
        let steps = 32;

        let mut sho_solver = Solver::new(method, sho_ode, y0, y0_prime);
        sho_solver.run(h, steps);
        let (ts, ys, ys_prime) = sho_solver.get_results_f64();

        for ((t, y), y_prime) in ts.iter().zip(ys).zip(ys_prime) {
            println!("({:.2}, {:.6}, {:.6})", t, y, y_prime);
        }

        assert!(true);
    }

    #[test]
    pub fn dho_test() {
        let method = RK4Method;

        let dho_ode = DampedHarmonicOscillator { k: 1.0, b: 10.0 };

        let y0 = 0.0;
        let y0_prime = 1.0;
        let h = 0.1;
        let steps = 32;

        let mut dho_solver = Solver::new(method, dho_ode, y0, y0_prime);
        dho_solver.run(h, steps);
        let (ts, ys, ys_prime) = dho_solver.get_results_f64();

        for ((t, y), y_prime) in ts.iter().zip(ys).zip(ys_prime) {
            println!("({:.2}, {:.6}, {:.6})", t, y, y_prime);
        }

        assert!(true);
    }
}
