use crate::math_module::core::{Scalar};
use crate::math_module::integrate::{DiffEqSystem, EulerMethod, Integrator, RK4Method};

struct SHO {
    m: Scalar,
    k: Scalar,
}

impl DiffEqSystem for SHO {
    type Vector = Scalar;
    fn derivative(&self, _t: Scalar, y: Self::Vector, _y_prime: Self::Vector) -> Self::Vector {
        -(self.k / self.m) * y
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn sho_test() {
        type Method = RK4Method;

        let sho_ode = SHO{m: Scalar(1.0), k: Scalar(1.0)};

        let t0 = Scalar(0.0);
        let y0 = Scalar(0.0);
        let y0_prime = Scalar(1.0);
        let h = Scalar(0.1);
        let step = 32;

        let mut t = t0;
        let mut y = y0;
        let mut y_prime = y0_prime;

        let mut ts = vec![t0];
        let mut ys = vec![y0];
        let mut ys_prime = vec![y0_prime];

        let mut result = vec![(t0, y0, y0_prime)];

        for _ in 0..step {
            t = t + h;
            (y, y_prime) = Method::step(&sho_ode, t, y, y_prime, h);
            ts.push(t);
            ys.push(y);
            ys_prime.push(y_prime);
            result.push((t, y, y_prime));
        }
        // dbg!((ts, ys, ys_prime));

        for (t, y, y_prime) in result {
            println!("({}, {}, {})", t, y, y_prime);
        }

        assert!(false);
    }
}
