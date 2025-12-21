use super::core::{LinearSpace, Scalar};

pub trait DiffEqSystem {
    type Vector: LinearSpace + Copy;

    fn derivative(&self, t: Scalar, y: Self::Vector, u_prime: Self::Vector) -> Self::Vector;
}

pub trait Integrator {
    fn step<S>(
        system: &S,
        t: Scalar,
        y: S::Vector,
        y_prime: S::Vector,
        h: Scalar,
    ) -> (S::Vector, S::Vector)
    where
        S: DiffEqSystem;
}

pub struct EulerMethod;
pub struct RK4Method;

impl Integrator for EulerMethod {
    fn step<S>(
        system: &S,
        t: Scalar,
        y: S::Vector,
        y_prime: S::Vector,
        h: Scalar,
    ) -> (S::Vector, S::Vector)
    where
        S: DiffEqSystem,
    {
        (
            y + y_prime * h,
            y_prime + system.derivative(t, y, y_prime) * h,
        )
    }
}

impl Integrator for RK4Method {
    fn step<S>(
        system: &S,
        t: Scalar,
        y: S::Vector,
        y_prime: S::Vector,
        h: Scalar,
    ) -> (S::Vector, S::Vector)
    where
        S: DiffEqSystem,
    {
        let k11 = y_prime;
        let k12 = system.derivative(t, y, y_prime);
        let k21 = y_prime + k12 * (h / 2.0);
        let k22 = system.derivative(t + h / 2.0, y + k11 * (h / 2.0), y_prime + k12 * (h / 2.0));
        let k31 = y_prime + k22 * (h / 2.0);
        let k32 = system.derivative(t + h / 2.0, y + k21 * (h / 2.0), y_prime + k22 * (h / 2.0));
        let k41 = y_prime + k32 * h;
        let k42 = system.derivative(t + h, y + k31 * h, y_prime + k32 * h);
        (
            y + (k11 + k21 * 2.0 + k31 * 2.0 + k41) * (h / 6.0),
            y_prime + (k12 + k22 * 2.0 + k32 * 2.0 + k42) * (h / 6.0),
        )
    }
}
