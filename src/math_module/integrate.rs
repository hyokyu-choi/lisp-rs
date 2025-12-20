// use std::process::Output;

use super::core::{Scalar, LinearSpace};

pub trait DiffEqSystem {
    type Vector: LinearSpace + Copy;

    fn derivative(&self, t: Scalar, u: Self::Vector) -> Self::Vector;
}

pub trait Integrator {
    fn step<S>(&self, system: &S, t: Scalar, u: S::Vector, h: Scalar) -> S::Vector
    where S: DiffEqSystem;
}

pub struct EulerMethod;
pub struct RK4Method;

impl Integrator for EulerMethod {
    fn step<S>(&self, system: &S, t: Scalar, u: S::Vector, h: Scalar) -> S::Vector
        where S: DiffEqSystem {
        u + system.derivative(t, u) * h
    }
}

impl Integrator for RK4Method {
    fn step<S>(&self, system: &S, t: Scalar, u: S::Vector, h: Scalar) -> S::Vector
        where
        S: DiffEqSystem,
    {
        let k1 = system.derivative(t, u);
        let k2 = system.derivative(t + h / 2.0, u + k1 * (h / 2.0));
        let k3 = system.derivative(t + h / 2.0, u + k2 * (h / 2.0));
        let k4 = system.derivative(t + h, u + k3 * h);
        u + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0)
    }
}
