use super::core::{LinearSpace, Scalar};

pub trait System {
    type Vector: LinearSpace + Copy;

    fn derivative(&self, t: Scalar, y: Self::Vector, u_prime: Self::Vector) -> Self::Vector;
}

pub trait Integrator {
    fn step<S>(
        &mut self,
        system: &S,
        t: Scalar,
        y: S::Vector,
        y_prime: S::Vector,
        h: Scalar,
    ) -> (S::Vector, S::Vector)
    where
        S: System;
}

pub struct Solver<I, S>
where
    I: Integrator,
    S: System,
{
    integrator: I,
    system: S,
    t: Scalar,
    y: S::Vector,
    y_prime: S::Vector,
    results: Results<S>,
}

#[derive(Clone)]
pub struct Results<S: System> {
    ts: Vec<Scalar>,
    ys: Vec<S::Vector>,
    ys_prime: Vec<S::Vector>,
}

impl<S: System> Results<S> {
    pub fn new(t0: Scalar, y0: S::Vector, y0_prime: S::Vector) -> Self {
        Self {
            ts: vec![t0],
            ys: vec![y0],
            ys_prime: vec![y0_prime],
        }
    }

    pub fn update(&mut self, t: Scalar, y: S::Vector, y_prime: S::Vector) {
        self.ts.push(t);
        self.ys.push(y);
        self.ys_prime.push(y_prime);
    }

    pub fn get(&self, index: usize) -> (Scalar, S::Vector, S::Vector) {
        (self.ts[index], self.ys[index], self.ys_prime[index])
    }

    pub fn get_f64(&self, index: usize) -> (f64, <S::Vector as LinearSpace>::Value, <S::Vector as LinearSpace>::Value) {
        (self.ts[index].get_value(), self.ys[index].get_value(), self.ys_prime[index].get_value())
    }

    pub fn get_ts(&self) -> Vec<Scalar> {
        self.ts.clone()
    }

    pub fn get_ts_f64(&self) -> Vec<<Scalar as LinearSpace>::Value> {
        self.ts.iter().map(|t| t.get_value()).collect()
    }

    pub fn get_ys(&self) -> Vec<S::Vector> {
        self.ys.clone()
    }

    pub fn get_ys_f64(&self) -> Vec<<S::Vector as LinearSpace>::Value> {
        self.ys.iter().map(|y| y.get_value()).collect()
    }

    pub fn get_ys_prime(&self) -> Vec<S::Vector> {
        self.ys_prime.clone()
    }

    pub fn get_ys_prime_f64(&self) -> Vec<<S::Vector as LinearSpace>::Value> {
        self.ys_prime.iter().map(|y| y.get_value()).collect()
    }
}

impl<I, S> Solver<I, S>
where
    I: Integrator,
    S: System,
{
    pub fn new(integrator: I, system: S, y0: S::Vector, y0_prime: S::Vector) -> Self {
        Self {
            integrator: integrator,
            system: system,
            t: Scalar::new(0.0),
            y: y0,
            y_prime: y0_prime,
            results: Results::new(Scalar::new(0.0), y0, y0_prime),
        }
    }

    pub fn run(&mut self, h: Scalar, steps: usize) {
        for _ in 0..steps {
            let (y, y_prime) = self
                .integrator
                .step(&self.system, self.t, self.y, self.y_prime, h);
            self.t = self.t + h;
            self.y = y;
            self.y_prime = y_prime;
            self.update();
        }
    }

    fn update(&mut self) {
        self.results.update(self.t, self.y, self.y_prime);
    }

    // TODO: implement
    pub fn get_results(&self) -> &Results<S> {
        &self.results
    }

    pub fn get_ts_f64(&self) -> Vec<f64> {
        self.results.get_ts_f64()
    }

    pub fn get_ys_f64(&self) -> Vec<<S::Vector as LinearSpace>::Value> {
        self.results.get_ys_f64()
    }

    pub fn get_ys_prime_f64(&self) -> Vec<<S::Vector as LinearSpace>::Value> {
        self.results.get_ys_prime_f64()
    }
}

pub struct EulerMethod;
pub struct RK4Method;

impl Integrator for EulerMethod {
    fn step<S>(
        &mut self,
        system: &S,
        t: Scalar,
        y: S::Vector,
        y_prime: S::Vector,
        h: Scalar,
    ) -> (S::Vector, S::Vector)
    where
        S: System,
    {
        (
            y + h * y_prime,
            y_prime + h * system.derivative(t, y, y_prime),
        )
    }
}

impl Integrator for RK4Method {
    fn step<S>(
        &mut self,
        system: &S,
        t: Scalar,
        y: S::Vector,
        y_prime: S::Vector,
        h: Scalar,
    ) -> (S::Vector, S::Vector)
    where
        S: System,
    {
        let k11 = y_prime;
        let k12 = system.derivative(t, y, y_prime);
        let k21 = y_prime + h * (k12 / 2.0);
        let k22 = system.derivative(t + h / 2.0, y + h * (k12 / 2.0), y_prime + h * (k12 / 2.0));
        let k31 = y_prime + k22 * (h / 2.0);
        let k32 = system.derivative(t + h / 2.0, y + h * (k21 / 2.0), y_prime + h * (k22 / 2.0));
        let k41 = y_prime + k32 * h;
        let k42 = system.derivative(t + h, y + h * k31, y_prime + h * k32);
        (
            y + (k11 + k21 * 2.0 + k31 * 2.0 + k41) * (h / 6.0),
            y_prime + (k12 + k22 * 2.0 + k32 * 2.0 + k42) * (h / 6.0),
        )
    }
}
