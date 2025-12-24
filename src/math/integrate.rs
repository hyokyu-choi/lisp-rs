use crate::math::core::LinearSpace;

pub trait System {
    type Vector: LinearSpace + Copy;

    fn derivative(&self, t: f64, y: Self::Vector, y_prime: Self::Vector) -> Self::Vector;
}

pub trait Integrator {
    fn step<S>(
        &mut self,
        system: &S,
        t: f64,
        y: S::Vector,
        y_prime: S::Vector,
        h: f64,
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
    t: f64,
    y: S::Vector,
    y_prime: S::Vector,
    results: Results<S>,
}

#[derive(Clone)]
pub struct Results<S: System> {
    ts: Vec<f64>,
    ys: Vec<S::Vector>,
    ys_prime: Vec<S::Vector>,
}

impl<S: System> Results<S> {
    pub fn new(t0: f64, y0: S::Vector, y0_prime: S::Vector) -> Self {
        Self {
            ts: vec![t0],
            ys: vec![y0],
            ys_prime: vec![y0_prime],
        }
    }

    pub fn update(&mut self, t: f64, y: S::Vector, y_prime: S::Vector) {
        self.ts.push(t);
        self.ys.push(y);
        self.ys_prime.push(y_prime);
    }

    pub fn get(&self, index: usize) -> (f64, S::Vector, S::Vector) {
        (self.ts[index], self.ys[index], self.ys_prime[index])
    }

    pub fn get_f64(
        &self,
        index: usize,
    ) -> (
        f64,
        <S::Vector as LinearSpace>::Data,
        <S::Vector as LinearSpace>::Data,
    ) {
        (
            self.ts[index],
            self.ys[index].get_data(),
            self.ys_prime[index].get_data(),
        )
    }

    pub fn get_ts(&self) -> Vec<f64> {
        self.ts.clone()
    }

    pub fn get_ts_f64(&self) -> Vec<f64> {
        self.ts.clone()
    }

    pub fn get_ys(&self) -> Vec<S::Vector> {
        self.ys.clone()
    }

    pub fn get_ys_f64(&self) -> Vec<<S::Vector as LinearSpace>::Data> {
        self.ys.iter().map(|y| y.get_data()).collect()
    }

    pub fn get_ys_prime(&self) -> Vec<S::Vector> {
        self.ys_prime.clone()
    }

    pub fn get_ys_prime_f64(&self) -> Vec<<S::Vector as LinearSpace>::Data> {
        self.ys_prime.iter().map(|y| y.get_data()).collect()
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
            t: 0.0,
            y: y0,
            y_prime: y0_prime,
            results: Results::new(0.0, y0, y0_prime),
        }
    }

    pub fn run(&mut self, h: f64, steps: usize) {
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

    pub fn get_results(&self) -> (Vec<f64>, Vec<S::Vector>, Vec<S::Vector>) {
        (self.get_ts(), self.get_ys(), self.get_ys_prime())
    }

    pub fn get_results_f64(
        &self,
    ) -> (
        Vec<<f64 as LinearSpace>::Data>,
        Vec<<S::Vector as LinearSpace>::Data>,
        Vec<<S::Vector as LinearSpace>::Data>,
    ) {
        (
            self.get_ts_f64(),
            self.get_ys_f64(),
            self.get_ys_prime_f64(),
        )
    }

    fn get_ts(&self) -> Vec<f64> {
        self.results.get_ts()
    }

    fn get_ts_f64(&self) -> Vec<f64> {
        self.results.get_ts_f64()
    }

    fn get_ys(&self) -> Vec<S::Vector> {
        self.results.get_ys()
    }

    fn get_ys_f64(&self) -> Vec<<S::Vector as LinearSpace>::Data> {
        self.results.get_ys_f64()
    }

    fn get_ys_prime(&self) -> Vec<S::Vector> {
        self.results.get_ys_prime()
    }

    fn get_ys_prime_f64(&self) -> Vec<<S::Vector as LinearSpace>::Data> {
        self.results.get_ys_prime_f64()
    }

    fn get_current(&self) -> (f64, S::Vector, S::Vector) {
        (self.t, self.y, self.y_prime)
    }
}

pub struct EulerMethod;
pub struct RK4Method;

/// Leapfrog Integration (Velocity Verlet)
///
/// Leapfrog integration is a method for 2nd-order ODEss without first derivatives.
/// $$ y'' = f(t, y, _y') $$
/// Stable for oscillatory motion.
/// **System::derivative must not use y_prime in return.**
pub struct LeapfrogMethod;

impl Integrator for EulerMethod {
    fn step<S>(
        &mut self,
        system: &S,
        t: f64,
        y: S::Vector,
        y_prime: S::Vector,
        h: f64,
    ) -> (S::Vector, S::Vector)
    where
        S: System,
    {
        (
            y + y_prime * h,
            y_prime + system.derivative(t, y, y_prime) * h,
        )
    }
}

impl Integrator for RK4Method {
    fn step<S>(
        &mut self,
        system: &S,
        t: f64,
        y: S::Vector,
        y_prime: S::Vector,
        h: f64,
    ) -> (S::Vector, S::Vector)
    where
        S: System,
    {
        let k11 = y_prime;
        let k12 = system.derivative(t, y, y_prime);
        let k21 = y_prime + (k12 / 2.0) * h;
        let k22 = system.derivative(t + h / 2.0, y + (k12 / 2.0) * h, y_prime + (k12 / 2.0) * h);
        let k31 = y_prime + k22 * (h / 2.0);
        let k32 = system.derivative(t + h / 2.0, y + (k21 / 2.0) * h, y_prime + (k22 / 2.0) * h);
        let k41 = y_prime + k32 * h;
        let k42 = system.derivative(t + h, y + k31 * h, y_prime + k32 * h);
        (
            y + (k11 + k21 * 2.0 + k31 * 2.0 + k41) * (h / 6.0),
            y_prime + (k12 + k22 * 2.0 + k32 * 2.0 + k42) * (h / 6.0),
        )
    }
}

impl Integrator for LeapfrogMethod {
    fn step<S>(
        &mut self,
        system: &S,
        t: f64,
        y: S::Vector,
        y_prime: S::Vector,
        h: f64,
    ) -> (S::Vector, S::Vector)
    where
        S: System,
    {
        let y_prime = y_prime + system.derivative(t, y, y_prime) * (h / 2.0);
        let y = y + y_prime * h;
        (
            y,
            y_prime + system.derivative(t + h, y, y_prime) * (h / 2.0),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-2;

    struct TestHarmonicOscillator;

    /// y'' + y = 0
    impl System for TestHarmonicOscillator {
        type Vector = f64;

        fn derivative(&self, _t: f64, y: Self::Vector, _y_prime: Self::Vector) -> Self::Vector {
            -y
        }
    }

    #[test]
    fn test_euler_method() {
        let system = TestHarmonicOscillator;
        let method = EulerMethod;
        let y0 = 1.0;
        let y0_prime = 0.0;
        let h = 0.01;
        let steps = 100;

        let mut test_solver = Solver::new(method, system, y0, y0_prime);
        test_solver.run(h, steps);
        let (_, y, _) = test_solver.get_current();

        let exact_y = 1.0f64.cos();
        let error = (y - exact_y).abs();

        assert!(error < EPS, "Euler method error too large: {}", error);
    }

    #[test]
    fn test_rk4_method() {
        let system = TestHarmonicOscillator;
        let integrator = RK4Method;
        let y0 = 1.0;
        let y0_prime = 0.0;
        let h = 0.01;
        let steps = 100;

        let mut test_solver = Solver::new(integrator, system, y0, y0_prime);
        test_solver.run(h, steps);
        let (_, y, _) = test_solver.get_current();

        let exact_y = 1.0f64.cos();
        let error = (y - exact_y).abs();

        assert!(error < EPS, "RK4 method error too large: {}", error);
    }

    #[test]
    fn test_rk4_better_than_euler() {
        let y0 = 1.0;
        let y0_prime = 0.0;
        let h = 0.1;
        let steps = 10;

        let system = TestHarmonicOscillator;
        let integrator = EulerMethod;
        let mut euler_test_solver = Solver::new(integrator, system, y0, y0_prime);
        euler_test_solver.run(h, steps);
        let (_, y_euler, _) = euler_test_solver.get_current();

        let system = TestHarmonicOscillator;
        let rk4 = RK4Method;
        let mut rk4_test_solver = Solver::new(rk4, system, y0, y0_prime);
        rk4_test_solver.run(h, steps);
        let (_, y_rk4, _) = rk4_test_solver.get_current();

        let exact_y = 1.0f64.cos();
        let error_euler = (y_euler - exact_y).abs();
        let error_rk4 = (y_rk4 - exact_y).abs();

        assert!(
            error_rk4 < error_euler,
            "RK4 method error are large than Euler method. RK4: {}  Euler: {}",
            error_rk4,
            error_euler
        )
    }

    #[test]
    fn test_leapfrog_energy_conservation() {
        let system = TestHarmonicOscillator;
        let integrator = LeapfrogMethod;

        let y0 = 1.0;
        let y0_prime = 0.0;
        let h = 0.01;
        let steps = 1000;

        let initial_energy = 0.5 * y0_prime * y0_prime + 0.5 * y0 * y0;

        let mut test_solver = Solver::new(integrator, system, y0, y0_prime);
        test_solver.run(h, steps);
        let (_, y, y_prime) = test_solver.get_current();

        let final_energy = 0.5 * y_prime * y_prime + 0.5 * y * y;

        let error = (final_energy - initial_energy).abs();
        println!(
            "Initial Energy: {}, Final Energy: {}, Error: {}",
            initial_energy, final_energy, error
        );

        assert!(error < 1e-5, "LeapFrogMethod: Energy is not conserved");
    }

    #[test]
    fn test_leapfrog_analytic_solution() {
        let integrator = LeapfrogMethod;
        let system = TestHarmonicOscillator;

        let y0 = 1.0;
        let y0_prime = 0.0;
        let h = 0.01;
        let steps = 100;

        let mut test_solver = Solver::new(integrator, system, y0, y0_prime);
        test_solver.run(h, steps);
        let (_, y, _) = test_solver.get_current();

        let exact_y = 1.0f64.cos();
        let error = (y - exact_y).abs();

        assert!(error < 1e-5, "Leapfrog method error too large: {}", error);
    }
}
