mod utils{
    pub mod plot;
}

mod physics_module {
    pub mod harmonic_oscillator;
}

mod math_module {
    pub mod core;
    pub mod integrate;
}

use math_module::core::Scalar;
use math_module::integrate::{RK4Method, Solver};
use physics_module::harmonic_oscillator::{DHO};
use utils::plot::plot_one;

fn main() {
    let dho_ode = DHO {
        omega_square: Scalar(7.0),
        b: Scalar(1.0),
    };

    let y0 = Scalar(1.5);
    let y0_prime = Scalar(0.0);
    let h = Scalar(0.01);
    let steps = 1000;

    let mut sho_solver = Solver::new(RK4Method, dho_ode, y0, y0_prime);
    sho_solver.run(h, steps);
    let ts = sho_solver.get_ts_f64();
    let ys = sho_solver.get_ys_f64().iter().map(|x| x[0]).collect();
    let ys_prime = sho_solver.get_ys_prime_f64().iter().map(|x| x[0]).collect();

    plot_one("Damped Harmonic Oscillation", ts, vec![ys, ys_prime], vec!["position", "velocity"], "plotters-doc-data/1.png");
}
