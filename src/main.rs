mod utils {
    pub mod plot;
}

mod physics {
    pub mod harmonic_oscillator;
}

mod math {
    pub mod core;
    pub mod integrate;
}

use math::core::Scalar;
use math::integrate::{RK4Method, Solver};
use physics::harmonic_oscillator::DampedHarmonicOscillator;
use utils::plot::plot_one;

fn main() {
    let dho_ode: DampedHarmonicOscillator = DampedHarmonicOscillator {
        omega_square: Scalar(15.0),
        b: Scalar(1.0),
    };

    let y0: Scalar = Scalar(1.5);
    let y0_prime: Scalar = Scalar(0.0);
    let h: Scalar = Scalar(0.01);
    let steps: usize = 1000;

    let mut sho_solver: Solver<RK4Method, DampedHarmonicOscillator> = Solver::new(RK4Method, dho_ode, y0, y0_prime);
    sho_solver.run(h, steps);
    let (ts, ys, ys_prime): (Vec<f64>, Vec<f64>, Vec<f64>) = sho_solver.get_results_f64();

    let _ = plot_one(
        "Damped Harmonic Oscillation",
        ts,
        vec![ys, ys_prime],
        vec!["position", "velocity"],
        "plotters-doc-data/1.png",
    );
}
