pub mod numerical_integration;

use numerical_integration::{integrate_step, euler_method, runge_kutta_4th, print_results};

fn main() {
    print_title("Euler Method");
    let (ts, ys) = integrate_step(myfunc, euler_method, 1.0, 0.0, 1.0, 0.1);
    print_results(ts, ys);

    println!();
    print_title("RK4");
    let (ts, ys) = integrate_step(myfunc, runge_kutta_4th, 1.0, 0.0, 1.0, 0.1);
    print_results(ts, ys);
}

fn myfunc(_t: f64, y: f64) -> f64 {
    // y' - y = 0
    // y = Ae^t
    y
}

fn print_title(title: &str) {
    println!("====== {title} ======")
}
