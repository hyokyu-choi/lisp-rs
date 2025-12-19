mod math_module {
    pub mod numerical_integration;
    pub mod numerical_integration_second_order;
    pub mod vector;
}

use std::f64::consts::PI;
// use math_module::numerical_integration::{integrate_step, euler_method, runge_kutta_4th, print_results};
use math_module::numerical_integration_second_order::{
    euler_method, integrate_step, print_results, runge_kutta_4th,
};
use math_module::vector::Vector2;

fn main() {
    print_title("Euler Method");
    let (ts, ys) = integrate_step(
        myfunc,
        euler_method,
        Vector2::new(0.0, 1.0),
        0.0,
        PI / 2.0,
        0.1,
    );
    print_results(ts, ys);

    println!();
    print_title("RK4");
    let (ts, ys) = integrate_step(
        myfunc,
        runge_kutta_4th,
        Vector2::new(0.0, 1.0),
        0.0,
        PI / 2.0,
        0.1,
    );
    print_results(ts, ys);
}

fn myfunc(_t: f64, u: Vector2) -> f64 {
    // y'' + y = 0
    // y = A*sin(t) + B*cos(t)
    // y' = A*cos(t) - B*sin(t)
    // dy/dt = u1
    // dy'/dt = -u0
    -u.x()
}

fn print_title(title: &str) {
    println!("========= {title} =========")
}
