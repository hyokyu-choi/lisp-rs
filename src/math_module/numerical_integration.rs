use super::vector::Vector2;
use crate::physics_module::constants::PI;
use crate::utils::print_title;

pub type IntegrationMethod<T> = fn(DerivativeEquation<T>, f64, T, f64) -> T;
pub type DerivativeEquation<T> = fn(f64, T) -> f64;

/// Euler Method
///
/// 초기값 t0, y0
/// dy/dt = f(t, y)
///
/// y_i+1 = y_i + h*f(t, y)
///
/// euler_method(f, t, y, h) -> y_i+1
pub fn euler_method(f: DerivativeEquation<f64>, t: f64, y: f64, h: f64) -> f64 {
    y + h * f(t, y)
}

/// Euler Method at 2nd order ODE
///
/// 초기값 t0, y0, y'0
/// d^2y/dt^2 = f(t, y, y')
///
/// du/dt = {u1, f(t, y, y')} = {u1, f(t, u)}
/// u_i+1 = u_i + h*f(t, u)
///
/// euler_method(f, t, y, h) -> y_i+1
pub fn euler_method_2(f: DerivativeEquation<Vector2>, t: f64, u: Vector2, h: f64) -> Vector2 {
    u + h * Vector2::new(u.y(), f(t, u))
}

/// 4th Runge-Kutta
///
/// 초기값 t0, y0
/// dy/dt = f(t, y)
///
/// k1 = f(t_i, y_i)
/// k2 = f(t_i + h/2, y_i + h*k1/2)
/// k3 = f(t_i + h/2, y_i + h*k2/2)
/// k4 = f(t_i + h, y_i + h*k3)
/// y_i+1 = y_i + h/6 * (k1 + 2*k2 + 2*k3 + k4)
///
/// RK4(f, t, y, h) -> y_i+1
pub fn runge_kutta_4th(f: DerivativeEquation<f64>, t: f64, y: f64, h: f64) -> f64 {
    let k1 = f(t, y);
    let k2 = f(t + h / 2.0, y + h * k1 / 2.0);
    let k3 = f(t + h / 2.0, y + h * k2 / 2.0);
    let k4 = f(t + h, y + h * k3);
    y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

/// 4th Runge-Kutta at 2nd order ODE
///
/// 초기값 t0, y0, y'0
/// d^2y/dt^2 = f(t, y, y') = f(t, u)
///
/// k1 = f(t_i, u_i)
/// k2 = f(t_i + h/2, u_i + h*k1/2)
/// k3 = f(t_i + h/2, u_i + h*k2/2)
/// k4 = f(t_i + h, u_i + h*k3)
/// u_i+1 = u_i + h/6 * (k1 + 2*k2 + 2*k3 + k4)
///
/// RK4(f, t, u, h) -> u_i+1
pub fn runge_kutta_4th_2(f: DerivativeEquation<Vector2>, t: f64, u: Vector2, h: f64) -> Vector2 {
    let k1 = Vector2::new(u.y(), f(t, u));
    let k2 = Vector2::new(u.y(), f(t + h / 2.0, u + h * k1 / 2.0));
    let k3 = Vector2::new(u.y(), f(t + h / 2.0, u + h * k2 / 2.0));
    let k4 = Vector2::new(u.y(), f(t + h, u + h * k3));
    u + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

pub fn print_results(ts: Vec<f64>, ys: Vec<f64>, h: f64) {
    println!("   t :            y");
    for i in 0..ts.len() {
        if i as i64 % (0.5/h) as i64 == 0 {
            println!("{:.2} : {:.10}", ts[i], ys[i])
        }
    }
}

pub fn print_results_2(ts: Vec<f64>, us: Vec<Vector2>, h: f64) {
    println!("   t :            u");
    for i in 0..ts.len() {
        if i as i64 % (0.5/h) as i64 == 0 {
            println!("{:.2} : {}", ts[i], us[i])
        }
    }
}

pub fn integrate_step(
    f: DerivativeEquation<f64>,
    integrate_method: IntegrationMethod<f64>,
    y0: f64,
    t_start: f64,
    t_end: f64,
    h: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut t = t_start; // 초기값
    let mut y = y0; // 초기값
    let mut ts = vec![t];
    let mut ys = vec![y];
    while t <= t_end {
        y = integrate_method(f, t, y, h);
        t += h;
        ts.push(t);
        ys.push(y);
    }
    (ts, ys)
}

pub fn integrate_step_2(
    f: DerivativeEquation<Vector2>,
    integrate_method: IntegrationMethod<Vector2>,
    u0: Vector2,
    t_start: f64,
    t_end: f64,
    h: f64,
) -> (Vec<f64>, Vec<Vector2>) {
    let mut t = t_start; // 초기값
    let mut u = u0; // 초기값
    let mut ts = vec![t];
    let mut us = vec![u];
    while t <= t_end {
        u = integrate_method(f, t, u, h);
        t += h;
        ts.push(t);
        us.push(u);
    }
    (ts, us)
}

/// y' - y = 0
/// y = Ae^t
/// dy/dt = y
pub fn test(h: f64) {
    let test_func = |_t: f64, y: f64| -> f64 { y };
    print_title("Euler Method");
    let (ts, ys) = integrate_step(test_func, euler_method, 1.0, 0.0, 2.0, h);
    print_results(ts, ys, h);

    println!();
    print_title("RK4");
    let (ts, ys) = integrate_step(test_func, runge_kutta_4th, 1.0, 0.0, 2.0, h);
    print_results(ts, ys, h);
}

/// y'' + y = 0
/// dy'/dt = -u0
/// y = A*sin(t) + B*cos(t)
/// y' = A*cos(t) - B*sin(t)
pub fn test_2(h: f64) {
    let test_func_2 = |_t: f64, u: Vector2| -> f64 { -u.x() };

    print_title("Euler Method");
    let (ts, ys) = integrate_step_2(
        test_func_2,
        euler_method_2,
        Vector2::new(0.0, 1.0),
        0.0,
        PI / 2.0,
        h,
    );
    print_results_2(ts, ys, h);

    println!();
    print_title("RK4");
    let (ts, ys) = integrate_step_2(
        test_func_2,
        runge_kutta_4th_2,
        Vector2::new(0.0, 1.0),
        0.0,
        PI / 2.0,
        h,
    );
    print_results_2(ts, ys, h);
}
