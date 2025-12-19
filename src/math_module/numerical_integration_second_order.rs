use super::vector::Vector2;

pub type IntegrationMethod = fn(DerivativeEquation, f64, Vector2, f64) -> Vector2;
pub type DerivativeEquation = fn(f64, Vector2) -> f64;

/// Euler Method
///
/// 초기값 t0, y0, y'0
/// d^2y/dt^2 = f(t, y, y')
///
/// u_vec = {y, y'}
/// du_vec/dt = {u1, f(t, u0, u1)} = {u1, f(t, u_vec)}
/// (update) u_i+1 = u_i + h*(du_vec/dt)
/// y_i+1 = y_i + h*y'_i
/// y'+i+1 = y'_i + h*f(t, u_vec)
///
/// euler_method(f, t, y, h) -> y_i+1
pub fn euler_method(f: DerivativeEquation, t: f64, u: Vector2, h: f64) -> Vector2 {
    u + h * Vector2::new(u.y(), f(t, u))
}

/// 4th Runge-Kutta
/// 초기값 t0, y0, y'0
/// d^2y/dt^2 = f(t, y, y') = f(t, u)
/// k1 = f(t_i, u_i)
/// k2 = f(t_i + h/2, u_i + h*k1/2)
/// k3 = f(t_i + h/2, u_i + h*k2/2)
/// k4 = f(t_i + h, u_i + h*k3)
/// u_i+1 = u_i + h/6 * (k1 + 2*k2 + 2*k3 + k4)
///
/// RK4(f, t, u, h) -> u_i+1
pub fn runge_kutta_4th(f: DerivativeEquation, t: f64, u: Vector2, h: f64) -> Vector2 {
    let k1 = Vector2::new(u.y(), f(t, u));
    let k2 = Vector2::new(u.y(), f(t + h / 2.0, u + h * k1 / 2.0));
    let k3 = Vector2::new(u.y(), f(t + h / 2.0, u + h * k2 / 2.0));
    let k4 = Vector2::new(u.y(), f(t + h, u + h * k3));
    u + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

pub fn print_results(ts: Vec<f64>, us: Vec<Vector2>) {
    println!("   t :            u");
    for (t, u) in ts.iter().zip(us.iter()) {
        println!("{t:.2} : {u}")
    }
}

pub fn integrate_step(
    f: DerivativeEquation,
    integrate_method: IntegrationMethod,
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
