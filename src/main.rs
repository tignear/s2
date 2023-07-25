use std::collections::HashMap;
use std::f64::consts::PI;
pub mod consts;
pub mod der;
use consts::*;
struct ComputeEnv {
    F0: f64,
    n: f64,
    cache_x: HashMap<i32, f64>,
    cache_dxdt: HashMap<i32, f64>,
}
impl ComputeEnv {
    fn x(&mut self, step: i32) -> f64 {
        let r = match self.cache_x.get(&step) {
            Some(x) => *x,
            None => {
                let t = f64::from(step) * consts::dt;
                if t < -0.3 {
                    return 0.0;
                }
                let x =
                    self.x(step - 1) + self.dxdt(step - 1) * dt + self.dxdt2(step - 1) * dt * dt;
                self.cache_x.insert(step, x);
                // println!("x@{step}:{x}");
                x
            }
        };

        r
    }

    fn dxdt2(&mut self, step: i32) -> f64 {
        let t = f64::from(step) * dt;
        self.F0 / m * self.phi_n(t) - G * A / (h_b * m) * self.x(step)
    }

    fn dxdt3(&mut self, step: i32) -> f64 {
        let t = f64::from(step) * dt;
        self.F0 / m * self.dphin_ndt(t) - G * A / (h_b * m) * self.dxdt(step)
    }

    fn dxdt(&mut self, step: i32) -> f64 {
        let t = f64::from(step) * dt;
        let r = match self.cache_dxdt.get(&step) {
            Some(dxdt) => *dxdt,
            None => {
                if t < -0.3 {
                    return 0.0;
                }

                let dxdt = self.dxdt(step - 1)
                    + self.dxdt2(step - 1) * dt
                    + self.dxdt3(step - 1) * dt * dt;
                self.cache_dxdt.insert(step, dxdt);
                // println!("dx@{step}:{dxdt}");
                dxdt
            }
        };

        r
    }
    fn phi_n(&self, t: f64) -> f64 {
        f64::sqrt(self.n / PI) * f64::exp(-self.n * t * t)
    }
    fn dphin_ndt(&self, t: f64) -> f64 {
        2.0 * t * (-self.n) * f64::sqrt(self.n / PI) * f64::exp(-self.n * t * t)
    }
}
fn compute_f0(omega: f64, n: f64) -> (f64, f64) {
    let mut F0_min = 0.0;
    let mut F0_max = 1000.0;
    let mut F0 = 0.0;
    let mut d = 0.0;
    const target_step: i32 = 300000;
    for it in 0..100 {
        eprintln!("{omega}:{it}");
        F0 = F0_max / 2.0 + F0_min / 2.0;
        let mut env = ComputeEnv {
            F0,
            n,
            cache_dxdt: HashMap::new(),
            cache_x: HashMap::new(),
        };
        for step in -target_step..target_step {
            env.dxdt(step);
        }
        let left = env.dxdt(target_step);

        let right = l * (omega - l * F0 / I);
        if left > right {
            F0_max = F0;
        } else {
            F0_min = F0;
        }
        d = (left - right).abs();
        if d < 1.0 {
            break;
        }
    }
    return (F0, d);
}

fn sub_main() {
    //2.0 * PI * 33.0;
    let n = 100.;
    for s in 0..6 {
        let omega = PI * f64::from(s) * 11.0;
        let (F0, d) = compute_f0(omega, n);
        let xs = der::solve_rk4(0.01, 10000, n, F0);
        println!(
            "{},{},{},{:?},{}",
            omega,
            F0,
            xs.iter().map(|(t, x)| *x).fold(0.0 / 0.0, f64::max),
            xs,
            d
        );
    }
}
fn main() {
    let xs = der::solve_rk4(0.3, 10000, 100.0, 100000.);
    xs.iter().for_each(|(t,x)|println!("{t}\t{x}"));
    return;
    const STACK_SIZE: usize = 512 * 1024 * 1024;
    std::thread::Builder::new()
        .stack_size(STACK_SIZE)
        .spawn(sub_main)
        .unwrap()
        .join()
        .unwrap();
}
