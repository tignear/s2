use std::f64::consts::PI;

pub const dt: f64 = 0.0000001;
pub const G: f64 = 79E9;
pub const d_d: f64 = 0.15;
pub const A: f64 = PI / 4.0 * d_d * d_d;
pub const m: f64 = 6.0;
pub const h_b: f64 = 0.05;
pub const l_a: f64 = 0.6;
pub const l_b: f64 = 0.93;
pub const d_B: f64 = 0.5;
pub const l: f64 = d_B / 2.0 + (l_a + l_b) - 0.05;
pub const M_B: f64 = 70.0;
pub const d_a: f64 = 0.1;
pub const M_a: f64 = 4.0;
pub const I_Ga: f64 = (d_a * d_a / 16.0 + l_a * l_a / 12.0) * M_a;
pub const I_GB: f64 = 1.0 / 8.0 * M_B * d_B * d_B;
pub const M_b: f64 = 1.0;
pub const I_Gb: f64 = 1.0 / 12.0 * M_b * (h_b * h_b + l_b * l_b);
pub const I: f64 = I_GB
    + (I_Ga + 1.0 / 4.0 * M_a * (d_B + l_a) * (d_B + l_a))
    + I_Gb
    + 1.0 / 4.0 * M_b * (d_B + 2.0 * l_a + l_b) * (d_B + 2.0 * l_a + l_b);