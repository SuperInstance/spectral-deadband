//! # spectral-deadband
//!
//! Deadband as spectral gap. Spin as time-as-rotation. Fractal conservation. Analog dial computation.
//! Pure Rust, zero dependencies.

pub mod deadband;
pub mod spin;
pub mod fractal;
pub mod dial;

pub use deadband::{Deadband, Action, spectral_deadband};
pub use spin::{Spin, SpinSystem};
pub use fractal::{FractalScale, fibonacci_graph, complete_graph};
pub use dial::{AnalogDial, DialBank};
