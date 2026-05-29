//! Module 1: Deadband — the spectral gap as a region of inaction.

/// Action decided by the thermostat metaphor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Action {
    Heat,
    Cool,
    Off,
}

/// A deadband is a region where the system does nothing.
///
/// In spectral terms: the gap between eigenvalues.
/// In physical terms: the mercury band where the thermostat doesn't switch.
/// In algorithmic terms: the tolerance where no action is needed.
#[derive(Debug, Clone, Copy)]
pub struct Deadband {
    pub lower: f64,
    pub upper: f64,
    /// Physical analog: how much force needed to escape the band.
    pub gravity: f64,
}

impl Deadband {
    pub fn new(center: f64, width: f64) -> Self {
        let half = width / 2.0;
        Deadband {
            lower: center - half,
            upper: center + half,
            gravity: width, // default gravity proportional to width
        }
    }

    pub fn with_gravity(mut self, gravity: f64) -> Self {
        self.gravity = gravity;
        self
    }

    /// Is the value inside the deadband?
    pub fn contains(&self, value: f64) -> bool {
        value >= self.lower && value <= self.upper
    }

    /// Width of the deadband.
    pub fn width(&self) -> f64 {
        self.upper - self.lower
    }

    /// Center of the deadband.
    pub fn center(&self) -> f64 {
        (self.lower + self.upper) / 2.0
    }

    /// Deadband ratio: width / range = spectral gap analog.
    pub fn deadband_ratio(&self, range: f64) -> f64 {
        if range.abs() < 1e-15 { return 0.0; }
        self.width() / range
    }

    /// The thermostat: returns None if in deadband, Some(Action) if outside.
    pub fn evaluate(&self, value: f64) -> Option<Action> {
        if self.contains(value) {
            None
        } else if value < self.lower {
            Some(Action::Heat)
        } else {
            Some(Action::Cool)
        }
    }
}

/// Compute a spectral deadband from a sorted list of eigenvalues.
///
/// The deadband width IS the spectral gap (distance between two eigenvalues).
/// The center IS the midpoint. The gravity IS related to eigenvalue magnitude.
pub fn spectral_deadband(eigenvalues: &[f64]) -> Deadband {
    if eigenvalues.len() < 2 {
        return Deadband::new(0.0, 0.0);
    }

    // Find the largest gap between consecutive eigenvalues
    let mut sorted: Vec<f64> = eigenvalues.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut max_gap = 0.0_f64;
    let mut gap_center = 0.0_f64;
    let mut gap_idx = 0usize;

    for i in 0..sorted.len() - 1 {
        let gap = sorted[i + 1] - sorted[i];
        if gap > max_gap {
            max_gap = gap;
            gap_center = (sorted[i] + sorted[i + 1]) / 2.0;
            gap_idx = i;
        }
    }

    // Gravity proportional to the eigenvalue magnitudes flanking the gap
    let gravity = (sorted[gap_idx].abs() + sorted[gap_idx + 1].abs()) / 2.0;

    Deadband {
        lower: gap_center - max_gap / 2.0,
        upper: gap_center + max_gap / 2.0,
        gravity,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deadband_contains() {
        let db = Deadband::new(5.0, 2.0); // [4.0, 6.0]
        assert!(db.contains(5.0));
        assert!(db.contains(4.0));
        assert!(db.contains(6.0));
        assert!(!db.contains(3.9));
        assert!(!db.contains(6.1));
    }

    #[test]
    fn test_deadband_evaluate_inside() {
        let db = Deadband::new(5.0, 2.0);
        assert_eq!(db.evaluate(5.0), None);
    }

    #[test]
    fn test_deadband_evaluate_below() {
        let db = Deadband::new(5.0, 2.0);
        assert_eq!(db.evaluate(3.0), Some(Action::Heat));
    }

    #[test]
    fn test_deadband_evaluate_above() {
        let db = Deadband::new(5.0, 2.0);
        assert_eq!(db.evaluate(7.0), Some(Action::Cool));
    }

    #[test]
    fn test_spectral_deadband_simple() {
        // Eigenvalues with a clear gap between 3.0 and 7.0
        let evals = vec![1.0, 2.0, 3.0, 7.0, 8.0];
        let db = spectral_deadband(&evals);
        // Gap is 4.0 between 3.0 and 7.0, center = 5.0
        assert!((db.center() - 5.0).abs() < 1e-10);
        assert!((db.width() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_spectral_deadband_single() {
        let db = spectral_deadband(&[1.0]);
        assert_eq!(db.width(), 0.0);
    }

    #[test]
    fn test_deadband_ratio() {
        let db = Deadband::new(5.0, 2.0);
        assert!((db.deadband_ratio(10.0) - 0.2).abs() < 1e-10);
    }
}
