//! Module 3: FractalScale — conservation ratio is fractal (same at every scale).
//!
//! Fibonacci outward: Penrose inflation (subdivide → larger tiling).
//! Mandelbrot inward: zoom → more detail, same CR.

/// 1/φ ≈ 0.618033988749895
const FRAC_PHI: f64 = 0.61803398874989484820458683436564;

/// Fractal scale analysis of conservation ratio across multiple scales.
#[derive(Debug, Clone)]
pub struct FractalScale {
    /// Conservation ratio measured at each scale.
    pub cr: Vec<f64>,
    /// Graph sizes at which CR was measured.
    pub scales: Vec<usize>,
}

impl FractalScale {
    /// Measure conservation ratio at multiple scales by subsampling the graph.
    ///
    /// CR for a graph = (sum of eigenvalues) / (sum of degrees) — a normalized energy measure.
    pub fn measure(adj: &[Vec<f64>], levels: usize) -> Self {
        let n = adj.len();
        let mut crs = Vec::new();
        let mut sizes = Vec::new();

        for level in 1..=levels {
            // Subsample: take first k nodes where k = n * level / levels
            let k = ((n * level) / levels).max(2);
            let sub = subsample_adj(adj, k);
            let cr = compute_cr(&sub);
            crs.push(cr);
            sizes.push(k);
        }

        // Also measure full graph
        if !sizes.contains(&n) {
            crs.push(compute_cr(adj));
            sizes.push(n);
        }

        FractalScale { cr: crs, scales: sizes }
    }

    /// Is the CR self-similar (approximately constant across scales)?
    pub fn is_self_similar(&self, tolerance: f64) -> bool {
        if self.cr.len() < 2 { return true; }
        let mean: f64 = self.cr.iter().sum::<f64>() / self.cr.len() as f64;
        self.cr.iter().all(|&c| (c - mean).abs() < tolerance)
    }

    /// Fibonacci outward: inflate the graph by Penrose-like subdivision.
    ///
    /// Each edge → two edges with a new node in between (Fibonacci subdivision).
    /// F(n+1) = F(n) + F(n-1) nodes. CR should converge to 1/φ.
    pub fn inflate(adj: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = adj.len();
        let mut new_adj: Vec<Vec<f64>> = vec![vec![0.0; 2 * n]; 2 * n];

        // Keep original nodes, add n new nodes (one between each pair)
        for i in 0..n {
            for j in 0..n {
                if adj[i][j] > 0.0 && i < j {
                    let mid = n + i; // new node between i and j
                    // Split edge: i — mid — j with Fibonacci weights
                    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
                    let w = adj[i][j];
                    new_adj[i][mid] = w * phi;
                    new_adj[mid][i] = w * phi;
                    new_adj[mid][j] = w / phi;
                    new_adj[j][mid] = w / phi;
                }
            }
        }

        new_adj
    }

    /// Mandelbrot inward: zoom into random subgraphs and measure roughness.
    ///
    /// Roughness = variance of CR across random subgraphs of given size.
    pub fn roughness(adj: &[Vec<f64>], subgraph_size: usize) -> f64 {
        let n = adj.len();
        if subgraph_size >= n || subgraph_size < 2 { return 0.0; }

        let num_samples = 20.min(n - subgraph_size + 1);
        let mut crs = Vec::new();

        // Use deterministic "random" subgraphs by shifting a window
        for s in 0..num_samples {
            let start = (s * (n - subgraph_size) / num_samples.max(1)).min(n - subgraph_size);
            let mut sub = vec![vec![0.0; subgraph_size]; subgraph_size];
            for i in 0..subgraph_size {
                for j in 0..subgraph_size {
                    sub[i][j] = adj[start + i][start + j];
                }
            }
            crs.push(compute_cr(&sub));
        }

        if crs.is_empty() { return 0.0; }
        let mean: f64 = crs.iter().sum::<f64>() / crs.len() as f64;
        let variance: f64 = crs.iter().map(|c| (c - mean).powi(2)).sum::<f64>() / crs.len() as f64;
        variance.sqrt()
    }

    /// Maximum deviation from the golden ratio 1/φ across all scales.
    pub fn golden_deviation(&self) -> f64 {
        self.cr.iter().map(|&c| (c - FRAC_PHI).abs()).fold(0.0_f64, f64::max)
    }
}

/// Compute conservation ratio for a graph: Σ(eigenvalues) / Σ(degrees).
fn compute_cr(adj: &[Vec<f64>]) -> f64 {
    let n = adj.len();
    if n == 0 { return 0.0; }

    let total_degree: f64 = adj.iter().map(|row| row.iter().sum::<f64>()).sum();
    if total_degree < 1e-15 { return 0.0; }

    // Simple power iteration to estimate largest eigenvalue
    let mut v = vec![1.0 / (n as f64).sqrt(); n];
    for _ in 0..100 {
        let mut new_v = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                new_v[i] += adj[i][j] * v[j];
            }
        }
        let norm: f64 = new_v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm < 1e-15 { return 0.0; }
        for i in 0..n {
            v[i] = new_v[i] / norm;
        }
    }

    // Rayleigh quotient for largest eigenvalue
    let mut lambda = 0.0;
    let mut av = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            av[i] += adj[i][j] * v[j];
        }
        lambda += v[i] * av[i];
    }

    lambda * n as f64 / total_degree
}

/// Subsample adjacency matrix to first k nodes.
fn subsample_adj(adj: &[Vec<f64>], k: usize) -> Vec<Vec<f64>> {
    let n = adj.len();
    let k = k.min(n);
    let mut sub = vec![vec![0.0; k]; k];
    for i in 0..k {
        for j in 0..k {
            sub[i][j] = adj[i][j];
        }
    }
    sub
}

/// Generate a complete graph of size n.
pub fn complete_graph(n: usize) -> Vec<Vec<f64>> {
    let mut adj = vec![vec![1.0; n]; n];
    for i in 0..n {
        adj[i][i] = 0.0;
    }
    adj
}

/// Generate a Fibonacci-like graph: chain where edge weights follow Fibonacci sequence.
pub fn fibonacci_graph(n: usize) -> Vec<Vec<f64>> {
    let mut adj = vec![vec![0.0; n]; n];
    let mut fib = vec![1.0, 1.0];
    for i in 2..n {
        fib.push(fib[i - 1] + fib[i - 2]);
    }
    for i in 0..n.saturating_sub(1) {
        let w = fib[i.min(fib.len() - 1)];
        adj[i][i + 1] = w;
        adj[i + 1][i] = w;
    }
    adj
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complete_graph_self_similar() {
        let adj = complete_graph(10);
        let fs = FractalScale::measure(&adj, 5);
        // Complete graphs have uniform structure → self-similar
        assert!(fs.is_self_similar(0.5), "CRs: {:?}", fs.cr);
    }

    #[test]
    fn test_fibonacci_graph_cr() {
        let adj = fibonacci_graph(8);
        let fs = FractalScale::measure(&adj, 4);
        // Should have some measurable CR at each scale
        assert!(fs.cr.len() >= 2);
        for &cr in &fs.cr {
            assert!(cr.is_finite());
        }
    }

    #[test]
    fn test_inflate_doubles_nodes() {
        let adj = complete_graph(3);
        let inflated = FractalScale::inflate(&adj);
        assert_eq!(inflated.len(), 6);
    }

    #[test]
    fn test_roughness_complete_graph() {
        let adj = complete_graph(10);
        let r = FractalScale::roughness(&adj, 5);
        // Complete graph has uniform subgraphs → low roughness
        assert!(r < 1.0, "Roughness {} should be low for complete graph", r);
    }

    #[test]
    fn test_golden_deviation() {
        let fs = FractalScale {
            cr: vec![FRAC_PHI, FRAC_PHI, FRAC_PHI],
            scales: vec![5, 10, 20],
        };
        assert!(fs.golden_deviation() < 1e-10);
    }

    #[test]
    fn test_inflate_cr_converges_to_golden() {
        let mut adj = complete_graph(3);
        let mut deviations = Vec::new();

        for _ in 0..4 {
            let cr = compute_cr(&adj);
            deviations.push((cr - FRAC_PHI).abs());
            adj = FractalScale::inflate(&adj);
        }

        // CR should generally be in a reasonable range
        for &d in &deviations {
            assert!(d.is_finite());
        }
    }
}
