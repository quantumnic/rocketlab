//! ASCII plotting for trajectory visualization.
//!
//! Terminal-friendly plots for altitude vs time, velocity vs time, etc.

/// A simple ASCII plot.
pub struct AsciiPlot {
    width: usize,
    height: usize,
    title: String,
    x_label: String,
    y_label: String,
}

impl AsciiPlot {
    pub fn new(title: &str, x_label: &str, y_label: &str) -> Self {
        Self {
            width: 70,
            height: 20,
            title: title.to_string(),
            x_label: x_label.to_string(),
            y_label: y_label.to_string(),
        }
    }

    pub fn with_size(mut self, width: usize, height: usize) -> Self {
        self.width = width;
        self.height = height;
        self
    }

    /// Render a single series as ASCII art.
    pub fn plot(&self, data: &[(f64, f64)]) -> String {
        if data.is_empty() {
            return String::from("(no data)");
        }

        let (x_min, x_max, y_min, y_max) = bounds(data);

        // Avoid zero ranges
        let x_range = if (x_max - x_min).abs() < 1e-10 {
            1.0
        } else {
            x_max - x_min
        };
        let y_range = if (y_max - y_min).abs() < 1e-10 {
            1.0
        } else {
            y_max - y_min
        };

        // Build grid
        let mut grid = vec![vec![' '; self.width]; self.height];

        // Plot points
        for &(x, y) in data {
            let col = ((x - x_min) / x_range * (self.width - 1) as f64).round() as usize;
            let row = ((y - y_min) / y_range * (self.height - 1) as f64).round() as usize;
            let col = col.min(self.width - 1);
            let row = row.min(self.height - 1);
            // Invert row (top = high y)
            let row = self.height - 1 - row;
            grid[row][col] = '█';
        }

        // Render
        let mut output = String::new();
        output.push_str(&format!("  {}\n", self.title));
        output.push_str(&format!("  {:>10} ┐\n", format_number(y_max)));

        let y_label_row = self.height / 2;
        for (i, row) in grid.iter().enumerate() {
            if i == y_label_row {
                output.push_str(&format!("{:>10} │", &self.y_label));
            } else {
                output.push_str("           │");
            }
            let line: String = row.iter().collect();
            output.push_str(&line);
            output.push_str("│\n");
        }

        output.push_str(&format!("  {:>10} ┘", format_number(y_min)));
        let x_axis_pad = self.width.saturating_sub(20);
        output.push_str(&"─".repeat(self.width));
        output.push('\n');
        output.push_str(&format!(
            "           {:<10}{}  {}\n",
            format_number(x_min),
            " ".repeat(x_axis_pad),
            format_number(x_max),
        ));
        output.push_str(&format!(
            "           {:^width$}\n",
            self.x_label,
            width = self.width
        ));

        output
    }

    /// Render two series overlaid (using different characters).
    pub fn plot_dual(
        &self,
        data1: &[(f64, f64)],
        label1: &str,
        data2: &[(f64, f64)],
        label2: &str,
    ) -> String {
        if data1.is_empty() && data2.is_empty() {
            return String::from("(no data)");
        }

        let all: Vec<(f64, f64)> = data1.iter().chain(data2.iter()).copied().collect();
        let (x_min, x_max, y_min, y_max) = bounds(&all);

        let x_range = if (x_max - x_min).abs() < 1e-10 {
            1.0
        } else {
            x_max - x_min
        };
        let y_range = if (y_max - y_min).abs() < 1e-10 {
            1.0
        } else {
            y_max - y_min
        };

        let mut grid = vec![vec![' '; self.width]; self.height];

        for &(x, y) in data1 {
            let col = ((x - x_min) / x_range * (self.width - 1) as f64).round() as usize;
            let row = ((y - y_min) / y_range * (self.height - 1) as f64).round() as usize;
            let col = col.min(self.width - 1);
            let row = (self.height - 1).saturating_sub(row.min(self.height - 1));
            grid[row][col] = '█';
        }

        for &(x, y) in data2 {
            let col = ((x - x_min) / x_range * (self.width - 1) as f64).round() as usize;
            let row = ((y - y_min) / y_range * (self.height - 1) as f64).round() as usize;
            let col = col.min(self.width - 1);
            let row = (self.height - 1).saturating_sub(row.min(self.height - 1));
            if grid[row][col] == '█' {
                grid[row][col] = '▓'; // overlap
            } else {
                grid[row][col] = '░';
            }
        }

        let mut output = String::new();
        output.push_str(&format!("  {}\n", self.title));
        output.push_str(&format!("  █ {}  ░ {}\n", label1, label2));
        output.push_str(&format!("  {:>10} ┐\n", format_number(y_max)));

        for row in &grid {
            output.push_str("           │");
            let line: String = row.iter().collect();
            output.push_str(&line);
            output.push_str("│\n");
        }

        output.push_str(&format!("  {:>10} ┘", format_number(y_min)));
        output.push_str(&"─".repeat(self.width));
        output.push('\n');
        output.push_str(&format!(
            "           {:<10}{}  {}\n",
            format_number(x_min),
            " ".repeat(self.width.saturating_sub(20)),
            format_number(x_max),
        ));

        output
    }
}

fn bounds(data: &[(f64, f64)]) -> (f64, f64, f64, f64) {
    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;

    for &(x, y) in data {
        x_min = x_min.min(x);
        x_max = x_max.max(x);
        y_min = y_min.min(y);
        y_max = y_max.max(y);
    }

    (x_min, x_max, y_min, y_max)
}

fn format_number(n: f64) -> String {
    let abs = n.abs();
    if abs >= 1e6 {
        format!("{:.1}M", n / 1e6)
    } else if abs >= 1e3 {
        format!("{:.1}k", n / 1e3)
    } else if abs >= 1.0 {
        format!("{:.1}", n)
    } else if abs >= 0.001 {
        format!("{:.4}", n)
    } else {
        format!("{:.2e}", n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_plot() {
        let data: Vec<(f64, f64)> = (0..100)
            .map(|i| {
                let t = i as f64;
                (t, t * t)
            })
            .collect();

        let plot = AsciiPlot::new("Test Plot", "Time (s)", "Value")
            .with_size(40, 10)
            .plot(&data);

        assert!(!plot.is_empty());
        assert!(plot.contains("Test Plot"));
        assert!(plot.contains('█'));
    }

    #[test]
    fn test_empty_plot() {
        let plot = AsciiPlot::new("Empty", "x", "y").plot(&[]);
        assert_eq!(plot, "(no data)");
    }

    #[test]
    fn test_dual_plot() {
        let d1: Vec<(f64, f64)> = (0..50).map(|i| (i as f64, i as f64)).collect();
        let d2: Vec<(f64, f64)> = (0..50).map(|i| (i as f64, 50.0 - i as f64)).collect();

        let plot = AsciiPlot::new("Dual", "x", "y")
            .with_size(40, 10)
            .plot_dual(&d1, "Series 1", &d2, "Series 2");

        assert!(plot.contains("Series 1"));
        assert!(plot.contains("Series 2"));
    }

    #[test]
    fn test_format_number() {
        assert_eq!(format_number(1_500_000.0), "1.5M");
        assert_eq!(format_number(5_000.0), "5.0k");
        assert_eq!(format_number(42.5), "42.5");
    }
}
