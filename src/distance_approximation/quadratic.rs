use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum InvertedQuadraticFunction {
    Quadratic {
        // 1 / a
        multiplier: f64,
        // - b / (2a)
        delta: f64,
        // c / a
        offset: f64,
        right_arc: bool,
    },
    Linear {
        n: f64,
        m: f64,
    },
    ConstantOutput {
        // like a vertical line
        c: f64,
    },
}

impl InvertedQuadraticFunction {
    // left_arc, right_arc, middle
    pub fn with(a: f64, b: f64, c: f64) -> Option<(Self, Self, f64)> {
        if a.abs() > 1e-8 {
            let multiplier = a.recip();
            let delta = -b * 0.5 * multiplier;
            let offset = c * multiplier;
            Some((
                Self::Quadratic {
                    multiplier,
                    delta,
                    offset,
                    right_arc: false,
                },
                Self::Quadratic {
                    multiplier,
                    delta,
                    offset,
                    right_arc: true,
                },
                delta,
            ))
        } else if b.abs() > 1e-8 {
            let inverse = Self::Linear { m: b, n: c };
            Some((inverse.clone(), inverse, 0.0))
        } else {
            None
        }
    }

    pub fn constant_output(c: f64) -> Self {
        Self::ConstantOutput { c }
    }

    fn base_multiplier(&self) -> f64 {
        match self {
            Self::Quadratic {
                right_arc: true, ..
            } => 1.0,
            Self::Quadratic { .. } => -1.0,
            _ => 0.0,
        }
    }

    pub fn x_value_for(&self, y_value: f64) -> Option<f64> {
        match &self {
            Self::Quadratic {
                multiplier,
                delta,
                offset,
                ..
            } => {
                let radicand = y_value * multiplier + delta.powi(2) - offset;
                if radicand >= 0.0 {
                    Some((radicand).sqrt() * self.base_multiplier() + delta)
                } else if radicand.abs() < 1e-3 {
                    Some(*delta)
                } else {
                    None
                }
            }
            Self::Linear { n, m } => Some((y_value - n) / m),
            Self::ConstantOutput { c } => Some(*c),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn low_precision_eq(a: f64, b: f64) -> bool {
        // use own epsilon since this problem has really bad fp behavior
        const EPS: f64 = 1e-3;

        if (a - b).abs() < EPS {
            true
        } else if a.abs().max(b.abs()) >= 1.0 {
            (a - b).abs() / a < EPS || (a - b).abs() / b < EPS
        } else {
            false
        }
    }

    proptest! {
        #[test]
        fn middle_is_actually_the_middle(
            a in -1e6f64..1e6f64,
            b in -1e6f64..1e6f64,
            c in -1e6f64..1e6f64,
        ) {
            prop_assume!(a.abs() > 1e-8);

            let (left, right, middle) = InvertedQuadraticFunction::with(a, b, c).unwrap();
            let ym = a * middle* middle + b * middle + c;
            let left_x = left.x_value_for(ym).unwrap();
            let right_x = right.x_value_for(ym).unwrap();

            prop_assert!(
                low_precision_eq(left_x, middle) && low_precision_eq(middle, right_x),
                "All of left inverted middle {}, middle {}, right inverted middle {} should be equal",
                left_x,
                middle,
                right_x
            );
        }

        #[test]
        fn inverted_linear_function_is_actually_inverse(
            b in -1e6f64..1e6f64,
            c in -1e6f64..1e6f64,
            y in -1e6f64..1e6f64,
        ) {
            prop_assume!(b.abs() > 1e-8);

            let (left, right, _) = InvertedQuadraticFunction::with(0.0, b, c).unwrap();
            for arc in &[left, right] {
                let x_value = arc.x_value_for(y).unwrap();
                let y_value = b * x_value + c;
                prop_assert!(
                    low_precision_eq(y_value, y),
                    "given y {} should be equal to calculated y {}",
                    y,
                    y_value
                );
            }
        }

        #[test]
        fn inverted_quadratic_function_is_actually_inverse(
            a in -1e6f64..1e6f64,
            b in -1e6f64..1e6f64,
            c in -1e6f64..1e6f64,
            y in -1e6f64..1e6f64,
        ) {
            prop_assume!(a.abs() > 1e-8);

            let xm = - b/ (2.0 * a);
            let ym = a * xm * xm + b * xm + c;

            let (left, right, _) = InvertedQuadraticFunction::with(a, b, c).unwrap();
            if a.signum() == (y - ym).signum() {
                for arc in &[left, right] {
                    let x_value = arc.x_value_for(y).unwrap();
                    let y_value = a * x_value * x_value + b * x_value + c;
                    prop_assert!(
                        low_precision_eq(y_value, y),
                        "given y {} should be equal to calculated y {}",
                        y,
                        y_value
                    );
                }
            } else {
                for arc in &[left, right] {
                    prop_assert!(
                        arc.x_value_for(y).is_none(),
                        "Gave x_value when there should be none"
                    );
                }
            }
        }
    }
}
