#![no_std]

use core::ops::Range;
use num_traits::{Float, PrimInt};

pub trait EwmaFilter<T> {
    /// Push the next sample in the series into the filter.
    /// Returns exponentially weighted moving average
    fn push_sample(&mut self, new_value: T) -> T;

    /// Returns cached exponentially weighted moving average
    fn ewma_average(&self) -> T;

    /// returns the local minima and maxima
    fn local_range(&self) -> Range<T>;
}

/// Implements exponential weighted moving average of time series samples,
/// including exponentially fading minimum and maximum
pub struct FloatSeriesEwmaFilter<T> {
    /// number of samples that have been pushed through the filter
    sample_count: usize,
    /// recent minimum value (not global minimum)
    local_min: T,
    /// recent maximum value (not global maximum)
    local_max: T,
    /// exponentially weighted moving average
    average: T,
    /// weighting factor-- bigger alpha causes faster fade of old values
    alpha: T,
}

impl<T> FloatSeriesEwmaFilter<T>
where
    T: Float + core::ops::AddAssign,
{
    pub fn new(alpha: T) -> Self {
        Self {
            sample_count: 0,
            alpha,
            local_min: T::zero(),
            local_max: T::zero(),
            average: T::zero(),
        }
    }

    pub fn default() -> Self {
        Self::new(T::from(0.01).unwrap())
    }
}

impl<T> EwmaFilter<T> for FloatSeriesEwmaFilter<T>
where
    T: Float + core::ops::AddAssign,
{
    /// Returns exponentially weighted moving average
    fn push_sample(&mut self, new_value: T) -> T {
        if self.sample_count == 0 {
            //seed the EMWA with the initial value
            self.local_min = new_value;
            self.local_max = new_value;
            self.average = new_value;
        } else {
            self.average += self.alpha * (new_value - self.average);

            // extrema fade toward average
            if new_value > self.local_max {
                self.local_max = new_value;
            } else if new_value > self.average {
                self.local_max += self.alpha * (new_value - self.local_max);
            }

            if new_value < self.local_min {
                self.local_min = new_value;
            } else if new_value < self.average {
                self.local_min += self.alpha * (new_value - self.local_min);
            }
        }
        self.sample_count += 1;

        self.average
    }

    fn ewma_average(&self) -> T {
        self.average
    }

    fn local_range(&self) -> Range<T> {
        self.local_min..self.local_max
    }
}

pub struct IntSeriesEwmaFilter<T> {
    /// sample count
    sample_count: usize,

    /// recent minimum value (not global minimum)
    local_min: T,
    /// recent maximum value (not global maximum)
    local_max: T,
    /// exponentially weighted moving average
    average: T,
    /// weighting factor-- bigger alpha causes faster fade of old values
    alpha_numerator: T,
    alpha_denominator: T,
}

impl<T> IntSeriesEwmaFilter<T>
where
    T: PrimInt + core::ops::AddAssign,
{
    pub fn new(alpha_numerator: T, alpha_denominator: T) -> Self {
        Self {
            sample_count: 0,
            alpha_numerator,
            alpha_denominator,
            local_min: T::zero(),
            local_max: T::zero(),
            average: T::zero(),
        }
    }

    pub fn default() -> Self {
        Self::new(T::one(), T::from(100).unwrap())
    }
}

impl<T> EwmaFilter<T> for IntSeriesEwmaFilter<T>
where
    T: PrimInt + core::ops::AddAssign,
{
    /// Returns exponentially weighted moving average
    fn push_sample(&mut self, new_value: T) -> T {
        if self.sample_count == 0 {
            //seed the EMWA with the initial value
            self.local_min = new_value;
            self.local_max = new_value;
            self.average = new_value;
        } else {
            self.average +=
                (self.alpha_numerator * (new_value - self.average)) / self.alpha_denominator;

            // extrema fade toward average
            if new_value > self.local_max {
                self.local_max = new_value;
            } else if new_value > self.average {
                self.local_max +=
                    (self.alpha_numerator * (new_value - self.local_max)) / self.alpha_denominator;
            }

            if new_value < self.local_min {
                self.local_min = new_value;
            } else if new_value < self.average {
                self.local_min +=
                    (self.alpha_numerator * (new_value - self.local_min)) / self.alpha_denominator;
            }
        }
        self.sample_count += 1;

        self.average
    }

    fn ewma_average(&self) -> T {
        self.average
    }

    fn local_range(&self) -> Range<T> {
        self.local_min..self.local_max
    }
}

#[cfg(test)]
mod tests {
    use crate::{EwmaFilter, FloatSeriesEwmaFilter, IntSeriesEwmaFilter};
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn float_basic() {
        let mut tracko: FloatSeriesEwmaFilter<f32> = FloatSeriesEwmaFilter::default();
        for i in 0..1000 {
            tracko.push_sample(i as f32);
        }
        assert_approx_eq!(tracko.ewma_average(), 900.0, 1f32);

        let mut tracko: FloatSeriesEwmaFilter<f32> = FloatSeriesEwmaFilter::new(0.01);
        for i in 0..1000 {
            tracko.push_sample(i as f32);
        }
        assert_approx_eq!(tracko.ewma_average(), 900.0, 1f32);
        let range = tracko.local_range();
        assert_eq!(range.end, 999.0);
        assert_eq!(range.start, 0.0);
    }

    #[test]
    fn integer_basic() {
        let mut tracko: IntSeriesEwmaFilter<u32> = IntSeriesEwmaFilter::default();
        for i in 0..1000 {
            tracko.push_sample(i);
        }
        assert_eq!(tracko.ewma_average(), 900);

        let mut tracko: IntSeriesEwmaFilter<u32> = IntSeriesEwmaFilter::new(1, 100);
        for i in 0..1000 {
            tracko.push_sample(i);
        }
        assert_eq!(tracko.ewma_average(), 900);

        let range = tracko.local_range();
        assert_eq!(range.end, 999);
        assert_eq!(range.start, 0);
    }
}
