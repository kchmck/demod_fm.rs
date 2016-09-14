//! This crates provides a simple [FM
//! signal](https://en.wikipedia.org/wiki/Frequency_modulation) demodulator for use with
//! software radio. It demodulates using phase difference approximation as described
//! below.
//!
//! ## Theory
//!
//! Consider the classical equation \[1] for an FM signal:
//!
//! > s(t) = a(t) cos(ω<sub>c</sub>t + φ(t))
//!
//! with
//!
//! > φ(t) = ω<sub>∆</sub>∫x(τ)dτ
//!
//! where the integral is evaluated from 0 to t and x(t) is the modulating signal to be
//! recovered.
//!
//! Differentiating this gives
//!
//! > dφ(t) / dt = ω<sub>∆</sub>x(t)
//!
//! so
//!
//! > x(t) = ω<sub>∆</sub><sup>-1</sup> dφ(t) / dt
//!
//! Differentiation in continuous time is approximated by [finite backward
//! difference](https://en.wikipedia.org/wiki/Finite_difference) in discrete time, so
//!
//! > x(t) ≈ ω<sub>∆</sub><sup>-1</sup> (φ[t] - φ[t-1]) / T
//!
//! Assuming a "normalized" period of T = 1, this becomes
//!
//! > x(t) ≈ w<sub>∆</sub><sup>-1</sup> (φ[t] - φ[t-1])
//!
//! This requires the change in phase between the current and previous sampling instants,
//! which can be computed from the corresponding I/Q samples. Given an FM signal s(t),
//! the received I/Q sequence will have components
//!
//! > i(t) = a(t) cos φ(t)
//!
//! > q(t) = a(t) sin φ(t)
//!
//! with each sample represented as
//!
//! > p(t) = i(t) + *j* q(t)
//!
//! Evaluating the [complex
//! argument](http://mathworld.wolfram.com/ComplexArgument.html) of this gives
//!
//! > arg(p(t)) = arctan[q(t) / i(t)] = arctan tan φ(t) = φ(t)
//!
//! so
//!
//! > arg(p(t)) - arg(p(t-1)) = φ(t) - φ(t-1)
//!
//! Applying the complex identities [arg(uv) ≡ arg(u) + arg(v) (mod (-π,
//! π\])](https://en.wikipedia.org/wiki/Argument_(complex_analysis)#Identities) and
//! [arg(u<sup>*</sup>) = -arg(u)](http://mathworld.wolfram.com/ComplexArgument.html),
//!
//! > arg(p(t)p(t-1)<sup>*</sup>) = arg(p(t)) - arg(p(t-1)) = φ(t) - φ(t-1)
//!
//! Combining all these results leads to the equation calculated at each sample:
//!
//! > x[t] = ω<sub>∆</sub><sup>-1</sup> arg(p[t]p[t - 1]<sup>*</sup>)
//!
//! using angular frequency deviation ω<sub>∆</sub> = 2π f<sub>∆</sub> and the current and
//! previous complex samples.
//!
//! ## References
//!
//! 1. "FM demodulation using a digital radio and digital signal processing", J.M. Shima,
//! 1995.

extern crate num;

use std::f32::consts::PI;

use num::complex::Complex32;

/// Demodulates an FM signal using a phase difference approximation.
pub struct FmDemod {
    /// Reciprocol of angular frequency deviation, ω<sub>∆</sub><sup>-1</sup>
    gain: f32,
    /// Previous sample, p[t-1].
    prev: Complex32,
}

impl FmDemod {
    /// Create a new `FmDemod` with the given frequency deviation f<sub>∆</sub> (Hz) and
    /// sample rate f<sub>s</sub> (Hz).
    ///
    /// The deviation must satisfy the Nyquist limit, f<sub>∆</sub> ≤ f<sub>s</sub> / 2.
    pub fn new(deviation: u32, sample_rate: u32) -> FmDemod {
        assert!(deviation <= sample_rate / 2);

        FmDemod {
            gain: (2.0 * PI * deviation as f32 / sample_rate as f32).recip(),
            prev: Complex32::new(0.0, 0.0),
        }
    }

    /// Feed in an FM sample, producing the next sample in the demodulated signal.
    pub fn feed(&mut self, sample: Complex32) -> f32 {
        // Compute x[t].
        let next = (sample * self.prev.conj()).arg() * self.gain;
        self.prev = sample;

        next
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_fm_demod() {
        // Tests a binary NRZ payload signal.

        // Samples per second.
        let samprate = 48000;
        // Frequency deviation.
        let dev = 4000;
        // Angular frequency deviation.
        let angdev = 2.0 * PI * dev as f32 / samprate as f32;

        // Data symbols to encode.
        let data = [-1, 1, 1, -1, 1, -1];

        // Generate "received" I/Q.
        let mut accum = 0.0f32;
        let mut sig = vec![];

        for &sym in data.iter() {
            // Use 2 samples per symbol.
            for _ in 0..2 {
                // Compute Riemann sum integral approximation.
                accum += angdev * sym as f32;
                sig.push(Complex32::new(accum.cos(), accum.sin()));
            }
        }

        let mut d = FmDemod::new(dev, samprate);
        let mut sig = sig.into_iter();

        // Load first sample into demodulator.
        d.feed(sig.next().unwrap());

        assert_eq!(d.feed(sig.next().unwrap()), -1.0);
        assert_eq!(d.feed(sig.next().unwrap()), 1.0);
        assert_eq!(d.feed(sig.next().unwrap()), 1.0);
        assert_eq!(d.feed(sig.next().unwrap()), 1.0);
        assert_eq!(d.feed(sig.next().unwrap()), 1.0);
        assert_eq!(d.feed(sig.next().unwrap()), -1.0);
        assert_eq!(d.feed(sig.next().unwrap()), -1.0);
        assert_eq!(d.feed(sig.next().unwrap()), 1.0);
        assert_eq!(d.feed(sig.next().unwrap()), 1.0);
        assert_eq!(d.feed(sig.next().unwrap()), -1.0);
        assert_eq!(d.feed(sig.next().unwrap()), -1.0);
    }
}
