// Copyright (C) 2026 Sebastian Toivonen
// SPDX-License-Identifier: GPL-3.0-or-later

import SebbuFFT
import RealModule
import ComplexModule

public extension FFT {
    /// Computes the discrete Fourier transform of equally spaced complex samples using FFTW.
    ///
    /// For `N` samples, the transform follows this convention:
    ///
    /// ```text
    /// X[k] = sum(x[n] * exp(-i * 2π * n * k / N), n = 0 ..< N)
    /// ```
    ///
    /// The result is unnormalized. Frequencies are returned in unshifted order, with
    /// zero and positive frequencies followed by negative frequencies. Use
    /// ``fftShift(frequencies:spectrum:)`` when a centered frequency axis is more useful.
    ///
    /// - Parameters:
    ///   - input: The equally spaced complex input samples.
    ///   - spacing: The distance between adjacent samples. Use a finite, positive
    ///     value. A spacing measured in seconds produces angular frequencies in
    ///     radians per second.
    /// - Returns: The angular-frequency bins and corresponding complex spectrum. Each
    ///   output array has `input.count` elements; empty input produces two empty arrays.
    @inlinable
    nonisolated static func fftw(_ input: [Complex<Double>], spacing: Double = 1.0, measure: Bool = false) -> (frequencies: [Double], spectrum: [Complex<Double>]) {
        if input.isEmpty { return ([], []) }
        let plan = FFTWPlan(sampleSize: input.count, measure: measure)
        let (frequencies, spectrum) = plan.execute(input, spacing: spacing)
        return (frequencies, spectrum)
    }

    /// Computes the inverse discrete Fourier transform of complex frequency bins using FFTW.
    ///
    /// For `N` bins, the transform follows this convention:
    ///
    /// ```text
    /// x[n] = (1 / N) * sum(X[k] * exp(i * 2π * n * k / N), k = 0 ..< N)
    /// ```
    ///
    /// Supply bins in the unshifted order produced by ``fft(_:spacing:)``. If the bins
    /// are centered, call ``ifftShift(_:_:)`` first.
    ///
    /// - Parameters:
    ///   - input: The complex frequency bins in unshifted order.
    ///   - spacing: The angular-frequency distance between adjacent bins. Use a finite,
    ///     positive value. The returned sample positions are separated by
    ///     `2π / (N * spacing)`.
    /// - Returns: The reconstructed sample positions and normalized complex signal.
    ///   Each output array has `input.count` elements; empty input produces two empty arrays.
    @inlinable
    nonisolated static func ifftw(_ input: [Complex<Double>], spacing: Double = 1.0, measure: Bool = false) -> (t: [Double], signal: [Complex<Double>]) {
        let plan = iFFTWPlan(sampleSize: input.count, measure: measure)
        let (t, signal) = plan.execute(input, spacing: spacing)
        return (t, signal)
    }
}
