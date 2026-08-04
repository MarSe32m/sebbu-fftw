// Copyright (C) 2026 Sebastian Toivonen
// SPDX-License-Identifier: GPL-3.0-or-later

@_exported import SebbuFFT
import CFFTW
import RealModule
import ComplexModule
import Synchronization

public extension FFT {
    @usableFromInline
    nonisolated internal static let _fftwPlanMutex = Mutex<Void>(())

    /// A reusable FFTW-backed plan for a one-dimensional forward complex
    /// discrete Fourier transform.
    ///
    /// An instance transforms arrays containing exactly ``sampleSize`` complex
    /// samples. It uses FFTW's forward sign convention:
    ///
    /// `Y[k] = Σ x[n] exp(-2π i n k / N)`
    ///
    /// The forward transform is unnormalized. ``iFFTWPlan`` applies the
    /// corresponding `1/N` normalization, so applying the inverse transform to
    /// the resulting spectrum recovers the original samples up to numerical
    /// roundoff.
    ///
    /// The output spectrum uses unshifted FFT bin order: the zero and positive
    /// frequency bins come first, followed by the negative frequency bins. For
    /// even transform sizes, the Nyquist bin is represented as a negative
    /// frequency.
    ///
    /// ## Planning
    ///
    /// By default, FFTW selects a plan using its inexpensive `FFTW_ESTIMATE`
    /// heuristic. Pass `measure: true` to use `FFTW_MEASURE`, which takes longer
    /// during initialization but may produce a faster execution plan.
    ///
    /// ## Thread safety
    ///
    /// A plan may be shared between threads. Concurrent calls to
    /// ``execute(_:spacing:)`` are supported because each call uses an independent,
    /// FFTW-aligned execution buffer. FFTW plan creation and destruction are
    /// serialized internally.
    final class FFTWPlan: Plan, @unchecked Sendable {
        /// The number of complex samples transformed by this plan.
        ///
        /// Every array passed to ``execute(_:spacing:)`` must contain exactly this
        /// many elements.
        public let sampleSize: Int

        private let planningBuffer: UnsafeMutablePointer<fftw_complex>?
        private let plan: OpaquePointer?

        /// Creates a reusable forward-transform plan.
        ///
        /// - Parameters:
        ///   - sampleSize: The number of complex samples in each transform.
        ///   - measure: Whether FFTW should measure candidate algorithms when
        ///     constructing the plan. This can make initialization substantially
        ///     slower, but may improve execution performance.
        /// - Precondition: `sampleSize` is non-negative and representable by `CInt`.
        public init(sampleSize: Int, measure: Bool = false) {
            guard sampleSize >= 0, let fftwSize = CInt(exactly: sampleSize) else {
                preconditionFailure(
                    "sampleSize must be positive and representable by CInt"
                )
            }
            if sampleSize == 0 {
                self.sampleSize = 0
                self.planningBuffer = nil
                self.plan = nil
                return
            }

            self.sampleSize = sampleSize

            let flags = measure ? FFTW_MEASURE : FFTW_ESTIMATE

            let (planningBuffer, plan) =
                FFT._fftwPlanMutex.withLock { _ in
                    guard let buffer = fftw_alloc_complex(sampleSize) else {
                        fatalError("Could not allocate the FFTW planning buffer")
                    }

                    guard let plan = fftw_plan_dft_1d(
                        fftwSize,
                        buffer,
                        buffer,
                        FFTW_FORWARD,
                        flags
                    ) else {
                        fftw_free(buffer)
                        fatalError("Could not create the FFTW plan")
                    }

                    return (buffer, plan)
                }

            self.planningBuffer = planningBuffer
            self.plan = plan
        }

        /// Computes an unnormalized forward discrete Fourier transform.
        ///
        /// For input-sample spacing `Δx`, the returned angular-frequency spacing is
        ///
        /// `Δω = 2π / (sampleSize * Δx)`.
        ///
        /// The returned spectrum is in unshifted FFT bin order. Use
        /// ``FFT/fftShift(frequencies:spectrum:)`` when a frequency axis centered
        /// around zero is required.
        ///
        /// This method does not modify `input` and may be called concurrently on
        /// the same plan.
        ///
        /// - Parameters:
        ///   - input: Equally spaced complex samples. The array must contain
        ///     exactly ``sampleSize`` elements.
        ///   - spacing: The positive, finite distance `Δx` between adjacent input
        ///     samples.
        /// - Returns: A tuple containing:
        ///   - `x`: Angular frequencies in radians per input-axis unit, in
        ///     unshifted FFT bin order.
        ///   - `y`: The unnormalized complex spectrum at those frequencies.
        /// - Precondition: `input.count == sampleSize`.
        /// - Precondition: `spacing` is finite and greater than zero.
        public func execute(
            _ input: [Complex<Double>],
            spacing: Double
        ) -> (x: [Double], y: [Complex<Double>]) {
            precondition(
                input.count == sampleSize,
                "The signal must have sampleSize elements"
            )
            precondition(
                spacing.isFinite && spacing > 0,
                "Spacing must be finite and positive"
            )
            precondition(
                MemoryLayout<Complex<Double>>.size ==
                    MemoryLayout<fftw_complex>.size &&
                MemoryLayout<Complex<Double>>.stride ==
                    MemoryLayout<fftw_complex>.stride,
                "Complex<Double> and fftw_complex have incompatible layouts"
            )
            if input.isEmpty { return ([], []) }

            let count = input.count
            let byteCount =
                count * MemoryLayout<fftw_complex>.stride

            guard let buffer = fftw_alloc_complex(count) else {
                fatalError("Could not allocate the FFTW execution buffer")
            }
            defer { fftw_free(buffer) }

            let bufferBytes = UnsafeMutableRawBufferPointer(
                start: UnsafeMutableRawPointer(buffer),
                count: byteCount
            )

            input.withUnsafeBytes {
                bufferBytes.copyMemory(from: $0)
            }

            // Deliberately not locked. FFTW permits concurrent execution
            // of the same plan when each invocation uses separate arrays.
            fftw_execute_dft(plan, buffer, buffer)

            var output = input
            output.withUnsafeMutableBytes {
                $0.copyMemory(
                    from: UnsafeRawBufferPointer(bufferBytes)
                )
            }

            let scale =
                2.0 * Double.pi / (Double(count) * spacing)
            let firstNegativeIndex = (count + 1) / 2

            let angularFrequencies = (0..<count).map { index in
                let bin = index < firstNegativeIndex
                    ? index
                    : index - count

                return Double(bin) * scale
            }

            return (angularFrequencies, output)
        }

        deinit {
            FFT._fftwPlanMutex.withLock { _ in
                if let plan {
                    fftw_destroy_plan(plan)
                }
                if let planningBuffer {
                    fftw_free(planningBuffer)
                }
            }
        }
    }


    /// A reusable FFTW-backed plan for a one-dimensional inverse complex
    /// discrete Fourier transform.
    ///
    /// An instance transforms arrays containing exactly ``sampleSize`` complex
    /// frequency-domain samples. It uses FFTW's backward sign convention and
    /// applies `1/N` normalization:
    ///
    /// `x[n] = (1/N) Σ Y[k] exp(+2π i n k / N)`
    ///
    /// Consequently, applying this plan to the unmodified output of a compatible
    /// ``FFTWPlan`` recovers the original signal up to numerical roundoff.
    ///
    /// Input spectra must use unshifted FFT bin order: the zero and positive
    /// frequency bins first, followed by the negative frequency bins. A shifted
    /// spectrum must be returned to unshifted order before execution.
    ///
    /// ## Planning
    ///
    /// By default, FFTW selects a plan using its inexpensive `FFTW_ESTIMATE`
    /// heuristic. Pass `measure: true` to use `FFTW_MEASURE`, which takes longer
    /// during initialization but may produce a faster execution plan.
    ///
    /// ## Thread safety
    ///
    /// A plan may be shared between threads. Concurrent calls to
    /// ``execute(_:spacing:)`` are supported because each call uses an independent,
    /// FFTW-aligned execution buffer. FFTW plan creation and destruction are
    /// serialized internally.
    final class iFFTWPlan: Plan, @unchecked Sendable {
        /// The number of complex frequency-domain samples transformed by this plan.
        ///
        /// Every array passed to ``execute(_:spacing:)`` must contain exactly this
        /// many elements.
        public let sampleSize: Int

        private let planningBuffer: UnsafeMutablePointer<fftw_complex>?
        private let plan: OpaquePointer?

        /// Creates a reusable inverse-transform plan.
        ///
        /// - Parameters:
        ///   - sampleSize: The number of complex frequency bins in each transform.
        ///   - measure: Whether FFTW should measure candidate algorithms when
        ///     constructing the plan. This can make initialization substantially
        ///     slower, but may improve execution performance.
        /// - Precondition: `sampleSize` is non-negative and representable by `CInt`.
        public init(sampleSize: Int, measure: Bool = false) {
            guard sampleSize >= 0, let fftwSize = CInt(exactly: sampleSize) else {
                preconditionFailure(
                    "sampleSize must be positive and representable by CInt"
                )
            }
            if sampleSize == 0 {
                self.sampleSize = 0
                self.planningBuffer = nil
                self.plan = nil
                return
            }

            self.sampleSize = sampleSize

            let flags = measure ? FFTW_MEASURE : FFTW_ESTIMATE

            let (planningBuffer, plan) =
                FFT._fftwPlanMutex.withLock { _ in
                    guard let buffer = fftw_alloc_complex(sampleSize) else {
                        fatalError("Could not allocate the FFTW planning buffer")
                    }

                    guard let plan = fftw_plan_dft_1d(
                        fftwSize,
                        buffer,
                        buffer,
                        FFTW_BACKWARD,
                        flags
                    ) else {
                        fftw_free(buffer)
                        fatalError("Could not create the FFTW plan")
                    }

                    return (buffer, plan)
                }

            self.planningBuffer = planningBuffer
            self.plan = plan
        }

        /// Computes a normalized inverse discrete Fourier transform.
        ///
        /// `spacing` is the angular-frequency separation `Δω` between adjacent FFT
        /// bins—not the original sample spacing. The returned sample spacing is
        ///
        /// `Δt = 2π / (sampleSize * Δω)`.
        ///
        /// For example, if a forward transform used sample spacing `Δt`, pass
        /// `2π / (sampleSize * Δt)` as the inverse transform's `spacing`.
        ///
        /// The returned axis begins at zero and is not shifted or centered. This
        /// method does not modify `input` and may be called concurrently on the
        /// same plan.
        ///
        /// - Parameters:
        ///   - input: A complex spectrum in unshifted FFT bin order. The array must
        ///     contain exactly ``sampleSize`` elements.
        ///   - spacing: The positive, finite angular-frequency spacing `Δω`.
        /// - Returns: A tuple containing:
        ///   - `x`: Sample coordinates `0, Δt, 2Δt, ...`.
        ///   - `y`: The normalized inverse-transformed complex signal.
        /// - Precondition: `input.count == sampleSize`.
        /// - Precondition: `spacing` is finite and greater than zero.
        public func execute(
            _ input: [Complex<Double>],
            spacing: Double = 1.0
        ) -> (x: [Double], y: [Complex<Double>]) {
            precondition(
                input.count == sampleSize,
                "The signal must have sampleSize elements"
            )
            precondition(
                spacing.isFinite && spacing > 0,
                "Spacing must be finite and positive"
            )
            precondition(
                MemoryLayout<Complex<Double>>.size ==
                    MemoryLayout<fftw_complex>.size &&
                MemoryLayout<Complex<Double>>.stride ==
                    MemoryLayout<fftw_complex>.stride,
                "Complex<Double> and fftw_complex have incompatible layouts"
            )
            if input.isEmpty { return ([], []) }

            let count = input.count
            let byteCount =
                count * MemoryLayout<fftw_complex>.stride

            guard let buffer = fftw_alloc_complex(count) else {
                fatalError("Could not allocate the FFTW execution buffer")
            }
            defer { fftw_free(buffer) }

            let bufferBytes = UnsafeMutableRawBufferPointer(
                start: UnsafeMutableRawPointer(buffer),
                count: byteCount
            )

            input.withUnsafeBytes {
                bufferBytes.copyMemory(from: $0)
            }

            // New-array execution allows this plan to be used concurrently.
            fftw_execute_dft(plan, buffer, buffer)

            var signal = input
            signal.withUnsafeMutableBytes {
                $0.copyMemory(
                    from: UnsafeRawBufferPointer(bufferBytes)
                )
            }

            // FFTW's backward transform is unnormalized.
            let normalization = 1.0 / Double(count)

            for index in signal.indices {
                signal[index].real *= normalization
                signal[index].imaginary *= normalization
            }

            let timeStep =
                2.0 * Double.pi / (Double(count) * spacing)

            let times = (0..<count).map {
                Double($0) * timeStep
            }

            return (times, signal)
        }

        deinit {
            FFT._fftwPlanMutex.withLock { _ in
                if let plan {
                    fftw_destroy_plan(plan)
                }
                if let planningBuffer {
                    fftw_free(planningBuffer)
                }
            }
        }
    }
}