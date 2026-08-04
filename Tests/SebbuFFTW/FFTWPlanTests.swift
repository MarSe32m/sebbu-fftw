// Copyright (C) 2026 Sebastian Toivonen
// SPDX-License-Identifier: GPL-3.0-or-later

import SebbuFFTW
import SebbuFFT
import ComplexModule
import Testing

private let accelerateSampleCounts = [8, 10, 12, 15, 16, 20, 24, 30, 32]

@Suite("FFTW plans")
struct FFTWPlanTests {
    @Test("FFTW forward plan matches the independent DFT", arguments: accelerateSampleCounts)
    func forwardMatchesNaiveDFT(sampleCount: Int) {
        let input = deterministicSignal(count: sampleCount, seed: 10)
        let plan = FFT.FFTWPlan(sampleSize: sampleCount)
        expectForwardPlanMatchesOracle(plan, input: input, spacing: 0.125)
    }

    @Test("FFTW inverse plan matches the independent DFT", arguments: accelerateSampleCounts)
    func inverseMatchesNaiveDFT(sampleCount: Int) {
        let input = deterministicSignal(count: sampleCount, seed: 11)
        let plan = FFT.iFFTWPlan(sampleSize: sampleCount)
        expectInversePlanMatchesOracle(plan, input: input, spacing: 0.375)
    }

    @Test("A FFTW plan can be reused")
    func planReuse() {
        let sampleCount = 16
        let forwardPlan = FFT.FFTWPlan(sampleSize: sampleCount)
        let inversePlan = FFT.iFFTWPlan(sampleSize: sampleCount)

        for seed in 12..<18 {
            let input = deterministicSignal(count: sampleCount, seed: seed)
            expectForwardPlanMatchesOracle(
                forwardPlan,
                input: input,
                spacing: 0.25
            )
            expectInversePlanMatchesOracle(
                inversePlan,
                input: input,
                spacing: 0.5
            )
        }
    }
}

@Suite("One-shot FFT API")
struct OneShotFFTWAPITests {
    @Test("FFT.fftw matches the documented forward transform", arguments: [1, 5, 8, 11, 16])
    func forwardTransform(sampleCount: Int) {
        let input = deterministicSignal(count: sampleCount, seed: 18)
        let spacing = 0.125
        let result = FFT.fftw(input, spacing: spacing)

        expectApproximatelyEqual(
            result.frequencies,
            expectedForwardAxis(count: sampleCount, spacing: spacing)
        )
        expectApproximatelyEqual(result.spectrum, naiveDFT(input))
    }

    @Test("FFT.ifftw matches the documented inverse transform", arguments: [1, 5, 8, 11, 16])
    func inverseTransform(sampleCount: Int) {
        let input = deterministicSignal(count: sampleCount, seed: 19)
        let spacing = 0.375
        let result = FFT.ifftw(input, spacing: spacing)

        expectApproximatelyEqual(
            result.t,
            expectedInverseAxis(count: sampleCount, spacing: spacing)
        )
        expectApproximatelyEqual(result.signal, naiveDFT(input, inverse: true))
    }

    @Test("Empty one-shot transforms return empty outputs")
    func emptyInput() {
        let forward = FFT.fftw([])
        let inverse = FFT.ifftw([])

        #expect(forward.frequencies.isEmpty)
        #expect(forward.spectrum.isEmpty)
        #expect(inverse.t.isEmpty)
        #expect(inverse.signal.isEmpty)
    }
}