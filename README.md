# sebbu-fftw

`sebbu-fftw` adds an [FFTW](https://www.fftw.org/)-backed implementation of
one-dimensional complex Fourier transforms to
[`sebbu-fft`](https://github.com/MarSe32m/sebbu-fft). It provides convenient
one-shot transforms and reusable FFTW plans while preserving the transform,
axis and normalization conventions of the backend-neutral `SebbuFFT` API.

The package is intentionally separate from `sebbu-fft`. Applications that do
not need FFTW can continue to use the portable, MIT-licensed core package
without adding a native library or GPL-licensed dependency.

## Features

- Forward and inverse transforms for `[Complex<Double>]`
- Unnormalized forward transforms and `1 / N`-normalized inverse transforms
- Angular-frequency axes compatible with `sebbu-fft`
- Convenient one-shot `FFT.fftw` and `FFT.ifftw` functions
- Reusable `FFT.FFTWPlan` and `FFT.iFFTWPlan` plan types
- A choice between FFTW's `FFTW_ESTIMATE` and `FFTW_MEASURE` planners
- Concurrent execution of a shared plan using independent aligned buffers
- No dependency on Foundation in the wrapper API
- Re-export of `SebbuFFT`, including its plan protocol and shift operations

## Package structure

| Package | Role | License |
| --- | --- | --- |
| [`sebbu-fft`](https://github.com/MarSe32m/sebbu-fft) | Backend-neutral API and portable implementations | MIT |
| [`sebbu-cfftw`](https://github.com/MarSe32m/sebbu-cfftw) | C interface and packaged FFTW library | FFTW is GPL-2.0-or-later |
| `sebbu-fftw` | Swift integration and reusable FFTW-backed plans | GPL-3.0-or-later |

Supported platforms and architectures are determined by the artifacts provided
by `sebbu-cfftw`.

## Installation

Add the package in `Package.swift`:

```swift
dependencies: [
    .package(
        url: "https://github.com/MarSe32m/sebbu-fftw.git",
        from: "0.1.0"
    )
]
```

Then add the library product to your target:

```swift
.target(
    name: "MyTarget",
    dependencies: [
        .product(name: "SebbuFFTW", package: "sebbu-fftw")
    ]
)
```

## Quick start

```swift
import ComplexModule
import SebbuFFTW

let samples: [Complex<Double>] = [
    Complex(1, 0),
    Complex(0, 0),
    Complex(-1, 0),
    Complex(0, 0),
]

let sampleSpacing = 0.001 // seconds

let (angularFrequencies, spectrum) = FFT.fftw(
    samples,
    spacing: sampleSpacing
)
```

`angularFrequencies` is measured in radians per second because the input
spacing is measured in seconds. The bins are returned in standard unshifted
FFT order: zero and positive frequencies first, followed by negative
frequencies.

To center the zero-frequency bin for display or analysis, use the shift
operations re-exported from `SebbuFFT`:

```swift
var centeredFrequencies = angularFrequencies
var centeredSpectrum = spectrum

FFT.fftShift(
    frequencies: &centeredFrequencies,
    spectrum: &centeredSpectrum
)
```

For a round trip, pass the angular-frequency bin spacingâ€”not the original
sample spacingâ€”to the inverse transform:

```swift
let angularFrequencySpacing =
    2 * Double.pi / (Double(samples.count) * sampleSpacing)

let (samplePositions, reconstructed) = FFT.ifftw(
    spectrum,
    spacing: angularFrequencySpacing
)
```

Apart from floating-point error, `reconstructed` contains the original samples,
and `samplePositions` has the original sample spacing.

## Transform conventions

For `N` samples, the forward transform is

```text
X[k] = sum(x[n] * exp(-i * 2Ï€ * n * k / N), n = 0 ..< N)
```

and the inverse transform is

```text
x[n] = (1 / N) * sum(X[k] * exp(i * 2Ï€ * n * k / N), k = 0 ..< N)
```

| Operation | `spacing` | Returned axis | Scaling |
| --- | --- | --- | --- |
| `FFT.fftw` | Sample spacing `Î”t` | `Ï‰â‚– = 2Ï€k / (NÎ”t)` | None |
| `FFT.ifftw` | Angular-frequency spacing `Î”Ï‰` | `tâ‚™ = 2Ï€n / (NÎ”Ï‰)` | `1 / N` |

Both operations use unshifted bin order. If a spectrum has been centered with
`FFT.fftShift`, restore it with `FFT.ifftShift` before applying the inverse
transform.

## Reusing plans

The one-shot functions create and destroy an FFTW plan on every call. For
repeated transforms of the same size, construct a plan once and reuse it:

```swift
let plan = FFT.FFTWPlan(
    sampleSize: samples.count,
    measure: true
)

for block in blocks {
    precondition(block.count == samples.count)

    let (frequencies, spectrum) = plan.execute(
        block,
        spacing: sampleSpacing
    )

    // Consume this block's spectrum.
}
```

Use `FFT.iFFTWPlan(sampleSize:measure:)` for repeated inverse transforms. Every
array passed to a plan must have the same number of elements as the
`sampleSize` supplied during initialization.

## Planning modes

By default, plans use `FFTW_ESTIMATE`. This makes plan creation inexpensive and
is appropriate for one-shot transforms or applications where startup time
matters.

Passing `measure: true` selects `FFTW_MEASURE`. FFTW then benchmarks candidate
algorithms when the plan is created. Planning may take substantially longer,
but the resulting plan can execute faster. This option is most useful when a
plan will be reused many times; repeatedly using `measure: true` with the
one-shot functions pays the planning cost on every call.

## Thread safety

A single `FFT.FFTWPlan` or `FFT.iFFTWPlan` instance may be shared between
threads. Each execution uses an independent FFTW-aligned buffer, so concurrent
calls do not mutate shared input or output storage. FFTW plan creation and
destruction are serialized internally.

## Selecting the FFTW backend

Importing `SebbuFFTW` does not change the backend selected by
`FFT.defaultFFTPlan` or the behavior of `FFT.fft`. Choose FFTW explicitly by
calling `FFT.fftw`, `FFT.ifftw`, `FFT.FFTWPlan` or `FFT.iFFTWPlan`.

This keeps backend selection visible at the call site and allows an application
to use the portable or Accelerate implementations from `sebbu-fft` alongside
FFTW plans.

## Numerical and API notes

- The current API supports one-dimensional transforms of `Complex<Double>`.
- `spacing` must be finite and greater than zero.
- One-shot transforms return empty arrays for empty input.
- A reusable plan requires every input array to match its configured size.
- Forward transforms are unnormalized; inverse transforms apply `1 / N`.
- Transform results and axes are newly allocated for each execution.
- Floating-point round trips are approximate.

## Acknowledgments

This package wraps FFTW, developed by Matteo Frigo and Steven G. Johnson. If
FFTW contributes to published work, its authors request citation of:

> Matteo Frigo and Steven G. Johnson, â€œThe Design and Implementation of FFTW3,â€
> *Proceedings of the IEEE* **93**(2), 216-231 (2005).
> [doi:10.1109/JPROC.2004.840301](https://doi.org/10.1109/JPROC.2004.840301)

## License

Copyright © 2026 Sebastian Toivonen.

The Swift source in `sebbu-fftw` is available under the GNU General Public
License, version 3 or any later version (`GPL-3.0-or-later`). See `LICENSE` for
the complete terms.

FFTW is distributed separately under the GNU General Public License, version 2
or any later version, and remains under the copyright of its respective
authors. Software distributed after linking with FFTW must comply with the
applicable FFTW license terms. Non-free FFTW licenses are also available from
MIT for uses that require different terms.

For FFTW's authoritative licensing information, see
[License and Copyright](https://www.fftw.org/doc/License-and-Copyright.html).