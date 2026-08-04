// swift-tools-version: 6.3

import PackageDescription

let package = Package(
    name: "sebbu-fftw",
    platforms: [
        .macOS(.v26)
    ],
    products: [
        .library(
            name: "SebbuFFTW",
            targets: ["SebbuFFTW"]
        ),
    ],
    dependencies: [
        .package(url: "https://github.com/apple/swift-numerics", from: "1.1.1"),
        .package(url: "https://github.com/MarSe32m/sebbu-cfftw", from: "3.3.11"),
        .package(url: "https://github.com/MarSe32m/sebbu-fft", from: "0.2.0")
    ],
    targets: [
        .target(
            name: "SebbuFFTW",
            dependencies: [
                .product(name: "Numerics", package: "swift-numerics"),
                .product(name: "CFFTW", package: "sebbu-cfftw"),
                .product(name: "SebbuFFT", package: "sebbu-fft")
            ]
        ),
        .testTarget(
            name: "SebbuFFTWTests",
            dependencies: ["SebbuFFTW"]
        ),
    ],
    swiftLanguageModes: [.v6]
)
