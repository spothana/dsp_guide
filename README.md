# DSP Study Guide — C Implementation

A reference implementation of the **common digital signal processing
algorithms**, organised by category, with every module annotated with
the problem it solves, its computational complexity, and the trade-offs
versus its alternatives.


## Project layout

```
dsp_study_guide/
├── CMakeLists.txt
├── include/
│   ├── dsp.h                   <- single master include
│   ├── common.h                <- complex type, math helpers
│   ├── transforms/             <- dft, fft, dct
│   ├── filtering/              <- fir, iir
│   ├── operations/             <- convolution, correlation
│   ├── spectral/               <- window functions
│   ├── sampling/               <- decimation, interpolation, resampling
│   └── wavelet/                <- discrete wavelet transform
├── src/                        <- implementations mirror include/
│   └── main.c                  <- annotated demo runner
├── tests/
│   └── test_all.c              <- 26 unit tests
└── docs/
    └── ALGORITHMS.md           <- complexity & trade-off cheat-sheet
```

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)

# Run the annotated demo
./dsp_demo

# Run the unit tests
ctest --output-on-failure
# or directly:
./test_all
```

Requires only a C11 compiler, CMake 3.16+, and libm. No external
dependencies.

## Algorithms covered

| Category | Algorithms |
|---|---|
| **Frequency-domain transforms** | DFT, FFT (radix-2 Cooley-Tukey), DCT-II/III |
| **Digital filtering** | FIR (windowed-sinc design), IIR (biquad, RBJ design) |
| **Signal operations** | Convolution (direct + FFT), cross/auto-correlation |
| **Spectral analysis** | Rectangular, Hamming, Hanning, Blackman windows |
| **Sample rate conversion** | Decimation, interpolation, rational resampling |
| **Multi-resolution analysis** | Discrete wavelet transform (Haar, Mallat pyramid) |

## Design notes

- **DFT** is the literal O(N²) definition — a correctness reference for
  the FFT, not for production use.
- **FFT** is an in-place iterative radix-2 transform. It requires
  power-of-two lengths; callers with other lengths zero-pad with
  `dsp_next_pow2`.
- **FIR** filters are designed by the windowed-sinc method with symmetric
  coefficients, giving exactly linear phase and guaranteed stability.
- **IIR** filters are second-order sections (biquads) with a Schur-Cohn
  stability check, since IIR stability is not automatic.
- **Convolution** is provided both directly (O(N·M)) and via the FFT
  (O(N log N)), demonstrating the time-multiplication / frequency-
  multiplication duality.
- **Wavelet** transform uses the Haar basis via Mallat's pyramid
  algorithm; the round-trip is numerically exact.

See `docs/ALGORITHMS.md` for the full complexity and trade-off table.

## Test coverage

26 tests covering: DFT/FFT agreement, FFT and DCT round-trips, FFT
power-of-two rejection, DCT energy compaction, FIR linear phase and DC
gain, FIR/IIR high-frequency attenuation, IIR stability detection,
direct/FFT convolution agreement, the convolution identity,
correlation-based delay estimation, auto-correlation peak, window
endpoint/symmetry properties, resampling output lengths, the wavelet
round-trip, and wavelet rejection of non-power-of-two lengths.
