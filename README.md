# DSP Guide — C Implementation

A reference implementation of the **common digital signal processing
algorithms**, organised by category, with every module annotated with
the problem it solves, its computational complexity, and the trade-offs
versus its alternatives.

## Project layout

```
dsp_guide/
├── CMakeLists.txt
├── include/
│   ├── dsp.h                   <- single master include
│   ├── common.h                <- complex type, math helpers
│   ├── transforms/             <- dft, fft, dct
│   ├── filtering/              <- fir, iir
│   ├── operations/             <- convolution, correlation
│   ├── spectral/               <- window functions
│   ├── sampling/               <- decimation, interpolation, resampling
│   ├── wavelet/                <- discrete wavelet transform
│   ├── coding/                 <- error detection, FEC, equalization
│   └── modulation/             <- QAM, OFDM, channel, coded OFDM
├── src/                        <- implementations mirror include/
│   └── main.c                  <- annotated demo runner
├── tests/
│   └── test_all.c              <- 79 unit tests
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
| **Error control coding** | CRC-32, parity, checksum; Hamming(7,4), Reed-Solomon, convolutional + Viterbi, LDPC; block & convolutional interleaving; LMS equalizer |
| **Modulation** | QPSK / 16-QAM / 64-QAM with soft demapping; multipath + AWGN channel; OFDM (IFFT/CP/FFT, per-subcarrier equalization); end-to-end coded OFDM |

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
- **Error detection** covers parity, the 16-bit Internet checksum, and
  table-driven CRC-32 (Ethernet polynomial).
- **Forward error correction** implements Hamming(7,4) (syndrome
  decoding), Reed-Solomon over GF(2^8) (Berlekamp-Massey, Chien search,
  Forney algorithm), a rate-1/2 convolutional code with a Viterbi
  decoder supporting both hard- and soft-decision metrics, and LDPC
  codes with two decoders. Turbo codes are documented but not coded.
- **LDPC** is built on a sparse parity-check matrix (a Tanner graph),
  with a regular-code generator and a fixed-matrix constructor. It
  provides three decoders: hard-decision bit-flipping, exact
  sum-product (belief propagation), and min-sum - the transcendental-
  free approximation used in real 5G/Wi-Fi 6 hardware. A BER-sweep
  helper measures any decoder over a simulated AWGN channel.
- **Interleaving** provides block (R×C matrix) and convolutional
  (staggered delay lines) interleavers. These add no redundancy; they
  scatter burst errors so the FEC stage sees only correctable counts.
  The demo shows a 6-symbol burst that destroys one Reed-Solomon
  codeword without interleaving but is fully survivable with it.
- **Channel equalization** is an LMS adaptive FIR filter, the DSP front-
  end stage that cleans channel distortion before the decoder runs.
- **Modulation** provides Gray-coded QAM (QPSK/16/64) with hard and soft
  demapping, a multipath + AWGN channel model, and OFDM built directly
  on the FFT (IFFT modulator, cyclic prefix, per-subcarrier equalizer).
  The coded-OFDM transceiver wires the FFT, QAM, channel, and LDPC
  modules into one end-to-end chain and measures raw vs coded BER.

See `docs/ALGORITHMS.md` for the full complexity and trade-off table.

## Test coverage

79 tests covering: DFT/FFT agreement, FFT and DCT round-trips, FFT
power-of-two rejection, DCT energy compaction, FIR linear phase and DC
gain, FIR/IIR high-frequency attenuation, IIR stability detection,
direct/FFT convolution agreement, the convolution identity,
correlation-based delay estimation, auto-correlation peak, window
endpoint/symmetry properties, resampling output lengths, the wavelet
round-trip, wavelet rejection of non-power-of-two lengths, parity and
checksum detection, the CRC-32 standard test vector and burst
detection, Hamming single-bit correction, Reed-Solomon burst
correction, hard- and soft-decision Viterbi decoding, block and convolutional
interleaver round-trips, interleaving-aided burst recovery, LDPC
construction and syndrome checks, LDPC bit-flipping,
sum-product and min-sum decoding, the BER-sweep noise monotonicity
and the soft-beats-hard decoder ranking, LMS equalizer convergence,
QAM round-trips and unit-energy normalization, the multipath/AWGN
channel, OFDM round-trips and per-subcarrier equalization, and the
end-to-end coded-OFDM chain.
