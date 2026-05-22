# DSP Algorithms — Complexity & Trade-off Cheat-Sheet

Every algorithm in this guide, with its cost and the decision rule for
when to reach for it. `N` is signal length; `M` is filter/kernel length.

## Frequency-domain transforms

| Algorithm | Complexity | Output | Use when |
|---|---|---|---|
| DFT  | O(N²)      | N complex bins | Reference / teaching only — too slow for real N |
| FFT  | O(N log N) | N complex bins | Any practical spectral analysis; identical result to DFT |
| DCT  | O(N²) here | N real coeffs  | Compression — energy compacts into few coefficients |

**DFT vs FFT.** Same mathematical result. The FFT exploits symmetry
(`W_N^{k+N/2} = -W_N^k`) and periodicity to divide and conquer. For
N = 1024 that is roughly a 100× speedup. Never compute the DFT directly
in production.

**DFT vs DCT.** Both move a signal to the frequency domain. The DCT uses
only cosines, stays real for real input, and handles block edges better
(implicit symmetric extension instead of periodic wraparound). Its strong
energy compaction is why JPEG and MPEG use it. If a question is about
image or video compression, the answer is DCT.

## Digital filtering

| Filter | Cost per sample | Stable? | Phase | Use when |
|---|---|---|---|---|
| FIR | O(M), large M | Always | Linear (if symmetric) | Phase/waveform shape must be preserved |
| IIR | O(order), small | Not guaranteed | Nonlinear | Efficiency matters; sharp cutoff, few coefficients |

**FIR vs IIR.** FIR guarantees stability and linear phase but needs many
taps for a sharp cutoff. IIR achieves the same magnitude response with a
handful of coefficients but risks instability (poles must stay inside the
unit circle) and distorts phase. Emphasis on phase → FIR. Emphasis on
efficiency or analog-filter equivalence → IIR.

## Signal operations

| Operation | Complexity | Use when |
|---|---|---|
| Direct convolution | O(N·M) | Short kernels |
| FFT convolution | O(N log N) | Long signals/kernels |
| Cross-correlation | O(N·M) | Delay/similarity between two signals |
| Auto-correlation | O(N²) | Periodicity, pitch, power spectrum |

**Convolution vs correlation.** Both slide one signal across another.
Convolution time-reverses one signal (this is filtering); correlation
does not (this is similarity measurement). In the frequency domain
convolution is `X·H`, correlation is `X·conj(Y)`.

## Spectral analysis — windowing

| Window | Peak side lobe | Main lobe | Use when |
|---|---|---|---|
| Rectangular | −13 dB | Narrowest | Resolving close, equal-amplitude tones |
| Hamming | −43 dB | Moderate | Good general-purpose default |
| Hanning | −31 dB | Moderate | General purpose, better roll-off |
| Blackman | −58 dB | Widest | Detecting a weak tone near a strong one |

The fundamental trade-off is **main-lobe width vs side-lobe level**: you
cannot have both fine frequency resolution and low leakage.

## Sample rate conversion

| Operation | Effect | Filter placement |
|---|---|---|
| Decimation (÷M) | Lower rate | Low-pass *before* discarding samples |
| Interpolation (×L) | Higher rate | Low-pass *after* zero-stuffing |
| Rational (×L/M) | Resample by L/M | Upsample, filter, downsample |

Mnemonic: the low-pass filter always sits on the **high-rate side** of
the operation. Decimation risks aliasing; interpolation risks imaging.

## Multi-resolution analysis

| Transform | Complexity | Use when |
|---|---|---|
| Wavelet (DWT, Haar) | O(N) | Non-stationary signals; transient localisation |

**FFT vs wavelet.** The FFT gives excellent frequency resolution but no
time localisation — it tells you *what* frequencies exist, not *when*.
The wavelet transform trades some frequency precision for time
localisation, with adaptive resolution (fine timing at high frequency,
fine frequency at low). Stationary signal → FFT. Time-varying spectrum
(speech, ECG, seismic) → wavelets.
