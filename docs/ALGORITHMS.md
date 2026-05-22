# DSP Algorithms — Complexity & Trade-off Cheat-Sheet

Every algorithm in this guide, with its cost and the decision rule for
when to reach for it. `N` is signal length; `M` is filter/kernel length.

## Frequency-domain transforms

| Algorithm | Complexity | Output | Use when |
|---|---|---|---|
| DFT  | O(N²)      | N complex bins | Reference / teaching only — too slow for real N |
| FFT  | O(N log N) | N complex bins | Any practical spectral analysis; identical result to DFT |
| DCT  | O(N²) here | N real coeffs  | Compression — energy compacts into few coefficients |
| Hilbert | O(N log N) | N real samples | Building the analytic signal — envelope, instantaneous frequency |

**DFT vs FFT.** Same mathematical result. The FFT exploits symmetry
(`W_N^{k+N/2} = -W_N^k`) and periodicity to divide and conquer. For
N = 1024 that is roughly a 100× speedup. Never compute the DFT directly
in production.

**DFT vs DCT.** Both move a signal to the frequency domain. The DCT uses
only cosines, stays real for real input, and handles block edges better
(implicit symmetric extension instead of periodic wraparound). Its strong
energy compaction is why JPEG and MPEG use it. If a question is about
image or video compression, the answer is DCT.

**Hilbert transform.** Shifts every frequency component by −90°. On its
own that is a curiosity; its value is the **analytic signal**
`z(t) = x(t) + j·H{x(t)}` — a complex signal with only positive
frequencies. Two quantities fall straight out of it: the **envelope**
`|z(t)|` (instantaneous amplitude — AM demodulation, onset detection)
and the **instantaneous frequency**, the time derivative of `arg z(t)`
(FM demodulation, chirp analysis). It is computed in the frequency
domain — FFT, keep DC and Nyquist, double the positive-frequency bins,
zero the negative ones, inverse FFT — which is exact, with no filter
design. The Wigner-Ville distribution also uses the analytic signal to
halve its cross terms.

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

## Advanced spectral estimation

The FFT periodogram (which the windows above prepare a signal for) has a
resolution fixed by the record length — two tones closer than about 1/N
blur into one peak. Advanced estimators do better by assuming structure.

| Method | Type | Resolution | Cost | Best when |
|---|---|---|---|---|
| AR | Parametric | High | Low–medium | Short records, sharp resonant peaks |
| ARMA | Parametric | High | High | Spectra with both peaks and notches |
| MUSIC | Subspace | Very high | High (grid search) | Closely spaced tones, direction finding |
| ESPRIT | Subspace | Very high | Medium | Fast estimation, no spectral search |

**Parametric methods** fit a model and read the spectrum off the model
parameters. **AR** treats the signal as an all-pole filter driven by
white noise; its coefficients come from the Yule-Walker equations
(solved efficiently by the Levinson-Durbin recursion) or from Burg's
forward-backward method, which is often more accurate for short records.
**ARMA** adds a moving-average (zeros) part for spectra that have notches
as well as peaks; its parameters are harder to fit, here via the
modified Yule-Walker method (AR part from the high-lag autocorrelations,
MA part from the residual).

**Subspace methods** split the data's covariance matrix, via an
eigen-decomposition, into a signal subspace (the large eigenvalues) and a
noise subspace (the small ones). **MUSIC** forms a pseudospectrum whose
peaks — where a steering vector is orthogonal to the noise subspace —
mark the tone frequencies. **ESPRIT** instead exploits the shift-
invariance of the signal subspace to solve for the frequencies directly,
with no spectral search. Both reach "super-resolution": they separate
tones far closer than the FFT's 1/N bin width.

The eigen-decomposition here uses the Jacobi method — repeated plane
rotations that zero the largest off-diagonal element — which suits the
small symmetric covariance matrices these estimators build. ESPRIT
recovers each tone from the eigenvalues of a shift-rotation matrix; its
accuracy is good but somewhat more sensitive to the covariance-matrix
size than MUSIC or the AR methods.

These methods are the workhorses of radar and sonar, array signal
processing (direction-of-arrival estimation), tone detection, and speech
analysis — anywhere the FFT's resolution falls short.

## Time-frequency analysis

Both the FFT and the estimators above produce one spectrum for the whole
signal. But speech, music, and radar chirps are *non-stationary* — their
frequency content changes over time. Time-frequency analysis recovers
the missing time axis: it shows *when* each frequency occurs.

| Method | Resolution | Cross terms | Invertible | Best for |
|---|---|---|---|---|
| STFT | Fixed (window-set) | None | Yes | Speech, audio, general purpose |
| Filter bank | Flexible tiling | None | Yes | Compression, multirate, channelizers |
| Wigner-Ville | Very high | Strong | Yes (with care) | Chirps, transients, high-resolution work |
| Wavelet (DWT) | Multi-scale | Low | Yes | Transients, denoising (see wavelet/) |

The **STFT** slides a short window along the signal and FFTs each piece;
the squared magnitude is the spectrogram. It is simple, invertible by
overlap-add, and real-time capable. Its one limitation is a *fixed*
resolution — the window length sets a single time/frequency trade-off
for all frequencies (the time-frequency uncertainty principle): a short
window localises time well but frequency poorly, and a long window does
the reverse.

A **filter bank** splits the signal into frequency subbands with a set
of bandpass filters, processes each band, then reconstructs. The
two-channel quadrature-mirror-filter (QMF) bank here decimates each of
its low and high subbands by two, so the total sample count is unchanged
(critical sampling). The QMF mirror relation between the analysis and
synthesis filters makes the aliasing from that decimation cancel,
giving near-perfect reconstruction. This subband structure underlies
audio codecs (MP3, AAC) and communications channelizers. (The guide's
wavelet transform is itself a wavelet filter bank.)

The **Wigner-Ville distribution** is quadratic, not a linear transform:
it is the Fourier transform of the signal's instantaneous
autocorrelation. It achieves the sharpest time-frequency localisation of
any method here, with no window at all — but for a multi-component
signal it shows **cross terms**: oscillating artefacts midway between
every pair of true components. That is its defining drawback. The
pseudo-WVD applies a smoothing window along the lag axis to suppress
those cross terms, trading some resolution for a cleaner picture. The
implementation works on the analytic signal, which halves the cross
terms and removes the frequency aliasing of a real input; its frequency
axis spans [0, 0.5) across N bins.

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

## Error control coding

Reliable digital communication: detect, or correct, the bit errors that
noise and channel distortion introduce.

### Error detection

| Technique | Redundancy | Catches | Use when |
|---|---|---|---|
| Parity | 1 bit | Any odd number of bit flips | Cheapest possible check |
| Checksum | 16 bits | General corruption | TCP/IP packets; weak vs multi-bit |
| CRC-32 | 32 bits | Burst errors | Ethernet, storage, wireless frames |

Detection only flags an error — it enables a retransmission request
(ARQ). It cannot repair anything.

### Forward error correction

| Code | Granularity | Corrects | Decoder | Use when |
|---|---|---|---|---|
| Hamming(7,4) | bit | 1-bit error / 2-bit detect | Syndrome lookup | ECC memory, satellite |
| Reed-Solomon | symbol | t symbol errors, n−k=2t | Berlekamp-Massey | CD/DVD, QR codes, DVB |
| Convolutional | bit stream | Scattered errors | Viterbi (trellis) | GSM, deep space |
| LDPC | bit | Near Shannon limit | Iterative belief propagation | 5G, Wi-Fi 6 |
| Turbo | bit | Near Shannon limit | Iterative (two decoders) | 4G LTE, deep space |

FEC adds redundancy up front so the receiver repairs errors with no
round-trip — essential for broadcast, storage, and real-time links.
Reed-Solomon's symbol-level math is what makes it strong against bursts:
a damaged symbol counts once however many of its bits flipped. This
guide implements Hamming, Reed-Solomon, the convolutional/Viterbi pair,
and LDPC; Turbo codes are described for completeness but not coded.

### LDPC codes

| Decoder | Input | Strength | Cost |
|---|---|---|---|
| Bit-flipping | Hard bits | Weak — corrects few errors | Integer, one pass cheap |
| Sum-product (BP) | Soft LLRs | Near Shannon limit | Iterative, tanh/atanh per edge |
| Min-sum | Soft LLRs | Close to sum-product | Iterative, no transcendentals |

An LDPC code is a *sparse* parity-check matrix H: a codeword c satisfies
H·c = 0 over GF(2). H is best seen as a bipartite **Tanner graph** —
variable nodes (code bits) on one side, check nodes (parity equations)
on the other, an edge per 1 in H.

Decoding passes messages along the graph edges. Bit-flipping passes hard
votes: each failing check is unhappy with its bits, and the most-suspect
bit is flipped each pass. Sum-product passes *probabilities* as
log-likelihood ratios — variable and check nodes refine each other's
belief every iteration. **Min-sum** keeps the same schedule but replaces
the exact tanh check-node rule with a minimum-of-magnitudes
approximation: no transcendental functions, fixed-point friendly, and —
with a normalising scale factor around 0.75 — within a hair of
sum-product. That is why min-sum is the decoder in real 5G and Wi-Fi 6
silicon. Sum-product's use of the demodulator's soft confidence is what
puts LDPC within a fraction of a decibel of the Shannon limit.

Decoders are characterised by a **BER sweep**: transmit many codewords
over a simulated AWGN channel at a range of noise levels and plot the
bit-error rate. The curve drops steeply once the noise clears a
threshold — the LDPC "waterfall". The sweep also makes the decoder
ranking concrete: soft decoders (sum-product, min-sum) sit far below
hard bit-flipping, and min-sum tracks sum-product closely.

LDPC vs Reed-Solomon: RS gives a *hard guarantee* of t corrected symbols
via exact algebra. LDPC offers *no such guarantee* — it is probabilistic
— but on a noisy channel it corrects far more, far closer to the limit.
The price is iterative decoding: more computation and some latency, the
same trade-off the guide meets in FIR-vs-IIR and direct-vs-FFT.

**Detection vs correction.** Detection is cheap but needs a feedback
channel to be useful (request a resend). Correction spends more
bandwidth on redundancy but needs no round-trip. Real systems often
combine them (hybrid ARQ): FEC handles common errors, a CRC catches the
rare residue and triggers a resend.

### Interleaving

| Type | Streaming | Latency / memory | Used in |
|---|---|---|---|
| Block | Needs a full R×C block | Higher | DVB-T, storage |
| Convolutional | Continuous | ~half the block cost | DVB, ADSL |

Interleaving adds **no redundancy and corrects nothing on its own** — it
only rearranges symbols. The transmitter scrambles the symbol order; the
receiver unscrambles it. A burst that hits the scrambled stream is
spread thin once de-interleaved, so each FEC codeword sees only a few
errors — within its correction limit. This is what lets a code rated for
`t` errors survive bursts far longer than `t`.

The block interleaver fills an R×C matrix row by row and reads it column
by column; map one codeword per row and a burst is dealt round-robin
across codewords. The convolutional interleaver uses a bank of
staggered delay lines, streams continuously with no block boundary, and
reaches the same spreading at roughly half the latency and memory.
Interleaver depth trades burst protection against latency.

  transmit:  encode → interleave → channel
  receive:   channel → deinterleave → decode

### DSP's role in error handling

| Function | What it does |
|---|---|
| Channel equalization | Adaptive filter (LMS/RLS) inverts channel distortion before decoding |
| Soft-decision decoding | Decoder uses analog confidence, not hard 0/1 — ~2 dB coding gain |
| Interleaving | Scatters burst errors across codewords so each is correctable |
| Carrier/timing recovery | Tracks frequency and symbol phase for accurate sampling |

The decoder is only as good as its input. The DSP front end cleans the
waveform — equalizing channel distortion, recovering carrier and timing
— so that the FEC stage sees errors it can actually correct. Soft-
decision decoding closes the loop by passing the demodulator's
confidence straight into the decoder instead of discarding it.

## Adaptive filters

A fixed filter has constant coefficients; an adaptive filter adjusts
them sample by sample to minimise the error between its output and a
desired reference. One filter then copes with an unknown or changing
system.

| Algorithm | Cost/sample | Convergence | Step tuning |
|---|---|---|---|
| LMS | O(L) | Slow | Hard — sensitive to mu and input power |
| NLMS | O(L) | Medium | Easy — step normalised by input power |
| RLS | O(L²) | Fast | Forgetting factor lambda |

All three share the skeleton: input `x(n)`, desired `d(n)`, error
`e(n) = d(n) − y(n)`, and a rule that nudges the weights to shrink
`e(n)`. They differ only in that rule. **LMS** takes a stochastic-
gradient step `w += mu·e·x` — cheap and robust but slow, and its
effective rate drifts with the input level. **NLMS** divides the step
by the input power `‖x‖²`, which removes that sensitivity and makes
`mu` easy to pick; it is the usual choice for echo and noise
cancellation. **RLS** recursively solves the full least-squares problem
over all past samples (exponentially weighted by the forgetting factor
`lambda`), maintaining an inverse-correlation matrix. It converges far
faster and tracks change better, at `O(L²)` cost — the same
fast-but-expensive trade-off seen in FIR-vs-IIR and direct-vs-FFT.

Adaptive filters power equalization (invert an unknown channel), system
identification (learn an unknown impulse response), noise and echo
cancellation (subtract noise correlated with a reference), prediction,
and beamforming.

## Kalman filtering and state estimation

Adaptive filters *learn an unknown system*. The Kalman filter solves a
different problem: it *estimates the hidden state* of a *known* dynamic
system from noisy measurements. It keeps a state estimate and a
covariance (how uncertain that estimate is), and runs a two-phase cycle.

| Phase | Effect on the state | Effect on uncertainty |
|---|---|---|
| Predict | Advance by the motion model `x ← F x` | Grows: `P ← F P Fᵀ + Q` |
| Update | Correct toward the measurement | Shrinks: `P ← (I − K H) P` |

The **Kalman gain** `K = P Hᵀ S⁻¹` is the heart of it: it weighs the
model against the measurement by their relative uncertainty. A precise
measurement pulls the estimate hard; a noisy one barely moves it. For a
linear system with Gaussian noise the Kalman filter is the optimal
(minimum-mean-square-error) estimator.

**Sensor fusion** is the headline use. Several noisy sensors measure the
same quantity; the update step naturally combines them, each pulling the
estimate in proportion to its precision, so the fused estimate beats any
single sensor. The static, single-sample form is inverse-variance
weighting: fuse readings with weights `1/variance`, and the fused
variance drops below every individual sensor's. A GPS + IMU navigation
filter is the classic running example.

The **Extended Kalman Filter (EKF)** handles *nonlinear* systems —
radar range/bearing, GPS — by linearising the model at each step:
the caller supplies the nonlinear transition and measurement functions
plus their Jacobians, and the EKF propagates the covariance through
those Jacobians. It is the standard nonlinear estimator, though (unlike
the linear filter) no longer provably optimal.

Connection to the adaptive filters above: **RLS is a Kalman filter** for
a particular state model — the "state" being the filter weights, assumed
constant. The two share the recursive predict/correct structure; the
Kalman filter is the general form.

## Modulation — QAM and OFDM

How bits become a transmitted waveform, and how a receiver gets them
back. This is the physical-layer counterpart to the error-control codes.

### QAM

| Order | Bits/symbol | Robustness |
|---|---|---|
| QPSK | 2 | Highest — points far apart |
| 16-QAM | 4 | Medium |
| 64-QAM | 6 | Lowest — needs a clean channel |

QAM packs a group of bits into one point of a 2-D (I/Q) constellation.
Higher orders carry more bits per symbol but pack points closer
together, so they need a higher SNR — the rate-vs-robustness trade-off
behind adaptive modulation in Wi-Fi and 5G. Gray coding labels adjacent
points so a slip to a neighbour costs a single bit. The soft demapper
emits a per-bit LLR instead of a hard bit, which is what a soft FEC
decoder (LDPC, Viterbi) needs.

### OFDM

OFDM splits data across many narrow subcarriers sent in parallel. Its
defining trick is that it *is* an FFT:

- **modulator** = load one QAM symbol per subcarrier, then IFFT
- **demodulator** = FFT, then read each subcarrier back
- **cyclic prefix** = copy each block's tail to its front; this turns
  the channel's linear convolution into circular convolution, so
  multipath becomes one complex gain per subcarrier

Equalization is therefore a single complex divide per subcarrier —
vastly cheaper than the long time-domain equalizer a single-carrier
system needs for the same multipath. This is why OFDM is the modulation
of Wi-Fi, LTE, 5G, DVB-T, and ADSL.

One implementation note: the IFFT's 1/N normalisation leaves the
time-domain signal at 1/N the power of the constellation. The modulator
multiplies by sqrt(N) (and the demodulator divides it back) so a
time-domain noise level corresponds to a meaningful per-symbol SNR —
the standard unitary OFDM convention.

### Coded OFDM — the full chain

Real systems wrap OFDM around an FEC code:

```
TX:  bits -> FEC encode -> QAM map -> OFDM (IFFT+CP) -> channel
RX:  channel -> OFDM (CP-strip+FFT) -> equalize -> QAM soft-demap
              -> FEC decode -> bits
```

OFDM tames multipath; the code cleans up the residual noise errors,
including the bits on deeply faded subcarriers. The QAM soft demapper
feeds per-bit LLRs (scaled by the channel's noise level) into the LDPC
decoder — the soft-information handoff that makes coded OFDM perform
close to channel capacity. This FEC + OFDM pairing is the data-carrying
core of every modern broadband wireless standard.

## Communications front end — pulse shaping & synchronization

The receiver-side analog refinements that sit before the demapper.

| Function | Algorithm | Role |
|---|---|---|
| Pulse shaping | Root-raised-cosine FIR | Bandlimits symbols, ISI-free at sample instants |
| Carrier recovery | Decision-directed PLL | Cancels frequency/phase offset that spins the constellation |
| Timing recovery | Gardner detector + interpolator | Locks the sampling instant to the symbol centre |

A symbol stream cannot be sent as bare impulses — sharp edges need
infinite bandwidth. The **root-raised-cosine** pulse is bandlimited yet
zero at every symbol instant but its own, so neighbours do not interfere
(the Nyquist ISI-free property). Splitting the raised cosine into a
transmit RRC and a receive RRC also makes the receive filter the matched
filter that maximises SNR. The roll-off factor trades bandwidth against
pulse-tail decay.

**Carrier** and **timing recovery** are feedback loops — the same
measure-error-and-nudge structure as the LMS equalizer. The carrier PLL
is decision-directed: it slices each sample to the nearest constellation
point and uses the angle between them as the phase error, with a
second-order loop filter that also tracks a constant frequency offset.
The Gardner timing detector uses the mid-symbol sample to estimate the
fractional timing offset, which an interpolating resampler then corrects.
Loop gain sets the speed/stability trade-off: small gain tracks smoothly,
large gain locks fast but jitters.

## Image processing — 2-D DSP

Every transform and filter here is **2-D**, and almost all are
*separable*: a 1-D operation applied to rows, then to columns. The 1-D
machinery from earlier modules does the real work.

### 2-D transforms

| Transform | Built from | Use |
|---|---|---|
| 2-D FFT | 1-D FFT on rows + columns | Spatial-frequency analysis, fast 2-D convolution |
| 2-D DCT | 1-D DCT on rows + columns | The JPEG transform — energy compaction for compression |
| 2-D Haar wavelet | 1-D Haar on rows + columns | The JPEG 2000 transform — multi-resolution analysis |

The 8x8 block DCT is the exact transform JPEG applies to each pixel
block; its energy compaction packs most of a block into a few
low-frequency coefficients. The 2-D wavelet instead splits an image into
four quadrants — a half-size approximation plus horizontal, vertical, and
diagonal detail — which is why JPEG 2000 uses it: the multi-resolution
structure scales and compresses better than fixed blocks.

### Spatial filtering

| Filter | Kind | Effect |
|---|---|---|
| Gaussian / box blur | Linear convolution | Smooths noise and detail |
| Sharpen | Linear convolution | Boosts local contrast |
| Sobel | Linear convolution | Gradient edge detector |
| Laplacian | Linear convolution | Orientation-free edge/detail operator |
| Median | Non-linear | Removes salt-and-pepper noise, keeps edges |

Blurring, sharpening, and edge detection are all **2-D convolution** —
slide a small kernel, take a weighted neighbourhood sum — the image-
domain counterpart of 1-D filtering. The **median filter** is the
exception: it replaces each pixel with its neighbourhood median, cannot
be written as a kernel, and removes impulse noise that a linear blur only
smears around.

### Point operators

Histogram equalization remaps pixels through the cumulative histogram so
the output uses the full tonal range — a contrast fix that depends only
on each pixel's own value. Otsu's method picks the threshold that best
splits the histogram into two classes by maximising the between-class
variance, binarizing an image with no manual tuning.

## Array signal processing

A single sensor measures a signal but not its direction. An array of
sensors can: a wave arriving at an angle reaches the sensors at slightly
different times, and that pattern of phase shifts encodes the direction.

| Method | Type | Resolution | Needs #sources | Cost |
|---|---|---|---|---|
| Delay-and-sum | Beamformer | Beamwidth ~ 1/aperture | No | Low |
| MVDR / Capon | Beamformer | High (sharp nulls) | No | Medium (matrix inverse) |
| Spatial MUSIC | Subspace DOA | Super-resolution | Yes | High (eigen + search) |
| Spatial ESPRIT | Subspace DOA | Super-resolution | Yes | Medium (eigen, no search) |

The model is a uniform linear array (ULA) of M sensors spaced d
wavelengths apart. A plane wave from angle theta gives sensor k a phase
advance of 2*pi*d*k*sin(theta); the vector of those phases is the
**steering vector** a(theta). With several time snapshots the M x M
covariance matrix R captures the spatial structure, and every method
works from R.

**Delay-and-sum** steers the array across all angles and measures the
power a(theta)^H R a(theta) — a spatial filter whose peaks mark the
sources. Simple, but its beamwidth is set by the array aperture, so it
blurs close sources. **MVDR (Capon)** keeps unit gain on the look
direction while minimising total output power, 1/(a^H R^-1 a), placing
sharp nulls on interferers and resolving far better.

**MUSIC** and **ESPRIT** are the same algorithms the guide uses for
temporal frequency estimation (see the advanced spectral estimation
section) — a direction of arrival and a temporal frequency are both just
the phase slope of a complex exponential. The spatial versions
eigen-decompose the *complex Hermitian* covariance into a signal and a
noise subspace. MUSIC scans a pseudospectrum whose peaks are the arrival
angles; ESPRIT exploits the ULA's shift-invariance (two overlapping
subarrays differ only by a phase rotation) to solve for the angles
directly, with no search. Both reach super-resolution: they separate
sources closer than the conventional beamwidth.

The complex eigen-decomposition uses a Hermitian Jacobi method (the
complex analogue of the real Jacobi solver used for the temporal
methods). ESPRIT recovers each angle from the eigenvalues of a
shift-rotation matrix, found by a shifted QR iteration with deflation —
a plain QR converges poorly because those eigenvalues all lie on the
unit circle.

These are the methods behind radar and sonar direction finding, smart
antennas in wireless systems, and acoustic source localisation.
