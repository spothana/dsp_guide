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
