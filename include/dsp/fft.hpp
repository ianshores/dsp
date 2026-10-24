// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Ian Shores
#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>
#if !defined(DSP_FFT_SCALAR) && (defined(__SSE2__) || defined(_M_X64))
#include <xmmintrin.h>
#define DSP_FFT_SSE 1
#ifdef DSP_FFT_FMA
#ifndef __FMA__
#error "FMA requires FMA compiler support and a compatible target CPU"
#endif
#include <immintrin.h>
#endif
#endif
#if defined(DSP_FFT_SSE) && !defined(DSP_DISABLE_AVX) &&                               \
    (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#define DSP_HAS_AVX_RECONSTRUCTION 1
#endif

namespace dsp {
// Even/odd real samples form an N/2-point complex FFT. Forward DIF and
// inverse DIT radix-4 stages share a mixed-radix permutation; no reorder pass.
// Spectrum contract: separate real/imaginary arrays, N/2+1 slots. Logical bin k
// lives at slot(k); DC is slot 0, Nyquist slot N/2, both with zero imaginary
// part. All ordinary slots use ordinary complex multiplication. IRs, history,
// accumulator, and inverse share this representation. The inverse has unit
// roundtrip gain. All plans/scratch belong to this instance, prepared off-thread.
// Inputs/outputs of the transform must not alias its private scratch buffers.
class RealFFT {
    friend struct FFTProbe;
    static std::size_t checked_size(std::size_t n) {
        if (n < 32 || (n & (n - 1)) != 0)
            throw std::invalid_argument("RealFFT size must be a power of two >= 32");
        return n;
    }
    using Complex = std::complex<float>;
    std::size_t n_, m_;
#ifdef DSP_HAS_AVX_RECONSTRUCTION
    const bool avx_ = __builtin_cpu_supports("avx") != 0;
#endif
    std::vector<float> real_root_r_, real_root_i_;
    struct PairRun {
        std::size_t a, b, count;
    };
    std::vector<PairRun> pair_runs_;
    std::vector<float> r_, i_;
    Complex get(std::size_t k) const noexcept { return Complex(r_[k], i_[k]); }
    void set(std::size_t k, Complex x) noexcept {
        r_[k] = x.real();
        i_[k] = x.imag();
    }
    struct Stage {
        std::size_t len, offset;
    };
    std::vector<Stage> stages_;
    std::vector<std::size_t> order_;
    std::vector<float> radix_r_, radix_i_;
#ifdef DSP_FFT_SSE
    struct Four {
        __m128 v;
        Four(__m128 x) : v(x) {}
        Four(float x) : v(_mm_set1_ps(x)) {}
        friend Four operator+(Four a, Four b) { return _mm_add_ps(a.v, b.v); }
        friend Four operator-(Four a, Four b) { return _mm_sub_ps(a.v, b.v); }
        friend Four operator*(Four a, Four b) { return _mm_mul_ps(a.v, b.v); }
    };
#endif
    template <class T> static T load(const float *p) noexcept {
#ifdef DSP_FFT_SSE
        if constexpr (!std::is_same<T, float>::value)
            return _mm_loadu_ps(p);
        else
#endif
            return *p;
    }
    template <class T> static void store(float *p, T x) noexcept {
#ifdef DSP_FFT_SSE
        if constexpr (!std::is_same<T, float>::value)
            _mm_storeu_ps(p, x.v);
        else
#endif
            *p = x;
    }
    template <class T> static T reverse_load(const float *p) noexcept {
#ifdef DSP_FFT_SSE
        if constexpr (!std::is_same<T, float>::value) {
            const auto x = _mm_loadu_ps(p - 3);
            return _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 1, 2, 3));
        } else
#endif
            return *p;
    }
    template <class T> static void reverse_store(float *p, T x) noexcept {
#ifdef DSP_FFT_SSE
        if constexpr (!std::is_same<T, float>::value)
            _mm_storeu_ps(p - 3, _mm_shuffle_ps(x.v, x.v, _MM_SHUFFLE(0, 1, 2, 3)));
        else
#endif
            *p = x;
    }
    template <bool packed> static std::size_t spectrum_index(std::size_t k) noexcept {
        if constexpr (packed)
            return (k / 8) * 16 + k % 8;
        else
            return k;
    }
    template <bool inverse, bool packed, class T>
    void reconstruct_pair(const float *sr, const float *si, float *dr, float *di, std::size_t a,
                          std::size_t b) noexcept {
        const auto sa = spectrum_index < inverse && packed > (a),
                   sb = spectrum_index < inverse && packed > (b);
        const auto da = spectrum_index < !inverse && packed > (a),
                   db = spectrum_index < !inverse && packed > (b);
        const T xr = load<T>(sr + sa), xi = load<T>(si + sa);
        const T yr = reverse_load<T>(sr + sb), yi = T(0.f) - reverse_load<T>(si + sb);
        const T ur = xr + yr, ui = xi + yi, vr = xr - yr, vi = xi - yi;
        const T wr = load<T>(real_root_r_.data() + a), wi = load<T>(real_root_i_.data() + a);
        const T tr = inverse ? wi * vr - wr * vi : wr * vi + wi * vr;
        const T ti = inverse ? wr * vr + wi * vi : wi * vi - wr * vr;
        store(dr + da, T(.5f) * (ur + tr));
        store(di + da, T(.5f) * (ui + ti));
        reverse_store(dr + db, T(.5f) * (ur - tr));
        reverse_store(di + db, T(.5f) * (ti - ui));
    }
#ifdef DSP_HAS_AVX_RECONSTRUCTION
    __attribute__((target("avx"),
                   always_inline)) static inline __m256 reverse_eight(__m256 x) noexcept {
        x = _mm256_permute_ps(x, _MM_SHUFFLE(0, 1, 2, 3));
        return _mm256_permute2f128_ps(x, x, 1);
    }
    template <bool inverse, bool packed>
    __attribute__((target("avx"))) void reconstruct_avx(const float *sr, const float *si, float *dr,
                                                        float *di) noexcept {
        const auto half = _mm256_set1_ps(.5f), zero = _mm256_setzero_ps();
        for (const auto &run : pair_runs_) {
            std::size_t j = 0;
            for (; j + 8 <= run.count; j += 8) {
                const auto a = run.a + j, b = run.b - j;
                const auto sa = spectrum_index < inverse && packed > (a),
                           sb = spectrum_index < inverse && packed > (b);
                const auto da = spectrum_index < !inverse && packed > (a),
                           db = spectrum_index < !inverse && packed > (b);
                const auto xr = _mm256_loadu_ps(sr + sa), xi = _mm256_loadu_ps(si + sa);
                const auto yr = reverse_eight(_mm256_loadu_ps(sr + sb - 7));
                const auto yi = _mm256_sub_ps(zero, reverse_eight(_mm256_loadu_ps(si + sb - 7)));
                const auto ur = _mm256_add_ps(xr, yr), ui = _mm256_add_ps(xi, yi);
                const auto vr = _mm256_sub_ps(xr, yr), vi = _mm256_sub_ps(xi, yi);
                const auto wr = _mm256_loadu_ps(real_root_r_.data() + a),
                           wi = _mm256_loadu_ps(real_root_i_.data() + a);
                const auto wrvi = _mm256_mul_ps(wr, vi), wivr = _mm256_mul_ps(wi, vr);
                const auto wrvr = _mm256_mul_ps(wr, vr), wivi = _mm256_mul_ps(wi, vi);
                const auto tr = inverse ? _mm256_sub_ps(wivr, wrvi) : _mm256_add_ps(wrvi, wivr);
                const auto ti = inverse ? _mm256_add_ps(wrvr, wivi) : _mm256_sub_ps(wivi, wrvr);
                _mm256_storeu_ps(dr + da, _mm256_mul_ps(half, _mm256_add_ps(ur, tr)));
                _mm256_storeu_ps(di + da, _mm256_mul_ps(half, _mm256_add_ps(ui, ti)));
                _mm256_storeu_ps(dr + db - 7,
                                 reverse_eight(_mm256_mul_ps(half, _mm256_sub_ps(ur, tr))));
                _mm256_storeu_ps(di + db - 7,
                                 reverse_eight(_mm256_mul_ps(half, _mm256_sub_ps(ti, ui))));
            }
            for (; j + 4 <= run.count; j += 4)
                reconstruct_pair<inverse, packed, Four>(sr, si, dr, di, run.a + j, run.b - j);
            for (; j < run.count; ++j)
                reconstruct_pair<inverse, packed, float>(sr, si, dr, di, run.a + j, run.b - j);
        }
    }
#endif
    template <bool inverse, bool packed>
    void reconstruct(const float *sr, const float *si, float *dr, float *di) noexcept {
#ifdef DSP_HAS_AVX_RECONSTRUCTION
        if (avx_) {
            reconstruct_avx<inverse, packed>(sr, si, dr, di);
            return;
        }
#endif
        for (const auto &run : pair_runs_) {
            std::size_t j = 0;
#ifdef DSP_FFT_SSE
            for (; j + 4 <= run.count; j += 4)
                reconstruct_pair<inverse, packed, Four>(sr, si, dr, di, run.a + j, run.b - j);
#endif
            for (; j < run.count; ++j)
                reconstruct_pair<inverse, packed, float>(sr, si, dr, di, run.a + j, run.b - j);
        }
    }
    template <bool inverse, bool twiddled, class T, class Roots>
    static void kernel4(T &ar, T &ai, T &br, T &bi, T &cr, T &ci, T &dr, T &di, std::size_t tw,
                        std::size_t q, Roots root) noexcept {
        auto multiply = [&](T &x, T &y, std::size_t w) {
            const T wr = root(w, false), wi = root(w, true) * T(inverse ? -1.f : 1.f);
#ifdef DSP_FFT_FMA
            if constexpr (!std::is_same<T, float>::value) {
                const T re = _mm_fmsub_ps(x.v, wr.v, _mm_mul_ps(y.v, wi.v));
                y = _mm_fmadd_ps(x.v, wi.v, _mm_mul_ps(y.v, wr.v));
                x = re;
            } else {
                const T re = std::fma(x, wr, -y * wi);
                y = std::fma(x, wi, y * wr);
                x = re;
            }
#else
            const T re = x * wr - y * wi;
            y = x * wi + y * wr;
            x = re;
#endif
        };
        if constexpr (inverse && twiddled) {
            multiply(br, bi, tw);
            multiply(cr, ci, tw + q);
            multiply(dr, di, tw + 2 * q);
        }
        const T ur = ar + cr, ui = ai + ci, vr = ar - cr, vi = ai - ci;
        const T sr = br + dr, si = bi + di, tr = br - dr, ti = bi - di;
        ar = ur + sr;
        ai = ui + si;
        cr = ur - sr;
        ci = ui - si;
        if constexpr (inverse) {
            br = vr - ti;
            bi = vi + tr;
            dr = vr + ti;
            di = vi - tr;
        } else {
            br = vr + ti;
            bi = vi - tr;
            dr = vr - ti;
            di = vi + tr;
        }
        if constexpr (!inverse && twiddled) {
            multiply(br, bi, tw);
            multiply(cr, ci, tw + q);
            multiply(dr, di, tw + 2 * q);
        }
    }
    template <bool inverse, class T>
    void butterfly4(std::size_t a, std::size_t q, std::size_t tw) noexcept {
        T ar = load<T>(&r_[a]), ai = load<T>(&i_[a]);
        T br = load<T>(&r_[a + q]), bi = load<T>(&i_[a + q]);
        T cr = load<T>(&r_[a + 2 * q]), ci = load<T>(&i_[a + 2 * q]);
        T dr = load<T>(&r_[a + 3 * q]), di = load<T>(&i_[a + 3 * q]);
        auto root = [&](std::size_t w, bool imaginary) {
            return load<T>((imaginary ? radix_i_ : radix_r_).data() + w);
        };
        kernel4<inverse, true>(ar, ai, br, bi, cr, ci, dr, di, tw, q, root);
        store(&r_[a], ar);
        store(&i_[a], ai);
        store(&r_[a + q], br);
        store(&i_[a + q], bi);
        store(&r_[a + 2 * q], cr);
        store(&i_[a + 2 * q], ci);
        store(&r_[a + 3 * q], dr);
        store(&i_[a + 3 * q], di);
    }
#ifdef DSP_FFT_SSE
    // Small terminal stages vectorize across independent groups. Q=1 uses
    // four groups, Q=2 uses two. All supported sizes contain complete batches.
    template <bool inverse, int Q> void terminal_stage(const Stage &stage) noexcept {
        for (std::size_t base = 0; base < m_; base += 16) {
            auto read = [&](const std::vector<float> &v, __m128 &a, __m128 &b, __m128 &c,
                            __m128 &d) {
                a = _mm_loadu_ps(v.data() + base);
                b = _mm_loadu_ps(v.data() + base + 4);
                c = _mm_loadu_ps(v.data() + base + 8);
                d = _mm_loadu_ps(v.data() + base + 12);
                if constexpr (Q == 1) {
                    _MM_TRANSPOSE4_PS(a, b, c, d);
                } else {
                    const auto aa = _mm_movelh_ps(a, c), bb = _mm_movehl_ps(c, a);
                    const auto cc = _mm_movelh_ps(b, d), dd = _mm_movehl_ps(d, b);
                    a = aa;
                    b = bb;
                    c = cc;
                    d = dd;
                }
            };
            __m128 ar0, br0, cr0, dr0, ai0, bi0, ci0, di0;
            read(r_, ar0, br0, cr0, dr0);
            read(i_, ai0, bi0, ci0, di0);
            Four ar(ar0), br(br0), cr(cr0), dr(dr0), ai(ai0), bi(bi0), ci(ci0), di(di0);
            auto root = [&](std::size_t w, bool imaginary) -> Four {
                const auto &v = imaginary ? radix_i_ : radix_r_;
                if constexpr (Q == 2)
                    return _mm_set_ps(v[w + 1], v[w], v[w + 1], v[w]);
                else
                    return Four(imaginary ? 0.f : 1.f);
            };
            kernel4<inverse, Q != 1>(ar, ai, br, bi, cr, ci, dr, di, stage.offset, Q, root);
            auto write = [&](std::vector<float> &v, __m128 a, __m128 b, __m128 c, __m128 d) {
                if constexpr (Q == 1) {
                    _MM_TRANSPOSE4_PS(a, b, c, d);
                } else {
                    const auto aa = _mm_movelh_ps(a, b), bb = _mm_movelh_ps(c, d);
                    const auto cc = _mm_movehl_ps(b, a), dd = _mm_movehl_ps(d, c);
                    a = aa;
                    b = bb;
                    c = cc;
                    d = dd;
                }
                _mm_storeu_ps(v.data() + base, a);
                _mm_storeu_ps(v.data() + base + 4, b);
                _mm_storeu_ps(v.data() + base + 8, c);
                _mm_storeu_ps(v.data() + base + 12, d);
            };
            write(r_, ar.v, br.v, cr.v, dr.v);
            write(i_, ai.v, bi.v, ci.v, di.v);
        }
    }
#endif
#ifdef DSP_HAS_AVX_RECONSTRUCTION
    __attribute__((target("avx"), always_inline)) static inline void
    multiply_eight(__m256 &x, __m256 &y, __m256 wr, __m256 wi) noexcept {
        const auto re = _mm256_sub_ps(_mm256_mul_ps(x, wr), _mm256_mul_ps(y, wi));
        y = _mm256_add_ps(_mm256_mul_ps(x, wi), _mm256_mul_ps(y, wr));
        x = re;
    }
    template <bool inverse>
    __attribute__((target("avx"))) void radix_stage_avx(const Stage &stage) noexcept {
        const auto q = stage.len / 4;
        for (std::size_t base = 0; base < m_; base += stage.len)
            for (std::size_t j = 0; j < q; j += 8) {
                const auto a = base + j, tw = stage.offset + j;
                auto ar = _mm256_loadu_ps(r_.data() + a), ai = _mm256_loadu_ps(i_.data() + a);
                auto br = _mm256_loadu_ps(r_.data() + a + q),
                     bi = _mm256_loadu_ps(i_.data() + a + q);
                auto cr = _mm256_loadu_ps(r_.data() + a + 2 * q),
                     ci = _mm256_loadu_ps(i_.data() + a + 2 * q);
                auto dr = _mm256_loadu_ps(r_.data() + a + 3 * q),
                     di = _mm256_loadu_ps(i_.data() + a + 3 * q);
                if constexpr (inverse) {
                    const auto sign = _mm256_set1_ps(-0.f);
                    multiply_eight(br, bi, _mm256_loadu_ps(radix_r_.data() + tw),
                                   _mm256_xor_ps(_mm256_loadu_ps(radix_i_.data() + tw), sign));
                    multiply_eight(cr, ci, _mm256_loadu_ps(radix_r_.data() + tw + q),
                                   _mm256_xor_ps(_mm256_loadu_ps(radix_i_.data() + tw + q), sign));
                    multiply_eight(
                        dr, di, _mm256_loadu_ps(radix_r_.data() + tw + 2 * q),
                        _mm256_xor_ps(_mm256_loadu_ps(radix_i_.data() + tw + 2 * q), sign));
                }
                const auto ur = _mm256_add_ps(ar, cr), ui = _mm256_add_ps(ai, ci);
                const auto vr = _mm256_sub_ps(ar, cr), vi = _mm256_sub_ps(ai, ci);
                const auto sr = _mm256_add_ps(br, dr), si = _mm256_add_ps(bi, di);
                const auto tr = _mm256_sub_ps(br, dr), ti = _mm256_sub_ps(bi, di);
                ar = _mm256_add_ps(ur, sr);
                ai = _mm256_add_ps(ui, si);
                cr = _mm256_sub_ps(ur, sr);
                ci = _mm256_sub_ps(ui, si);
                if constexpr (inverse) {
                    br = _mm256_sub_ps(vr, ti);
                    bi = _mm256_add_ps(vi, tr);
                    dr = _mm256_add_ps(vr, ti);
                    di = _mm256_sub_ps(vi, tr);
                } else {
                    br = _mm256_add_ps(vr, ti);
                    bi = _mm256_sub_ps(vi, tr);
                    dr = _mm256_sub_ps(vr, ti);
                    di = _mm256_add_ps(vi, tr);
                    multiply_eight(br, bi, _mm256_loadu_ps(radix_r_.data() + tw),
                                   _mm256_loadu_ps(radix_i_.data() + tw));
                    multiply_eight(cr, ci, _mm256_loadu_ps(radix_r_.data() + tw + q),
                                   _mm256_loadu_ps(radix_i_.data() + tw + q));
                    multiply_eight(dr, di, _mm256_loadu_ps(radix_r_.data() + tw + 2 * q),
                                   _mm256_loadu_ps(radix_i_.data() + tw + 2 * q));
                }
                _mm256_storeu_ps(r_.data() + a, ar);
                _mm256_storeu_ps(i_.data() + a, ai);
                _mm256_storeu_ps(r_.data() + a + q, br);
                _mm256_storeu_ps(i_.data() + a + q, bi);
                _mm256_storeu_ps(r_.data() + a + 2 * q, cr);
                _mm256_storeu_ps(i_.data() + a + 2 * q, ci);
                _mm256_storeu_ps(r_.data() + a + 3 * q, dr);
                _mm256_storeu_ps(i_.data() + a + 3 * q, di);
            }
    }
#endif
    bool use_fixed_schedule() const noexcept {
#ifdef DSP_HAS_AVX_RECONSTRUCTION
        return !avx_;
#else
        return true;
#endif
    }
    template <bool inverse> void radix_stage(const Stage &stage) noexcept {
        const auto q = stage.len / 4;
#ifdef DSP_HAS_AVX_RECONSTRUCTION
        if (avx_ && q >= 8) {
            radix_stage_avx<inverse>(stage);
            return;
        }
#endif
#ifdef DSP_FFT_SSE
        if (q == 1) {
            terminal_stage<inverse, 1>(stage);
            return;
        }
        if (q == 2) {
            terminal_stage<inverse, 2>(stage);
            return;
        }
#endif
        for (std::size_t base = 0; base < m_; base += stage.len) {
            std::size_t j = 0;
#ifdef DSP_FFT_SSE
            for (; j + 4 <= q; j += 4)
                butterfly4<inverse, Four>(base + j, q, stage.offset + j);
#endif
            for (; j < q; ++j)
                butterfly4<inverse, float>(base + j, q, stage.offset + j);
        }
    }
    void radix2_last() noexcept {
        for (std::size_t k = 0; k < m_; k += 2) {
            const float ar = r_[k], ai = i_[k], br = r_[k + 1], bi = i_[k + 1];
            r_[k] = ar + br;
            i_[k] = ai + bi;
            r_[k + 1] = ar - br;
            i_[k + 1] = ai - bi;
        }
    }
    template <class T> void half_first(const float *time, std::size_t j) noexcept {
        const auto q = m_ / 4;
        auto input_pair = [&](std::size_t k, T &real, T &imag) {
#ifdef DSP_FFT_SSE
            if constexpr (!std::is_same<T, float>::value) {
                const auto a = _mm_loadu_ps(time + 2 * k), b = _mm_loadu_ps(time + 2 * k + 4);
                real = _mm_shuffle_ps(a, b, _MM_SHUFFLE(2, 0, 2, 0));
                imag = _mm_shuffle_ps(a, b, _MM_SHUFFLE(3, 1, 3, 1));
            } else
#endif
            {
                real = T(time[2 * k]);
                imag = T(time[2 * k + 1]);
            }
        };
        T ar(0.f), ai(0.f), br(0.f), bi(0.f), cr(0.f), ci(0.f), dr(0.f), di(0.f);
        input_pair(j, ar, ai);
        input_pair(j + q, br, bi);
        auto root = [&](std::size_t w, bool imaginary) {
            return load<T>((imaginary ? radix_i_ : radix_r_).data() + w);
        };
        kernel4<false, true>(ar, ai, br, bi, cr, ci, dr, di, j, q, root);
        store(r_.data() + j, ar);
        store(i_.data() + j, ai);
        store(r_.data() + j + q, br);
        store(i_.data() + j + q, bi);
        store(r_.data() + j + 2 * q, cr);
        store(i_.data() + j + 2 * q, ci);
        store(r_.data() + j + 3 * q, dr);
        store(i_.data() + j + 3 * q, di);
    }
    void complex_half_forward(const float *time) noexcept {
        // Exactly half the time samples are zero. Fuse input packing with the
        // first DIF stage, without loading or storing the known-zero half.
        std::size_t j = 0;
#ifdef DSP_FFT_SSE
        for (; j + 4 <= m_ / 4; j += 4)
            half_first<Four>(time, j);
#endif
        for (; j < m_ / 4; ++j)
            half_first<float>(time, j);
        if (use_fixed_schedule() && m_ == 256)
            fixed_stages<false, 256, 64, 192>();
        else if (use_fixed_schedule() && m_ == 2048)
            fixed_stages<false, 2048, 512, 1536>();
        else {
            for (std::size_t k = 1; k < stages_.size(); ++k)
                radix_stage<false>(stages_[k]);
            if ((m_ >> (2 * stages_.size())) == 2)
                radix2_last();
        }
    }
    template <bool inverse, std::size_t M, std::size_t Len, std::size_t Offset>
    void fixed_stages() noexcept {
        if constexpr (Len >= 4) {
            if constexpr (inverse)
                fixed_stages<inverse, M, Len / 4, Offset + 3 * Len / 4>();
            constexpr auto q = Len / 4;
#ifdef DSP_FFT_SSE
            if constexpr (q == 1)
                terminal_stage<inverse, 1>({Len, Offset});
            else if constexpr (q == 2)
                terminal_stage<inverse, 2>({Len, Offset});
            else
#endif
                for (std::size_t base = 0; base < M; base += Len) {
                    std::size_t j = 0;
#ifdef DSP_FFT_SSE
                    for (; j + 4 <= q; j += 4)
                        butterfly4<inverse, Four>(base + j, q, Offset + j);
#endif
                    for (; j < q; ++j)
                        butterfly4<inverse, float>(base + j, q, Offset + j);
                }
            if constexpr (!inverse)
                fixed_stages<inverse, M, Len / 4, Offset + 3 * Len / 4>();
        } else if constexpr (Len == 2)
            radix2_last();
    }
    void complex_fft(bool inverse) noexcept {
        if (use_fixed_schedule() && m_ == 256) {
            if (inverse)
                fixed_stages<true, 256, 256, 0>();
            else
                fixed_stages<false, 256, 256, 0>();
            return;
        }
        if (use_fixed_schedule() && m_ == 2048) {
            if (inverse)
                fixed_stages<true, 2048, 2048, 0>();
            else
                fixed_stages<false, 2048, 2048, 0>();
            return;
        }
        const bool odd = ((m_ >> (2 * stages_.size())) == 2);
        if (inverse) {
            if (odd)
                radix2_last();
            for (auto s = stages_.rbegin(); s != stages_.rend(); ++s)
                radix_stage<true>(*s);
        } else {
            for (const auto &s : stages_)
                radix_stage<false>(s);
            if (odd)
                radix2_last();
        }
    }

  public:
    explicit RealFFT(std::size_t n)
        : n_(checked_size(n)), m_(n_ / 2), real_root_r_(m_), real_root_i_(m_), r_(m_), i_(m_) {
        order_.resize(m_);
        for (std::size_t k = 0; k < m_; ++k) {
            std::size_t v = k, len = m_, index = 0;
            while (len >= 4) {
                len /= 4;
                index += (v % 4) * len;
                v /= 4;
            }
            if (len == 2)
                index += v % 2;
            order_[k] = index;
        }
        for (std::size_t k = 0; k < m_; ++k) {
            const double angle = -6.2831853071795864769 * double(k) / n_;
            real_root_r_[slot(k)] = float(std::cos(angle));
            real_root_i_[slot(k)] = float(std::sin(angle));
        }
        // In digit-reversed order, conjugate partners form reverse-contiguous
        // runs. Recurse into the low quarter for natural bins divisible by four.
        for (std::size_t len = m_; len >= 4; len /= 4) {
            const auto q = len / 4;
            pair_runs_.push_back({q, len - 1, q});
            if (q / 2)
                pair_runs_.push_back({2 * q, 3 * q - 1, q / 2});
        }
        for (std::size_t len = m_; len >= 4; len /= 4) {
            const auto q = len / 4, offset = radix_r_.size();
            stages_.push_back({len, offset});
            for (std::size_t power = 1; power <= 3; ++power)
                for (std::size_t j = 0; j < q; ++j) {
                    const double angle = -6.2831853071795864769 * double(j * power) / len;
                    radix_r_.push_back(float(std::cos(angle)));
                    radix_i_.push_back(float(std::sin(angle)));
                }
        }
    }
    std::size_t slot(std::size_t k) const noexcept { return k == m_ ? m_ : order_[k]; }
    std::size_t memory_bytes() const noexcept {
        return sizeof(*this) + (real_root_r_.capacity() + real_root_i_.capacity()) * sizeof(float) +
               pair_runs_.capacity() * sizeof(PairRun) +
               (r_.capacity() + i_.capacity() + radix_r_.capacity() + radix_i_.capacity()) *
                   sizeof(float) +
               order_.capacity() * sizeof(std::size_t) + stages_.capacity() * sizeof(Stage);
    }
    template <bool half_input, bool packed = false>
    void forward_impl(const float *time, float *real, float *imag) noexcept {
        if constexpr (half_input)
            complex_half_forward(time);
        else {
            for (std::size_t j = 0; j < m_; ++j)
                set(j, Complex(time[2 * j], time[2 * j + 1]));
            complex_fft(false);
        }
        real[0] = get(0).real() + get(0).imag();
        imag[0] = 0;
        real[spectrum_index<packed>(m_)] = get(0).real() - get(0).imag();
        imag[spectrum_index<packed>(m_)] = 0;
        reconstruct<false, packed>(r_.data(), i_.data(), real, imag);
        const auto half = get(slot(m_ / 2));
        real[spectrum_index<packed>(slot(m_ / 2))] = half.real();
        imag[spectrum_index<packed>(slot(m_ / 2))] = -half.imag();
    }
    void forward(const float *time, float *real, float *imag) noexcept {
        forward_impl<false>(time, real, imag);
    }
    // Reads N/2 samples; the other N/2 samples are implicitly zero.
    void forward_half(const float *time, float *real, float *imag) noexcept {
        forward_impl<true>(time, real, imag);
    }
    template <bool packed = false>
    void inverse_impl(const float *real, const float *imag, float *time) noexcept {
        set(0, Complex(.5f * (real[0] + real[spectrum_index<packed>(m_)]),
                       .5f * (real[0] - real[spectrum_index<packed>(m_)])));
        reconstruct<true, packed>(real, imag, r_.data(), i_.data());
        set(slot(m_ / 2), Complex(real[spectrum_index<packed>(slot(m_ / 2))],
                                  -imag[spectrum_index<packed>(slot(m_ / 2))]));
        complex_fft(true);
        const float scale = 1.f / float(m_);
        for (std::size_t j = 0; j < m_; ++j) {
            time[2 * j] = get(j).real() * scale;
            time[2 * j + 1] = get(j).imag() * scale;
        }
    }
    void inverse(const float *real, const float *imag, float *time) noexcept {
        inverse_impl<false>(real, imag, time);
    }
    // Eight real values followed by eight imaginary values per physical bin group.
    // Requires N+16 floats, including separate DC/Nyquist and zero padding.
    void forward_half_packed(const float *time, float *spectrum) noexcept {
        forward_impl<true, true>(time, spectrum, spectrum + 8);
    }
    void inverse_packed(const float *spectrum, float *time) noexcept {
        inverse_impl<true>(spectrum, spectrum + 8, time);
    }
};
} // namespace dsp
