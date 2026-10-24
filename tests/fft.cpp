// SPDX-License-Identifier: MIT
#include "dsp/fft.hpp"
#include <complex>
#include <iostream>
#include <random>
#include <stdexcept>
int main() {
    try {
        for (std::size_t bad : {0u, 31u, 33u, 1000u}) {
            try {
                shores_dsp::RealFFT plan(bad);
                throw std::runtime_error("invalid FFT size was accepted");
            } catch (const std::invalid_argument &) {
            }
        }
        double worst = 0, rms_worst = 0, spectrum_worst = 0;
        std::mt19937 gen(171);
        std::uniform_real_distribution<float> noise(-1.f, 1.f);
        int cases = 0;
        for (std::size_t n = 32; n <= 16384; n *= 2) {
            shores_dsp::RealFFT plan(n);
            std::vector<float> in(n), out(n), real(n / 2 + 1), imag(n / 2 + 1);
            std::vector<bool> seen(n / 2 + 1, false);
            for (std::size_t k = 0; k <= n / 2; ++k) {
                const auto slot = plan.slot(k);
                if (slot > n / 2 || seen[slot])
                    throw std::runtime_error("spectrum mapping is not a permutation");
                seen[slot] = true;
            }
            for (int type = 0; type < 6; ++type) {
                for (std::size_t i = 0; i < n; ++i) {
                    if (type == 0)
                        in[i] = i == 0 ? 1.f : 0.f;
                    if (type == 1)
                        in[i] = 1.f;
                    if (type == 2)
                        in[i] = i & 1 ? -1.f : 1.f;
                    if (type == 3)
                        in[i] = float(std::cos(6.283185307179586 * 7 * i / n));
                    if (type == 4)
                        in[i] = noise(gen);
                    if (type == 5)
                        in[i] = float(noise(gen) * std::exp(-30. * i / n));
                }
#ifdef DSP_TEST_PACKED
                std::fill(in.begin() + n / 2, in.end(), 0.f);
#ifdef DSP_TEST_PACKED8
                constexpr std::size_t width = 8;
#else
                constexpr std::size_t width = 4;
#endif
                std::vector<float> packed(n + 2 * width, 0.f);
                plan.forward_half_packed(in.data(), packed.data());
                for (std::size_t k = 0; k <= n / 2; ++k) {
                    real[k] = packed[(k / width) * 2 * width + k % width];
                    imag[k] = packed[(k / width) * 2 * width + k % width + width];
                }
#elif defined(DSP_TEST_HALF_FORWARD)
                std::fill(in.begin() + n / 2, in.end(), 0.f);
                plan.forward_half(in.data(), real.data(), imag.data());
#else
                plan.forward(in.data(), real.data(), imag.data());
#endif
                // Independent double-precision full complex reference, checking
                // every bin (in addition to the direct DFT samples below).
                std::vector<std::complex<double>> reference(in.begin(), in.end());
                for (std::size_t k = 1, j = 0; k < n; ++k) {
                    auto bit = n / 2;
                    for (; j & bit; bit >>= 1)
                        j ^= bit;
                    j ^= bit;
                    if (k < j)
                        std::swap(reference[k], reference[j]);
                }
                for (std::size_t len = 2; len <= n; len *= 2) {
                    for (std::size_t base = 0; base < n; base += len)
                        for (std::size_t k = 0; k < len / 2; ++k) {
                            const double a = -6.2831853071795864769 * k / len;
                            const auto v = reference[base + k + len / 2] *
                                           std::complex<double>(std::cos(a), std::sin(a));
                            const auto u = reference[base + k];
                            reference[base + k] = u + v;
                            reference[base + k + len / 2] = u - v;
                        }
                }
                for (std::size_t k = 0; k <= n / 2; ++k) {
                    const auto slot = plan.slot(k);
                    const double e =
                        std::abs(reference[k] - std::complex<double>(real[slot], imag[slot]));
                    spectrum_worst = std::max(spectrum_worst, e);
                    if (e > 3.e-6 * n)
                        throw std::runtime_error("FFT differs from independent complex reference");
                }
                if (imag[plan.slot(0)] != 0 || imag[plan.slot(n / 2)] != 0)
                    throw std::runtime_error("DC/Nyquist are not real");
                for (std::size_t k = 0; k <= n / 2; ++k) {
                    if (n > 256 && k != 0 && k != 1 && k != 2 && k != 7 && k != n / 7 &&
                        k != n / 4 && k != n / 2 - 1 && k != n / 2)
                        continue;
                    std::complex<double> ref{};
                    for (std::size_t j = 0; j < n; ++j) {
                        const double a = -6.2831853071795864769 * double(j) * k / n;
                        ref += double(in[j]) * std::complex<double>(std::cos(a), std::sin(a));
                    }
                    const auto slot = plan.slot(k);
                    double err = std::abs(ref - std::complex<double>(real[slot], imag[slot]));
                    spectrum_worst = std::max(spectrum_worst, err);
                    if (err > 3.e-6 * n)
                        throw std::runtime_error("FFT spectrum differs from double DFT");
                }
#ifdef DSP_TEST_PACKED
                plan.inverse_packed(packed.data(), out.data());
#else
                plan.inverse(real.data(), imag.data(), out.data());
#endif
                double sq = 0;
                for (std::size_t j = 0; j < n; ++j) {
                    double e = double(out[j]) - in[j];
                    worst = std::max(worst, std::abs(e));
                    sq += e * e;
                }
                rms_worst = std::max(rms_worst, std::sqrt(sq / n));
                ++cases;
            }
        }
        if (worst > 3.e-6)
            throw std::runtime_error("real FFT roundtrip error exceeds tolerance");
        std::cout
            << "PASS: " << cases
            << " FFT cases, N=32..16384; impulse/DC/Nyquist/tone/random/decay; max roundtrip error "
            << worst << ", worst RMS " << rms_worst << ", max sampled spectrum error "
            << spectrum_worst << '\n';
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
