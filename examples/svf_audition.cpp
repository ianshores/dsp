// SPDX-License-Identifier: MIT
// Original synthetic audition material: no recordings, speech assets or IRs.
#include "dsp/biquad.hpp"
#include "dsp/svf.hpp"
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

using namespace shores_dsp;
constexpr std::size_t rate = 48000, frames = rate * 16;
constexpr double tau = 6.2831853071795864769;

static void little(std::ofstream &file, std::uint32_t value, int bytes) {
    for (int i = 0; i < bytes; ++i) file.put(char((value >> (8 * i)) & 255u));
}
static void wav(const std::filesystem::path &path, const std::vector<float> &samples,
                int channels, double gain) {
    std::ofstream file(path, std::ios::binary);
    file.exceptions(std::ios::badbit | std::ios::failbit);
    const auto size = std::uint32_t(samples.size() * 2);
    file.write("RIFF", 4); little(file, 36 + size, 4); file.write("WAVEfmt ", 8);
    little(file, 16, 4); little(file, 1, 2); little(file, std::uint32_t(channels), 2);
    little(file, rate, 4); little(file, rate * std::uint32_t(channels) * 2, 4);
    little(file, std::uint32_t(channels) * 2, 2); little(file, 16, 2);
    file.write("data", 4); little(file, size, 4);
    for (float sample : samples) {
        const auto value = std::int16_t(std::lround(sample * gain * 32767.));
        little(file, std::uint16_t(value), 2);
    }
}

int main(int argc, char **argv) {
    try {
        const std::filesystem::path directory = argc > 1 ? argv[1] : "build/audition";
        std::filesystem::create_directories(directory);
        const auto dry_path = directory / "01-dry-synthetic.wav";
        const auto bq_path = directory / "02-biquad.wav";
        const auto svf_path = directory / "03-svf.wav";
        const auto stereo_path = directory / "04-left-biquad-right-svf.wav";
        for (const auto &path : {dry_path, bq_path, svf_path, stereo_path, directory / "README.txt"})
            if (std::filesystem::exists(path)) throw std::runtime_error("audition target exists; choose a new output directory");
        std::vector<float> dry(frames), biquad(frames), svf(frames), stereo(frames * 2);
        std::vector<SvfCoefficients> svf_events(frames);
        std::vector<BiquadCoefficients> bq_events(frames);
        std::uint32_t rng = 0x50395u;
        // Harmonic, vowel-like bursts interleaved with noise/fricative pulses.
        // Repeated motifs and broad spectrum make filtering easy to compare.
        for (std::size_t i = 0; i < frames; ++i) {
            const double t = double(i) / rate, within = std::fmod(t, 0.32);
            const double envelope = std::sin(3.14159265358979323846 * std::min(1., within / 0.26));
            double voice = 0;
            for (int harmonic = 1; harmonic <= 60; ++harmonic) {
                const double hz = 140 * harmonic;
                const double weight = 0.03 + std::exp(-std::pow((hz - 650) / 350, 2)) +
                                      0.5 * std::exp(-std::pow((hz - 1800) / 500, 2));
                // Integrated phase for the varying fundamental, no discontinuous
                // multiplication of instantaneous frequency by absolute time.
                const double cycles = 130 * t - 20 / (tau * 1.3) * std::cos(tau * 1.3 * t);
                voice += weight * std::sin(tau * harmonic * cycles) / harmonic;
            }
            rng = 1664525u * rng + 1013904223u;
            const double noise = double(rng) / 4294967296.0 - 0.5;
            const double fricative = std::fmod(t, 0.96) > 0.64 ? 0.35 : 0.04;
            const double fade = std::min({1., t / 0.02, (16. - t) / 0.02});
            dry[i] = float((0.8 * voice + fricative * noise) * envelope * fade);
            const int scene = int(i / (4 * rate));
            double hz = 3600, q = 0.7071067811865476;
            auto smode = SvfType::LowPass;
            auto bmode = BiquadType::LowPass;
            if (scene == 1) {
                hz = 300 * std::pow(8000. / 300, 0.5 + 0.5 * std::sin(tau * 0.5 * t));
                q = 0.8 + 3.2 * (0.5 + 0.5 * std::sin(tau * 0.37 * t));
            } else if (scene == 2) {
                hz = 300 * std::pow(20., 0.5 + 0.5 * std::sin(tau * 80 * t));
                q = 0.7 + 7.3 * (0.5 + 0.5 * std::sin(tau * 61 * t));
                smode = SvfType::BandPass; bmode = BiquadType::BandPass;
            } else if (scene == 3) {
                hz = 600 * std::pow(8000. / 600, 0.5 + 0.5 * std::sin(tau * 2000 * t));
                q = 0.7 + 3.3 * (0.5 + 0.5 * std::sin(tau * 1703 * t));
            }
            svf_events[i] = design_svf(rate, {smode, hz, q});
            bq_events[i] = design_biquad(rate, {bmode, hz, q});
        }
        biquad = svf = dry;
        // First scene: same 300 Hz high-pass + 3.6 kHz low-pass radio band.
        Svf high_svf(design_svf(rate, {SvfType::HighPass, 300}));
        Biquad high_bq(design_biquad(rate, {BiquadType::HighPass, 300}));
        high_svf.process(svf.data(), svf.data(), 4 * rate);
        high_bq.process(biquad.data(), biquad.data(), 4 * rate);
        Svf sf;
        Biquad bf;
        for (std::size_t i = 0; i < frames; ++i) {
            if (i % (4 * rate) == 0) bf.reset();
            bf.set_coefficients(bq_events[i]);
            bf.process(biquad.data() + i, biquad.data() + i, 1);
        }
        for (std::size_t scene = 0; scene < 4; ++scene) {
            sf.reset();
            const auto at = scene * 4 * rate;
            sf.process_modulated(svf.data() + at, svf.data() + at, svf_events.data() + at, 4 * rate);
        }
        double peak = 0, difference_energy = 0;
        for (std::size_t i = 0; i < frames; ++i) {
            for (double value : {double(dry[i]), double(biquad[i]), double(svf[i])}) {
                if (!std::isfinite(value)) throw std::runtime_error("nonfinite audition output");
                peak = std::max(peak, std::abs(value));
            }
            const double delta = double(biquad[i]) - svf[i];
            difference_energy += delta * delta;
            stereo[2 * i] = biquad[i]; stereo[2 * i + 1] = svf[i];
        }
        const double gain = peak > 0 ? std::min(1., 0.8 / peak) : 1.;
        wav(dry_path, dry, 1, gain); wav(bq_path, biquad, 1, gain);
        wav(svf_path, svf, 1, gain); wav(stereo_path, stereo, 2, gain);
        std::ofstream guide(directory / "README.txt");
        guide << "Original synthetic harmonic/fricative bursts, 16 seconds, 48 kHz PCM16.\n"
                 "Use a player of your choice at a comfortable level. Compare 02 and 03\n"
                 "from the same timestamp; 01 is dry. 04 has biquad left and SVF right\n"
                 "for soloing channels, not a blind listening test.\n"
                 "0-4 s: static HP 300 Hz + LP 3.6 kHz, Q=1/sqrt(2).\n"
                 "4-8 s: slow 0.5 Hz LP sweep, 300-8000 Hz, Q 0.8-4.\n"
                 "8-12 s: 80 Hz BP modulation, 300-6000 Hz, Q 0.7-8; unity center BP.\n"
                 "12-16 s: 2 kHz LP modulation, 600-8000 Hz, Q 0.7-4.\n"
                 "Both engines use identical cutoff/Q and per-sample prepared updates.\n"
                 "Every scene resets state; boundaries can click. No coefficient ramps,\n"
                 "extra smoothing, limiter, saturation, speech recording, IR or plugin.\n"
                 "One gain shared by all files preserves relative levels; no individual\n"
                 "normalization. This is material for Ian to judge, not an engine choice.\n"
              << "raw overall peak=" << peak << "; shared gain=" << gain
              << "; raw A/B difference RMS=" << std::sqrt(difference_energy / frames) << '\n';
        std::cout << "wrote audition to " << directory << "; raw peak=" << peak
                  << "; shared gain=" << gain << "; raw difference RMS="
                  << std::sqrt(difference_energy / frames) << '\n';
    } catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
