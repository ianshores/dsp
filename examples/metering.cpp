// SPDX-License-Identifier: MIT
// Minimal synthetic composition; no audio hardware, host adapter or plugin UI.
#include "dsp/metering.hpp"
#include <iostream>
#include <iomanip>

int main(int argc,char**argv){
    try {
    using namespace shores_dsp;
    const double fs=argc>1?std::stod(argv[1]):48000;
    constexpr double pi=3.14159265358979323846;
    std::array<PeakDetector,2> peaks;
    std::array<EnvelopeFollower,2> envelopes;
    std::array<PeakHold,2> holds;
    std::array<MeanSquareAccumulator,2> rms;
    std::array<KWeightingFilter,2> k;
    std::array<TruePeakOversampler,2> tp;
    ChannelWeighting channels{ChannelRole::Left,ChannelRole::Right};
    SlidingWindow momentary(fs,.4),short_term(fs,3);LoudnessGate integrated(100,fs);
    for(std::size_t c=0;c<2;++c){k[c].configure(fs);tp[c].configure(fs);envelopes[c].configure(fs,.01,.3);holds[c]=PeakHold(fs,1,20);}
    const auto frames=std::size_t(std::llround(fs*5));
    // All constructors/allocations precede this simulated callback loop. In a
    // host, each frame's two inputs come from caller buffers rather than sin().
    for(std::size_t at=0;at<frames;at+=256){
        const auto n=std::min<std::size_t>(256,frames-at);
        for(std::size_t i=0;i<n;++i){
            const float x=float(std::pow(10.,-23./20)*std::sin(2*pi*1000*double(at+i)/fs));
            std::array<double,2> e{};
            for(std::size_t c=0;c<2;++c){
                peaks[c].process(x);envelopes[c].process(x);holds[c].process(std::abs(x));rms[c].process(x);tp[c].process(x);
                const double y=k[c].process(x);e[c]=y*y;
            }
            const double energy=channels.process(e.data());momentary.process(energy);short_term.process(energy);integrated.process(energy);
        }
    }
    // End-of-stream and analysis/UI-boundary work, outside callback:
    for(auto &f:tp)f.flush();
    const auto result=integrated.result();
    std::cout<<std::setprecision(12)<<"synthetic stereo: rate="<<fs<<" TP_factor="<<tp[0].factor()<<" peak="<<peaks[0].value()<<" envelope="<<envelopes[0].value()
             <<" hold="<<holds[0].value()<<" RMS="<<rms[0].rms()<<" true_peak_dBTP="<<tp[0].dbtp()
             <<" M="<<energy_to_lufs(momentary.mean())<<" S="<<energy_to_lufs(short_term.mean())
             <<" I="<<result.lufs<<" LUFS; history_complete="<<result.complete<<'\n';
    return result.complete?0:1;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
