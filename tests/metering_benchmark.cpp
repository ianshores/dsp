// SPDX-License-Identifier: MIT
#include "dsp/metering.hpp"
#include <chrono>
#include <cstdio>
#include <random>

using namespace shores_dsp;
constexpr std::size_t frames=48000*20,block=256;
constexpr int trials=9;
template<class T> static void observe(const T&x){
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(&x) : "memory");
#else
    (void)x;
#endif
}
static double ns(std::chrono::steady_clock::time_point start){
    return std::chrono::duration<double,std::nano>(std::chrono::steady_clock::now()-start).count();
}
static void row(const char*name,const char*mode,int channels,int trial,const char*unit,double cost,double checksum){
    if(!std::isfinite(checksum))throw std::runtime_error("nonfinite benchmark checksum");
    std::printf("%s,%s,%d,%d,%s,%.6f,%.12g\n",name,mode,channels,trial,unit,cost,checksum);
}
template<class Factory,class Process> static void run(const char*name,int channels,Factory factory,Process process,
        const std::vector<float>&input,const std::vector<double>&energies){
    for(int trial=0;trial<trials;++trial){
        const auto start_setup=std::chrono::steady_clock::now();auto state=factory();observe(state);
        const double setup=ns(start_setup);double checksum=0;
        const auto start=std::chrono::steady_clock::now();
        for(std::size_t at=0;at<frames;at+=block)
            checksum+=process(state,input.data()+at,energies.data()+at,std::min(block,frames-at));
        observe(state);const double processing=ns(start)/frames;
        row(name,"setup",channels,trial,"instance",setup,checksum);
        row(name,"process",channels,trial,"frame",processing,checksum);
    }
}

// Benchmark-only composition, not an additional library module. Includes every
// requested primitive; Q/IIR is the existing Biquad through the K cascade.
struct Composition {
    std::array<PeakDetector,6> peak;
    std::array<EnvelopeFollower,6> env;
    std::array<PeakHold,6> hold;
    std::array<MeanSquareAccumulator,6> square;
    std::array<KWeightingFilter,6> k;
    std::array<TruePeakOversampler,6> tp;
    SlidingWindow m{48000,.4},s{48000,3};LoudnessGate gate;
    ChannelWeighting channels;
    explicit Composition(int n):channels(n==2?ChannelWeighting{ChannelRole::Left,ChannelRole::Right}:
        ChannelWeighting{ChannelRole::Left,ChannelRole::Right,ChannelRole::Centre,ChannelRole::Lfe,ChannelRole::LeftSurround,ChannelRole::RightSurround}){}
    double process(const float*x,std::size_t n){
        for(std::size_t i=0;i<n;++i){
            std::array<double,6> e{};
            for(std::size_t c=0;c<channels.channels();++c){
                const float v=c%2?x[i]:-x[i];peak[c].process(v);env[c].process(v);hold[c].process(std::abs(v));square[c].process(v);tp[c].process(v);
                const double y=k[c].process(v);e[c]=y*y;
            }
            const double energy=channels.process(e.data());m.process(energy);s.process(energy);gate.process(energy);
        }
        return m.mean()+s.mean()+tp[0].value();
    }
};

int main(){
    std::vector<float> input(frames);std::vector<double> e(frames);
    std::mt19937 rng(41);for(std::size_t i=0;i<frames;++i){input[i]=float(double(rng())/4294967296.-.5);e[i]=double(input[i])*input[i];}
    std::puts("primitive,mode,channels,trial,unit,ns_per_unit,checksum");
    run("peak",1,[]{return PeakDetector{};},[](auto&f,const float*x,const double*,std::size_t n){return f.process(x,n);},input,e);
    run("envelope",1,[]{return EnvelopeFollower{};},[](auto&f,const float*x,const double*,std::size_t n){for(std::size_t i=0;i<n;++i)f.process(x[i]);return f.value();},input,e);
    run("hold",1,[]{return PeakHold{};},[](auto&f,const float*x,const double*,std::size_t n){for(std::size_t i=0;i<n;++i)f.process(std::abs(x[i]));return f.value();},input,e);
    run("mean_square",1,[]{return MeanSquareAccumulator{};},[](auto&f,const float*x,const double*,std::size_t n){f.process(x,n);return f.mean_square();},input,e);
    run("window_400ms",1,[]{return SlidingWindow{};},[](auto&f,const float*,const double*y,std::size_t n){for(std::size_t i=0;i<n;++i)f.process(y[i]);return f.mean();},input,e);
    const auto bq=design_biquad(48000,{BiquadType::HighPass,300});
    run("biquad",1,[&]{return BiquadFilter{bq};},[](auto&f,const float*x,const double*,std::size_t n){std::array<float,256>y{};f.process(x,y.data(),n);return double(y[n-1]);},input,e);
    run("k_weighting",1,[]{return KWeightingFilter{};},[](auto&f,const float*x,const double*,std::size_t n){std::array<float,256>y{};f.process(x,y.data(),n);return double(y[n-1]);},input,e);
    run("channel_weight_5_1",6,[]{return ChannelWeighting{ChannelRole::Left,ChannelRole::Right,ChannelRole::Centre,ChannelRole::Lfe,ChannelRole::LeftSurround,ChannelRole::RightSurround};},
        [](auto&f,const float*,const double*y,std::size_t n){double v=0;for(std::size_t i=0;i<n;++i){const double a[]{y[i],y[i],y[i],y[i],y[i],y[i]};v=f.process(a);}return v;},input,e);
    run("gate",1,[]{return LoudnessGate{};},[](auto&f,const float*,const double*y,std::size_t n){f.process(y,n);return double(f.stored_blocks());},input,e);
    run("true_peak_4x",1,[]{return TruePeakOversampler{};},[](auto&f,const float*x,const double*,std::size_t n){return f.process(x,n);},input,e);
    for(int channels:{2,6})run("composition",channels,[&]{return Composition{channels};},
        [](auto&f,const float*x,const double*,std::size_t n){return f.process(x,n);},input,e);
    // Analysis queries are measured separately from the callback's O(1) gate
    // insertion. 20 seconds produces 197 overlapping gating blocks.
    LoudnessGate query;query.process(e.data(),e.size());
    for(int trial=0;trial<trials;++trial){
        const auto start=std::chrono::steady_clock::now();const auto r=query.result();observe(r);
        row("gate_197_blocks","query",1,trial,"query",ns(start),r.lufs);
    }
    EnvelopeFollower env;
    for(int trial=0;trial<trials;++trial){
        const auto start=std::chrono::steady_clock::now();
        for(std::size_t i=0;i<frames;++i){env.configure(48000,.001+e[i],.1+e[i]);observe(env);}
        row("envelope_ballistics","update",1,trial,"event",ns(start)/frames,env.value());
    }
}
