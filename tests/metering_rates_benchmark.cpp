// SPDX-License-Identifier: MIT
#include "dsp/metering.hpp"
#include <chrono>
#include <cstdio>
#include <random>
using namespace shores_dsp;
using Clock=std::chrono::steady_clock;
template<class T>static void observe(const T&v){
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(&v) : "memory");
#else
    (void)v;
#endif
}
static double elapsed(Clock::time_point t){return std::chrono::duration<double,std::nano>(Clock::now()-t).count();}
static void row(const char*name,const char*mode,double fs,std::size_t frames,int trial,double cost,double checksum){
    if(!std::isfinite(checksum))throw std::runtime_error("nonfinite benchmark checksum");
    std::printf("%s,%s,%.12g,1,256,%zu,%d,%.6f,%.12g\n",name,mode,fs,frames,trial,cost,checksum);
}
struct Window3s {
    SlidingWindow window;
    explicit Window3s(double fs):window(fs,3){}
    void configure(double fs){window.configure(fs,3);}
};
template<class T>static void prepare(T&){}
static void prepare(Window3s&f){
    for(std::size_t i=0;i<f.window.capacity();++i)f.window.process(.01);
}
template<class Factory,class Process>static void run(const char*name,double fs,Factory factory,Process process,const std::vector<float>&x){
    for(int t=0;t<9;++t){
        auto start=Clock::now();auto f=factory();observe(f);const auto setup=elapsed(start);double sum=0;
        prepare(f); // steady full-window removal, excluded from setup/process timers
        start=Clock::now();
        for(std::size_t at=0;at<x.size();at+=256){const auto n=std::min<std::size_t>(256,x.size()-at);sum+=process(f,x.data()+at,n);}
        observe(f);const auto cost=elapsed(start)/x.size();row(name,"setup_ns_instance",fs,x.size(),t,setup,sum);row(name,"process_ns_frame",fs,x.size(),t,cost,sum);
        // configure is a reset/new measurement epoch, setup only. Vary by a
        // fractional Hz to include both design and rate-dependent storage.
        start=Clock::now();for(int i=0;i<32;++i){f.configure(fs==768000?fs-.25*(i%2):fs+.25*(i%2));observe(f);}
        row(name,"reconfigure_ns_event",fs,x.size(),t,elapsed(start)/32,sum);
    }
}
int main(){
    std::puts("primitive,mode,sample_rate,channels,block,frames,trial,ns_per_unit,checksum");
    for(double fs:{1000.,4000.,44100.,48000.,48000.25,88200.,96000.,176400.,192000.,768000.}){
        const auto n=std::max<std::size_t>(8192,std::size_t(std::min(384000.,fs*2)));
        std::vector<float>x(n);std::mt19937 rng(41);for(auto &v:x)v=float(double(rng())/4294967296.-.5);
        if(fs>=4000)run("k_weighting",fs,[&]{return KWeightingFilter(fs);},[](auto&f,const float*x,std::size_t n){float y[256];f.process(x,y,n);return double(y[n-1]);},x);
        run("true_peak",fs,[&]{return TruePeakOversampler(fs);},[](auto&f,const float*x,std::size_t n){return f.process(x,n);},x);
        run("window_3s",fs,[&]{return Window3s(fs);},[](auto&f,const float*x,std::size_t n){for(std::size_t i=0;i<n;++i)f.window.process(double(x[i])*x[i]);return f.window.mean();},x);
        run("gate",fs,[&]{return LoudnessGate(100,fs);},[](auto&f,const float*x,std::size_t n){for(std::size_t i=0;i<n;++i)f.process(double(x[i])*x[i]);return double(f.stored_blocks());},x);
    }
}
