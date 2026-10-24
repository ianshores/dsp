// SPDX-License-Identifier: MIT
#include "dsp/metering.hpp"
#include <cstdlib>
#include <iostream>
#include <new>

// Keep instrumentation boundaries without imposing a vendor attribute on
// compilers that do not recognize it. Only GCC/Clang builds are validated here.
#if defined(__GNUC__) || defined(__clang__)
#define DSP_METERING_TEST_NOINLINE [[gnu::noinline]]
#elif defined(_MSC_VER)
#define DSP_METERING_TEST_NOINLINE __declspec(noinline)
#else
#define DSP_METERING_TEST_NOINLINE
#endif

static bool observing=false;
static std::size_t allocations=0;
DSP_METERING_TEST_NOINLINE void *operator new(std::size_t n) {
    if(observing)++allocations;
    if(auto*p=std::malloc(n?n:1))return p;
    throw std::bad_alloc();
}
DSP_METERING_TEST_NOINLINE void *operator new[](std::size_t n){return ::operator new(n);}
DSP_METERING_TEST_NOINLINE void operator delete(void*p) noexcept{std::free(p);}
DSP_METERING_TEST_NOINLINE void operator delete[](void*p) noexcept{std::free(p);}
DSP_METERING_TEST_NOINLINE void operator delete(void*p,std::size_t) noexcept{std::free(p);}
DSP_METERING_TEST_NOINLINE void operator delete[](void*p,std::size_t) noexcept{std::free(p);}

#undef DSP_METERING_TEST_NOINLINE

int main(){
    using namespace shores_dsp;
    PeakDetector peak;EnvelopeFollower env;PeakHold hold;MeanSquareAccumulator square;
    SlidingWindow window;KWeightingFilter k;LoudnessGate gate;TruePeakOversampler tp;
    ChannelWeighting channels{ChannelRole::Left,ChannelRole::Right};BiquadFilter biquad;
    KWeightingFilter fractional_k(44117.3);LoudnessGate fractional_gate(100,44117.3);
    TruePeakOversampler fractional_tp(44117.3),maximum_factor_tp(1000);SlidingWindow fractional_window(44117.3,3);
    double checksum=0;
    observing=true;
    for(int i=0;i<96000;++i){
        const float x=float((i%997-499)/1000.);
        checksum+=peak.process(x)+env.process(x)+hold.process(std::abs(x));square.process(x);
        const double y=k.process(x),energies[]{y*y,y*y};
        const double e=channels.process(energies);gate.process(e);checksum+=window.process(e)+tp.process(x);
        const double fy=fractional_k.process(x);fractional_gate.process(fy*fy);
        checksum+=fractional_tp.process(x)+fractional_window.process(fy*fy);
        if(i<2048)checksum+=maximum_factor_tp.process(x);
        float b=x;biquad.process(&b,&b,1);checksum+=b;
    }
    peak.reset();env.reset();hold.reset();square.reset();window.reset();k.reset();gate.reset();tp.reset();channels.reset();biquad.reset();
    fractional_k.reset();fractional_gate.reset();fractional_tp.reset();maximum_factor_tp.reset();fractional_window.reset();
    observing=false;
    if(allocations!=0||!std::isfinite(checksum)){std::cerr<<"callback/reset allocated or nonfinite\n";return 1;}
    std::cout<<"PASS: zero heap allocations during all primitive processing/reset paths (96000 frames); setup excluded\n";
}
