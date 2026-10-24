// SPDX-License-Identifier: MIT
#include "dsp/metering.hpp"
#include <cstdlib>
#include <iostream>
#include <new>

static bool observing=false;
static std::size_t allocations=0;
void *operator new(std::size_t n) {
    if(observing)++allocations;
    if(auto*p=std::malloc(n?n:1))return p;
    throw std::bad_alloc();
}
void *operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void*p) noexcept{std::free(p);}
void operator delete[](void*p) noexcept{std::free(p);}
void operator delete(void*p,std::size_t) noexcept{std::free(p);}
void operator delete[](void*p,std::size_t) noexcept{std::free(p);}

int main(){
    using namespace shores_dsp;
    PeakDetector peak;EnvelopeFollower env;PeakHold hold;MeanSquareAccumulator square;
    SlidingWindow window;KWeightingFilter k;LoudnessGate gate;TruePeakOversampler tp;
    ChannelWeighting channels{ChannelRole::Left,ChannelRole::Right};BiquadFilter biquad;
    double checksum=0;
    observing=true;
    for(int i=0;i<96000;++i){
        const float x=float((i%997-499)/1000.);
        checksum+=peak.process(x)+env.process(x)+hold.process(std::abs(x));square.process(x);
        const double y=k.process(x),energies[]{y*y,y*y};
        const double e=channels.process(energies);gate.process(e);checksum+=window.process(e)+tp.process(x);
        float b=x;biquad.process(&b,&b,1);checksum+=b;
    }
    peak.reset();env.reset();hold.reset();square.reset();window.reset();k.reset();gate.reset();tp.reset();channels.reset();biquad.reset();
    observing=false;
    if(allocations!=0||!std::isfinite(checksum)){std::cerr<<"callback/reset allocated or nonfinite\n";return 1;}
    std::cout<<"PASS: zero heap allocations during all primitive processing/reset paths (96000 frames); setup excluded\n";
}
