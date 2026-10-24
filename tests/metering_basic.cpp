// SPDX-License-Identifier: MIT
#include "dsp/metering/basic.hpp"
#include <array>
#include <iostream>
#include <random>
#include <type_traits>

using namespace shores_dsp;
static void require(bool value, const char *why) { if (!value) throw std::runtime_error(why); }
static void near(double a, double b, double tolerance, const char *why) {
    require(std::isfinite(a) && std::isfinite(b) && std::abs(a-b) <= tolerance, why);
}
template<class F> static void rejects(F f) {
    bool caught = false;
    try { f(); } catch (const std::invalid_argument &) { caught = true; }
    require(caught, "bad setup accepted");
}
int main() {
    try {
        static_assert(std::is_same_v<BiquadFilter, Biquad>, "must reuse existing biquad");
        PeakDetector peak;
        const float values[]{0, -0.25f, 0.125f, -2.f, 0.5f};
        near(peak.process(values, 5), 2, 0, "absolute sample peak");
        peak.process(nullptr, 0); peak.reset(); near(peak.value(), 0, 0, "peak reset");
        EnvelopeFollower env(1000, 0.01, 0.1);
        for (int n=0; n<10; ++n) env.process(-1);
        near(env.value(), 1-std::exp(-1.), 1.e-14, "attack 1/e time constant");
        const auto top = env.value();
        for (int n=0; n<100; ++n) env.process(0);
        near(env.value(), top*std::exp(-1.), 1.e-14, "release 1/e time constant");
        env.configure(48000, 0, 0); near(env.process(-0.75f), 0.75, 0, "instant attack");
        near(env.process(0), 0, 0, "instant release"); env.reset();
        EnvelopeFollower block_env, sample_env;
        std::array<float, 3073> audio{}; std::array<double, 3073> blocks{}, singles{};
        std::mt19937 rng(2718);
        for (auto &x:audio) x=float(double(rng())/4294967296.-0.5);
        for (std::size_t n=0;n<audio.size();++n) singles[n]=sample_env.process(audio[n]);
        std::size_t at=0;
        for (std::size_t block : {1u,17u,256u,7u,1024u,63u,1705u}) {
            block=std::min(block,audio.size()-at); block_env.process(nullptr,nullptr,0);
            block_env.process(audio.data()+at,blocks.data()+at,block); at+=block;
        }
        require(at==audio.size() && blocks==singles, "envelope variable block invariance");
        PeakHold hold(1000,0.003,20);
        hold.process(1);
        for(int i=0;i<3;++i) near(hold.process(0),1,0,"hold duration");
        for(int i=0;i<1000;++i) hold.process(0);
        near(hold.value(),0.1,1.e-14,"20 dB/second decay");
        hold.process(2); near(hold.value(),2,0,"peak refresh"); hold.reset();
        PeakHold instant(1000,0,20); instant.process(1);
        near(instant.process(0),std::pow(10.,-20./20000),1.e-15,"zero hold");
        MeanSquareAccumulator energy;
        long double sum=0;
        for(float x:audio) { energy.process(x); sum+=static_cast<long double>(x)*x; }
        near(energy.mean_square(),double(sum/audio.size()),2.e-16,"independent energy sum");
        near(energy.rms(),std::sqrt(double(sum/audio.size())),1.e-15,"RMS");
        energy.process(nullptr,0); energy.reset(); require(energy.count()==0,"energy reset");
        MeanSquareAccumulator limited(2);
        limited.process(1); limited.process(2); require(!limited.process(3),"energy capacity accepted");
        require(limited.overflowed() && limited.count()==2,"overflow must freeze prefix");
        near(limited.mean_square(),2.5,0,"prefix mean"); limited.reset();
        require(!limited.overflowed() && limited.valid(),"reset overflow");
        SlidingWindow window(1000,0.017);
        std::vector<double> history;
        for(int i=0;i<20000;++i) {
            const double x=double(rng())/4294967296.; history.push_back(x*x);
            window.process(x*x); long double reference=0;
            const auto first=history.size()>17?history.size()-17:0;
            for(std::size_t j=first;j<history.size();++j) reference+=history[j];
            near(window.mean(),double(reference/(history.size()-first)),2.e-14,"independent moving sum");
        }
        require(window.ready() && window.count()==17,"window readiness/capacity");
        window.reset(); require(!window.ready() && window.count()==0,"window reset");
        near(window.process(9),9,0,"window partial startup");
        window.reset(); for(int i=0;i<40;++i) window.process(0);
        near(window.mean(),0,0,"window silence after reset");
        window.process(1.e90);
        for(int i=0;i<17;++i) window.process(1);
        near(window.mean(),1,1.e-14,"large peak must leave moving window");
        const auto nan=std::numeric_limits<float>::quiet_NaN();
        const auto inf=std::numeric_limits<float>::infinity();
        for(float bad:{nan,inf,-inf}) {
            peak.process(bad); env.process(bad); hold.process(bad); energy.process(bad); window.process(bad);
        }
        require(!peak.valid()&&!env.valid()&&!hold.valid()&&!energy.valid()&&!window.valid(),"sticky invalid flags");
        window.reset(); window.process(-1); require(!window.valid(),"negative energy accepted");
        window.reset(); window.process(1.e101); require(!window.valid(),"excess energy accepted");
        peak.reset(); env.reset(); hold.reset(); energy.reset(); window.reset();
        require(peak.valid()&&env.valid()&&hold.valid()&&energy.valid()&&window.valid(),"reset status");
        for(double bad:{0.,999.,768001.,double(nan)}) rejects([&]{EnvelopeFollower x(bad);});
        rejects([]{EnvelopeFollower x(48000,-1,0);}); rejects([]{PeakHold x(48000,61);});
        rejects([]{PeakHold x(48000,1,-1);}); rejects([]{SlidingWindow x(48000,0);});
        rejects([]{SlidingWindow x(768000,10);}); rejects([]{MeanSquareAccumulator x(0);});
        std::cout<<"PASS: basic peak, attack/release, hold/decay, independent mean-square/window, "
                     "variable blocks, silence/reset, invalid inputs, setup/capacity limits\n";
    } catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 1; }
}
