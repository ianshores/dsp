// SPDX-License-Identifier: MIT
#include "dsp/metering.hpp"
#include <complex>
#include <iostream>
#include <random>

using namespace shores_dsp;
constexpr long double pi=3.1415926535897932384626433832795L;
static void require(bool ok,const char *s) { if(!ok) throw std::runtime_error(s); }
static void near(double a,double b,double t,const char*s) {
    require(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=t,s);
}
template<class F> static void rejects(F f) {
    bool caught=false; try {f();}catch(const std::invalid_argument&){caught=true;}
    require(caught,"bad setup accepted");
}
struct DF1 {
    long double b0,b1,b2,a1,a2,x1=0,x2=0,y1=0,y2=0;
    long double process(long double x) {
        const auto y=b0*x+b1*x1+b2*x2-a1*y1-a2*y2;
        x2=x1;x1=x;y2=y1;y1=y;return y;
    }
};
// Brute-force independent gate scans raw per-frame energies for each 19200
// frame block; it does not reuse production's four-hop decomposition/sums.
static double reference_gate(const std::vector<double>& e) {
    std::vector<long double> blocks;
    for(std::size_t start=0;start+19200<=e.size();start+=4800) {
        long double sum=0;for(std::size_t i=start;i<start+19200;++i)sum+=e[i];
        blocks.push_back(sum/19200);
    }
    long double sum=0; std::size_t n=0;
    const auto absolute=std::pow(10.L,(-70.L+0.691L)/10);
    for(auto x:blocks)if(x>absolute){sum+=x;++n;}
    if(!n)return -std::numeric_limits<double>::infinity();
    const auto threshold=std::max(absolute,sum/n/10);
    sum=0;n=0;for(auto x:blocks)if(x>threshold){sum+=x;++n;}
    return n?double(-0.691L+10*std::log10(sum/n)):-std::numeric_limits<double>::infinity();
}

static void k_weighting() {
    KWeightingFilter filter;
    DF1 shelf{1.53512485958697L,-2.69169618940638L,1.19839281085285L,-1.69065929318241L,0.73248077421585L};
    DF1 high{1,-2,1,-1.99004745483398L,0.99007225036621L};
    std::mt19937 rng(17);double worst=0;
    std::vector<float> input(100000),output(input.size());
    for(auto &x:input)x=float(double(rng())/4294967296.-0.5);
    filter.process(input.data(),output.data(),input.size());
    for(std::size_t i=0;i<input.size();++i) {
        const double y=double(high.process(shelf.process(input[i])));
        worst=std::max(worst,std::abs(y-output[i]));
        near(output[i],y,5.e-6,"K weighting independent DF1");
    }
    filter.reset(); auto in_place=input;
    for(std::size_t at=0;at<input.size();) {
        const auto n=std::min<std::size_t>(at%257+1,input.size()-at);
        filter.process(nullptr,nullptr,0);filter.process(in_place.data()+at,in_place.data()+at,n);at+=n;
    }
    require(in_place==output,"K weighting variable/in-place/reset");
    filter.process(std::numeric_limits<float>::infinity());require(!filter.valid(),"invalid K input");
    require(std::isfinite(filter.process(1)),"nonfinite K state poisoned");filter.reset();
    require(std::isfinite(filter.process(std::numeric_limits<float>::max()))&&!filter.valid(),"K intermediate overflow sanitation");
    filter.reset();
    near(filter.process(0),0,0,"K silence/reset");
    require(filter.valid(),"K reset status");
    std::cout<<"K reference maximum absolute error="<<worst<<'\n';
}

static void gate_and_layout() {
    using R=ChannelRole;
    ChannelWeighting stereo{R::Left,R::Right}, surround{R::Left,R::Right,R::Centre,R::Lfe,R::LeftSurround,R::RightSurround};
    double energy[]{1,2,3,1.e90,4,5};
    near(stereo.process(energy),3,0,"stereo sum");
    near(surround.process(energy),6+1.41*9,1.e-14,"surround sum and LFE exclusion");
    ChannelWeighting mono{R::Centre}, lfe{R::Lfe};
    near(mono.process(energy),1,0,"mono no dual-mono doubling");near(lfe.process(energy),0,0,"LFE zero");
    ChannelWeighting permuted{R::RightSurround,R::Left};double e2[]{4,5};
    near(permuted.process(e2),1.41*4+5,1.e-14,"explicit channel order");
    double invalid_energy[]{std::numeric_limits<double>::quiet_NaN(),1};
    near(stereo.process(invalid_energy),1,0,"channel invalid energy sanitation");
    require(!stereo.valid(),"channel sticky invalid status");stereo.reset();require(stereo.valid(),"channel status reset");
    rejects([]{ChannelWeighting x{R::Left,R::Left};});rejects([]{ChannelWeighting x{}; /* default mono is valid */
        ChannelWeighting empty(std::initializer_list<R>{});});
    rejects([]{ChannelWeighting x{static_cast<R>(99)};});
    std::vector<double> e(48000*8+321);
    for(std::size_t i=0;i<e.size();++i) {
        const int s=int(i/48000);
        e[i]=s==0?0:s==1?1.e-10:s==2?1.e-6:s<6?0.01:0.00001;
        e[i]*=1+0.2*std::sin(double(i)*0.17);
    }
    LoudnessGate full(100),chunks(100);
    full.process(e.data(),e.size());
    for(std::size_t at=0;at<e.size();) {
        const auto n=std::min<std::size_t>(1+at%997,e.size()-at);
        chunks.process(nullptr,0);chunks.process(e.data()+at,n);at+=n;
    }
    const auto result=full.result();
    near(result.lufs,reference_gate(e),1.e-12,"gate independent full-block integration");
    near(chunks.result().lufs,result.lufs,0,"gate block invariance");
    require(result.complete&&result.stored_blocks==(e.size()-19200)/4800+1,"gate block/hop/tail semantics");
    require(result.included_blocks<result.absolute_blocks,"relative gate not exercised");
    full.reset();for(int i=0;i<19199;++i)full.process(1);
    require(full.stored_blocks()==0,"partial 400 ms block used");
    full.process(1);require(full.stored_blocks()==1,"first block off by one");
    for(int i=0;i<4799;++i)full.process(1);
    require(full.stored_blocks()==1,"hop early");full.process(1);require(full.stored_blocks()==2,"hop late");
    LoudnessGate capacity(1);for(int i=0;i<24000;++i)capacity.process(1);
    require(capacity.overflowed()&&!capacity.result().complete&&capacity.stored_blocks()==1,"gate overflow freeze/status");
    capacity.reset();require(!capacity.overflowed()&&capacity.result().complete,"gate reset overflow");
    for(int i=0;i<19200;++i)capacity.process(0);
    require(std::isinf(capacity.result().lufs)&&capacity.result().lufs<0,"gated silence");
    capacity.process(std::numeric_limits<double>::quiet_NaN());require(!capacity.result().complete,"gate invalid flag");
    for(double bad:{44100.,96000.,0.,std::numeric_limits<double>::quiet_NaN()}) {
        rejects([&]{KWeightingFilter x(bad);});rejects([&]{LoudnessGate x(10,bad);});
    }
    rejects([]{LoudnessGate x(0);});rejects([]{LoudnessGate x(1000001);});
    std::cout<<"gate independent result="<<result.lufs<<" LUFS; blocks="<<result.stored_blocks<<'\n';
}

// Independently synthesized Tech 3341 Table 1 cases 1-6 (not downloaded EBU
// fixtures). Unrounded analytic tones; 1-5 stereo, 6 explicit five-channel sum.
static void synthetic_ebu_loudness() {
    struct Segment { double seconds,db; };
    const std::vector<std::vector<Segment>> cases{
        {{20,-23}},{{20,-33}},{{10,-36},{60,-23},{10,-36}},
        {{10,-72},{10,-36},{60,-23},{10,-36},{10,-72}},
        {{20,-26},{20.1,-20},{20,-26}},{{20,-23}}};
    for(std::size_t test=0;test<cases.size();++test) {
        std::array<KWeightingFilter,5> filters;
        ChannelWeighting weights=test==5?ChannelWeighting{ChannelRole::Left,ChannelRole::Right,ChannelRole::Centre,ChannelRole::LeftSurround,ChannelRole::RightSurround}:
                                             ChannelWeighting{ChannelRole::Left,ChannelRole::Right};
        LoudnessGate gate;SlidingWindow momentary(48000,.4),short_term(48000,3);
        std::size_t at=0;
        for(auto segment:cases[test]) {
            const auto count=std::size_t(std::llround(segment.seconds*48000));
            for(std::size_t n=0;n<count;++n,++at) {
                std::array<double,5> energies{};
                const auto tone=std::sin(2*double(pi)*1000*double(at)/48000);
                for(std::size_t c=0;c<weights.channels();++c) {
                    const double db=test==5?(c<2?-28:c==2?-24:-30):segment.db;
                    const float x=float(std::pow(10.,db/20)*tone);
                    const double y=filters[c].process(x);energies[c]=y*y;
                }
                const auto e=weights.process(energies.data());gate.process(e);momentary.process(e);short_term.process(e);
            }
        }
        const auto expected=test==1?-33.:-23.; const auto got=gate.result();
        require(got.complete,"unexpected synthetic capacity/invalidity");
        near(got.lufs,expected,.1,"synthetic EBU integrated tolerance");
        if(test<2) {
            near(energy_to_lufs(momentary.mean()),expected,.1,"synthetic EBU M");
            near(energy_to_lufs(short_term.mean()),expected,.1,"synthetic EBU S");
        }
        std::cout<<"synthetic EBU case "<<test+1<<" integrated="<<got.lufs<<" expected="<<expected<<" tolerance=0.1 LU\n";
    }
}

int main() {
    try { k_weighting();gate_and_layout();synthetic_ebu_loudness();
        std::cout<<"PASS: exact 48k K cascade, independent gate, explicit layouts/LFE, synthetic loudness cases; no blanket compliance claim\n";
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
