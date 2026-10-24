// SPDX-License-Identifier: MIT
#include "dsp/metering/true_peak.hpp"
#include <iostream>
#include <random>

using namespace shores_dsp;
constexpr double pi=3.1415926535897932384626433832795;
static void require(bool ok,const char*s){if(!ok)throw std::runtime_error(s);}
static void near(double a,double b,double t,const char*s){require(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=t,s);}
// Independently transcribed full-rate impulse response, ITU BS.1770-5 Annex 2.
// Reference below explicitly inserts zeros and convolves in long double; no
// production history/phase loop or production coefficients are accessed.
static const std::array<long double,48> h{
 .001708984375L,-.0291748046875L,-.0189208984375L,-.00830078125L,
 .010986328125L,.029296875L,.0330810546875L,.014892578125L,
 -.0196533203125L,-.0517578125L,-.0582275390625L,-.026611328125L,
 .033203125L,.089111328125L,.1015625L,.047607421875L,
 -.0594482421875L,-.16650390625L,-.2003173828125L,-.102294921875L,
 .1373291015625L,.465087890625L,.77978515625L,.97216796875L,
 .97216796875L,.77978515625L,.465087890625L,.1373291015625L,
 -.102294921875L,-.2003173828125L,-.16650390625L,-.0594482421875L,
 .047607421875L,.1015625L,.089111328125L,.033203125L,
 -.026611328125L,-.0582275390625L,-.0517578125L,-.0196533203125L,
 .014892578125L,.0330810546875L,.029296875L,.010986328125L,
 -.00830078125L,-.0189208984375L,-.0291748046875L,.001708984375L};
static void convolution_reference() {
    std::mt19937 rng(33);std::vector<float>x(4096);
    for(auto &v:x)v=float(double(rng())/4294967296.-.5);
    x.back()=2; // forces end-of-stream tail to matter
    std::vector<long double> zero_stuffed(x.size()*4+44,0),reference(zero_stuffed.size(),0);
    for(std::size_t i=0;i<x.size();++i)zero_stuffed[4*i]=x[i];
    for(std::size_t n=0;n<reference.size();++n)
        for(std::size_t k=0;k<h.size()&&k<=n;++k)reference[n]+=h[k]*zero_stuffed[n-k];
    TruePeakOversampler tp,blocks;
    double worst=0,expected_peak=0;
    for(std::size_t n=0;n<x.size()+11;++n){
        long double expected=0;
        for(int phase=0;phase<4;++phase)expected=std::max(expected,std::abs(reference[4*n+std::size_t(phase)]));
        const double actual=tp.process(n<x.size()?x[n]:0);
        worst=std::max(worst,std::abs(actual-double(expected)));
        near(actual,double(expected),2.e-15,"independent zero-stuffed FIR convolution");
        expected_peak=std::max(expected_peak,double(expected));
    }
    std::size_t at=0;
    while(at<x.size()){
        const auto n=std::min<std::size_t>(1+at%257,x.size()-at);
        blocks.process(nullptr,0);blocks.process(x.data()+at,n);at+=n;
    }
    near(blocks.flush(),expected_peak,2.e-15,"flush/variable blocks");
    near(blocks.process(0),0,0,"flush did not drain delayed tail");
    tp.reset();require(tp.valid()&&tp.value()==0&&std::isinf(tp.dbtp()),"true peak reset/silence");
    for(float bad:{std::numeric_limits<float>::quiet_NaN(),std::numeric_limits<float>::infinity()})
        require(std::isfinite(tp.process(bad)),"nonfinite true peak state");
    require(!tp.valid(),"true peak invalid flag");tp.reset();
    std::cout<<"true-peak direct-convolution maximum error="<<worst<<'\n';
}
static double tone(double normalized_frequency,double phase,double amplitude) {
    TruePeakOversampler tp;
    for(int n=0;n<48000;++n) {
        const double fade=std::min({1.,double(n)/480,double(47999-n)/480});
        tp.process(float(amplitude*fade*std::sin(2*pi*normalized_frequency*n+phase)));
    }
    tp.flush();return tp.dbtp();
}
int main(){
    try{
        convolution_reference();
        const std::array<double,5> f{.25,.25,1./6,1./8,.25};
        const std::array<double,5> phase{0,45,60,67.5,45};
        for(int i=0;i<5;++i){
            const double got=tone(f[std::size_t(i)],phase[std::size_t(i)]*pi/180,i==4?1.41:.5);
            const double nominal=i==4?3.:-6.;
            require(got>=nominal-.4&&got<=nominal+.2,"synthetic EBU 15-19 tolerance");
            std::cout<<"synthetic EBU case "<<15+i<<" dBTP="<<got<<" expected="<<nominal<<" +0.2/-0.4\n";
        }
        // Quantify finite-band tone accuracy independently of the published FIR
        // recurrence. Includes frequencies outside the minimum EBU tone cases.
        double minimum=100,maximum=-100;
        for(double f0:{.01,.1,.125,1./6,.25,.35,.4,20000./48000})
            for(int p=0;p<8;++p){
                const double error=tone(f0,double(p)*pi/16,.5)-20*std::log10(.5);
                minimum=std::min(minimum,error);maximum=std::max(maximum,error);
            }
        require(minimum>-.8&&maximum<.3,"true peak tone accuracy regression bound");
        std::cout<<"tested 480 Hz..20 kHz tone/phase error range="<<minimum<<".."<<maximum<<" dB; not universal tolerance\n";
        bool rejected=false;try{TruePeakOversampler bad(96000);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"unsupported true peak rate accepted");
        std::cout<<"PASS: bandlimited 4x FIR, independent convolution, synthetic EBU 15-19, flush, silence/reset/invalid and phase sweep; official transient fixtures pending\n";
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
