// SPDX-License-Identifier: MIT
#include "dsp/metering/rate.hpp" // must compile without transitive includes
#include "dsp/metering.hpp"
#include <complex>
#include <iostream>
#include <random>

using namespace shores_dsp;
constexpr long double pi=3.1415926535897932384626433832795L;
static void require(bool v,const char*s){if(!v)throw std::runtime_error(s);}
static void near(double a,double b,double t,const char*s){if(!(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=t)){std::cerr<<s<<": "<<a<<" vs "<<b<<" tolerance "<<t<<'\n';throw std::runtime_error(s);}}
template<class F>static void rejects(F f){bool caught=false;try{f();}catch(const std::invalid_argument&){caught=true;}require(caught,"invalid setup accepted");}
static const std::vector<double> rates{4000,4000.25,4001.25,8000,22050,44100,48000,88200,96000,176400,192000,44117.3,48000.25,96000.123,349525.5,768000};

// Independent physical-Hz parameterization: De Man shelf constants and RLB
// frequency/Q. RLB numerator scaled to the published 48 kHz Nyquist gain,
// unlike the raw [1,-2,1] numerator convention at all sample rates.
static std::complex<long double> reference_k(double fs,double hz){
    using C=std::complex<long double>;
    const C z=std::exp(C(0,-2*pi*hz/fs));
    const auto section=[&](long double f,long double q,long double vh,long double vb,long double dc){
        const auto k=std::tan(pi*f/fs),d=1+k/q+k*k;
        const auto b0=(vh+vb*k/q+dc*k*k)/d;
        const auto b1=2*(dc*k*k-vh)/d,b2=(vh-vb*k/q+dc*k*k)/d;
        return (b0+b1*z+b2*z*z)/(1.L+2*(k*k-1)/d*z+(1-k/q+k*k)/d*z*z);
    };
    const auto vh=std::pow(10.L,3.99984385397L/20),vb=std::pow(vh,.499666774155L);
    const long double rlb_gain=4.L/(1.L+1.99004745483398L+.99007225036621L);
    return section(1681.9744509555319L,.7071752369554193L,vh,vb,1)*
           section(38.13547087613982L,.5003270373253953L,rlb_gain,0,0);
}
static std::complex<long double> response(BiquadCoefficients c,double fs,double hz){
    const auto z=std::exp(std::complex<long double>(0,-2*pi*hz/fs));
    return (static_cast<long double>(c.b0)+static_cast<long double>(c.b1)*z+static_cast<long double>(c.b2)*z*z)/
           (1.L+static_cast<long double>(c.a1)*z+static_cast<long double>(c.a2)*z*z);
}
static void k_rates(){
    double worst_db=0,worst_energy=0;
    for(double fs:rates){
        const auto c=design_k_weighting(fs);
        for(double hz:{10.,25.,38.13547087613982,100.,500.,997.,1681.9744509555319,2000.,10000.,20000.}){
            if(hz>=.49*fs)continue;
            const auto expected=reference_k(fs,hz);
            const auto got=response(c.shelf,fs,hz)*response(c.highpass,fs,hz);
            const double error=double(20*std::log10(std::abs(got/expected)));
            worst_db=std::max(worst_db,std::abs(error));near(error,0,2.e-6,"independent physical-Hz K response");
        }
        KWeightingFilter filter(fs),chunks(fs);
        const std::size_t count=std::size_t(std::ceil(fs*.15));
        std::vector<float> input(count),output(count),other(count);
        for(std::size_t i=0;i<count;++i)input[i]=float(.2*std::sin(2*pi*997*i/fs));
        filter.process(input.data(),output.data(),count);
        for(std::size_t at=0;at<count;){const auto n=std::min<std::size_t>(1+at%257,count-at);chunks.process(nullptr,nullptr,0);chunks.process(input.data()+at,other.data()+at,n);at+=n;}
        require(other==output,"K variable blocks");
        // Steady tone energy uses independent complex response, after 100 ms.
        long double measured=0,expected=0;
        const auto h=reference_k(fs,997);const auto start=std::size_t(std::ceil(fs*.1));
        for(std::size_t i=start;i<count;++i){measured+=static_cast<long double>(output[i])*output[i];
            const auto y=.2L*std::imag(h*std::exp(std::complex<long double>(0,2*pi*997*i/fs)));expected+=y*y;}
        const auto error=double(10*std::log10(measured/expected));worst_energy=std::max(worst_energy,std::abs(error));near(error,0,.00002,"independent K tone energy");
        filter.reset();near(filter.process(0),0,0,"K reset silence");
        filter.process(std::numeric_limits<float>::quiet_NaN());require(!filter.valid(),"K invalid input");
        filter.configure(48000);require(filter.valid()&&filter.sample_rate()==48000,"K setup epoch");
        rejects([&]{filter.configure(3999.99);});require(filter.sample_rate()==48000,"K failed setup mutated rate");
    }
    std::cout<<"K independent response max_abs_dB="<<worst_db<<" tone_energy_max_abs_dB="<<worst_energy<<'\n';
}
static double scan_gate(const std::vector<double>&e,double fs,std::size_t &count){
    const auto width=std::size_t(std::floor(static_cast<long double>(fs)*2/5+.5L));
    std::vector<long double> values;count=0;
    for(std::size_t j=0;;++j){
        const auto start=std::size_t(std::floor(static_cast<long double>(j)*fs/10+.5L));
        if(start+width>e.size())break;
        long double sum=0;for(std::size_t i=start;i<start+width;++i)sum+=e[i];
        values.push_back(sum/width);++count;
    }
    const auto absolute=std::pow(10.L,(-70.L+.691L)/10);long double sum=0;std::size_t n=0;
    for(auto v:values)if(v>absolute){sum+=v;++n;}
    if(!n)return -std::numeric_limits<double>::infinity();
    const auto threshold=std::max(absolute,sum/n/10);sum=0;n=0;
    for(auto v:values)if(v>threshold){sum+=v;++n;}
    return double(-.691L+10*std::log10(sum/n));
}
static void timing(){
    double worst=0;
    for(double fs:rates){
        metering_detail::HopClock clock(fs);std::uint64_t elapsed=0;
        // Absolute independent long-double schedules, one million hops: checks
        // remainder carry without waiting for hours of audio or global overflow.
        for(std::uint64_t j=1;j<=1000000;++j){elapsed+=clock.next();
            const long double ideal=static_cast<long double>(j)*fs/10;
            require(std::abs(static_cast<long double>(elapsed)-ideal)<=.50001L,"clock accumulated drift");}
        clock.reset();require(clock.next()==std::size_t(std::floor(static_cast<long double>(fs)/10+.5L)),"clock reset");
        EnvelopeFollower envelope(fs,.01,.1);double reference_envelope=0;
        const auto alpha=std::exp(-1/(fs*.01));
        for(int i=0;i<1000;++i){reference_envelope=alpha*reference_envelope+(1-alpha);near(envelope.process(1),reference_envelope,2.e-13,"fractional envelope timing");}
        SlidingWindow window(fs,.4),short_term(fs,3);
        require(short_term.capacity()==std::size_t(std::floor(static_cast<long double>(fs)*3+.5L)),"three second capacity at rate");
        const auto width=std::size_t(std::floor(static_cast<long double>(fs)*2/5+.5L));
        require(window.capacity()==width,"400 ms rounding");
        LoudnessGate gate(30,fs),chunks(30,fs);
        std::vector<double> e(std::size_t(std::ceil(fs*1.731)));
        for(std::size_t i=0;i<e.size();++i){const double t=double(i)/fs;e[i]=(t<.4?1.e-10:t<.7?.000001:t<1.2?.01:.0001)*(1+.3*std::sin(i*.17));}
        std::size_t expected_boundaries=0;
        for(std::size_t i=0;i<e.size();++i){if(gate.process(e[i])){
            const auto boundary=width+std::size_t(std::floor(static_cast<long double>(expected_boundaries)*fs/10+.5L));
            require(i+1==boundary,"fractional gate boundary");++expected_boundaries;}}
        for(std::size_t at=0;at<e.size();){const auto n=std::min<std::size_t>(1+at%503,e.size()-at);chunks.process(nullptr,0);chunks.process(e.data()+at,n);at+=n;}
        std::size_t expected_count=0;const auto reference=scan_gate(e,fs,expected_count);const auto got=gate.result();
        near(got.lufs,reference,2.e-12,"fractional gate independent scan");
        worst=std::max(worst,std::abs(got.lufs-reference));near(chunks.result().lufs,got.lufs,0,"fractional gate block independence");
        require(got.stored_blocks==expected_count&&expected_boundaries==expected_count,"gate count/tail");
        gate.reset();for(std::size_t i=0;i<width;++i)gate.process(0);require(std::isinf(gate.result().lufs),"gate reset silence");
        gate.process(std::numeric_limits<double>::infinity());require(!gate.valid(),"gate invalid input");
        gate.configure(48000.25);require(gate.valid()&&!gate.stored_blocks(),"gate setup resets");
        rejects([&]{gate.configure(0);});require(gate.sample_rate()==48000.25,"failed gate setup changed epoch");
        window.configure(fs,3);require(!window.count(),"window setup reset");rejects([&]{window.configure(fs,0);});require(window.sample_rate()==fs,"failed window setup changed epoch");
        // Channel layout mapping is independent of rate and remains explicit.
        ChannelWeighting w{ChannelRole::RightSurround,ChannelRole::Lfe,ChannelRole::Centre};const double energies[]{.1,1.e90,.2};near(w.process(energies),.341,1.e-15,"layout at rate");
    }
    std::cout<<"fractional schedule: 1000000 hops/rate, <=0.50001 sample phase error; gate_max_abs_LU="<<worst<<'\n';
}
static void true_peak_rates(){
    double min_error=100,max_error=-100,multitone_worst=0,pulse_worst=0;
    auto tp_rates=rates;tp_rates.insert(tp_rates.end(),{1000,1000.25,192000.125,std::nextafter(48000.,0.),std::nextafter(48000.,100000.)});
    for(double fs:tp_rates){
        TruePeakOversampler tp(fs),chunks(fs);
        const auto legacy=fs==48000;
        const auto factor=legacy?4:std::max<std::size_t>(2,std::size_t(std::ceil(384000.L/fs)));
        require(tp.factor()==factor&&tp.taps_per_phase()==(legacy?12:64),"rate-aware FIR size");
        near(tp.latency_frames(),legacy?5.875:32,0,"FIR delay metadata");require(tp.flush_frames()==(legacy?11:63),"FIR flush metadata");
        // Closed-form bandlimited tones, not resampled official fixtures. The
        // fade keeps file-edge transients away from the steady-state maximum.
        for(double norm:{.01,.1,.25,.4,.45}){
            if(norm*fs>20000)continue;
            for(double phase:{0.,.37,1.13}){
                tp.reset();const std::size_t n=legacy?4096:1024;
                for(std::size_t i=0;i<n;++i){const auto fade=std::min({1.,double(i)/128,double(n-1-i)/128});tp.process(float(.5*fade*std::sin(2*pi*norm*i+phase)));}
                tp.flush();const auto error=tp.dbtp()-20*std::log10(.5);
                if(!legacy){min_error=std::min(min_error,error);max_error=std::max(max_error,error);near(error,0,.3,"generated FIR analytic tone peak");}
            }
        }
        if(!legacy){
            // Independent dense continuous multisine reference. All harmonics
            // fit the stated validation band. No production FIR data is used.
            const int top=std::max(1,int(std::floor(std::min(.4,20000/fs)*128)));
            const std::array<int,4> harmonic{1,std::max(1,top/3),std::max(1,2*top/3),top};
            const auto signal=[&](long double t){long double y=0;for(int j=0;j<4;++j)y+=.2L*std::sin(2*pi*harmonic[std::size_t(j)]*t/128+.37L*j);return y;};
            long double reference=0;for(int i=0;i<32768;++i)reference=std::max(reference,std::abs(signal(static_cast<long double>(i)/256)));
            tp.reset();for(int i=0;i<1024;++i){const auto fade=std::min({1.,double(i)/128,double(1023-i)/128});tp.process(float(fade*signal(i)));}tp.flush();
            const auto error=20*std::log10(tp.value()/double(reference));multitone_worst=std::max(multitone_worst,std::abs(error));near(error,0,.3,"dense analytic multisine TP reference");
            // Fractionally shifted bandlimited sinc pulse; truncation at +/-80
            // input frames is an explicit generated finite-file test, not an
            // official transient fixture or an infinite-sinc proof.
            tp.reset();for(int i=0;i<162;++i){const auto t=2*std::min(.4,20000/fs)*(i-80.37L);tp.process(float(std::sin(pi*t)/(pi*t)));}tp.flush();
            const auto pulse_error=20*std::log10(tp.value());pulse_worst=std::max(pulse_worst,std::abs(pulse_error));near(pulse_error,0,.3,"generated shifted sinc pulse TP");
        }
        // Arbitrary signed input, variable blocks, tail containing final impulse.
        std::vector<float>x(193);for(std::size_t i=0;i<x.size();++i)x[i]=float(.2*std::sin(i*.73));x.back()=2;
        tp.reset();for(float v:x)tp.process(v);tp.flush();
        for(std::size_t at=0;at<x.size();){const auto n=std::min<std::size_t>(at%17+1,x.size()-at);chunks.process(nullptr,0);chunks.process(x.data()+at,n);at+=n;}
        near(chunks.flush(),tp.value(),0,"rate TP block independence/flush");near(chunks.process(0),0,1.e-30,"TP tail not drained");
        // Delayed impulse phase-zero peak at frame 32; no input peak shortcut.
        if(!legacy){tp.reset();for(std::size_t i=0;i<32;++i)require(tp.process(i==0?1:0)<1,"TP delay early impulse peak");near(tp.process(0),1,1.e-14,"TP impulse delay");}
        tp.reset();near(tp.process(0),0,0,"TP reset silence");tp.process(std::numeric_limits<float>::quiet_NaN());require(!tp.valid(),"TP invalid input");
        tp.configure(44100.25);require(tp.valid()&&tp.value()==0&&tp.sample_rate()==44100.25,"TP setup epoch");
        tp.process(1);const auto retained=tp.value();
        rejects([&]{tp.configure(768001);});require(tp.sample_rate()==44100.25&&tp.value()==retained,"TP failed setup changed epoch");
    }
    // Rate-independent kernel checked against a separate continuous-time sum
    // with phase DC normalization, including the entire zero-extended tail.
    double worst=0;
    for(double fs:{1000.,44100.,48000.25,192000.,768000.}){
        TruePeakOversampler tp(fs);const auto f=tp.factor();const std::vector<float>x{.3f,-.7f,.1f,1.2f,-.2f};
        for(std::size_t n=0;n<x.size()+63;++n){long double expected=0;
            for(std::size_t p=0;p<f;++p){long double y=0,dc=0;
                for(int k=0;k<64;++k){const auto t=static_cast<long double>(k)-32+static_cast<long double>(p)/f;
                    const auto h=(t==0?1.L:std::sin(pi*t)/(pi*t))*(.42L+.5L*std::cos(pi*t/32)+.08L*std::cos(2*pi*t/32));dc+=h;
                    if(n>=std::size_t(k)&&n-std::size_t(k)<x.size())y+=h*x[n-std::size_t(k)];}
                expected=std::max(expected,std::abs(y/dc));}
            const auto actual=tp.process(n<x.size()?x[n]:0);worst=std::max(worst,std::abs(actual-double(expected)));near(actual,double(expected),2.e-15,"generated FIR independent convolution");}
    }
    std::cout<<"generated TP analytic tone error_dB="<<min_error<<".."<<max_error<<" multisine_max_abs_dB="<<multitone_worst<<" sinc_pulse_max_abs_dB="<<pulse_worst<<" independent_convolution_max_abs="<<worst<<"; tested cases only\n";
}
static void controls(){
    for(double fs:{-1.,0.,999.999,768000.001,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()}){
        rejects([&]{KWeightingFilter k(fs);});rejects([&]{TruePeakOversampler tp(fs);});rejects([&]{LoudnessGate g(1,fs);});rejects([&]{SlidingWindow w(fs,.4);});
        rejects([&]{EnvelopeFollower f(fs);});rejects([&]{PeakHold f(fs);});}
    for(double fs:{1000.,1000.25,3999.999})rejects([&]{KWeightingFilter k(fs);});
    for(double bad:{-.1,0.,10.1,std::numeric_limits<double>::quiet_NaN()})rejects([&]{SlidingWindow w(48000,bad);});
    // Exact rational half-ties and their immediate binary64 neighbours.
    require(SlidingWindow(4001.25,.4).capacity()==1601,"400ms half-up tie");
    require(SlidingWindow(std::nextafter(4001.25,0.),.4).capacity()==1600,"400ms below half tie");
    require(SlidingWindow(std::nextafter(4001.25,5000.),.4).capacity()==1601,"400ms above half tie");
    require(SlidingWindow(1000.5,1).capacity()==1001,"generic window half-up tie");
    require(SlidingWindow(std::nextafter(1000.5,0.),1).capacity()==1000,"generic window below half tie");
    SlidingWindow maximum(768000,3);require(maximum.capacity()==2304000,"maximum three-second allocation");rejects([]{SlidingWindow w(768000,3.00001);});
    LoudnessGate bounded(1,48000.25);for(int i=0;i<25000;++i)bounded.process(1);require(bounded.overflowed()&&!bounded.result().complete,"fractional history overflow");
}
int main(){try{k_rates();timing();true_peak_rates();controls();std::cout<<"PASS arbitrary-rate numeric/energy/timing/TP/setup/layout limits; no blanket conformance claim\n";}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
