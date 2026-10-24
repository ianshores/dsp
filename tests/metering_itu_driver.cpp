// SPDX-License-Identifier: MIT
// Offline fixture adapter only: explicit PCM16LE stdin, no parsing in library.
#include "dsp/metering.hpp"
#include <iomanip>
#include <iostream>
#include <string>

using namespace shores_dsp;
static void number(double x){if(std::isfinite(x))std::cout<<x;else std::cout<<"null";}
int main(int argc,char**argv){
    try{
        if(argc!=2)throw std::invalid_argument("usage: dsp_metering_itu_driver CHANNELS < PCM16LE");
        const int count=std::stoi(argv[1]);
        using R=ChannelRole;
        ChannelWeighting channels=count==1?ChannelWeighting{R::Centre}:count==2?ChannelWeighting{R::Left,R::Right}:
            count==5?ChannelWeighting{R::Left,R::Right,R::Centre,R::LeftSurround,R::RightSurround}:
            count==6?ChannelWeighting{R::Left,R::Right,R::Centre,R::Lfe,R::LeftSurround,R::RightSurround}:
            throw std::invalid_argument("fixture layout must be mono, stereo, 5.0 or ordered 5.1");
        std::array<KWeightingFilter,6> filters;
        SlidingWindow momentary(48000,.4),short_term(48000,3);LoudnessGate gate;
        std::array<unsigned char,256*6*2> input{};
        for(;;){
            std::cin.read(reinterpret_cast<char*>(input.data()),std::streamsize(256*count*2));
            const auto bytes=std::size_t(std::cin.gcount());
            if(bytes%std::size_t(count*2))throw std::runtime_error("truncated PCM frame");
            for(std::size_t at=0;at<bytes;at+=std::size_t(count*2)){
                std::array<double,6> energies{};
                for(int c=0;c<count;++c){
                    const auto p=at+std::size_t(c*2);int x=int(input[p])+(int(input[p+1])<<8);
                    if(x>=32768)x-=65536;
                    const double y=filters[std::size_t(c)].process(float(x)/32768.f);energies[std::size_t(c)]=y*y;
                }
                const auto e=channels.process(energies.data());gate.process(e);momentary.process(e);short_term.process(e);
            }
            if(bytes==0||std::cin.eof())break;
            if(!std::cin)throw std::runtime_error("PCM input failed");
        }
        const auto result=gate.result();bool valid=result.complete&&channels.valid()&&momentary.valid()&&short_term.valid();
        for(const auto &f:filters)valid=valid&&f.valid();
        std::cout<<std::setprecision(17)<<"{\"integrated_lufs\":";number(result.lufs);
        std::cout<<",\"momentary_lufs\":";number(momentary.ready()?energy_to_lufs(momentary.mean()):-INFINITY);
        std::cout<<",\"short_term_lufs\":";number(short_term.ready()?energy_to_lufs(short_term.mean()):-INFINITY);
        std::cout<<",\"blocks\":"<<result.stored_blocks<<",\"complete\":"<<(valid?"true":"false")<<"}\n";
        return valid?0:1;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
