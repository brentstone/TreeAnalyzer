
#include "DataFormats/interface/FatJet.h"

namespace TAna {
//--------------------------------------------------------------------------------------------------
void FatJet::addFatJetInfo(const float bbt, const float tau1, const float tau2, const float tau3){
    _bbt = bbt;
    _tau1 =tau1;
    _tau2 =tau2;
    _tau3 = tau3;
}
//--------------------------------------------------------------------------------------------------
void FatJet::addSubJet(const SubJet& sj) {_sjs.push_back(sj);}
//--------------------------------------------------------------------------------------------------
float     FatJet::bbt()       const{return _bbt;}
float     FatJet::tau1()      const{return _tau1;}
float     FatJet::tau2()      const{return _tau2;}
float     FatJet::tau3()      const{return _tau3;}
float     FatJet::tau2otau1() const{return _tau1 == 0 ? 99 : _tau2/_tau1;}
float     FatJet::tau3otau1() const{return _tau1 == 0 ? 99 : _tau3/_tau1;}
float     FatJet::tau3otau2() const{return _tau2 == 0 ? 99 : _tau3/_tau2;}
//--------------------------------------------------------------------------------------------------
size     FatJet::nSubJets()  const{return _sjs.size();}
float     FatJet::minSJCSV()  const{
    if(!_sjs.size()) return 0;
    float minSJCSV = 1;
    for(const auto& sj : _sjs) { minSJCSV = std::min(minSJCSV,sj.csv());}
    return std::max(minSJCSV,float(0.0));
}
float     FatJet::maxSJCSV()  const{
    float maxSJCSV = 0;
    for(const auto& sj : _sjs) { maxSJCSV = std::max(maxSJCSV,sj.csv());}
    return maxSJCSV;
}
MomentumF FatJet::sdMom()     const{
    MomentumF sd;
    for(const auto& sj : _sjs) { sd.p4() += sj.p4(); }
    return sd;
}
MomentumF FatJet::rawSdMom()  const{
    MomentumF sd;
    for(const auto& sj : _sjs) { sd.p4() += sj.rawMom().p4(); }
    return sd;
}
//--------------------------------------------------------------------------------------------------
const SubJet&     FatJet::subJet(const size idx)  const {
    if(idx >= nSubJets())
        throw std::out_of_range("FatJet::subJet() -> Not a valid SubJet idx!)");
    return _sjs[idx];
}
SubJet&           FatJet::subJet(const size idx) {
    if(idx >= nSubJets())
        throw std::out_of_range("FatJet::subJet() -> Not a valid SubJet idx!)");
    return _sjs[idx];
}
}
