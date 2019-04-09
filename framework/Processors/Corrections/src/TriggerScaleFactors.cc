
#include "Processors/Corrections/interface/TriggerScaleFactors.h"

namespace TAna {
TriggerScaleFactors::TriggerScaleFactors(const std::string& dataDir, const std::string& triggerSFFile, bool verbose ){
    TFile * file = TObjectHelper::getFile(dataDir+triggerSFFile,"read",verbose);
    electronSFs.reset(new  TObjectHelper::Hist1DContainer(file,"electronSF",verbose) );
    muonSFs.reset(new  TObjectHelper::Hist1DContainer(file,"muonSF",verbose) );
    delete file;
}

void TriggerScaleFactors::setParameters(const std::string& dataDir, const LeptonParameters& lepParam, bool verbose) {
	TFile *file = TObjectHelper::getFile(dataDir+lepParam.leptonCorrSFFile,"read",verbose);
    electronSFs.reset(new  TObjectHelper::Hist1DContainer(file,"electronSF",verbose) );
    muonSFs.reset(new  TObjectHelper::Hist1DContainer(file,"muonSF",verbose) );
    delete file;
}

float TriggerScaleFactors::getElectronTriggerSF(const float ht) const{
    return electronSFs->getBinContentByValue(ht).val();
}
float TriggerScaleFactors::getMuonTriggerSF(const float ht) const{
    return muonSFs->getBinContentByValue(ht).val();
}
float TriggerScaleFactors::getLeptonTriggerSF(const float ht, const bool leadingLepIsMuon) const {
    if(leadingLepIsMuon) return getMuonTriggerSF(ht) ;
    return getElectronTriggerSF(ht);
}
}



