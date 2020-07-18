#include "../interface/BTagScaleFactors.h"
#include "DataFormats/interface/Jet.h"
#include "TreeReaders/interface/JetReader.h"
#include "DataFormats/interface/FatJet.h"
#include "Processors/Corrections/interface/BTagCalibrationStandalone.h"
#include "Configuration/interface/ReaderConstants.h"


namespace TAna {
using namespace CorrHelp;
using namespace BTagging;

//_____________________________________________________________________________
BTagScaleFactors::BTagScaleFactors(const std::string& dataDir, BTagging::BTAGWP maxWP,
        std::string heavySFName,std::string lightSFName) :
        dataDir(dataDir), maxWP(maxWP), heavySFName(heavySFName), lightSFName(lightSFName) {}
//_____________________________________________________________________________
BTagScaleFactors::~BTagScaleFactors() {}
//_____________________________________________________________________________
void BTagScaleFactors::setParameters(const std::string& sfFile, const std::string& effFile,
        bool verbose){
    efficiencies.clear();
    calibReaders.clear();
    calib.reset(new BTagCalibration("CSVv2", dataDir+sfFile));

    auto makeNewReader = [&] (const BTagEntry::OperatingPoint op){
        if(op == BTagEntry::OP_RESHAPING) {
            calibReaders.emplace_back(new BTagCalibrationReader(op, "central",{"up_hf", "down_hf"}));
        } else {
            calibReaders.emplace_back(new BTagCalibrationReader(op, "central",{"up", "down"}));
        }
        calibReaders.back()->load(*calib, BTagEntry::FLAV_B,heavySFName);
        calibReaders.back()->load(*calib, BTagEntry::FLAV_C,heavySFName);
        calibReaders.back()->load(*calib, BTagEntry::FLAV_UDSG,lightSFName);
    };

    if(doReshaping) {
        makeNewReader(BTagEntry::OP_RESHAPING );
        return;
    }

    if(maxWP >= BTAG_L) makeNewReader(BTagEntry::OP_LOOSE );
    if(maxWP >= BTAG_M) makeNewReader(BTagEntry::OP_MEDIUM);
    if(maxWP >= BTAG_T) makeNewReader(BTagEntry::OP_TIGHT );

    TFile * efile = TObjectHelper::getFile(dataDir+effFile,"read",verbose);
    std::vector<TString> taggers {"loose","med","tight"};
    std::vector<TString> flvs {"b","c","l"};
    for(unsigned int iF = FLV_B; iF <= FLV_L; ++iF){
        efficiencies.push_back(std::vector< std::unique_ptr<TObjectHelper::Hist2DContainer> >());
        //the tagger is displaced by 1 for the inclusive entry
        for(unsigned int iT = BTAG_L; iT <= maxWP; ++iT){
            efficiencies[iF].emplace_back(new  TObjectHelper::Hist2DContainer(efile,
                    TString::Format("%s_%s",flvs[iF].Data(),taggers[iT-1].Data()).Data(),verbose));
        }
    }
    efile->Close();
}
//_____________________________________________________________________________
float BTagScaleFactors::getJetEff(const BTagging::FLAVOR flv, const float pt, const float eta,
        const BTagging::BTAGWP wp) const {
    return efficiencies[flv][wp -1]->getBinContentByValue(pt, eta < 0 ? -eta : eta).val();
}
//_____________________________________________________________________________
float BTagScaleFactors::getJetSF(const BTagging::FLAVOR flv, const float pt, const float eta, const float disc,
        const BTagging::BTAGWP wp, CorrHelp::CORRTYPE corrT) const{
    //wp[0] = inclusive
    // BTagEntry::JetFlavor set to have same order as BTagging::Flavor
	return calibReaders[wp-1]->eval_auto_bounds(systNames[corrT],BTagEntry::JetFlavor(flv),eta,pt,disc);
}
//____________________________________________________________________________
float BTagScaleFactors::getJetCorr(const BaseRecoJet* jet,  CorrHelp::CORRTYPE lightT,
        CorrHelp::CORRTYPE heavyT) const {
    const float pt = jet->pt();
    const float eta = jet->eta();
    const auto flv = jetFlavor(*jet);
    const float disc = (jet->*btagCorrGetBTagVal)();
    const CorrHelp::CORRTYPE corrT =  flv == FLV_L ? lightT : heavyT;
    if(corrT == NONE) return 1.0;

    float lE = 1;
    float hE = 1;
    float lSF = 1;
    float hSF = 1;

//    printf("\nInput Jet:\nflv = %d\npt = %f\neta = %f\ndisc = %f",flv,pt,eta,disc);

    assignVals(flv,pt,eta,disc,corrT,lE,hE,lSF,hSF);
//    printf("\nlE = %f\nhE = %f\nlSF = %F\nhSF = %f\n",lE,hE,lSF,hSF);

    const float origEff = lE - hE;
    const float newEff = lSF*lE - hSF*hE;
    if(origEff == 0 || newEff == 0 || origEff > 1.0 || newEff > 1.0){
        TString errStr = TString::Format(
                "BTagScaleFactors::getJetCorr() -> Bad efficiencies: corr(pt,eta,flv):lowEff,highEff,lowSF,highSF,origEff,newEff -> %u(%f,%f,%u)%f,%f,%f,%f,%f,%f",
                corrT,pt,eta,flv,lE,hE,lSF,hSF,origEff,newEff);
        throw std::invalid_argument(errStr.Data());
    }

//    std::cout<<"sf = "<<newEff/origEff<<std::endl<<std::endl;
    return newEff/origEff;
}
//_____________________________________________________________________________
void BTagScaleFactors::setReshaping() {
	heavySFName = "iterativefit";
	lightSFName = "iterativefit";
	doReshaping = true;
}
//_____________________________________________________________________________
//_____________________________________________________________________________
JetBTagScaleFactors::JetBTagScaleFactors(const std::string& dataDir) :
                BTagScaleFactors(dataDir,BTAG_T,"comb","incl") {}
//_____________________________________________________________________________
JetBTagScaleFactors::~JetBTagScaleFactors() {}
//_____________________________________________________________________________
void JetBTagScaleFactors::setParameters(const JetParameters& parameters, bool verbose) {
    btagCorrWP = parameters.jetBtagCorrWP;
    btagCorrGetBTagVal = parameters.jetBtagCorrGetBTagVal;
    BTagScaleFactors::setParameters(parameters.jetBtagCorrSFFile,parameters.jetBtagCorrEffFile,
            verbose);
}
//_____________________________________________________________________________
void JetBTagScaleFactors::assignVals(const BTagging::FLAVOR flv,const float pt, const float eta,
        const float disc, const CorrHelp::CORRTYPE corrT,
        float& lE,float& hE,float& lSF,float& hSF) const {
    auto gE = [&](const BTagging::BTAGWP wp)->float{return getJetEff(flv,pt,eta,wp);};
    auto gS = [&](const BTagging::BTAGWP wp)->float{return getJetSF(flv,pt,eta,disc,wp,corrT);};

    if(disc <  btagCorrWP[BTAG_L]){
        lE = 1.0;
        hE = gE(BTAG_L);
        lSF = 1.0;
        hSF = gS(BTAG_L);
    } else if(disc <  btagCorrWP[BTAG_M]){
        lE = gE(BTAG_L);
        hE = gE(BTAG_M);
        lSF = gS(BTAG_L);
        hSF = gS(BTAG_M);
    } else if(disc <  btagCorrWP[BTAG_T]){
        lE = gE(BTAG_M);
        hE = gE(BTAG_T);
        lSF = gS(BTAG_M);
        hSF = gS(BTAG_T);
    } else {
        lE = gE(BTAG_T);
        hE = 0;
        lSF = gS(BTAG_T);
        hSF = 0;
    }
}
//_____________________________________________________________________________
float JetBTagScaleFactors::getSF(const std::vector<const Jet*>& jets, CorrHelp::CORRTYPE lightT,
        CorrHelp::CORRTYPE heavyT) const{
	float SF = 1.0;
    if(doReshaping) {
        for(const auto* j : jets) SF *= getReshapingCorr(j,lightT,heavyT);
    } else {
        for(const auto* j : jets) SF *= getJetCorr(j,lightT,heavyT);
    }
    return SF;
}
//_____________________________________________________________________________
float BTagScaleFactors::getReshapingCorr(const BaseRecoJet* j, CorrHelp::CORRTYPE lightT,
		CorrHelp::CORRTYPE heavyT) const {
    const float pt = j->pt();
    const float eta = j->eta();
    const auto flv = jetFlavor(*j);
    const float disc = (j->*btagCorrGetBTagVal)();
    const CorrHelp::CORRTYPE corrT =  flv == FLV_L ? lightT : heavyT;

//    printf("Input Jet:\npt = %f\neta = %f\nflv = %d\ndisc = %f\n",pt,eta,flv,disc);
	return calibReaders[0]->eval_auto_bounds(systNames_Reshaping[corrT],BTagEntry::JetFlavor(flv),eta,pt,disc);
}
//____________________________________________________________________________
//_____________________________________________________________________________

SubJetBTagScaleFactors::SubJetBTagScaleFactors(const std::string& dataDir) :
                BTagScaleFactors(dataDir,BTAG_M,"lt","incl") {}
//_____________________________________________________________________________
SubJetBTagScaleFactors::~SubJetBTagScaleFactors() {}
//_____________________________________________________________________________
void SubJetBTagScaleFactors::setParameters(const JetParameters& parameters, bool verbose) {
    btagCorrWP = parameters.sjBtagCorrWP;
    btagCorrGetBTagVal = parameters.sjBtagCorrGetBTagVal;
    BTagScaleFactors::setParameters(parameters.sjBtagCorrSFFile,parameters.sjBtagCorrEffFile,
            verbose);
}
//_____________________________________________________________________________
void SubJetBTagScaleFactors::assignVals(const BTagging::FLAVOR flv,const float pt, const float eta,
        const float disc, const CorrHelp::CORRTYPE corrT,
        float& lE,float& hE,float& lSF,float& hSF) const {
    auto gE = [&](const BTagging::BTAGWP wp)->float{return getJetEff(flv,pt,eta,wp);};
    auto gS = [&](const BTagging::BTAGWP wp)->float{return getJetSF(flv,pt,eta,disc,wp,corrT);};

    if(disc <  btagCorrWP[BTAG_L]){
        lE = 1.0;
        hE = gE(BTAG_L);
        lSF = 1.0;
        hSF = gS(BTAG_L);
    } else if(disc <  btagCorrWP[BTAG_M]){
        lE = gE(BTAG_L);
        hE = gE(BTAG_M);
        lSF = gS(BTAG_L);
        hSF = gS(BTAG_M);
    } else {
        lE = gE(BTAG_M);
        hE = 0;
        lSF = gS(BTAG_M);
        hSF = 0;

    }
}
//_____________________________________________________________________________
float SubJetBTagScaleFactors::getSF(const JetParameters& parameters,
        const std::vector<const FatJet*>& fatJets, CorrHelp::CORRTYPE lightT,
        CorrHelp::CORRTYPE heavyT) const{
    std::vector<const BaseRecoJet*> subjets;
    for(const auto* fj : fatJets){
        if(fj){
//        	printf("jet: pt = %f, eta = %f\n",fj->pt(),fj->eta());
            for(const auto& sj : fj->subJets()){
//            	printf("sj (%d): pt = %f, eta = %f, dcsv = %f\n",jetFlavor(sj),sj.pt(),sj.eta(),sj.deep_csv());
                if(sj.absEta() < parameters.maxBTagJetETA
                        && sj.pt() >= parameters.minBtagJetPT)
                    subjets.push_back(static_cast<const BaseRecoJet*>(&sj) );
            }
//            std::cout<<std::endl;
        }
    }

    float SF = 1.0;
    for(const auto* j : subjets) SF *= getJetCorr(j,lightT,heavyT);
    return SF;
}

}

