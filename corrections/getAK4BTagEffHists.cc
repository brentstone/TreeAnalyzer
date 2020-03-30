
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "Processors/EventSelection/interface/EventSelection.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "Processors/Variables/interface/BTagging.h"

using namespace TAna;
using namespace FillerConstants;
using namespace std;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){

//        turnOffCorr(CORR_TRIG);
//        turnOffCorr(CORR_PU  );
//        turnOffCorr(CORR_LEP );
        turnOffCorr(CORR_SJBTAG);
        turnOffCorr(CORR_AK4BTAG);
//        turnOffCorr(CORR_SDMASS);
//        turnOffCorr(CORR_JER);
    }

    bool isDeepFlavTagged(const Jet* j, BTagging::BTAGWP bwp) {
    	if(j->deep_flavor() >= parameters.jets.DeepFlavor_WP[bwp]) return true;
    	return false;
    }

    void testJets(TString pref, TString id) {
        auto jets = PhysicsUtilities::selObjsMom(reader_jet->jets,20,2.4);

        for(const auto* j : jets){
            if(!j->passTightID()) continue;
            auto flvI = BTagging::jetFlavor(*j);
            TString flvS = "l";
            if(flvI == BTagging::FLV_B) flvS = "b";
            else if(flvI == BTagging::FLV_C) flvS = "c";

            const float pt = j->pt();
            const float absETA = j->absEta();

            auto fill = [&](const TString& label) {
                TString genS = isSignal() ? "sig" : "bkg";

                plotter.getOrMake2DPre(pref+"_"+id+"_"+flvS, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(pt,absETA,weight);
                plotter.getOrMake2DPre(genS+"_"+id+"_"+flvS, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(pt,absETA,weight);

                plotter.getOrMake2DPre(pref+"_"+id+"_"+flvS+"_fulleta", label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,48,-2.4,2.4)->Fill(pt,j->eta(),weight);
                plotter.getOrMake2DPre(genS+"_"+id+"_"+flvS+"_fulleta", label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,48,-2.4,2.4)->Fill(pt,j->eta(),weight);

                if(!isSignal() && mcProc != FillerConstants::QCD) {
                    plotter.getOrMake2DPre("bkg_noQCD_"+id+"_"+flvS, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(pt,absETA,weight);
                }
            };

            fill("incl");
            if(isDeepFlavTagged(j,BTagging::BTAG_L)) fill("loose");
            if(isDeepFlavTagged(j,BTagging::BTAG_M)) fill("med");
            if(isDeepFlavTagged(j,BTagging::BTAG_T)) fill("tight");
        }
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(isRealData()) return false;
        if(!passEventFilters) return false;
        if(!passTriggerPreselection && !passTriggerPreselection2l) return false;

        TString prefix = smpName;

        testJets(prefix,"noHbbReq");

        if(!hbbCand) return false;
        testJets(prefix,"reqHbb");

        if(hbbCand->sdMom().mass() < 10) return false;
        testJets(prefix,"goodMSD");

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getAK4BTagEffHists(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
