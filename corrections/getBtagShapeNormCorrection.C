
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Processors/Variables/interface/LeptonSelection.h"
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
//        turnOffCorr(CORR_AK4BTAG);
//        turnOffCorr(CORR_SDMASS);
//        turnOffCorr(CORR_JER);
    }

    void getCorr(TString pref, TString id) {
        float btagN = getAK4BTagWeights();

        auto plt = [&]() {
            plotter.getOrMake1DPre(pref+"_"+id+"_nobWT", "ht",";H_{T}[GeV]",1000,0,5000)->Fill(ht,weight/btagN);
            plotter.getOrMake1DPre(pref+"_"+id+"_wbWT", "ht",";H_{T}[GeV]",1000,0,5000)->Fill(ht,weight);
        };

        plt();
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(isRealData()) return false;
        if(!passEventFilters) return false;
        TString prefix = smpName;

        getCorr(prefix,"base");

        if(!passTriggerPreselection && !passTriggerPreselection2l) return false;
        getCorr(prefix,"passTrig");

        if(!hbbCand) return false;
        getCorr(prefix,"reqHbb");

        if(hbbCand->sdMom().mass() < 10) return false;
        getCorr(prefix,"goodMSD");

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getBtagShapeNormCorrection(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
