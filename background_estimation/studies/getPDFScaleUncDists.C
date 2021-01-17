
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "TreeReaders/interface/EventReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Configuration/interface/FillerConstants.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/EventSelection/interface/EventSelection.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Corrections/interface/TriggerScaleFactors.h"
#include "Processors/Corrections/interface/BTagScaleFactors.h"
#include "Processors/Corrections/interface/FatJetBTagScaleFactors.h"

#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "Processors/Corrections/interface/JetAndMETCorrections.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "AnalysisSupport/Utilities/interface/Types.h"


using namespace TAna;
using namespace CorrHelp;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed) {
    }

    bool pass1lSR() {
    	if(lepChan != SINGLELEP || !hbbCand || !wjjCand) return false;
    	if(hh.mass() < 700) return false;
    	if(hbbMass < 30 || hbbMass > 210) return false;
    	if(hbbTag < 0.8) return false;
    	if(nMedBTags_HbbV) return false;
    	if(wjjCand->tau2otau1() > 0.75) return false;
    	if(hWW.pt() / hh.mass() < 0.3) return false;
    	if(hwwLi > 11) return false;
    	return true;
    }

    bool pass2lSR() {
    	if(lepChan != DILEP || !hbbCand) return false;
    	if(hh.mass() < 700) return false;
    	if(hbbMass < 30 || hbbMass > 210) return false;
    	if(hbbTag < 0.8) return false;
    	if(nMedBTags_HbbV) return false;
    	if(llMass < 6 || llMass > 75) return false;
    	if(llDR > 1.0) return false;
    	if(llMetDphi > TMath::PiOver2()) return false;
    	if(reader_event->met.pt() < 85) return false;
    	return true;
    }

    void makePlots(TString sn, std::vector<float>& w_scale, std::vector<float>& w_pdf) {

    	auto plt = [&](TString wtS, double sf = 1.0) {
    		plotter.getOrMake1DPre(sn+"_"+wtS,"yield",";yield",1,-1,1)->Fill(0.0,weight*sf);
    	};

    	plt("nom");
    	for(unsigned i=0; i<w_scale.size(); i++) plt(TString::Format("s_%u",i),w_scale[i]);
    	for(unsigned i=0; i<w_pdf.size(); i++) plt(TString::Format("p_%u",i),w_pdf[i]);
    }

    bool runEvent() override {
    	if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
    	if(!passEventFilters) return false;

    	std::vector<float> w_scale, w_pdf;

        for(unsigned int i = 1; i < 9; ++i){
            w_scale.push_back(reader_event->genWeights[i]);
        }
        for(unsigned int i = 473; i < 574; ++i) {//Number 473 is the nominal
            w_pdf.push_back(reader_event->genWeights[i]);
        }

        makePlots(smpName,w_scale,w_pdf);
        if(pass1lSR()) makePlots(smpName+"_1l",w_scale,w_pdf);
        if(pass2lSR()) makePlots(smpName+"_2l",w_scale,w_pdf);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};



#endif

void getPDFScaleUncDists(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();
    a.write(outFileName);
}
