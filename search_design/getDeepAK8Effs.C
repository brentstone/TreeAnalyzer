
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DefaultSearchRegionAnalyzer.h"
#include "Configuration/interface/FillerConstants.h"
#include "Processors/EventSelection/interface/EventSelection.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/JetKinematics.h"
#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Variables/interface/DileptonSelection.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Variables/interface/Hww2lSolver.h"

#include "Processors/Variables/interface/SignalHelper.h"

using namespace TAna;
using namespace std;
using namespace FillerConstants;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
: DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){

        turnOffCorr(CORR_AK8BTAG);
    }

    void getEffDists(TString sn, const FatJet *fj) {

		plotter.getOrMake1DPre(sn,"inclB_mx",";M_{X} [GeV]",100,700,4600)->Fill(signal_mass,weight);

		bool passLoose = fj->deep_MDZHbb() >= 0.8;
		bool passTight = fj->deep_MDZHbb() >= 0.97;

		if(passLoose) plotter.getOrMake1DPre(sn,"passL_mx",";M_{X} [GeV]",100,700,4600)->Fill(signal_mass,weight);
		if(passTight) plotter.getOrMake1DPre(sn,"passT_mx",";M_{X} [GeV]",100,700,4600)->Fill(signal_mass,weight);

    }

    bool hasGoodHbb(SignalHelper& sh, bool isTrue1l) {
    	const FatJet *hbb = sh.recoHbb;

    	if(hbb->pt() < 200) return false;
    	if(hbb->absEta() > 2.4) return false;

    	for(const auto& sj : hbb->subJets()) {
    		if(sj.pt() < 20 || sj.absEta() > 2.4) return false;
    	}

    	if(isTrue1l) {
    		if(PhysicsUtilities::deltaR(*sh.genHbb,*diHiggsEvt.hww) < 1.6) return false;
    		if(PhysicsUtilities::absDeltaPhi(*sh.genHbb,*sh.genlep1) < 2.0) return false;
    	} else {
    		const MomentumF gendilep = sh.genlep1->p4() + sh.genlep2->p4();
    		if(PhysicsUtilities::absDeltaPhi(*sh.genHbb,gendilep) < 2.0) return false;
    	}

    	return true;
    }

    void filterEvents(TString sn, SignalHelper& sh, bool isTrue1l) {
    	if(isTrue1l && !sh.hasMatchedSingleLep()) return;
    	else if(!isTrue1l && !sh.hasMatchedDileps()) return;

    	getEffDists(sn+"_onlyMatchHbb",sh.recoHbb);

    	// filter on typical Hbb kinematic requirements
    	if(!hasGoodHbb(sh,isTrue1l)) return;

    	getEffDists(sn+"_goodHbb",sh.recoHbb);

    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(ht < 400) return false;

	    TString sn = smpName;
	    if(!isSignal()) return false;
	    if(diHiggsEvt.type != DiHiggsEvent::DILEP && diHiggsEvt.type < DiHiggsEvent::TAU_MU) return false;

	    bool isTrue1l = (diHiggsEvt.type != DiHiggsEvent::DILEP);
	    TString ls = isTrue1l ? "true1l" : "true2l";
		plotter.getOrMake1DPre(sn,ls+"_mx",";M_{X} [GeV]",100,700,4600)->Fill(signal_mass,weight);


        SignalHelper sigInfo(diHiggsEvt,reader_muon,reader_electron);
        sigInfo.minElRecoPt = 10;
        sigInfo.minMuRecoPt = 10;
        sigInfo.maxMuRecoEta = 2.4;
		sigInfo.maxElRecoEta = 2.5;
		sigInfo.setRecoLeptons(0.1);
		sigInfo.setRecoHbb(reader_fatjet->jets,0.1);

		if(!sigInfo.hasMatchedHbb()) return false;
		filterEvents(sn,sigInfo,isTrue1l);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getDeepAK8Effs(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
