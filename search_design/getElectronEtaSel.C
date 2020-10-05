
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
#include "Processors/Variables/interface/HiggsSolver.h"

using namespace TAna;
using namespace std;
using namespace FillerConstants;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
: DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){}


    bool pass1lSel(float tau21, float liRatio, float nAk4bJets, float ptom, float dak8) {
    	if(liRatio > 11) return false;
    	if(nAk4bJets) return false;
    	if(tau21 > 0.75) return false;
    	if(ptom < 0.3) return false;
    	if(dak8 < 0.8) return false;

    	return true;
    }

    void testSingleElectron(TString sn, std::vector<float> etas) {

    	auto go = [&](float eta) {
    		LeptonParameters params = parameters.leptons;
    		params.el_maxETA = eta;

    		std::vector<const Lepton*> leps1 = LeptonProcessor::getLeptons(params,*reader_muon,*reader_electron);
    		if(!leps1.size()) return;
    	    if(!EventSelection::passTriggerPreselection(parameters.event,*reader_event,ht,leps1)) return;

    	    const Lepton* lep = leps1.front();

    	    const FatJet *hbbTemp = 0;
    	    const FatJet *wjjTemp = 0;
    	    if(reader_fatjet && reader_fatjet_noLep) {
    	        fjProc->loadFatJets(parameters.fatJets,*reader_fatjet,*reader_fatjet_noLep,lep);
    	        hbbTemp = fjProc->getHBBCand();
    	        wjjTemp = fjProc->getWjjCand();
    	    } else return;

    	    if(!hbbTemp || ! wjjTemp) return;
    	    if(hbbCand && hbbTemp->eta() != hbbCand->eta()) std::cout<<"chose different Hbb than hbbCand";
    	    if(wjjCand && wjjTemp->eta() != wjjCand->eta()) std::cout<<"chose different Wjj than wjjCand";

    	    float dak8 = hbbTemp->deep_MDZHbb();
    	    float mbb = isCorrOn(CORR_SDMASS) ? hbbFJSFProc->getCorrSDMass(hbbTemp) : hbbTemp->sdMom().mass();

    	    std::vector<const Jet*> jets_veto = PhysicsUtilities::selObjsD(jets,
    	            [&](const Jet* j){return  PhysicsUtilities::deltaR2(*j,*hbbTemp ) >= 1.2*1.2;});

    	    int nAk4Jets = PhysicsUtilities::selObjsMomD(jets_veto,
    	            parameters.jets.minBtagJetPT,parameters.jets.maxBTagJetETA,
    	            [&](const Jet* j){return BTagging::passJetBTagWP(parameters.jets,*j);} ).size();

    	    HSolverLiInfo info;
    	    const double qqSDMass = wjjTemp->sdMom().mass();
    	    float liRatio  = hSolverLi->minimize(lep->p4(),reader_event->met.p4(),
    	            wjjTemp->p4(), qqSDMass, info);

//    	    MomentumF nu = info.neutrino;
//    	    MomentumF lepW = info.wlnu;
//    	    MomentumF hadW = info.wqqjet;
    	    MomentumF HWW = info.hWW;

    	    MomentumF HH = hbbTemp->p4() + HWW.p4();
    	    float mhh = HH.mass();

    	    if(!pass1lSel(wjjTemp->tau2otau1(),liRatio,nAk4Jets,HWW.pt()/mhh,dak8)) return;
    	    if(mbb < 30 || mbb > 210) return;
    	    if(mhh < 700 || mhh > 4000) return;

    	    TString etaS = TString::Format("_eta%f",eta); etaS.ReplaceAll("0",""); etaS.ReplaceAll(".","p");
    	    if(isSignal()) {
        		plotter.getOrMake1DPre(sn+etaS,"mx",";M_{X}",100,700,4600)->Fill(signal_mass,weight);
    	    } else {
        		plotter.getOrMake1DPre(sn+etaS,"mhh",";M_{HH}",100,700,4000)->Fill(mhh,weight);
    	    }
    	};

    	for(const auto& eta : etas) {
    		go(eta);
    	}
    }

    bool runEvent() override {
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(ht < 400) return false;
        if(lepChan == DILEP) return false;
	    if(selectedDileptons.size() >= 2) return false;

	    TString sn = smpName;
        if (mcProc == FillerConstants::TTBAR && smDecayEvt.nLepsTT >= 0 && smDecayEvt.nLepsTT <= 2) {
        	sn += TString::Format("%d",smDecayEvt.nLepsTT);
        }
        if(isSignal()) sn = sn.BeginsWith("rad") ? "radion" : "blkgrav";

        std::vector<float> etas = {2.5,2.1,1.479};
        testSingleElectron(sn,etas);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getElectronEtaSel(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
