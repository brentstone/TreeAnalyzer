
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
#include "Processors/Variables/interface/SignalHelper.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "TreeReaders/interface/JetReader.h"

#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "Processors/Variables/interface/FatJetSelection.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Variables/interface/Hww2lSolver.h"
#include "Processors/Corrections/interface/TriggerScaleFactors.h"


using namespace TAna;
using namespace std;
using namespace FillerConstants;

class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed)
: DefaultSearchRegionAnalyzer(fileName,treeName,treeInt, randSeed){

//        turnOffCorr(CORR_TRIG);
//        turnOffCorr(CORR_PU  );
//        turnOffCorr(CORR_LEP );
//        turnOffCorr(CORR_SJBTAG);
//        turnOffCorr(CORR_AK4BTAG);
//        turnOffCorr(CORR_SDMASS);
//        turnOffCorr(CORR_JER);
//        turnOnCorr(CORR_HEM1516);
    }

	void plt(TString sn, TString id) {
		plotter.getOrMake1DPre(sn,id+"_mx",";m_{X} [GeV]",100,700,4600)->Fill(signal_mass,weight);
	}

    bool passLepID(const Lepton *lep, bool is1l) {
    	if(is1l) {
    		if(lep->isMuon()) return ((const Muon*)lep)->passMedID();
    		else return ((const Electron*)lep)->passMVA90ID_noIso();
    	} else {
    		if(lep->isMuon()) return ((const Muon*)lep)->passLooseID();
    		else return ((const Electron*)lep)->passMedID_noIso();
    	}
    }

    bool passLepIso(const Lepton *lep, bool is1l) {
    	return (lep->miniIso() <= (is1l ? 0.2 : 0.15));
    }

    void getEfficiency1(TString sn, SignalHelper& sh, bool isMuon) {

    	TString chS = isMuon ? "m" : "e";

    	plt(sn,chS+"_base");

    	if(!sh.recolep1) return;
    	plt(sn,chS+"_reco");

    	bool passID1 = passLepID(sh.recolep1,true);
    	bool passIso1 = passLepIso(sh.recolep1,true);

    	if(passID1) {
    		plt(sn,chS+"_recoID");
    		if(passIso1) plt(sn,chS+"_recoIDIso");
    	}

    }

    void getEfficiency2(TString sn, SignalHelper& sh, bool isMuon1, bool isMuon2) {

    	TString ch1 = isMuon1 ? "m" : "e";
    	TString ch2 = isMuon2 ? "m" : "e";

    	plt(sn,ch1+"_base");
    	plt(sn,ch2+"_base");
    	plt(sn,ch1+ch2+"_base");

    	bool passID1 = sh.recolep1 ? passLepID(sh.recolep1,false) : false;
    	bool passID2 = sh.recolep2 ? passLepID(sh.recolep2,false) : false;

    	bool passIso1 = sh.recolep1 ? passLepIso(sh.recolep1,false) : false;
    	bool passIso2 = sh.recolep2 ? passLepIso(sh.recolep2,false) : false;

    	if(sh.recolep1) {
    		plt(sn,ch1+"_reco");
    		if(passID1) {
        		plt(sn,ch1+"_recoID");
        		if(passIso1) plt(sn,ch1+"_recoIDIso");
    		}
    	}
    	if(sh.recolep2) {
    		plt(sn,ch2+"_reco");
    		if(passID2) {
        		plt(sn,ch2+"_recoID");
        		if(passIso2) plt(sn,ch2+"_recoIDIso");
    		}
    	}
    	if(sh.recolep1 && sh.recolep2) {
    		plt(sn,ch1+ch2+"_reco");
    		if(passID1 && passID2) {
    			plt(sn,ch1+ch2+"_recoID");
    			if(passIso1 && passIso2) plt(sn,ch1+ch2+"_recoIDIso");
    		}
    	}

    }

    bool runEvent() override {
//    	std::cout<<std::endl<<"EVENT "<<reader_event->event.val()<<std::endl;

        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        if(!isSignal()) return false;

    	SignalHelper sh(diHiggsEvt,reader_muon,reader_electron);
    	if(!sh.isTrue1l && !sh.isTrue2l) return false;

        TString prefix = smpName.Contains("radion") ? "radion" : "blkgrav";

    	if(sh.isTrue1l && !sh.genlep1) return false;
    	if(sh.isTrue2l && (!sh.genlep1 || !sh.genlep2)) return false;

    	sh.setRecoLeptons(0.1);
    	sh.setRecoHbb(reader_fatjet->jets,0.15);
    	TString lepStr = sh.getLepStr();


        if(sh.isTrue2l) {
        	prefix += "_2l";
        	bool isMuon1 = lepStr.BeginsWith("m");
        	bool isMuon2 = lepStr.EndsWith("m");
        	bool goodGenPt2 = sh.genlep2->pt() > 10;
        	bool goodGenPt1 = (sh.genlep1->pt() > (isMuon1 ? 27 : 30)) || (sh.genlep2->pt() > (isMuon2 ? 27 : 30));
        	bool goodGenEta1 = sh.genlep1->absEta() < (isMuon1 ? 2.4 : 1.479);
        	bool goodGenEta2 = sh.genlep2->absEta() < (isMuon1 ? 2.4 : 1.479);

        	getEfficiency2(prefix,sh,isMuon1,isMuon2);
        	if(goodGenPt1 && goodGenPt2) {
        		getEfficiency2(prefix+"_genpt",sh,isMuon1,isMuon2);
        		if(goodGenEta1 && goodGenEta2) getEfficiency2(prefix+"_genpteta",sh,isMuon1,isMuon2);
        	}
        } else {
        	prefix += "_1l";
        	bool isMuon = (lepStr == "m");

        	bool goodGenPt = (sh.genlep1->pt() > (isMuon ? 27 : 30));
        	bool goodGenEta = sh.genlep1->absEta() < (isMuon ? 2.4 : 1.479);

        	getEfficiency1(prefix,sh,isMuon);
        	if(goodGenPt) {
            	getEfficiency1(prefix+"_genpt",sh,isMuon);
            	if(goodGenEta) getEfficiency1(prefix+"_genpteta",sh,isMuon);
        	}
        }


        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    int counter = 0;

};

#endif

void getSignalLeptonEff(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
