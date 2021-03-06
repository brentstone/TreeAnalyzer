#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/DileptonSearchRegionAnalyzer.h"
#include "TreeAnalyzer/interface/BaseTreeCopier.h"
#include "Configuration/interface/FillerConstants.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"

#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/JetKinematics.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "DataFormats/interface/GenParticle.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Variables/interface/LeptonSelection.h"

#include "TSystem.h"
using namespace TAna;
using namespace std;

class Analyzer : public DileptonSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DileptonSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){}

    void loadVariables() override {
        reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event);
        reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiJet",isRealData());       load(reader_fatjet);
        reader_jet     =std::make_shared<JetReader>     ("ak4PuppiJet",isRealData(),false); load(reader_jet);
        reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
        reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon);

        if(!isRealData()){
            reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");          load(reader_genpart   );
        }
    }

    bool passMMeID(const Lepton* lep1, const Lepton* lep2) {
    	bool pass = false;
    	if (((const Electron*)lep1)->passMedID_noISO() && ((const Electron*)lep2)->passMedID_noISO()) pass = true;
    	return pass;
    }

    void plot(TString sn, const Lepton* lep1, const Lepton* lep2, const FatJet* hbb) {
    	const MomentumF dilepmom = lep1->p4() + lep2->p4();
    	if (dilepmom.mass() < 12 || dilepmom.mass() > 75) return; // mass cut on 12 < Mll < 75

    	const MomentumF bbllmom = dilepmom.p4() + hbb->p4();
    	const MomentumF metLL = dilepmom.p4() + reader_event->met.p4();
    	double dPhi_metLL = PhysicsUtilities::deltaPhi(reader_event->met,dilepmom);
    	double dPhi_ll = PhysicsUtilities::deltaPhi(*lep1,*lep2);
    	double dR_ll = PhysicsUtilities::deltaR(*lep1,*lep2);

    	plotter.getOrMake1DPre(sn,"met",";MET",50,0,800)->Fill(reader_event->met.Et(),weight);
    	plotter.getOrMake1DPre(sn,"dPhi_metLL",";#Delta#Phi_{ll,miss}",50,0,3.14)->Fill(fabs(dPhi_metLL),weight);
    	plotter.getOrMake1DPre(sn,"met_o_ptbb","",50,0,2)->Fill(reader_event->met.Et()/hbb->pt(),weight);
    	plotter.getOrMake1DPre(sn,"met_o_dPhill",";GeV",100,0,2000)->Fill(reader_event->met.Et()/dPhi_ll,weight);
    	plotter.getOrMake1DPre(sn,"met_o_dRll",";GeV",100,0,2000)->Fill(reader_event->met.Et()/dR_ll,weight);
    	plotter.getOrMake1DPre(sn,"met_o_metLL",";GeV",100,0,5)->Fill(reader_event->met.Et()/metLL.Et(),weight);
    	plotter.getOrMake1DPre(sn,"met_o_ptLL",";GeV",100,0,5)->Fill(reader_event->met.Et()/dilepmom.pt(),weight);

    	plotter.getOrMake1DPre(sn,"pt_o_m",";GeV",100,0,1)->Fill(hww.pt()/hh.mass(),weight);
    	plotter.getOrMake1DPre(sn,"met_o_m",";GeV",100,0,1)->Fill(reader_event->met.pt()/hh.mass(),weight);
    	plotter.getOrMake1DPre(sn,"met_o_ptHww",";GeV",100,0,1)->Fill(reader_event->met.pt()/hww.pt(),weight);

    	plotter.getOrMake1DPre(sn,"met_o_ht",";GeV",100,0,1)->Fill(reader_event->met.pt()/ht_puppiNoLepOverlap,weight);
    	plotter.getOrMake1DPre(sn,"met_o_rootHTwlep",";GeV",100,0,50)->Fill(reader_event->met.pt()/sqrt(ht_puppi),weight);
    	plotter.getOrMake1DPre(sn,"met_o_rootHTnolep",";GeV",100,0,50)->Fill(reader_event->met.pt()/sqrt(ht_puppiNoLepOverlap),weight);
    	plotter.getOrMake1DPre(sn,"met_o_rootPTbb",";GeV",100,0,50)->Fill(reader_event->met.pt()/sqrt(hbbCand->pt()),weight);


    	// 2D plots
    	plotter.getOrMake2DPre(sn,"met_x_ht",";met;ht",50,0,1000,50,0,2000)->Fill(reader_event->met.pt(),ht_puppi,weight);
    	plotter.getOrMake2DPre(sn,"met_x_ptww",";met;ptww",50,0,1000,50,0,2000)->Fill(reader_event->met.pt(),hww.pt(),weight);
    	plotter.getOrMake2DPre(sn,"met_x_ptbb",";met;ptbb",50,0,1000,50,0,2000)->Fill(reader_event->met.pt(),hbbCand->pt(),weight);

    }
    void debugEvt(TString sn, const Lepton* lep1, const Lepton* lep2, const FatJet* hbb) {
    	ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);
    	cout<<endl<<sn<<endl;
    	printf("hbbjet (%i): {E = %.2f, pt = %.2f, eta = %.2f, phi = %.2f}\n",hbbCSVCat,hbb->E(),hbb->pt(),hbb->eta(),hbb->phi());
    	printf("lep1 (%i): {E = %.2f, pt = %.2f, eta = %.2f, phi = %.2f}\n",lep1->q(),lep1->E(),lep1->pt(),lep1->eta(),lep1->phi());
    	printf("lep2 (%i): {E = %.2f, pt = %.2f, eta = %.2f, phi = %.2f}\n",lep2->q(),lep2->E(),lep2->pt(),lep2->eta(),lep2->phi());
    	cout<<endl;
    }

    bool runEvent() override {
        if(!DileptonSearchRegionAnalyzer::runEvent()) return false;
        if(reader_event->process == FillerConstants::SIGNAL && diHiggsEvt.type != DiHiggsEvent::DILEP) return false;

        // current selection
        if(!hbbCand) return false;
        if(selectedDileptons.size() != 2) return false;
        if(hbbCSVCat < BTagging::CSVSJ_MF) return false;
        if (hbbMass < 30 || hbbMass > 210) return false;
        if (ht_puppi < 400) return false;
        if (PhysicsUtilities::deltaR2(*selectedDileptons[0],*selectedDileptons[1]) > 1.6*1.6) return false;
        if (nMedBTags_HbbV != 0) return false;
        if (hh.mass() < 700) return false;

        // assuming the input dataset is a skim using the Dilepton + Hbb selection
        TString chan = dilepMap[dilepChan];
        TString sn = smpName+chan;
        if (dilepChan == ee && !passMMeID(selectedDileptons[0],selectedDileptons[1])) return false;

        plot(sn,selectedDileptons[0],selectedDileptons[1],hbbCand);
        if (!isSignal()) plot("bkg"+chan,selectedDileptons[0],selectedDileptons[1],hbbCand);

        const MomentumF dilepMom = selectedDileptons[0]->p4() + selectedDileptons[1]->p4();

        if (fabs(PhysicsUtilities::deltaPhi(reader_event->met,dilepMom)) < TMath::PiOver2()) {
            plot(sn+"wDPHIcut_",selectedDileptons[0],selectedDileptons[1],hbbCand);
            if (!isSignal()) plot("bkg"+chan+"wDPHIcut_",selectedDileptons[0],selectedDileptons[1],hbbCand);
        }
        if (hh.mass() > 2000) {
            plot(sn+"Mhh2000_",selectedDileptons[0],selectedDileptons[1],hbbCand);
            if (!isSignal()) plot("bkg"+chan+"Mhh2000_",selectedDileptons[0],selectedDileptons[1],hbbCand);

            if (fabs(PhysicsUtilities::deltaPhi(reader_event->met,dilepMom)) < TMath::PiOver2()) {
                plot(sn+"wDPHIcut_Mhh2000_",selectedDileptons[0],selectedDileptons[1],hbbCand);
                if (!isSignal()) plot("bkg"+chan+"wDPHIcut_Mhh2000_",selectedDileptons[0],selectedDileptons[1],hbbCand);
            }
        }

        return true;
    }
    HistGetter plotter;
    void write(TString fileName){
    	plotter.write(fileName);
    }
};

#endif

void getMetVars(std::string fileName, int treeInt, int randSeed, std::string outFileName){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.analyze();
    a.write(outFileName);
}
void getMetVars(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec, float numEvent){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
