
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
using namespace FillerConstants;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed) {
    }


    std::pair<double,double> METXYCorr_Met_MetPhi(double originalMet, double originalMet_phi){

      std::pair<double,double>  TheXYCorr_Met_MetPhi(originalMet,originalMet_phi);

      DataEra era = DataEra(*reader_event->dataEra);
      int npv = *reader_event->npv;
      if(npv>100) npv=100;

      double METxcorr(0.),METycorr(0.);

      if(isRealData()) {
//    	  std::cout<<"isrealData"<<std::endl;
//Current recommendation for 2016 and 2018 DATA
    	  DataRun datarun = DataRun(*reader_event->dataRun);

    	  if(datarun==RUN2016B) {
    	  	METxcorr = -(-0.0478335*npv -0.108032);
    	  	METycorr = -(0.125148*npv +0.355672);
    	  }
    	  if(datarun==RUN2016C) {
    	  	METxcorr = -(-0.0916985*npv +0.393247);
    	  	METycorr = -(0.151445*npv +0.114491);
    	  }
    	  if(datarun==RUN2016D) {
    	  	METxcorr = -(-0.0581169*npv +0.567316);
    	  	METycorr = -(0.147549*npv +0.403088);
    	  }
    	  if(datarun==RUN2016E) {
    	  	METxcorr = -(-0.065622*npv +0.536856);
    	  	METycorr = -(0.188532*npv +0.495346);
    	  }
    	  if(datarun==RUN2016F) {
    	  	METxcorr = -(-0.0313322*npv +0.39866);
    	  	METycorr = -(0.16081*npv +0.960177);
    	  }
    	  if(datarun==RUN2016G) {
    	  	METxcorr = -(0.040803*npv -0.290384);
    	  	METycorr = -(0.0961935*npv +0.666096);
    	  }
    	  if(datarun==RUN2016H) {
    	  	METxcorr = -(0.0330868*npv -0.209534);
    	  	METycorr = -(0.141513*npv +0.816732);
    	  }
// -----------------------------------------------------------------
    	  if(datarun==RUN2017B) {
    	  	METxcorr = -(-0.19563*npv +1.51859);
    	  	METycorr = -(0.306987*npv +-1.84713);
    	  }
    	  if(datarun==RUN2017C) {
    	  	METxcorr = -(-0.161661*npv +0.589933);
    	  	METycorr = -(0.233569*npv +-0.995546);
    	  }
    	  if(datarun==RUN2017D) {
    	  	METxcorr = -(-0.180911*npv +1.23553);
    	  	METycorr = -(0.240155*npv +-1.27449);
    	  }
    	  if(datarun==RUN2017E) {
    	  	METxcorr = -(-0.149494*npv +0.901305);
    	  	METycorr = -(0.178212*npv +-0.535537);
    	  }
    	  if(datarun==RUN2017F) {
    	  	METxcorr = -(-0.165154*npv +1.02018);
    	  	METycorr = -(0.253794*npv +0.75776);
    	  }
// -------------------------------------------------------------------
    	  if(datarun==RUN2018A) {
    	  	METxcorr = -(0.362865*npv -1.94505);
    	  	METycorr = -(0.0709085*npv -0.307365);
    	  }
    	  if(datarun==RUN2018B) {
    	  	METxcorr = -(0.492083*npv -2.93552);
    	  	METycorr = -(0.17874*npv -0.786844);
    	  }
    	  if(datarun==RUN2018C) {
    	  	METxcorr = -(0.521349*npv -1.44544);
    	  	METycorr = -(0.118956*npv -1.96434);
    	  }
    	  if(datarun==RUN2018D) {
    	  	METxcorr = -(0.531151*npv -1.37568);
    	  	METycorr = -(0.0884639*npv -1.57089);
    	  }
      } else {
//    	  std::cout<<"isMC"<<std::endl;
// ---------------------------- MC ----------------------- //
          if(era==ERA_2016) {
          	METxcorr = -(-0.195191*npv -0.170948);
          	METycorr = -(-0.0311891*npv +0.787627);
          }

          if(era==ERA_2017) {
          	METxcorr = -(-0.182569*npv +0.276542);
          	METycorr = -(0.155652*npv +-0.417633);
          }

          if(era==ERA_2018) {
          	METxcorr = -(0.296713*npv -0.141506);
          	METycorr = -(0.115685*npv +0.0128193);
          }
      }

      double CorrectedMET_x = originalMet *cos( originalMet_phi)+METxcorr;
      double CorrectedMET_y = originalMet *sin( originalMet_phi)+METycorr;

      double CorrectedMET = sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y);
      double CorrectedMETPhi;
      if(CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
      else if(CorrectedMET_x==0 && CorrectedMET_y<0 )CorrectedMETPhi = -TMath::Pi();
      else if(CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x);
      else if(CorrectedMET_x <0&& CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) + TMath::Pi();
      else if(CorrectedMET_x <0&& CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) - TMath::Pi();
      else CorrectedMETPhi =0;

      TheXYCorr_Met_MetPhi.first= CorrectedMET;
      TheXYCorr_Met_MetPhi.second= CorrectedMETPhi;

      return TheXYCorr_Met_MetPhi;

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

    void makePlots(TString sn) {

    	double metpt = reader_event->met.pt();
    	double metphi = reader_event->met.phi();

    	auto plt = [&](TString wtS) {
    		plotter.getOrMake1DPre(sn+"_"+wtS,"met",";met",100,0,1000)->Fill(metpt,weight);
    		plotter.getOrMake1DPre(sn+"_"+wtS,"metphi",";metphi",100,-3.14,3.14)->Fill(metphi,weight);
    	};

    	plt("nom");

    	auto corrs = METXYCorr_Met_MetPhi(metpt,metphi);

    	metpt = corrs.first;
    	metphi = corrs.second;

    	plt("corr");

    }

    bool hasHemLep(const Lepton* lep1, const Lepton* lep2) {
    	bool isMuon1 = lep1->isMuon();
    	float lep1ETA = lep1->eta();
    	float lep1Phi = lep1->phi();

    	bool isMuon2 = lep2->isMuon();
    	float lep2ETA = lep2->eta();
    	float lep2Phi = lep2->phi();

        if(!isMuon1 && lep1ETA <= -1.479 && lep1Phi >= -1.55 && lep1Phi <= -0.9) return true;
        if(!isMuon2 && lep2ETA <= -1.479 && lep2Phi >= -1.55 && lep2Phi <= -0.9) return true;
        return false;
    }


    bool runEvent() override {
    	if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
    	if(!passEventFilters) return false;
    	if(lepChan == NOCHANNEL || !hbbCand) return false;

        makePlots(smpName);

        if(lepChan == DILEP) {
        	if(hasHemLep(dilep1,dilep2)) {
        		if(isRealData()) return false;
        		else weight *= (21.08/59.74);
        	}
        }

        if(pass1lSR()) makePlots(smpName+"_1l");
        if(pass2lSR()) makePlots(smpName+"_2l");

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};


#endif

void getMetPhiCorrs(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
