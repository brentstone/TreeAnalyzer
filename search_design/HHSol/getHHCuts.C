#include <AnalysisSupport/Utilities/interface/Types.h>
#include <Math/GenVector/LorentzVector.h>
#include <TH1.h>
#include <TMath.h>
#include "DataFormats/interface/FatJet.h"
#include "DataFormats/interface/GenParticle.h"
#include "DataFormats/interface/Lepton.h"
#include "DataFormats/interface/Momentum.h"
#include "Processors/GenTools/interface/DiHiggsEvent.h"
#include "Configuration/interface/ReaderConstants.h"
#include <TString.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include "TreeAnalyzer/framework/Processors/Variables/interface/HiggsSolver.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#if !defined(__CINT__) || defined(__MAKECINT__)


#include "HHSolTreeAnalyzer.h"



using namespace TAna;
using namespace FillerConstants;


class Analyzer : public HHSolTreeAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed,int step) :
        HHSolTreeAnalyzer(fileName,treeName,treeInt, randSeed), step(step), HSolver(""){
        if(fileName.find("2016") !=  std::string::npos){
            parameters = ReaderConstants::set2016Parameters();
        } else if(fileName.find("2017") !=  std::string::npos){
            parameters = ReaderConstants::set2017Parameters();
        } else if(fileName.find("2018") !=  std::string::npos){
            parameters = ReaderConstants::set2018Parameters();
        } else {
            parameters.hww.liFileName = "ERROR";
        }
        std::cout << "Loaded: "<< parameters.hww.liFileName<<std::endl;
        HSolver.setParamters(parameters.hww);
    }





    bool runEvent() override {
        if(!HHSolTreeAnalyzer::runEvent()) return false;
        if((*bbJet_SDmass < 100 || *bbJet_SDmass>150)) return false;
        if(*nAK4Btags != 0) return false;
        if(*hbbCat < 4) return false;

        HSolverLiInfo hwwInfo;
        double ll = HSolver.minimize(lepton,met,qqJet,*qqJet_SDmass,hwwInfo);

        const double lHH = (hwwInfo.hWW + bbJet.p4()).mass();
        const double ptom = hwwInfo.hWW.pt() / lHH;
        auto bN = HSolverBasic::getInvisible(met, (qqJet.p4()+lepton.p4()));

        const double md = PhysicsUtilities::deltaR( bN.p4() + lepton.p4() ,qqJet)
        * (bN.p4() + lepton.p4() + qqJet.p4()).pt()/2.0;

//        if(lHH < 700) return false;
        if(ptom < 0.3) return false;


        auto mkDPlts = [&](const TString& prefix, float lV){

            plotter.getOrMake1DPre(prefix, "mhh",";#it{m}_{HH}",40,0,4000)
                    ->Fill(lHH,weight);

            plotter.getOrMake1DPre(prefix+"_mIncl", "likeli",";likeli",255,-1,50)
                    ->Fill(lV,weight);

            if(std::abs(lHH - 800) < 800*0.10)
            plotter.getOrMake1DPre(prefix+"_m800", "likeli",";likeli",255,-1,50)
                    ->Fill(lV,weight);
            if(std::abs(lHH - 1000) < 1000*0.10)
            plotter.getOrMake1DPre(prefix+"_m1000", "likeli",";likeli",255,-1,50)
                    ->Fill(lV,weight);
            if(std::abs(lHH - 2000) < 2000*0.10)
            plotter.getOrMake1DPre(prefix+"_m2000", "likeli",";likeli",255,-1,50)
                    ->Fill(lV,weight);
            if(std::abs(lHH - 2500) < 2500*0.10)
            plotter.getOrMake1DPre(prefix+"_m2500", "likeli",";likeli",255,-1,50)
                    ->Fill(lV,weight);
            if(std::abs(lHH - 3000) < 3000*0.10)
            plotter.getOrMake1DPre(prefix+"_m3000", "likeli",";likeli",255,-1,50)
                    ->Fill(lV,weight);

        };

        auto mkCutLEPPlts = [&](const TString& prefix, float ll){
            mkDPlts(prefix+"_emu", ll);
            if(*isMuon == 1)
                mkDPlts(prefix+"_mu", ll);
            else
                mkDPlts(prefix+"_e", ll);
        };


        auto mkCutBTAGPlts = [&](const TString& prefix, float ll){
            mkCutLEPPlts(prefix+"_LMT", ll);
            if(*hbbCat == 4)
                mkCutLEPPlts(prefix+"_L", ll);
            if(*hbbCat == 5)
                mkCutLEPPlts(prefix+"_M", ll);
            if(*hbbCat == 6)
                mkCutLEPPlts(prefix+"_T", ll);
            if(*hbbCat >= 5)
                mkCutLEPPlts(prefix+"_MT", ll);
        };

        auto mkCut1Plts = [&](const TString& prefix){
            mkCutBTAGPlts(prefix+"_tauIncl", ll);

            if(*qqJet_t2ot1 < 0.75)
                mkCutBTAGPlts(prefix+"_tau0p75",ll);


            if(*qqJet_t2ot1 < 0.55)
                mkCutBTAGPlts(prefix+"_tau0p55",ll);
            if(*qqJet_t2ot1 < 0.45)
                mkCutBTAGPlts(prefix+"_tau0p45",ll);
            if(*qqJet_t2ot1 < 0.40)
                mkCutBTAGPlts(prefix+"_tau0p40",ll);
            if(*qqJet_t2ot1 < 0.35)
                mkCutBTAGPlts(prefix+"_tau0p35",ll);

            if(*qqJet_t2ot1 >= 0.75) return;


            if(*qqJet_t2ot1 >= 0.55)
                mkCutBTAGPlts(prefix+"_tau0p55to0p75",ll);
            if(*qqJet_t2ot1 >= 0.45)
                mkCutBTAGPlts(prefix+"_tau0p45to0p75",ll);
            if(*qqJet_t2ot1 >= 0.40)
                mkCutBTAGPlts(prefix+"_tau0p40to0p75",ll);
            if(*qqJet_t2ot1 >= 0.35)
                mkCutBTAGPlts(prefix+"_tau0p35to0p75",ll);


            if(md < 125){
                if(*qqJet_t2ot1 >= 0.55 && *qqJet_t2ot1 < 0.75)
                    mkCutBTAGPlts(prefix+"_oLPMD",ll);
                if(*qqJet_t2ot1 < 0.55)
                    mkCutBTAGPlts(prefix+"_oHPMD",ll);
            }


            if((*qqJet_t2ot1 < 0.75&& ll<11) && (ll>=2.5 || *qqJet_t2ot1 >= 0.55) )
                mkCutBTAGPlts(prefix+"_nLP",ll);
            if(ll<2.5 && *qqJet_t2ot1 < 0.55)
                mkCutBTAGPlts(prefix+"_nHP",ll);

        };

        if(isSignal()) mkCut1Plts(smpName);
        else {
            if(*process <= FillerConstants::WJETS)
                mkCut1Plts("bkg");
            mkCut1Plts("bkgWQCD");
            mkCut1Plts(smpName);
        }

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    int step;

    HSolverLi HSolver;
    ParameterSet parameters;
};

#endif

void getHHCuts(int step, std::string fileName, std::string outFileName,
        float xSec=-1, float numEvent=-1){

    Analyzer a(fileName,"treeMaker/Events",1,1,step);
    std::string outN = "hSolTrees_getCuts";


    outN += "_"+ outFileName+".root";

    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outN);

}
