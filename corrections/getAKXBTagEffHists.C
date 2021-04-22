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

//    void loadVariables() override {
//        reader_event       =loadReader<EventReader>   ("event",isRealData());
//        reader_fatjet      =loadReader<FatJetReader>  ("ak8PuppiJet",isRealData(),true,true);
//        reader_fatjet_noLep=loadReader<FatJetReader>  ("ak8PuppiNoLepJet",isRealData(),false,true);
//        reader_jet         =loadReader<JetReader>     ("ak4Jet",isRealData());
//        reader_electron    =loadReader<ElectronReader>("electron");
//        reader_muon        =loadReader<MuonReader>    ("muon",isRealData());
//
//        if(!isRealData()){
//            reader_genpart =loadReader<GenParticleReader>   ("genParticle");
//        }
//
//        checkConfig();
//    }

    bool isDeepFlavTagged(const Jet* j, BTagging::BTAGWP bwp) {
    	if(j->deep_flavor() >= parameters.jets.DeepFlavor_WP[bwp]) return true;
    	return false;
    }

    void testAK8Jets(TString pref, TString id) {
    	if(!reader_fatjet) return;
        auto jets = PhysicsUtilities::selObjsMom(reader_fatjet->jets,50,2.4);
        if(lepChan == NOCHANNEL) return;

        if(!jets.size()) return;

        for(const auto* j : jets) {
        	if(j->nSubJets() != 2) continue;
        	if(selectedLepton && PhysicsUtilities::deltaR(*selectedLepton,*j) <= 0.8) return;
        	if(selectedDileptons.size()) {
        		if(PhysicsUtilities::deltaR(*selectedDileptons[0],*j) <= 0.8) return;
        		if(selectedDileptons.size() >= 2) {
            		if(PhysicsUtilities::deltaR(*selectedDileptons[1],*j) <= 0.8) return;
        		}
        	}

        	const auto& sj1 = j->subJet(0);
        	const auto& sj2 = j->subJet(1);

        	int flv1 = BTagging::jetFlavor(sj1);
        	int flv2 = BTagging::jetFlavor(sj2);

            TString flvS1 = "l";
            if(flv1 == BTagging::FLV_B) flvS1 = "b";
            else if(flv1 == BTagging::FLV_C) flvS1 = "c";

            TString flvS2 = "l";
            if(flv2 == BTagging::FLV_B) flvS2 = "b";
            else if(flv2 == BTagging::FLV_C) flvS2 = "c";

            const float pt1 = sj1.pt();
            const float eta1 = sj1.eta();
            const float pt2 = sj2.pt();
            const float eta2 = sj2.eta();

            const float ptJ = j->pt();
            const float etaJ = j->eta();

            auto fillSJ = [&](const TString& label, bool isSJ1) {
                TString genS = isSignal() ? "sig" : "bkg";
                TString suf, flvS;
                float pt, eta;
                if(isSJ1) {
                	suf = "_sj1";
                	flvS = flvS1;
                	pt = pt1;
                	eta = eta1;
                } else {
                	suf = "_sj2";
                	flvS = flvS2;
                	pt = pt2;
                	eta = eta2;
                }

                plotter.getOrMake2DPre(pref+"_"+id+"_"+flvS+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(pt,fabs(eta),weight);
                plotter.getOrMake2DPre(genS+"_"+id+"_"+flvS+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(pt,fabs(eta),weight);

                plotter.getOrMake2DPre(pref+"_"+id+"_"+flvS+"_fulleta"+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,48,-2.4,2.4)->Fill(pt,eta,weight);
                plotter.getOrMake2DPre(genS+"_"+id+"_"+flvS+"_fulleta"+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,48,-2.4,2.4)->Fill(pt,eta,weight);

                if(!isSignal() && mcProc != FillerConstants::QCD) {
                    plotter.getOrMake2DPre("bkg_noQCD_"+id+"_"+flvS+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(pt,fabs(eta),weight);
                    plotter.getOrMake2DPre("bkg_noQCD_"+id+"_"+flvS+"_fulleta"+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,48,-2.4,2.4)->Fill(pt,fabs(eta),weight);
                }
            };

            auto fillFJ = [&](const TString& label) {
                TString genS = isSignal() ? "sig" : "bkg";
                TString suf = "_fj";
                TString flvS = flvS1+flvS2;


                plotter.getOrMake2DPre(pref+"_"+id+"_"+flvS+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(ptJ,fabs(etaJ),weight);
                plotter.getOrMake2DPre(genS+"_"+id+"_"+flvS+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(ptJ,fabs(etaJ),weight);

                plotter.getOrMake2DPre(pref+"_"+id+"_"+flvS+"_fulleta"+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,48,-2.4,2.4)->Fill(ptJ,etaJ,weight);
                plotter.getOrMake2DPre(genS+"_"+id+"_"+flvS+"_fulleta"+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,48,-2.4,2.4)->Fill(ptJ,etaJ,weight);

                if(!isSignal() && mcProc != FillerConstants::QCD) {
                    plotter.getOrMake2DPre("bkg_noQCD_"+id+"_"+flvS+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(ptJ,fabs(etaJ),weight);
                    plotter.getOrMake2DPre("bkg_noQCD_"+id+"_"+flvS+"_fulleta"+suf, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,48,-2.4,2.4)->Fill(ptJ,fabs(etaJ),weight);
                }
            };

            // subjet
            if(pt1 >= 20 && fabs(eta1) <= 2.4) {
                fillSJ("incl",true);
                if(BTagging::passSubjetBTagLWP(parameters.jets,sj1)) fillSJ("loose",true);
                if(BTagging::passSubjetBTagMWP(parameters.jets,sj1)) fillSJ("med",true);
            }
            if(pt2 >= 20 && fabs(eta2) <= 2.4) {
                fillSJ("incl",false);
                if(BTagging::passSubjetBTagLWP(parameters.jets,sj2)) fillSJ("loose",false);
                if(BTagging::passSubjetBTagMWP(parameters.jets,sj2)) fillSJ("med",false);
            }

            // fatjet
            if(pt1 >= 20 && pt2 >= 20 && fabs(eta1) <= 2.4 && fabs(eta2) <= 2.4 && fabs(etaJ) <= 2.4) {
            	fillFJ("incl");
            	if(j->deep_MDZHbb() >= 0.8)  fillFJ("loose");
            	if(j->deep_MDZHbb() >= 0.97) fillFJ("tight");
            }

        }

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
                    plotter.getOrMake2DPre("bkg_noQCD_"+id+"_"+flvS+"_fulleta", label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,48,-2.4,2.4)->Fill(pt,absETA,weight);
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

        testJets(prefix,"ak4_noHbbReq");
        testAK8Jets(prefix,"ak8_noHbbReq");

        if(!hbbCand) return false;
        testJets(prefix,"ak4_reqHbb");
        testAK8Jets(prefix,"ak8_reqHbb");

        if(hbbCand->sdMom().mass() < 10) return false;
        testJets(prefix,"ak4_goodMSD");
        testAK8Jets(prefix,"ak8_goodMSD");

        return true;
    }


    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

};

#endif

void getAKXBTagEffHists(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
