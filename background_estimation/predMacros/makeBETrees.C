
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TreeAnalyzer/interface/BaseTreeCopier.h"
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


#include "TPRegexp.h"
using namespace TAna;
using namespace CorrHelp;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed) {
        addUncVariables = (treeType == TREE_OTHER);
        doPDFWeights = (treeInt > TREE_OTHER);
    }

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed, CORRTYPE jerUNC, CORRTYPE jesUNC,CORRTYPE metUNC,CORRTYPE hemUNC) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed){
    	JERProc->setCorrType(jerUNC);
        if(jesUNC == UP || jesUNC == DOWN ){
            JESUncProc ->setCorrType(jesUNC);
            turnOnCorr(CORR_JES);
        }
        if(metUNC == UP || metUNC == DOWN){
            METUncProc ->setCorrType(metUNC);
            turnOnCorr(CORR_MET);
        }
        if(hemUNC != NONE) {
            HEMIssueProc->setCorrType(hemUNC);
            turnOnCorr(CORR_HEM1516);
        }

    }



    virtual BaseEventAnalyzer * setupEventAnalyzer() override {return new CopierEventAnalyzer();}

    virtual void bookOutputVariables() override {

        if(isRealData()){
        	outTree->addSingle(dataset_,  "",  "dataset");
        	outTree->addSingle(dataRun_,  "",  "dataRun");
        } else {
        	outTree->addSingle(process_,  "",  "process");
        	outTree->addSingle(dhType_,  "",  "dhType");
        	outTree->addSingle(xsec_,  "",  "xsec");
        	outTree->addSingle(trig_N_,  "",  "trig_N");
        	outTree->addSingle(pu_N_,  "",  "pu_N");
        	outTree->addSingle(lep_N_,  "",  "lep_N");
        	outTree->addSingle(btag_N_,  "",  "btag_N");
//        	outTree->addSingle(sjbtag_N_, "", "sjbtag_N");
        	outTree->addSingle(fjbtag_N_, "", "fjbtag_N");
        	outTree->addSingle(topPt_N_, "",  "topPt_N");
        	outTree->addSingle(prefire_N_, "",  "prefire_N");
        	outTree->addSingle(hbbDecayType_,  "",  "hbbDecayType");
        	outTree->addSingle(eQuarksInHbb_,  "",  "eQuarksInHbb");
        	outTree->addSingle(nLepsTT_,  "",  "nLepsTT");
            outTree->addSingle(sampParam_,  "",  "sampParam");
        }

        outTree->addSingle(era_,"","era");
        outTree->addSingle(ht_,  "",  "ht");
        outTree->addSingle(met_,  "",  "met");
        outTree->addSingle(metphi_,  "",  "metphi");

        outTree->addSingle(event_, "", "event");
        outTree->addSingle(run_, "", "run");
        outTree->addSingle(lepChan_,  "",  "lepChan");
        outTree->addSingle(hasHEMLep_,  "",  "hasHEMLep");

        outTree->addSingle(isMuon1_, "", "isMuon1");
        outTree->addSingle(isMuon2_, "", "isMuon2");
        outTree->addSingle(lep1PT_, "", "lep1PT");
        outTree->addSingle(lep2PT_, "", "lep2PT");
        outTree->addSingle(lep1ETA_,  "",  "lep1ETA");
        outTree->addSingle(lep2ETA_,  "",  "lep2ETA");
        outTree->addSingle(lep1Phi_,  "",  "lep1Phi");
        outTree->addSingle(lep2Phi_,  "",  "lep2Phi");

        outTree->addSingle(dilepPT_, "", "dilepPT");
        outTree->addSingle(dilepMass_, "", "dilepMass");
        outTree->addSingle(dilepDR_, "", "dilepDR");
        outTree->addSingle(llMetDphi_, "", "llMetDphi");

        outTree->addSingle(hbbMass_,  "",  "hbbMass");
        outTree->addSingle(hbbMassUp_,  "",  "hbbMassUp");
        outTree->addSingle(hbbMassDown_,  "",  "hbbMassDown");

        outTree->addSingle(hbbPT_,  "",  "hbbPT");
        outTree->addSingle(hbbCSVCat_,  "",  "hbbCSVCat");
        outTree->addSingle(hbbTag_,  "",  "hbbTag");

        outTree->addSingle(hhMass_,  "",  "hhMass");
        outTree->addSingle(hhMassBasic_,  "",  "hhMassBasic");
        outTree->addSingle(hwwPT_,  "",  "hwwPT");
        outTree->addSingle(hwwLi_,  "",  "hwwLi");
        outTree->addSingle(wjjTau2o1_,  "",  "wjjTau2o1");
        outTree->addSingle(wjjMass_,  "",  "wjjMass");
        outTree->addSingle(wlnuMass_,  "",  "wlnuMass");
        outTree->addSingle(wlnuPT_,  "",  "wlnuPT");
        outTree->addSingle(wwDM_,  "",  "wwDM");

        outTree->addSingle(nAK4Btags_,  "",  "nAK4Btags");

        if(addUncVariables){
        	outTree->addSingle(w_muIDUp_,  "",  "w_muIDUp");
        	outTree->addSingle(w_muIDDown_,  "",  "w_muIDDown");
        	outTree->addSingle(w_elRecoUp_,  "",  "w_elRecoUp");
        	outTree->addSingle(w_elRecoDown_,  "",  "w_elRecoDown");
        	outTree->addSingle(w_elIDUp_,  "",  "w_elIDUp");
        	outTree->addSingle(w_elIDDown_,  "",  "w_elIDDown");
        	outTree->addSingle(w_bAk4_realUp_,  "",  "w_bAk4_realUp");
        	outTree->addSingle(w_bAk4_realDown_,  "",  "w_bAk4_realDown");
        	outTree->addSingle(w_bAk4_fakeUp_,  "",  "w_bAk4_fakeUp");
        	outTree->addSingle(w_bAk4_fakeDown_,  "",  "w_bAk4_fakeDown");
        	outTree->addSingle(w_bAk8_Up_,  "",  "w_bAk8_Up");
        	outTree->addSingle(w_bAk8_Down_,  "",  "w_bAk8_Down");
        	outTree->addSingle(w_puUp_,  "",  "w_puUp");
        	outTree->addSingle(w_puDown_,  "",  "w_puDown");
        	outTree->addSingle(w_prefireUp_,  "",  "w_prefireUp");
        	outTree->addSingle(w_prefireDown_,  "",  "w_prefireDown");

//        	outTree->addSingle(w_muISOUp_,  "",  "w_muISOUp");
//        	outTree->addSingle(w_elISOUp_,  "",  "w_elISOUp");
//        	outTree->addSingle(w_muISODown_,  "",  "w_muISODown");
//        	outTree->addSingle(w_elISODown_,  "",  "w_elISODown");

        	if(doPDFWeights) {
            	outTree->addVector(w_pdf_, "", "w_pdf_N", "w_pdf");
            	outTree->addVector(w_scale_, "", "w_scale_N", "w_scale");
        	}
        }

    }

    bool runEvent() override {
        bool passPre1 = true;
        bool passPre2 = true;
        if(!DefaultSearchRegionAnalyzer::runEvent() || !passEventFilters) {
        	passPre1 = false;
        	passPre2 = false;
        }

        if(lepChan != SINGLELEP || !hbbCand || !wjjCand)  passPre1 = false;
        if(lepChan != DILEP || !hbbCand)                  passPre2 = false;

        if(passPre2)       lepChan_ = DILEP;
        else if (passPre1) lepChan_ = SINGLELEP;
        else               lepChan_ = NOCHANNEL;

        if(!doPDFWeights && lepChan_ == NOCHANNEL) return false;

        switch(FillerConstants::DataEra(*reader_event->dataEra)){
        case FillerConstants::ERA_2018:
            era_ = 2018;
            break;
        case FillerConstants::ERA_2017:
            era_ = 2017;
            break;
        case FillerConstants::ERA_2016:
            era_ = 2016;
            break;
        default:
            era_ = 0;
        }

        ht_        = ht;
        met_       = reader_event->met.pt();
        metphi_    = reader_event->met.phi();
        event_     = *reader_event->event;
        run_       = *reader_event->run;
        sampParam_ = *reader_event->sampParam;

        if(isRealData()){
            dataset_ = *reader_event->dataset;
            dataRun_ = *reader_event->dataRun;
        }

        if(lepChan == SINGLELEP){
            isMuon1_  = selectedLepton->isMuon();
            lep1PT_   = selectedLepton->pt();
            lep1ETA_  = selectedLepton->eta();
            lep1Phi_  = selectedLepton->phi();

            if(wjjCand) {
                hwwLi_  = hwwLi;
                hwwPT_ = hWW.pt();
                wlnuMass_ = wlnu.mass();
                wlnuPT_    = wlnu.pt();
                wjjTau2o1_ = wjjCand->tau2otau1();
                wjjMass_   = wjjCand->sdMom().mass();
                wjjPT_     = wqq.pt();
                wwDM_      = wwDM;

                if(hbbCand) {
                    hhMass_  = hh.mass();
                    hhMassBasic_   = hh_basic.mass();
                } else {
                    hhMass_  = 0;
                    hhMassBasic_   = 0;
                }
            } else {
                hwwLi_  = 0;
                hwwPT_ = 0;
                wlnuMass_ = 0;
                wlnuPT_    = 0;
                wjjTau2o1_ = 0;
                wjjMass_   = 0;
                wjjPT_     = 0;
                wwDM_      = 0;

                hhMass_  = 0;
                hhMassBasic_   = 0;
            }

            isMuon2_ = 0;
            lep2PT_  = 0;
            lep2ETA_ = 0;
            lep2Phi_ = 0;
            dilepPT_ = 0;
            dilepMass_ = 0;
            dilepDR_   = 0;
            llMetDphi_ = 0;

        } else if (lepChan == DILEP) {
            isMuon1_ = dilep1->isMuon();
            isMuon2_ = dilep2->isMuon();
            lep1PT_  = dilep1->pt();
            lep2PT_  = dilep2->pt();
            lep1ETA_  = dilep1->eta();
            lep2ETA_  = dilep2->eta();
            lep1Phi_  = dilep1->phi();
            lep2Phi_  = dilep2->phi();

            dilepPT_   = (dilep1->p4()+dilep2->p4()).pt();
            dilepMass_ = llMass;
            dilepDR_   = llDR;
            llMetDphi_ = llMetDphi;
            hwwPT_ = hWW.pt();

            if(FillerConstants::DataEra(*reader_event->dataEra) == FillerConstants::ERA_2018) {
                if(!isMuon1_ && lep1ETA_ <= -1.479 && lep1Phi_ >= -1.55 && lep1Phi_ <= -0.9) hasHEMLep_ = true;
                if(!isMuon2_ && lep2ETA_ <= -1.479 && lep2Phi_ >= -1.55 && lep2Phi_ <= -0.9) hasHEMLep_ = true;
            }

            if(hbbCand) {
                hhMass_  = hh.mass();
            } else {
                hhMass_  = 0;
            }

            hwwLi_   = 0;
            wlnuMass_ = 0;
            wlnuPT_   = 0;
            wjjTau2o1_ = 0;
            wjjMass_   = 0;
            wjjPT_     = 0;
            hhMassBasic_ = 0;
            wwDM_ = 0;
        } else {
            isMuon1_ = 0;
            isMuon2_ = 0;
            lep1PT_  = 0;
            lep2PT_  = 0;
            lep1ETA_  = 0;
            lep2ETA_  = 0;
            lep1Phi_  = 0;
            lep2Phi_  = 0;
        }

        if(hbbCand) {
            hbbMass_ = hbbMass;
            hbbMassUp_ = hbbFJSFProc->getCorrSDMassUp(hbbCand);
            hbbMassDown_ = hbbFJSFProc->getCorrSDMassDown(hbbCand);

            hbbPT_ = hbbCand->pt();
            hbbCSVCat_ = hbbCSVCat;
            hbbTag_ = hbbTag;
        } else {
            hbbMass_ = 0;
            hbbMassUp_ = 0;
            hbbMassDown_ = 0;
            hbbPT_ = 0;
            hbbCSVCat_ = 0;
            hbbTag_ = 0;
        }
        nAK4Btags_ = std::min(nMedBTags_HbbV,250);

        if(!isRealData()) {
            fillWeights();
            fillGenVariables();
        }

        if(addUncVariables) fillUncertaintyVariables(doPDFWeights);

        return true;
    }

    void fillWeights() {
        xsec_    = getXSecWeight();
        pu_N_    = getPUWeight();
        lep_N_   = getLeptonWeight();
        fjbtag_N_= getFJBTagWeights();
        btag_N_  = getAK4BTagWeights();
//        sjbtag_N_= getSJBTagWeights();
        topPt_N_ = getTopPTWeight();
        trig_N_  = getTriggerWeight();
        prefire_N_ = (FillerConstants::DataEra(*reader_event->dataEra) == FillerConstants::ERA_2018) ? 1.0 : *reader_event->prefweight;
    }

    void fillGenVariables() {
        process_ = *reader_event->process;
        dhType_  = diHiggsEvt.type;

        auto getDecayType = [&](const FatJet* fjet) {
            double maxDR2 = 0.8*0.8;

            int topDecayType = 0; // NONE b wj wjb wjj wjjb bb wjbb wjjbb
            int maxQuarksFromTop = 0;
            int totQuarksFromTops = 0;
            int WDecayType  = 0; // NONE b wj wjb wjj wjjb bb wjbb wjjbb
            int maxQuarksFromW = 0;
            int totQuarksFromWs = 0;
            int numB = 0;

            for(const auto& d : smDecayEvt.topDecays) {
                if(d.type == TopDecay::BAD) continue;
                if(d.type != TopDecay::HAD){
                    if (PhysicsUtilities::deltaR2(*d.b,*fjet) < maxDR2) {
                        totQuarksFromTops++;
                        numB++;
                        if (maxQuarksFromTop == 0) maxQuarksFromTop = 1;
                    }
                } else {
                    if (!d.b) continue;
                    bool passB = (PhysicsUtilities::deltaR2(*d.b,*fjet) < maxDR2);
                    int nW = (PhysicsUtilities::deltaR2(*d.W_decay.dau1,*fjet) < maxDR2) +
                            (PhysicsUtilities::deltaR2(*d.W_decay.dau2,*fjet) < maxDR2);
                    int nT = nW + passB;
                    totQuarksFromTops += nT;
                    numB += passB;
                    if (nT > maxQuarksFromTop) maxQuarksFromTop = nT;
                }
            }

            if (totQuarksFromTops == 1){
                if(numB == 1) topDecayType = 1;
                else topDecayType = 2;
            } else if(totQuarksFromTops == 2){
                if (numB == 1) topDecayType = 3;
                else if (numB == 2) topDecayType = 6;
                else topDecayType = 4;
            } else if(totQuarksFromTops == 3) {
                if (numB == 1) topDecayType = 5;
                else topDecayType = 7;
            } else if (totQuarksFromTops == 4) topDecayType = 8;

            for (const auto& d : smDecayEvt.bosonDecays) {
                if(d.type != BosonDecay::Z_HAD && d.type != BosonDecay::W_HAD ) continue;
                int nW = (PhysicsUtilities::deltaR2(*d.dau1,*fjet) < maxDR2) +
                        (PhysicsUtilities::deltaR2(*d.dau2,*fjet) < maxDR2);
                if (d.dau1->absPdgId() == ParticleInfo::p_b &&
                        PhysicsUtilities::deltaR2(*d.dau1,*fjet) < maxDR2) numB++;
                if (d.dau2->absPdgId() == ParticleInfo::p_b &&
                        PhysicsUtilities::deltaR2(*d.dau2,*fjet) < maxDR2) numB++;

                totQuarksFromWs += nW;
                if (nW > maxQuarksFromW) maxQuarksFromW = nW;
            }

            if (totQuarksFromWs == 1) WDecayType = 2;
            else if (totQuarksFromWs == 2) WDecayType = 4;

            size8 decayType = 0;
            size8 nExtraQuarks = 0;
            if(maxQuarksFromTop >= maxQuarksFromW){
                decayType = topDecayType;
                nExtraQuarks = (totQuarksFromTops - maxQuarksFromTop) + totQuarksFromWs;
            } else {
                decayType = WDecayType;
                nExtraQuarks = (totQuarksFromWs - maxQuarksFromW) + totQuarksFromTops;
            }

            return std::make_pair(decayType,nExtraQuarks);
        };

        if (hbbCand) {
            std::pair<size8,size8> beVars = getDecayType(hbbCand);
            hbbDecayType_ = beVars.first;
            eQuarksInHbb_ = beVars.second;
        }

        if (smDecayEvt.nLepsTT != -1) nLepsTT_ = smDecayEvt.nLepsTT;
    }

    void fillUncertaintyVariables(bool doPDFWeights) {
        auto lepProc = &*(lepChan==DILEP ? dileptonSFProc : leptonSFProc);
        const float nomMu = lepProc->getMuonSF();
        const float nomEl = lepProc->getElectronSF();

        w_muIDUp_     = float(lepProc->getMuonSF(NONE,UP)*nomEl);
//            w_muISOUp_    = float(leptonSFProc->getMuonSF(NONE,NOMINAL,UP)*nomEl);
        w_elRecoUp_   = float(lepProc->getElectronSF(UP,NOMINAL)*nomMu);
        w_elIDUp_     = float(lepProc->getElectronSF(NOMINAL,UP)*nomMu);
//            w_elISOUp_    = float(leptonSFProc->getElectronSF(NOMINAL,NOMINAL,UP)*nomMu);
        w_bAk4_realUp_   = float(ak4btagSFProc->getSF(jets_HbbV,NOMINAL,UP));
        w_bAk4_realDown_ = float(ak4btagSFProc->getSF(jets_HbbV,NOMINAL,DOWN));
        w_bAk4_fakeUp_   = float(ak4btagSFProc->getSF(jets_HbbV,UP,NOMINAL));
        w_bAk4_fakeDown_ = float(ak4btagSFProc->getSF(jets_HbbV,DOWN,NOMINAL));

        w_bAk8_Up_       = float(fjbtagSFProc->getSF(parameters.fatJets,hbbCand,mcProc == FillerConstants::TTBAR,UP));
        w_bAk8_Down_     = float(fjbtagSFProc->getSF(parameters.fatJets,hbbCand,mcProc == FillerConstants::TTBAR,DOWN));

        w_puUp_       = float(puSFProc->getCorrection(*reader_event->nTruePUInts,CorrHelp::UP));

        w_muIDDown_   = float(lepProc->getMuonSF(NONE,DOWN)*nomEl);
//            w_muISODown_  = float(leptonSFProc->getMuonSF(NONE,NOMINAL,DOWN)*nomEl);
        w_elRecoDown_ = float(lepProc->getElectronSF(DOWN,NOMINAL)*nomMu);
        w_elIDDown_   = float(lepProc->getElectronSF(NOMINAL,DOWN)*nomMu);
//            w_elISODown_  = float(leptonSFProc->getElectronSF(NOMINAL,NOMINAL,DOWN)*nomMu);
        w_puDown_     = float(puSFProc->getCorrection(*reader_event->nTruePUInts,CorrHelp::DOWN));

        w_prefireUp_ = FillerConstants::DataEra(*reader_event->dataEra) == FillerConstants::ERA_2018 ? 1.0 : *reader_event->prefweightup;
        w_prefireDown_ = FillerConstants::DataEra(*reader_event->dataEra) == FillerConstants::ERA_2018 ? 1.0 : *reader_event->prefweightdown;

        if(doPDFWeights) {
            for(unsigned int i = 1; i < 9; ++i){
    //            if(i == 5 || i ==7)continue; //told to ignore
                w_scale_->push_back(reader_event->genWeights[i]);
            }
            for(unsigned int i = 473; i < 574; ++i) //Number 473 is the nominal
                w_pdf_->push_back(reader_event->genWeights[i]);
        }

    }

    //Event information and weights
    size8 process_    = 0;
    size8 dhType_     = 0;
    size8 dataset_    = 0;
    size8 dataRun_    = 0;
    float xsec_       = 0;
    float trig_N_     = 0;
    float pu_N_       = 0;
    float lep_N_      = 0;
    float btag_N_     = 0;
//    float sjbtag_N_   = 0;
    float fjbtag_N_   = 0;
    float topPt_N_    = 0;
    float prefire_N_  = 0;

    size8 lepChan_    = 0;
    size64 event_     = 0;
    size   run_       = 0;
    size16 sampParam_ = 0;
    size16 era_       = 0;
    size8 hasHEMLep_  = 0;

    //SR variables
    float ht_        = 0;
    float met_       = 0;
    float metphi_    = 0;

    size8 isMuon1_    = 0;
    size8 isMuon2_    = 0;
    float lep1PT_     = 0;
    float lep2PT_     = 0;
    float lep1ETA_    = 0;
    float lep2ETA_    = 0;
    float lep1Phi_    = 0;
    float lep2Phi_    = 0;

    float dilepPT_    = 0;
    float dilepMass_  = 0;
    float dilepDR_    = 0;
    float llMetDphi_  = 0;

    float hbbMass_   = 0;
    float hbbMassUp_   = 0;
    float hbbMassDown_   = 0;

    float hbbPT_     = 0;
    size8 hbbCSVCat_ = 0;
    float hbbTag_ = 0;

    float hhMass_    = 0;
    float hwwPT_     = 0;
    float hhMassBasic_ = 0;
    float hwwLi_   = 0;

    float wjjTau2o1_ = 0;
    float wjjMass_   = 0;
    float wjjPT_     = 0;
    float wlnuMass_  = 0;
    float wlnuPT_    = 0;
    float wwDM_      = 0;

    size8 nAK4Btags_ = 0;

    //BE extra variables
    size8 hbbDecayType_   = 0;
    size8 eQuarksInHbb_   = 0;

    size8 nLepsTT_ = 250;

    //systematic variables
    float w_muIDUp_      = 0;
    float w_elRecoUp_    = 0;
    float w_elIDUp_      = 0;
    float w_bAk4_realUp_    = 0;
    float w_bAk4_fakeUp_    = 0;
    float w_bAk4_realDown_  = 0;
    float w_bAk4_fakeDown_  = 0;
    float w_bAk8_Up_   = 0;
    float w_bAk8_Down_ = 0;
    float w_puUp_        = 0;
    float w_muIDDown_    = 0;
    float w_elRecoDown_  = 0;
    float w_elIDDown_    = 0;
    float w_puDown_      = 0;
    float w_prefireUp_   = 0;
    float w_prefireDown_ = 0;
//    float w_muISOUp_     = 0;
//    float w_elISOUp_     = 0;
//    float w_muISODown_   = 0;
//    float w_elISODown_   = 0;

    ASTypes::spv_float w_pdf_     = ASTypes::make_spv_float();
    ASTypes::spv_float w_scale_     = ASTypes::make_spv_float();

    bool addUncVariables = false;
    bool doPDFWeights = false;

};




void doOne(const std::string& fileName, const int treeInt,int randSeed, const std::string& outFileName, float xSec = -1, float numEvent = -1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    if(xSec > 0) a.setSampleInfo(xSec,numEvent);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();
}

void doOneVar(const std::string& fileName, const int treeInt,int randSeed, const std::string& outFileName, const CORRTYPE jerType, CORRTYPE jesUNC, CORRTYPE metUNC, CORRTYPE hemUNC, float xSec = -1, float numEvent = -1){
	Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed,jerType,jesUNC,metUNC,hemUNC);
    if(xSec > 0) a.setSampleInfo(xSec,numEvent);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();
}

#endif

void makeBETrees(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    doOne(fileName,treeInt,randSeed,outFileName,xSec,numEvent);
    if(treeInt==2 && (outFileName.find("Radion") != std::string::npos || outFileName.find("Bulk") != std::string::npos)) {
        size_t lastindex = outFileName.find_last_of(".");
        std::string extLessName = outFileName.substr(0, lastindex);
        doOneVar(fileName,treeInt,randSeed+1,extLessName+"_JERUp.root"  ,UP     ,NONE,NONE,NONE,xSec,numEvent);
        doOneVar(fileName,treeInt,randSeed+2,extLessName+"_JERDown.root",DOWN   ,NONE,NONE,NONE,xSec,numEvent);
        doOneVar(fileName,treeInt,randSeed+3,extLessName+"_JESDown.root",NOMINAL,DOWN,NONE,NONE,xSec,numEvent);
        doOneVar(fileName,treeInt,randSeed+4,extLessName+"_JESUp.root"  ,NOMINAL,UP  ,NONE,NONE,xSec,numEvent);
        doOneVar(fileName,treeInt,randSeed+5,extLessName+"_METDown.root",NOMINAL,NONE,DOWN,NONE,xSec,numEvent);
        doOneVar(fileName,treeInt,randSeed+6,extLessName+"_METUp.root"  ,NOMINAL,NONE,UP  ,NONE,xSec,numEvent);
        doOneVar(fileName,treeInt,randSeed+7,extLessName+"_HEM.root"    ,NOMINAL,NONE,NONE,DOWN,xSec,numEvent);
    }

}
