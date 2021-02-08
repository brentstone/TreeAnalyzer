
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
#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "Processors/Corrections/interface/JetAndMETCorrections.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "Processors/Corrections/interface/FatJetBTagScaleFactors.h"
#include "AnalysisSupport/Utilities/interface/Types.h"


#include "TPRegexp.h"
using namespace TAna;
using namespace CorrHelp;
class Analyzer : public DefaultSearchRegionAnalyzer {
public:

    Analyzer(std::string fileName, std::string treeName, int treeInt, int randSeed) : DefaultSearchRegionAnalyzer(fileName,treeName,treeInt,randSeed) {
        fjbtagSFProc_new.reset(new FatJetBTagScaleFactors (dataDirectory));
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
        	outTree->addSingle(hbbDecayType_,  "",  "hbbDecayType");
        	outTree->addSingle(nLepsTT_,  "",  "nLepsTT");
            outTree->addSingle(sampParam_,  "",  "sampParam");
        }

        outTree->addSingle(era_,"","era");
    	outTree->addSingle(ht_,  "",  "ht");
    	outTree->addSingle(met_,  "",  "met");
        outTree->addSingle(event_, "", "event");
        outTree->addSingle(run_, "", "run");
        outTree->addSingle(lepChan_,  "",  "lepChan");

        outTree->addSingle(isMuon1_, "", "isMuon1");
        outTree->addSingle(isMuon2_, "", "isMuon2");
        outTree->addSingle(lep1PT_, "", "lep1PT");
        outTree->addSingle(lep2PT_, "", "lep2PT");
    	outTree->addSingle(lep1ETA_,  "",  "lep1ETA");
    	outTree->addSingle(lep2ETA_,  "",  "lep2ETA");

        outTree->addSingle(dilepMass_, "", "dilepMass");
        outTree->addSingle(dilepDR_, "", "dilepDR");
        outTree->addSingle(llMetDphi_, "", "llMetDphi");

    	outTree->addSingle(hbbMass_,  "",  "hbbMass");
        outTree->addSingle(hbbPT_,  "",  "hbbPT");

    	outTree->addSingle(hbbCSVCat_,  "",  "hbbCSVCat");
    	outTree->addSingle(hbbTag_,  "",  "hbbTag");

    	outTree->addSingle(hhMass_,  "",  "hhMass");
    	outTree->addSingle(hwwPT_,  "",  "hwwPT");
        outTree->addSingle(hwwLi_,  "",  "hwwLi");
    	outTree->addSingle(wjjTau2o1_,  "",  "wjjTau2o1");

    	outTree->addSingle(nAK4Btags_,  "",  "nAK4Btags");
    	outTree->addSingle(nJets_hbbV,"","nJets_hbbV");

    	outTree->addSingle(w_ak8_central,  "",  "w_ak8_central");
    	outTree->addSingle(w_ak8_central_new,  "",  "w_ak8_central_new");
    	outTree->addSingle(w_ak8_up,  "",  "w_ak8_up");
    	outTree->addSingle(w_ak8_up_new,  "",  "w_ak8_up_new");
    	outTree->addSingle(w_ak8_down,  "",  "w_ak8_down");
    	outTree->addSingle(w_ak8_down_new,  "",  "w_ak8_down_new");
    }

    float getTTBAR_dak8_SF(const FatJet* hbb, CORRTYPE corr = NOMINAL) {
    	float sf = 1.0, norm = 1.0;
    	float pt = hbb->pt();
    	FillerConstants::DataEra yr = FillerConstants::DataEra(*reader_event->dataEra);

    	if(yr == FillerConstants::ERA_2016) {
    		if(pt <= 600) {
    			sf = 1.039;
    			norm = 0.72;

    			if(corr == UP) {
    				sf +=0.061;
    				norm += 0.05;
    			} else if(corr == DOWN) {
    				sf -= 0.058;
    				norm -= 0.05;
    			}

    		} else if(pt <= 800) {
    			sf = 1.035;
    			norm = 0.65;

    			if(corr == UP) {
    				sf += 0.105;
    				norm += 0.06;
    			} else if(corr == DOWN) {
    				sf -= 0.098;
    				norm -= 0.06;
    			}

    		} else {
    			sf = 1.301;
    			norm = 0.52;

    			if(corr == UP) {
    				sf += 0.325;
    				norm += 0.07;
    			} else if(corr == DOWN) {
    				sf -= 0.266;
    				norm -= 0.07;
    			}
    		}

    	} else if(yr == FillerConstants::ERA_2017) {

    		if(pt <= 600) {
    			sf = 0.91;
    			norm = 0.85;

    			if(corr == UP) {
    				sf += 0.05;
    				norm += 0.06;
    			} else if(corr == DOWN) {
    				sf -= 0.05;
    				norm -= 0.06;
    			}

    		} else if(pt <= 800) {
    			sf = 0.93;
    			norm = 0.87;

    			if(corr == UP) {
    				sf += 0.11;
    				norm += 0.08;
    			} else if(corr == DOWN) {
    				sf -= 0.09;
    				norm -= 0.08;
    			}

    		} else {
    			sf = 1.07;
    			norm = 0.74;

    			if(corr == UP) {
    				sf += 0.28;
    				norm += 0.09;
    			} else if(corr == DOWN) {
    				sf -= 0.25;
    				norm -= 0.09;
    			}
    		}

    	} else if(yr == FillerConstants::ERA_2018) {

    		if(pt <= 600) {
    			sf = 0.89;
    			norm = 0.83;

    			if(corr == UP) {
    				sf += 0.04;
    				norm += 0.06;
    			} else if(corr == DOWN) {
    				sf -= 0.05;
    				norm -= 0.06;
    			}

    		} else if(pt <= 800) {
    			sf = 0.94;
    			norm = 0.89;

    			if(corr == UP) {
    				sf += 0.08;
    				norm += 0.08;
    			} else if(corr == DOWN) {
    				sf -= 0.08;
    				norm -= 0.08;
    			}

    		} else {
    			sf = 1.05;
    			norm = 0.86;

    			if(corr == UP) {
    				sf += 0.21;
    				norm += 0.09;
    			} else if(corr == DOWN) {
    				sf -= 0.19;
    				norm -= 0.09;
    			}
    		}

    	} else {
    		std::cout<<"getTTBAR_dak8_SF fcked up"<<std::endl;
    	}



    	return sf*norm;
    }

    bool runEvent() override {
        bool passPre1 = true;
        bool passPre2 = true;
        if(!DefaultSearchRegionAnalyzer::runEvent() || !passEventFilters) {
        	passPre1 = false;
        	passPre2 = false;
        }

        if(!setupNewFJProc){
        	setupNewFJProc = true;
        	newParams = parameters.fatJets;
        	newParams.fatJetBtagSFFile = "corrections/btagging/deepak8_hbbSF_new.csv";
        	fjbtagSFProc_new->setParameters(newParams,int(*reader_event->dataEra));
        }

        if(lepChan != SINGLELEP || !hbbCand || !wjjCand)  passPre1 = false;
        if(lepChan != DILEP || !hbbCand)                  passPre2 = false;

        if(passPre2)       lepChan_ = DILEP;
        else if (passPre1) lepChan_ = SINGLELEP;
        else               lepChan_ = NOCHANNEL;

        if(lepChan_ == NOCHANNEL) return false;

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

            if(wjjCand) {
                hwwLi_  = hwwLi;
                hwwPT_ = hWW.pt();
                wjjTau2o1_ = wjjCand->tau2otau1();

                if(hbbCand) {
                    hhMass_  = hh.mass();
                } else {
                    hhMass_  = 0;
                }
            } else {
                hwwLi_  = 0;
                hwwPT_ = 0;
                wjjTau2o1_ = 0;
                hhMass_  = 0;
            }

            isMuon2_ = 0;
            lep2PT_  = 0;
            lep2ETA_ = 0;
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
        	dilepMass_ = llMass;
        	dilepDR_   = llDR;
        	llMetDphi_ = llMetDphi;
            hwwPT_ = hWW.pt();

            if(hbbCand) {
                hhMass_  = hh.mass();
            } else {
                hhMass_  = 0;
            }

            hwwLi_   = 0;
            wjjTau2o1_ = 0;
        } else {
            isMuon1_ = 0;
            isMuon2_ = 0;
            lep1PT_  = 0;
            lep2PT_  = 0;
            lep1ETA_  = 0;
            lep2ETA_  = 0;
        }

        if(hbbCand) {
            hbbMass_ = hbbMass;
            hbbCSVCat_ = hbbCSVCat;
            hbbTag_ = hbbTag;
            hbbPT_ = hbbCand->pt();
        } else {
            hbbMass_ = 0;
            hbbCSVCat_ = 0;
            hbbTag_ = 0;
        }
        nAK4Btags_   = std::min(nMedBTags_HbbV,250);

        if(!isRealData()) {
            fillWeights();
            fillGenVariables();

            bool isTTBAR = (FillerConstants::MCProcess(*reader_event->process) == FillerConstants::TTBAR);

            w_ak8_central  = float(fjbtagSFProc->getSF(parameters.fatJets,{hbbCand}));
            w_ak8_up       = float(fjbtagSFProc->getSF(parameters.fatJets,{hbbCand},UP));
            w_ak8_down     = float(fjbtagSFProc->getSF(parameters.fatJets,{hbbCand},DOWN));

            if(isTTBAR) {
            	w_ak8_central_new = getTTBAR_dak8_SF(hbbCand,NOMINAL);
            	w_ak8_up_new      = getTTBAR_dak8_SF(hbbCand,UP);
            	w_ak8_down_new    = getTTBAR_dak8_SF(hbbCand,DOWN);
            } else {
                w_ak8_central_new  = float(fjbtagSFProc_new->getSF(newParams,{hbbCand}));
                w_ak8_up_new       = float(fjbtagSFProc_new->getSF(newParams,{hbbCand},UP));
                w_ak8_down_new     = float(fjbtagSFProc_new->getSF(newParams,{hbbCand},DOWN));
            }

            nJets_hbbV = jets_HbbV.size();
        }

        return true;
    }

    void fillWeights() {
        xsec_    = getXSecWeight();
        pu_N_    = getPUWeight();
        lep_N_   = getLeptonWeight();
        trig_N_  = getTriggerWeight();
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
        }

        if (smDecayEvt.nLepsTT != -1) nLepsTT_ = smDecayEvt.nLepsTT;
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
    float w_ak8_central   = 0;
    float w_ak8_central_new   = 0;

    size8 lepChan_    = 0;
    size64 event_     = 0;
    size   run_       = 0;
    size16 sampParam_ = 0;
    size16 era_ = 0;

    //SR variables
    float ht_        = 0;
    float met_       = 0;

    size8 isMuon1_    = 0;
    size8 isMuon2_    = 0;
    float lep1PT_     = 0;
    float lep2PT_     = 0;
    float lep1ETA_    = 0;
    float lep2ETA_    = 0;

    float dilepMass_  = 0;
    float dilepDR_    = 0;
    float llMetDphi_  = 0;

    float hbbMass_   = 0;
    float hbbPT_     = 0;
    size8 hbbCSVCat_ = 0;
    float hbbTag_ = 0;

    float hhMass_    = 0;
    float hwwPT_     = 0;
    float hwwLi_   = 0;

    float wjjTau2o1_ = 0;

    size8 nAK4Btags_ = 0;
    size8 nJets_hbbV = 0;

    //BE extra variables
    size8 hbbDecayType_   = 0;
    size8 nLepsTT_ = 250;

    //systematic variables
    float w_ak8_up    = 0;
    float w_ak8_down    = 0;
    float w_ak8_up_new  = 0;
    float w_ak8_down_new  = 0;

    bool setupNewFJProc = false;
    std::unique_ptr<FatJetBTagScaleFactors> fjbtagSFProc_new;
    FatJetParameters newParams;

};




void doOne(const std::string& fileName, const int treeInt,int randSeed,  const std::string& outFileName, float xSec = -1, float numEvent = -1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    if(xSec > 0) a.setSampleInfo(xSec,numEvent);
    a.initializeTreeCopy(outFileName,BaseTreeAnalyzer::COPY_NONE);
    a.analyze();
}

#endif

void makeBtagTrees(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    doOne(fileName,treeInt,randSeed,outFileName,xSec,numEvent);
}
