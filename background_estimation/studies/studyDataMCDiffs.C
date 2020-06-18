
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

    std::pair<int,int> getDecayType(const FatJet* fjet) {
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

        return std::make_pair(decayType,numB);
    }

    int getTrueSubjetFlavor() {

    	int flv1 = BTagging::jetFlavor(hbbCand->subJets()[0]);
    	int flv2 = BTagging::jetFlavor(hbbCand->subJets()[1]);

    	int comb = 0;
    	if(flv1 == BTagging::FLV_L) {
    		if(flv2 == BTagging::FLV_C) comb = 1;
    		else if(flv2 == BTagging::FLV_B) comb = 2;
    	} else if(flv1 == BTagging::FLV_C) {
    		if(flv2 == BTagging::FLV_L) comb = 1;
    		else if(flv2 == BTagging::FLV_C) comb = 3;
    		else if(flv2 == BTagging::FLV_B) comb = 4;
    	} else if(flv1 == BTagging::FLV_B) {
    		if(flv2 == BTagging::FLV_L) comb = 2;
    		else if(flv2 == BTagging::FLV_C) comb = 4;
    		else if(flv2 == BTagging::FLV_B) comb = 5;
    	}
        return comb;
    }

    bool pass1lCuts() {
    	if(wjjCand->tau2otau1() > 0.75) return false;
    	if(hWW.pt() / hh.mass() < 0.3) return false;
    	if(hwwLi > 11) return false;
    	return true;
    }

    bool pass2lCuts() {
    	if(llMass < 6 || llMass > 75) return false;
    	if(llDR > 1.0) return false;
    	if(llMetDphi > TMath::PiOver2()) return false;
    	if(reader_event->met.pt() < 85) return false;
    	return true;
    }

    void studyTrueHbbProperties(TString prefix) {
        if(!hbbCand) return;
    	if(hh.mass() < 700 || hh.mass() > 4000) return;
    	if(hbbMass < 30 || hbbMass > 210) return;
    	if(lepChan == NOCHANNEL) return;
    	if(lepChan == SINGLELEP && !wjjCand) return;

        TString chS = TString::Format("%dl",lepChan);

        std::pair<int,int> genStuff = getDecayType(hbbCand);
        int dtype = genStuff.first;
        int numB = genStuff.second;
        int trueFlav = getTrueSubjetFlavor();
//        if(trueFlav == 4) ParticleInfo::printGenInfo(reader_genpart->genParticles,-1);

        TString dStr = "tw";
        if(dtype==0) dStr = "qg";
        else if(dtype==4) dStr = "mw";
        else if(dtype==5) dStr = "mt";

        auto mkplots = [&](TString pre) {
        	if(!isRealData()) {
            	plotter.getOrMake2DPre(pre,"hbbCatxNumB",";hbbCat;numB",6,0.5,6.5,6,-0.5,5.5)->Fill(hbbCSVCat,numB,weight);
            	plotter.getOrMake2DPre(pre,"hbbCatxsjFlav",";hbbCat;sjFlav",6,0.5,6.5,6,-0.5,5.5)->Fill(hbbCSVCat,trueFlav,weight);
            	plotter.getOrMake2DPre(pre,"mbbxsjFlav",";mbb;sjFlav",30,30,210,6,-0.5,5.5)->Fill(hbbMass,trueFlav,weight);

            	std::vector<TString> sjFlavSs = {"ll","lc","lb","cc","cb","bb"};
            	plotter.getOrMake2DPre(pre,sjFlavSs[trueFlav]+"_mbbxhbbCat",";mbb;hbbCat",30,30,210,6,0.5,6.5)->Fill(hbbMass,hbbCSVCat,weight);
        	}
        	plotter.getOrMake2DPre(pre,"mbbxhbbCat",";mbb;hbbCat",30,30,210,6,0.5,6.5)->Fill(hbbMass,hbbCSVCat,weight);
        };

        auto plt = [&](TString str) {
        	mkplots("bkg_"+chS+"_"+str);
        	mkplots(prefix+"_"+chS+"_"+str);
        	mkplots(dStr+"_"+chS+"_"+str);
        	mkplots(prefix+"_"+dStr+"_"+chS+"_"+str);
        };

        auto go = [&](TString str) {
            if(!isRealData()) plt(str);
            else
            	if(hbbMass<=96 || hbbMass>=150) mkplots("data_"+chS+"_"+str);
        };

// -------- do selections --------------
        go("objSel");

        if(hbbCSVCat >= 4) go("hbbTagging");
        if(lepChan == SINGLELEP) {
        	if(pass1lCuts()) {
        		go("passCutsNoBT");
        		if(!nMedBTags_HbbV) {
        			go("passCutsWithBV");
        	        if(hbbCSVCat >= 4) go("fullSel");
        		}
        	}
        } else if(lepChan == DILEP) {
        	if(pass2lCuts()) {
        		go("passCutsNoBT");
        		if(!nMedBTags_HbbV) {
        			go("passCutsWithBV");
        	        if(hbbCSVCat >= 4) go("fullSel");
        		}
        	}
        }
    }

    bool isDeepFlavTagged(const Jet* j, BTagging::BTAGWP bwp) {
    	if(j->deep_flavor() >= parameters.jets.DeepFlavor_WP[bwp]) return true;
    	return false;
    }

    void testAK4Jets(TString pref, TString id) {
        auto jets = PhysicsUtilities::selObjsMom(reader_jet->jets,20,2.4);

        for(const auto* j : jets){
            if(!j->passTightID()) continue;
            auto flvI = BTagging::jetFlavor(*j);
            TString flvS = "l";
            if(flvI == BTagging::FLV_B) flvS = "b";
            else if(flvI == BTagging::FLV_C) flvS = "c";

            const float pt = j->pt();
            const float absETA = j->absEta();

            map<TString,BTagging::BTAGWP> wps = { {"incl",BTagging::BTAG_INCL}, {"loose",BTagging::BTAG_L},
            	{"med",BTagging::BTAG_M}, {"tight",BTagging::BTAG_T}};

            auto fill = [&](const TString& label) {
                plotter.getOrMake1DPre(pref+"_"+id,"yield",";cat",4,-0.5,3.5)->Fill(wps[label],weight);
                plotter.getOrMake2DPre(pref+"_"+id, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(pt,absETA,weight);

                if(!isRealData()) {
                    plotter.getOrMake1DPre(pref+"_"+id+"_"+flvS,"yield",";cat",4,-0.5,3.5)->Fill(wps[label],weight);
                    plotter.getOrMake2DPre(pref+"_"+id+"_"+flvS, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(pt,absETA,weight);
                }
            };

            if(!isDeepFlavTagged(j,BTagging::BTAG_L)) fill("none");
            else {
                if(!isDeepFlavTagged(j,BTagging::BTAG_M)) fill("loose");
                else {
                    if(!isDeepFlavTagged(j,BTagging::BTAG_T)) fill("med");
                    else fill("tight");
                }
            }

//            if(isDeepFlavTagged(j,BTagging::BTAG_L)) fill("loose");
//            if(isDeepFlavTagged(j,BTagging::BTAG_M)) fill("med");
//            if(isDeepFlavTagged(j,BTagging::BTAG_T)) fill("tight");
        }
    }

    void testAK8Jets(TString pref, TString id) {
    	if(!reader_fatjet) return;
        auto jets = PhysicsUtilities::selObjsMom(reader_fatjet->jets,50,2.4);
        if(lepChan != NOCHANNEL) return;

        if(!jets.size()) return;

        bool eventAlreadyPlotted = false;
        std::map<TString,int> cats = { {"ll",0}, {"lc",1}, {"lb",2},
        		{"cc",3},{"cb",4}, {"bb",5} };

        for(const auto* j : jets) {
        	if(j->nSubJets() != 2) continue;

        	const auto& sj1 = j->subJet(0);
        	const auto& sj2 = j->subJet(1);

        	int flv1 = -1, flv2 = -1;
        	TString flvS1 = "l", flvS2 = "l", flvSS = "ll";

        	if(!isRealData()) {
            	flv1 = BTagging::jetFlavor(sj1);
            	flv2 = BTagging::jetFlavor(sj2);

                if(flv1 == BTagging::FLV_B) {
                	flvS1 = "b";
                	if(flv2 == BTagging::FLV_B) {
                		flvS2 = "b";
                		flvSS = "bb";
                	} else if (flv2 == BTagging::FLV_C) {
                		flvS2 = "c";
                		flvSS = "cb";
                	} else {
                		flvSS = "lb";
                	}
                } else if (flv1 == BTagging::FLV_C) {
                	flvS1 = "c";
                	if(flv2 == BTagging::FLV_B) {
                		flvS2 = "b";
                		flvSS = "cb";
                	} else if (flv2 == BTagging::FLV_C) {
                		flvS2 = "c";
                		flvSS = "cc";
                	} else {
                		flvSS = "lc";
                	}
                } else {
                	if(flv2 == BTagging::FLV_B) {
                		flvS2 = "b";
                		flvSS = "lb";
                	} else if (flv2 == BTagging::FLV_C) {
                		flvS2 = "c";
                		flvSS = "lc";
                	}
                }
        	}

            const float pt1 = sj1.pt();
            const float eta1 = sj1.eta();
            const float pt2 = sj2.pt();
            const float eta2 = sj2.eta();

            auto fillForSJ = [&](const TString& hbbS, const TString& label, bool isSJ1) {
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

                plotter.getOrMake2DPre(pref+"_"+id+suf+hbbS, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(pt,fabs(eta),weight);
                plotter.getOrMake1DPre(pref+"_"+id+hbbS,"dcsv",";deep_csv",100,0,1)->Fill(sj1.deep_csv(),weight);
                plotter.getOrMake1DPre(pref+"_"+id+hbbS,"dcsv",";deep_csv",100,0,1)->Fill(sj2.deep_csv(),weight);

                if(!isRealData()) {
                    plotter.getOrMake2DPre(pref+"_"+id+"_"+flvS+suf+hbbS, label,";jet p_{T}[GeV];jet |#eta|",196,20,1000,24,0,2.4)->Fill(pt,fabs(eta),weight);
                    plotter.getOrMake1DPre(pref+"_"+id+"_"+flvS1+hbbS,"dcsv",";deep_csv",100,0,1)->Fill(sj1.deep_csv(),weight);
                    plotter.getOrMake1DPre(pref+"_"+id+"_"+flvS2+hbbS,"dcsv",";deep_csv",100,0,1)->Fill(sj2.deep_csv(),weight);
                }
            };

            auto fillFJ = [&](const TString& hbbS) {
            	int cat = BTagging::getCSVSJCat(parameters.jets,j->subJets());
                plotter.getOrMake1DPre(pref+"_"+id+hbbS,"yield",";cat",6,0.5,6.5)->Fill(cat,weight);
                plotter.getOrMake2DPre(pref+"_"+id+hbbS,"mj_x_yield",";jet mass;cat",100,0,300,6,0.5,6.5)->Fill(j->sdMom().mass(),cat,weight);

                if(!isRealData()) {
                    plotter.getOrMake1DPre(pref+"_"+id+"_"+flvSS+hbbS,"yield",";cat",6,0.5,6.5)->Fill(cat,weight);
                    plotter.getOrMake2DPre(pref+"_"+id+"_"+flvSS+hbbS,"mj_x_yield",";jet mass;cat",100,0,300,6,0.5,6.5)->Fill(j->sdMom().mass(),cat,weight);
                }
            };

            auto doSJ = [&](const TString& label, bool isSJ1) {
            	fillForSJ("",label,isSJ1);
//            	if(hbbCand && j->index() == hbbCand->index()) {
//                	fillForSJ("_hbb",label,isSJ1);
//            	}
            };

            auto doFJ = [&]() {
            	fillFJ("");
//            	if(hbbCand && j->index() == hbbCand->index()) {
//                	fillFJ("_hbb");
//            	}
            };

            bool isGoodSJ1 = (pt1 >= 20 && fabs(eta1) <= 2.4);
            bool isGoodSJ2 = (pt2 >= 20 && fabs(eta2) <= 2.4);

            if(isGoodSJ1) {
                doSJ("incl",true);
                if(BTagging::passSubjetBTagLWP(parameters.jets,sj1)) doSJ("loose",true);
                if(BTagging::passSubjetBTagMWP(parameters.jets,sj1)) doSJ("med",true);
            }
            if(isGoodSJ2) {
                doSJ("incl",false);
                if(BTagging::passSubjetBTagLWP(parameters.jets,sj2)) doSJ("loose",false);
                if(BTagging::passSubjetBTagMWP(parameters.jets,sj2)) doSJ("med",false);
            }

            if(isGoodSJ1 && isGoodSJ2) {
            	makeSomeAK8Plots(pref+"_"+id,j,eventAlreadyPlotted);
            	doFJ();
            }
        }
    }

    bool passElID(const Electron& el, int i, bool useIso) {
    	bool doesPass = false;
    	if(useIso) {
        	switch(i) {
        		case 0: doesPass = true; break;
        		case 1: doesPass = el.passLooseID(); break;
        		case 2: doesPass = el.passMedID(); break;
        		case 3: doesPass = el.passTightID(); break;
        		case 4: doesPass = el.passHEEPID(); break;
        		case 5: doesPass = el.passMVALooseID(); break;
        		case 6: doesPass = el.passMVA90ID(); break;
        		case 7: doesPass = el.passMVA80ID(); break;
        		default: break;
        	}
    	} else {
        	switch(i) {
    			case 0: doesPass = true; break;
        		case 1: doesPass = el.passLooseID_noIso(); break;
        		case 2: doesPass = el.passMedID_noIso(); break;
        		case 3: doesPass = el.passTightID_noIso(); break;
        		case 4: doesPass = el.passHEEPID_noIso(); break;
        		case 5: doesPass = el.passMVALooseID_noIso(); break;
        		case 6: doesPass = el.passMVA90ID_noIso(); break;
        		case 7: doesPass = el.passMVA80ID_noIso(); break;
        		default: break;
        	}
    	}
    	return doesPass;
    }

    bool passMuID(const Muon& mu, int i) {
    	bool doesPass = false;
        switch (i) {
        	case 0: doesPass = true; break;
        	case 1: doesPass = mu.passSoftID(); break;
        	case 2: doesPass = mu.passLooseID(); break;
        	case 3: doesPass = mu.passMedID(); break;
        	case 4: doesPass = mu.passTightID(); break;
        	case 5: doesPass = mu.passHighPT(); break;
        	default: break;
        }
        return doesPass;
    }

    void makeSomeAK8Plots(TString pre, const FatJet* j, bool& evtFilled) {

    	if(!evtFilled) {
        	plotter.getOrMake1DPre(pre+"_spec","met",";met",100,0,1000)->Fill(reader_event->met.pt(),weight);
        	plotter.getOrMake1DPre(pre+"_spec","ht",";ht",100,0,4000)->Fill(ht,weight);
        	plotter.getOrMake1DPre(pre+"_spec","met_o_ht",";ht",100,0,2)->Fill(reader_event->met.pt()/ht,weight);
        	plotter.getOrMake1DPre(pre+"_spec","met_o_rtht",";",200,0,100)->Fill(reader_event->met.pt()/sqrt(ht),weight);
        	plotter.getOrMake1DPre(pre+"_spec","met_o_rtrtht",";",200,0,400)->Fill(reader_event->met.pt()/pow(ht,0.25),weight);

        	plotter.getOrMake1DPre(pre+"_spec","maxjetpt",";pt",100,0,1000)->Fill(j->pt(),weight);
        	plotter.getOrMake1DPre(pre+"_spec","maxjeteta",";eta",100,-2.5,2.5)->Fill(j->eta(),weight);
        	plotter.getOrMake1DPre(pre+"_spec","maxtau21",";tau21",100,0,1)->Fill(j->tau2otau1(),weight);

        	for(const auto& l : reader_electron->electrons) {
        		if(l.pt() < 5 || l.absEta() > 2.5) continue;
            	plotter.getOrMake1DPre(pre+"_specEl","pt",";pt",100,0,1000)->Fill(l.pt(),weight);
            	plotter.getOrMake1DPre(pre+"_specEl","sip",";sip",200,0,10)->Fill(l.sip3D(),weight);
            	plotter.getOrMake1DPre(pre+"_specEl","dz",";dz",200,0,1)->Fill(fabs(l.dz()),weight);
            	plotter.getOrMake1DPre(pre+"_specEl","d0",";d0",200,0,1)->Fill(fabs(l.d0()),weight);
            	plotter.getOrMake1DPre(pre+"_specEl","eta",";eta",100,-2.5,2.5)->Fill(l.eta(),weight);
            	plotter.getOrMake1DPre(pre+"_specEl","miniIso",";miniIso",200,0,1)->Fill(l.miniIso(),weight);
            	plotter.getOrMake1DPre(pre+"_specEl","pfIso",";pfIso",200,0,1)->Fill(l.pfIso(),weight);

            	plotter.getOrMake2DPre(pre+"_specEl","pt_x_idnoIso",";pt;id",100,5,1000,8,-0.5,7.5)->Fill(l.pt(),0.0,weight);
            	plotter.getOrMake2DPre(pre+"_specEl","pt_x_pfIso",";pt;iso",100,5,1000,200,0,1)->Fill(l.pt(),l.pfIso(),weight);
            	plotter.getOrMake2DPre(pre+"_specEl","pt_x_miniIso",";pt;iso",100,5,1000,200,0,1)->Fill(l.pt(),l.miniIso(),weight);

            	for(int i=1; i<8; i++) {
            		if(passElID(l,i,true))
            			plotter.getOrMake2DPre(pre+"_specEl","pt_x_idwIso",";pt;id",100,5,1000,8,-0.5,7.5)->Fill(l.pt(),i,weight);
            		if(passElID(l,i,false))
            			plotter.getOrMake2DPre(pre+"_specEl","pt_x_idnoIso",";pt;id",100,5,1000,8,-0.5,7.5)->Fill(l.pt(),i,weight);
            	}
        	}

        	for(const auto& l : reader_muon->muons) {
        		if(l.pt() < 5 || l.absEta() > 2.4) continue;
            	plotter.getOrMake1DPre(pre+"_specMu","pt",";pt",100,0,1000)->Fill(l.pt(),weight);
            	plotter.getOrMake1DPre(pre+"_specMu","sip",";sip",200,0,10)->Fill(l.sip3D(),weight);
            	plotter.getOrMake1DPre(pre+"_specMu","dz",";dz",200,0,1)->Fill(fabs(l.dz()),weight);
            	plotter.getOrMake1DPre(pre+"_specMu","d0",";d0",200,0,1)->Fill(fabs(l.d0()),weight);
            	plotter.getOrMake1DPre(pre+"_specMu","eta",";eta",100,-2.5,2.5)->Fill(l.eta(),weight);
            	plotter.getOrMake1DPre(pre+"_specMu","miniIso",";miniIso",200,0,1)->Fill(l.miniIso(),weight);
            	plotter.getOrMake1DPre(pre+"_specMu","pfIso",";pfIso",200,0,1)->Fill(l.pfIso(),weight);

            	plotter.getOrMake2DPre(pre+"_specMu","pt_x_id",";pt;id",100,5,1000,8,-0.5,7.5)->Fill(l.pt(),0.0,weight);
            	plotter.getOrMake2DPre(pre+"_specMu","pt_x_pfIso",";pt;iso",100,5,1000,200,0,1)->Fill(l.pt(),l.pfIso(),weight);
            	plotter.getOrMake2DPre(pre+"_specMu","pt_x_miniIso",";pt;iso",100,5,1000,200,0,1)->Fill(l.pt(),l.miniIso(),weight);

            	for(int i=1; i<6; i++) {
            		if(passMuID(l,i))
            			plotter.getOrMake2DPre(pre+"_specMu","pt_x_id",";pt;id",100,5,1000,8,-0.5,5.5)->Fill(l.pt(),i,weight);
            	}
        	}
        	evtFilled = true;
    	}

    	plotter.getOrMake1DPre(pre+"_spec","jetpt",";pt",100,0,1000)->Fill(j->pt(),weight);
    	plotter.getOrMake1DPre(pre+"_spec","jeteta",";eta",100,-2.5,2.5)->Fill(j->eta(),weight);
    	plotter.getOrMake1DPre(pre+"_spec","tau21",";tau21",100,0,1)->Fill(j->tau2otau1(),weight);

    }

    void studyJetBtagging(TString pre, bool applIso, bool applMet) {
    	auto go = [&](TString suf) {
        	testAK4Jets(pre,"ak4"+suf);
        	testAK8Jets(pre,"ak8"+suf);
    	};

    	auto passIso = [&]() {
        	for(const auto& el : reader_electron->electrons) {
        		if(el.pt() < 5 || el.absEta() > 2.5) continue;
        		if(el.pfIso() < 0.01) return false;
        	}
        	for(const auto& mu : reader_muon->muons) {
        		if(mu.pt() < 5 || mu.absEta() > 2.4) continue;
        		if(mu.miniIso() < 0.01) return false;
        	}
    		return true;
    	};

    	if(!applIso && !applMet) {
    		go("");
    		return;
    	}

    	TString suf = "";
    	bool canPlot = true;
    	if(applIso) {
        	if(!passIso()) return;
        	suf = "_applIso";
        	if(applMet) {
        		if(reader_event->met.pt() > 150) return;
        		else suf += "Met";
        	}
    	} else if(applMet) {
    		if(reader_event->met.pt() > 150) return;
    		suf = "_applMet";
    	}

    	go(suf);
    }

    void mkJetRegPlots(TString pre) {
    	plotter.getOrMake1DPre(pre,"met",";met",100,0,1000)->Fill(reader_event->met.pt(),weight);
    	plotter.getOrMake1DPre(pre,"ht",";ht",100,0,4000)->Fill(ht,weight);
    	plotter.getOrMake1DPre(pre,"met_o_ht",";ht",100,0,2)->Fill(reader_event->met.pt()/ht,weight);
    	plotter.getOrMake1DPre(pre,"met_o_rtht",";ht",100,0,400)->Fill(reader_event->met.pt()/sqrt(ht),weight);

    	if(jets_HbbV.size()) {
        	plotter.getOrMake1DPre(pre,"minjetpt",";pt",100,0,1000)->Fill(jets_HbbV.back()->pt(),weight);
        	plotter.getOrMake1DPre(pre,"maxjetpt",";pt",100,0,1000)->Fill(jets_HbbV.front()->pt(),weight);
    	}

    	studyJetBtagging(pre,false,false);
    	studyJetBtagging(pre,true,false);
    	studyJetBtagging(pre,false,true);
    	studyJetBtagging(pre,true,true);
    }

    void isolateJetRegion(TString prefix) {
//    	if(ht < 400) return;
    	if(ht < 800) return;

    	if(!reader_jet->jets.size()) return;
    	if(!reader_fatjet->jets.size()) return;
//    	if(reader_event->met > 100) return;

//    	if(jets_HbbV.size() && jets_HbbV.front()->pt() < 180) return;

    	LeptonParameters pars1 = parameters.leptons;

    	pars1.el_maxDZ = 9999;
    	pars1.el_maxD0 = 9999;
    	pars1.el_maxSip3D = 9999;
    	pars1.mu_maxDZ = 9999;
    	pars1.mu_maxD0 = 9999;
    	pars1.mu_maxSip3D = 9999;
    	pars1.el_getID = &Electron::passLooseID_noIso;
//    	pars1.el_getISO = &Electron::pfIso;
    	pars1.el_maxISO = 0.4;
    	pars1.mu_getID = &Muon::passSoftID;
//    	pars1.mu_getISO = &Muon::pfIso;
    	pars1.mu_maxISO = 0.4;

    	auto leps = LeptonProcessor::getLeptons(pars1,*reader_muon,*reader_electron);
    	if(leps.size()) return;

    	mkJetRegPlots(prefix);
    }

    bool runEvent() override {
    	if(counter > 10) return false;
//    	counter++;
//    	if(reader_event->event.val() != 222026) return false;
        if(!DefaultSearchRegionAnalyzer::runEvent()) return false;
        if(!passEventFilters) return false;
        TString prefix = smpName;
        if (mcProc == FillerConstants::TTBAR && smDecayEvt.nLepsTT >= 0 && smDecayEvt.nLepsTT <= 2) {
        	prefix += TString::Format("%d",smDecayEvt.nLepsTT);
        	weight *= 0.870428;
        }

//        studyTrueHbbProperties(prefix);

        isolateJetRegion(prefix);

        return true;
    }

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;
    int counter = 0;

};

#endif

void studyDataMCDiffs(std::string fileName, int treeInt, int randSeed, std::string outFileName, float xSec=-1, float numEvent=-1){
    Analyzer a(fileName,"treeMaker/Events",treeInt,randSeed);
    a.setSampleInfo(xSec,numEvent);
    a.analyze();
    a.write(outFileName);
}
