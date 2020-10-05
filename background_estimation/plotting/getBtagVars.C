#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "Configuration/interface/FillerConstants.h"

using namespace std;

void getBtagVars() {
//	TString fArea = TString::Format("/Users/brentstone/Dropbox/Physics/HHbbWW/",year);
//	if(year == 0) fArea.ReplaceAll("0","Run2");
    HistGetter plotter;

    TFile *f1 = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/btagging/btagTrees2018.root");


    float hbbMass, hhMass, xsec, ht, trigN, puN, weight, hbbPt, hbbTag, lepN;
    float dilepMass, dilepDR, llMetDphi, met, tau21, hwwPt, hwwLi, lep1Eta, lep2Eta, lep1Pt, lep2Pt;
    float ak4_bR_up, ak4_bR_dn, ak4_bF_up, ak4_bF_dn, ak8_bR_up, ak8_bR_dn, ak8_bF_up, ak8_bF_dn;
    float ak8_b_nom, ak4_b_nom;
    UChar_t hbbCSVCat, hbbDecayType, nAK4Btags, lepChan, isMuon1, isMuon2,
		process, dataset, nLepsTT, nJets_hbbV;

//    float ttbarSF = 1.0;
//    if(year==2016) ttbarSF = 0.859681;
//    else if(year==2017) ttbarSF = 0.914884;
//    else if(year==2018) ttbarSF = 0.856187;
//    else if(year == 0)  ttbarSF = 0.870428;

    auto plots = [&](TString pref, TString bCat, TString fN, float sf) {
        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mbb",";M_{bb}",30,30,210)->Fill(hbbMass,weight*sf);
        plotter.getOrMake1DPre(pref+"_"+bCat,fN+"_mhh",";M_{HH}",132,700,4000)->Fill(hhMass,weight*sf);
    	plotter.getOrMake2DPre(pref+"_"+bCat,fN+"_hbbCatMbb",";hbbCat;mbb",6,0.5,6.5,30,30,210)
    	                		->Fill(hbbCSVCat,hbbMass,weight*sf);

        plotter.getOrMake2DPre(pref+"_"+bCat,fN+"_nJetsMbb",";N_{jets};M_{bb}",10,-0.5,9.5,30,30,210)
        		->Fill(nJets_hbbV,hbbMass,weight*sf);

    };

    auto pltVars = [&](TString pref, TString bCat, TString fN) {
    	std::vector<float> wts {ak4_bR_up,ak4_bR_dn,ak4_bF_up,ak4_bF_dn,
    		ak8_bR_up,ak8_bR_dn,ak8_bF_up,ak8_bF_dn,
			ak4_bF_up*ak8_bF_up, ak4_bF_dn*ak8_bF_dn, ak4_bR_up*ak8_bR_up,
			ak4_bR_dn*ak8_bR_dn};

    	std::vector<TString> wtIds = {"ak4_rU","ak4_rD","ak4_fU","ak4_fD","ak8_rU",
    			"ak8_rD","ak8_fU","ak8_fD","b_fU","b_fD","b_rU","b_rD"};

		plots(pref,bCat,fN+"_nom",1);

    	for(unsigned i=0; i<wts.size(); ++i) {
    		float corr = 0;
    		if(wtIds[i].Contains("ak4")) {
    			corr = wts[i]/ak4_b_nom;
    		} else if(wtIds[i].Contains("ak8")) {
    			corr = wts[i]/ak8_b_nom;
    		} else {
    			corr = wts[i]/(ak8_b_nom*ak4_b_nom);
    		}

    		plots(pref,bCat,fN+"_"+wtIds[i],corr);
    	}

    };

    auto mkPlots = [&](TString procName, TString fN, TString region, int nQuarks) {

        TString bName = "tw";
        if(nQuarks == 0) bName = "qg";
        else if (nQuarks==4) bName = "mw";
        else if (nQuarks==5) bName = "mt";

        pltVars(bName,region,fN);
    	pltVars(procName,region,fN);

        if(nQuarks==0 && process!=8) pltVars(bName+"_noQCD",region,fN);
        if(process!=8) pltVars("bkg_noQCD",region,fN);
    	pltVars("bkg" ,region,fN);

        TString ttS = TString::Format("%d",int(nLepsTT));
        if(process==2) {
            pltVars(procName+ttS,region,fN);
//            pltVars(procName+"_"+bName,region,fN);
//            pltVars(procName+ttS+"_"+bName,region,fN);
        }
    };

    auto passSR1 = [&]() {
    	if (lepChan != 1) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700 || hhMass > 4000) return false;
    	if (nAK4Btags != 0) return false;
//    	if (hbbTag < 0.5) return false;
    	if (hbbCSVCat < 4) return false;
    	if (tau21 > 0.75) return false;
    	if (hwwPt / hhMass < 0.3) return false;
    	if (hwwLi > 11) return false;

    	return true;
    };

    auto passSR2 = [&]() {
    	if (lepChan != 2) return false;
    	if (hbbMass < 30 || hbbMass > 210) return false;
    	if (hhMass < 700 || hhMass > 4000) return false;
    	if (nAK4Btags != 0) return false;
//    	if (hbbTag < 0.5) return false;
    	if (hbbCSVCat < 4) return false;
    	if (dilepMass < 6 || dilepMass > 75) return false;
    	if (dilepDR > 1.0) return false;
    	if (llMetDphi > TMath::PiOver2()) return false;
    	if (met < 85) return false;

    	return true;
    };


    TTree *t = (TTree*)f1->Get("treeMaker/Events");

    t->SetBranchAddress("hwwLi",&hwwLi);
    t->SetBranchAddress("hbbMass",&hbbMass);
    t->SetBranchAddress("hhMass",&hhMass);
    t->SetBranchAddress("hwwPT",&hwwPt);
    t->SetBranchAddress("wjjTau2o1",&tau21);
    t->SetBranchAddress("ht",&ht);
    t->SetBranchAddress("hbbCSVCat",&hbbCSVCat);
    t->SetBranchAddress("nAK4Btags",&nAK4Btags);
    t->SetBranchAddress("dilepMass",&dilepMass);
    t->SetBranchAddress("llMetDphi",&llMetDphi);
    t->SetBranchAddress("met",&met);
    t->SetBranchAddress("dilepDR",&dilepDR);
    t->SetBranchAddress("lepChan",&lepChan);
    t->SetBranchAddress("isMuon1",&isMuon1);
    t->SetBranchAddress("isMuon2",&isMuon2);
    t->SetBranchAddress("hbbTag",&hbbTag);
    t->SetBranchAddress("lep1PT",&lep1Pt);
    t->SetBranchAddress("lep2PT",&lep2Pt);
    t->SetBranchAddress("lep1ETA",&lep1Eta);
    t->SetBranchAddress("lep2ETA",&lep2Eta);

	t->SetBranchAddress("process",&process);
	t->SetBranchAddress("xsec",&xsec);
	t->SetBranchAddress("pu_N",&puN);
	t->SetBranchAddress("trig_N",&trigN);
	t->SetBranchAddress("lep_N",&lepN);
    t->SetBranchAddress("hbbDecayType",&hbbDecayType);
    t->SetBranchAddress("nLepsTT",&nLepsTT);
    t->SetBranchAddress("nJets_hbbV",&nJets_hbbV);

    t->SetBranchAddress("w_ak4_central",&ak4_b_nom);
    t->SetBranchAddress("w_ak4_realUp",&ak4_bR_up);
    t->SetBranchAddress("w_ak4_fakeUp",&ak4_bF_up);
    t->SetBranchAddress("w_ak4_realDown",&ak4_bR_dn);
    t->SetBranchAddress("w_ak4_fakeDown",&ak4_bF_dn);

    t->SetBranchAddress("w_ak8_central",&ak8_b_nom);
    t->SetBranchAddress("w_ak8_realUp",&ak8_bR_up);
    t->SetBranchAddress("w_ak8_fakeUp",&ak8_bF_up);
    t->SetBranchAddress("w_ak8_realDown",&ak8_bR_dn);
    t->SetBranchAddress("w_ak8_fakeDown",&ak8_bF_dn);

    printf("Will process %lld entries in total\n",t->GetEntries());
    for (unsigned int i=0; i<t->GetEntries(); i++) {
        if (i%100000 == 0) printf("processing evt %d\n",i);
//        std::cout<<int(process)<<std::endl;
//        if(process == 1) continue;
        t->GetEntry(i);
        TString pref = FillerConstants::MCProcessNames[process];
//        std::cout<<pref<<<<std::endl;
//        if(i>100) return;

        // set weights
        weight = xsec*trigN*puN*lepN*ak4_b_nom*ak8_b_nom;
        if(process == 2) weight *= 0.87; // ttbar scale factor

        bool is1l = passSR1();
        bool is2l = passSR2();

        if(is1l || is2l) {
        	mkPlots(pref,"","lepSel",hbbDecayType);
        	if(is1l) {
            	mkPlots(pref,"","SR1",hbbDecayType);
        	} else if(is2l) {
            	mkPlots(pref,"","SR2",hbbDecayType);
        	}
        }
    }

    plotter.write("debugBTAG.root");
}
