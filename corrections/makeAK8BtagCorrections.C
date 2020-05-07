#include "HistoPlotting/include/PlotTools.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

using namespace std;

void makeAK8BtagCorrections(int year) {

	if(year < 2016 && year > 2018) throw std::invalid_argument("year must be either 2016, 2017, or 2018");

    vector<double> etabins = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.4};
    int nETA = etabins.size()-1;
    vector<double> ptbins  = {20,25,30,35,40,50,60,75,100,125,150,175,200,250,300,350,400,500,600,800,1000};
    int nPT = ptbins.size()-1;

	TString path = "/Users/brentstone/Dropbox/Physics/HHbbWW/btagging/";
	TFile *fin = TFile::Open(path+TString::Format("btagEffs%d.root",year));

	TString idS = "ak8_reqHbb";
	TString titleS = ";jet p_{T}[GeV];jet |#eta|";

    TString denN = "incl";
    vector<TString> numNs {"loose","med"};
    vector<TString> flvs {"l","c","b"};

    vector<TH2*> effs;

    auto getHist = [&](TString id, TString wp, TString name) {
    	TH2 *h = (TH2*)fin->Get(id+"_sj1_"+wp);
    	if(!h) throw std::invalid_argument("no hist found");
    	h->Add((TH2*)fin->Get(id+"_sj2_"+wp),1);
    	return PlotTools::rebin2D(h,name,titleS,nPT,&ptbins[0],nETA,&etabins[0]);
    };

    for(const auto& fl : flvs) {
    	TH2 *hd = getHist("bkg_noQCD_"+idS+"_"+fl,"incl","den");
    	TH2 *hd1 = getHist("radion_m1000_"+idS+"_"+fl,"incl","den_rad1000");
    	TH2 *hd3 = getHist("radion_m3000_"+idS+"_"+fl,"incl","den_rad3000");

//    	TH2 *hd = (TH2*)fin->Get("bkg_noQCD_ak8_"+idS+"_"+fl+"_incl_sj1");
//    	if(!hd) throw std::invalid_argument("no den found");
//    	hd->Add((TH2*)fin->Get("bkg_noQCD_ak8_"+idS+"_"+fl+"_incl_sj2"),1);
//    	hd = PlotTools::rebin2D(hd,"den",titleS,nPT,&ptbins[0],nETA,&etabins[0]);

    	for(const auto& num : numNs) {
    		TH2 *hn = getHist("bkg_noQCD_"+idS+"_"+fl,num,fl+"_"+num);
    		TH2 *hn1 = getHist("radion_m1000_"+idS+"_"+fl,num,fl+"_"+num+"_rad1000");
    		TH2 *hn3 = getHist("radion_m3000_"+idS+"_"+fl,num,fl+"_"+num+"_rad3000");

//        	TH2 *hn = (TH2*)fin->Get("bkg_noQCD_ak8_"+idS+"_"+fl+"_"+num+"_sj1");
//        	if(!hn) throw std::invalid_argument("no num found");
//        	hn->Add((TH2*)fin->Get("bkg_noQCD_ak8_"+idS+"_"+fl+"_"+num+"_sj2"),1);
//        	hn = PlotTools::rebin2D(hn,fl+"_"+num,titleS,nPT,&ptbins[0],nETA,&etabins[0]);

        	hn->Divide(hn,hd,1,1,"b");
        	hn1->Divide(hn1,hd1,1,1,"b");
        	hn3->Divide(hn3,hd3,1,1,"b");

        	effs.push_back(hn);
        	effs.push_back(hn1);
        	effs.push_back(hn3);

    	}
    	delete hd;
    	delete hd1;
    	delete hd3;
    }

    TFile *fout = new TFile(TString::Format("%s/ak8_deepcsvEff_%d.root",path.Data(),year),"recreate");
    fout->cd();
    for(const auto& h : effs) h->Write();
    fout->Close();
    delete fout;
    return;
}
