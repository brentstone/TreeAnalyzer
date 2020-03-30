#include "HistoPlotting/include/PlotTools.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

using namespace std;

void makeAK4BtagCorrections(int year) {

	if(year < 2016 && year > 2018) throw std::invalid_argument("year must be either 2016, 2017, or 2018");

    vector<double> etabins = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.4};
    int nETA = etabins.size()-1;
    vector<double> ptbins  = {20,25,30,35,40,50,60,75,100,125,150,175,200,250,300,350,400,500,600,800,1000};
    int nPT = ptbins.size()-1;

	TString path = "/Users/brentstone/Dropbox/Physics/HHbbWW/btagging/";
	TFile *fin = TFile::Open(path+TString::Format("btagEffs_%d.root",year));

	TString idS = "reqHbb";
	TString titleS = ";jet p_{T}[GeV];jet |#eta|";

    TString denN = "incl";
    vector<TString> numNs {"loose","med","tight"};
    vector<TString> flvs {"l","c","b"};

    vector<TH2*> effs;
    for(const auto& fl : flvs) {
    	TH2 *hd = (TH2*)fin->Get("bkg_noQCD_"+idS+"_"+fl+"_incl");
    	if(!hd) throw std::invalid_argument("no den found");
    	hd = PlotTools::rebin2D(hd,"den",titleS,nPT,&ptbins[0],nETA,&etabins[0]);

    	for(const auto& num : numNs) {
        	TH2 *hn = (TH2*)fin->Get("bkg_noQCD_"+idS+"_"+fl+"_"+num);
        	if(!hn) throw std::invalid_argument("no num found");
        	hn = PlotTools::rebin2D(hn,fl+"_"+num,titleS,nPT,&ptbins[0],nETA,&etabins[0]);
        	hn->Divide(hn,hd,1,1,"b");
        	effs.push_back(hn);
    	}
    	delete hd;
    }

    TFile *fout = new TFile(TString::Format("%s/ak4_deepJetEff_%d.root",path.Data(),year),"recreate");
    fout->cd();
    for(const auto& h : effs) h->Write();
    fout->Close();
    delete fout;
    return;
}
