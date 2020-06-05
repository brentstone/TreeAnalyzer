
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include <map>

using namespace std;

TH1* getTotBkg(TFile *f, TString hS) {
	vector<TString> bkgs = {"ttbar0","ttbar1","ttbar2","singlet","zjets","wjets","qcd","diboson","ttx","hx"};

	TH1 *h = (TH1*)f->Get(bkgs[0]+"_"+hS);
	if(!h) std::cout<<"whatup"<<std::endl;
	h = (TH1*)h->Clone("bkg_"+hS);
	for(unsigned i=1; i<bkgs.size(); i++) {
		h->Add((TH1*)f->Get(bkgs[i]+"_"+hS),1);
	}
	return h;
}

void doAK8YieldComp(TFile *f) {
	vector<TString> tfs = {"ll","lc","lb","cc","cb","bb"};
	vector<TString> sels = {""};

	for(const auto& s : sels) {
		printf("slurm\n");
		Plotter *p = new Plotter();
		TH1 *hd = (TH1*)f->Get("data_ak8_applIsoMet_"+s+"yield");
		if(!hd) printf("no data hist\n");
		p->addHist(hd,"data");
		for(const auto& tf : tfs) {
			TH1 *hb = getTotBkg(f,"ak8_applIsoMet_"+tf+"_"+s+"yield");
			if(!hb) printf("no mc hist\n");
			p->addStackHist(hb,tf);
		}
		printf("churm\n");
		p->drawSplitRatio(-1,"s"+s,false,false,"s"+s);
		printf("flurm\n");
	}

}

void doAK4YieldComp(TFile *f) {
	vector<TString> tfs = {"l","c","b"};

	Plotter *p = new Plotter();
	TH1 *hd = (TH1*)f->Get("data_ak4_yield");
	p->addHist(hd,"data");
	for(const auto& tf : tfs) {
	    TH1 *hb = getTotBkg(f,"ak4_"+tf+"_yield");
	    p->addStackHist(hb,tf);
	}
	p->drawSplitRatio(-1,"ak4",false,false,"ak4");
}

void drawDataMCJetDiffs() {

	TFile *f = new TFile("all.root");
	doAK8YieldComp(f);
	printf("s");
	doAK4YieldComp(f);
	printf("t");
}
