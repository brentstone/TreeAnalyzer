
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include <map>

using namespace std;

TH1* getTotBkg(TFile *f, TString hS, bool doSlice, bool binInX=true, float low=0, float high=0) {
	vector<TString> bkgs = {"ttbar0","ttbar1","ttbar2","singlet","zjets","wjets","qcd","diboson","ttx","hx"};

	TH1 *h = 0;
//	printf("dbg1\n");
	if(!doSlice) {
//		printf("dbg1.1\n");
		h = (TH1*)f->Get(bkgs[0]+"_"+hS);
//		if(!h) std::cout<<"whatup"<<std::endl;
		h = (TH1*)h->Clone("bkg_"+hS);
		for(unsigned i=1; i<bkgs.size(); i++) {
//			printf("dbg1.2\n");
			h->Add((TH1*)f->Get(bkgs[i]+"_"+hS),1);
		}
	} else {
//		printf("dbg2.1\n");
		TH2 *hh = (TH2*)f->Get(bkgs[0]+"_"+hS);
//		std::cout<<bkgs[0]+"_"+hS<<std::endl;
//		if(!hh) std::cout<<"whatup"<<std::endl;
//		printf("dbg2.19\n");

		if(binInX) h = (TH1*)hh->ProjectionY("bkg",hh->GetXaxis()->FindFixBin(low),hh->GetXaxis()->FindFixBin(high));
		else       h = (TH1*)hh->ProjectionX("bkg",hh->GetYaxis()->FindFixBin(low),hh->GetYaxis()->FindFixBin(high));

//		printf("dbg2.2\n");

		delete hh;
//		printf("dbg2.3\n");

		for(unsigned i=1; i<bkgs.size(); i++) {
//			printf("dbg2.3.1\n");

			TH2 *hhh = (TH2*)f->Get(bkgs[i]+"_"+hS);
			h->Add((TH1*)hhh->ProjectionY(bkgs[i],hhh->GetXaxis()->FindFixBin(low),hhh->GetXaxis()->FindFixBin(high)),1);
		}
	}

	return h;
}

void doAK8YieldComp(TFile *f) {
	vector<TString> tfs = {"ll","lc","lb","cc","cb","bb"};
	vector<pair<float,float>> bins = { {0,250}, {30,210}, {30,60}, {60,90}, {30,90}, {150,210} };

	TH2 *hdd = (TH2*)f->Get("data_ak8_mj_x_yield");

//	printf("slurm\n");
//	Plotter *p = new Plotter();
//	TH1 *hd = (TH1*)f->Get("data_ak8_yield");
//	if(!hd) printf("no data hist\n");
//	p->addHist(hd,"data");
//	for(const auto& tf : tfs) {
//		printf("churm\n");
//		TH1 *hb = getTotBkg(f,"ak8_"+tf+"_yield",false);
//		if(!hb) printf("no mc hist\n");
//		p->addStackHist(hb,tf);
//	}
//	printf("churm2\n");
//	p->drawSplitRatio(-1,"ak8",false,false,"ak8");
//	printf("flurm\n");

	Plotter *pr = new Plotter();
	for(const auto& bin : bins) {
		Plotter *pp = new Plotter();
		TH1 *d = (TH1*)hdd->ProjectionY("data",hdd->GetXaxis()->FindFixBin(bin.first+0.01),
				hdd->GetXaxis()->FindFixBin(bin.second-0.01));
		pp->addHist(d,"data");
//		std::cout<<"data bins = "<<d->GetNbinsX()<<std::endl;
		for(const auto& tf : tfs) {
			TH1 *hb = getTotBkg(f,"ak8_"+tf+"_mj_x_yield",true,true,bin.first,bin.second);
//			std::cout<<tf<<" bins = "<<hb->GetNbinsX()<<std::endl;

			if(!hb) printf("no mc hist\n");
			pp->addStackHist(hb,tf);
		}
		TH1 *b = getTotBkg(f,"ak8_mj_x_yield",true,true,bin.first,bin.second);

//		std::cout<<"jerrold"<<std::endl;
//		pp->normalize();
		pp->drawSplitRatio(-1,TString::Format("ak8_%.0fto%.0f",bin.first,bin.second),false,false,
				TString::Format("ak8_%.0fto%.0f",bin.first,bin.second));
		d->Scale(1./d->Integral(0,10000));
		b->Scale(1./b->Integral(0,10000));
		d->Divide(b);
		pr->addHist(d,TString::Format("%.0fto%.0f",bin.first,bin.second));

//		std::cout<<"is a dangus"<<std::endl;
	}
	pr->draw(false,"normalized data/MC ratios");

}

void doAK4YieldComp(TFile *f) {
	vector<TString> tfs = {"l","c","b"};

	Plotter *p = new Plotter();
	TH1 *hd = (TH1*)f->Get("data_ak4_yield");
	p->addHist(hd,"data");
	for(const auto& tf : tfs) {
	    TH1 *hb = getTotBkg(f,"ak4_"+tf+"_yield",false);
	    p->addStackHist(hb,tf);
	}
	p->drawSplitRatio(-1,"ak4",false,false,"ak4");
}

void drawDataMCJetDiffs() {

	TFile *f = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/plots/studyDataMCDiffs/dataMCDiffs18.root");
	doAK8YieldComp(f);
	doAK4YieldComp(f);
}
