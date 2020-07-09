
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

	auto setSJBinLabels = [&](Plotter *p) {
		p->xAxis()->SetBinLabel(1,"NN");
		p->xAxis()->SetBinLabel(2,"LN");
		p->xAxis()->SetBinLabel(3,"LL");
		p->xAxis()->SetBinLabel(4,"MN");
		p->xAxis()->SetBinLabel(5,"ML");
		p->xAxis()->SetBinLabel(6,"MM");
	};

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
//		if(bin.first != 30) continue;
//		if(bin.second != 210) continue;
		std::cout<<"bin "<<bin.first<<" to "<<bin.second<<std::endl;
		std::cout<<std::endl<<"Data:"<<std::endl;

		Plotter *pp = new Plotter();
		TH1 *d = (TH1*)hdd->ProjectionY("data",hdd->GetXaxis()->FindFixBin(bin.first+0.01),
				hdd->GetXaxis()->FindFixBin(bin.second-0.01));
		float loo_d = d->GetBinContent(4);
		float med_d = d->GetBinContent(5);
		float tgt_d = d->GetBinContent(6);
		float nt_d = d->GetBinContent(1);
		float sr_d = d->Integral(4,6);
		float tot_d = d->Integral(1,6);

		std::cout<<"Loose = "<<loo_d<<std::endl;
		std::cout<<"Med = "<<med_d<<std::endl;
		std::cout<<"Tight = "<<tgt_d<<std::endl;
		std::cout<<"SR = "<<sr_d<<std::endl;
		std::cout<<"NonTopCR = "<<nt_d<<std::endl;

		pp->addHist(d,"data");
//		std::cout<<"data bins = "<<d->GetNbinsX()<<std::endl;
		for(const auto& tf : tfs) {
			TH1 *hb = getTotBkg(f,"ak8_"+tf+"_mj_x_yield",true,true,bin.first,bin.second);
//			std::cout<<tf<<" bins = "<<hb->GetNbinsX()<<std::endl;

			if(!hb) printf("no mc hist\n");
			pp->addStackHist(hb,tf);
		}
		TH1 *b = getTotBkg(f,"ak8_mj_x_yield",true,true,bin.first,bin.second);
		std::cout<<std::endl<<"MC:"<<std::endl;
		float loo_b = b->GetBinContent(4);
		float med_b = b->GetBinContent(5);
		float tgt_b = b->GetBinContent(6);
		float nt_b = b->GetBinContent(1);

		float sr_b = b->Integral(4,6);
		float tot_b = b->Integral(1,6);

		std::cout<<"Loose = "<<loo_b<<std::endl;
		std::cout<<"Med = "<<med_b<<std::endl;
		std::cout<<"Tight = "<<tgt_b<<std::endl;
		std::cout<<"SR = "<<sr_b<<std::endl;
		std::cout<<"NonTopCR = "<<nt_b<<std::endl;

		std::cout<<std::endl<<"SR data/MC = "<<((sr_d)/(sr_b))<<std::endl;
		std::cout<<"NonTopCR data/MC = "<<(nt_d/(nt_b))<<std::endl;
		std::cout<<"Norm CORR:"<<std::endl;
		std::cout<<"SR data/MC = "<<((sr_d/tot_d)/(sr_b/tot_b))<<std::endl;
		std::cout<<"NonTopCR data/MC = "<<((nt_d/tot_d)/(nt_b/tot_b))<<std::endl;

//		std::cout<<"jerrold"<<std::endl;
		pp->normalize();
		pp->setYTitle("a.u.");
		pp->setXTitle("sj b-tagging category");
		pp->setYTitleBot("data/bkg");
		pp->drawSplitRatio(-1,TString::Format("ak8_%.0fto%.0f",bin.first,bin.second),false,false,
				TString::Format("ak8_%.0fto%.0f",bin.first,bin.second));
//		setSJBinLabels(pp);

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

void drawDataMCJetDiffs(int yr) {

	TString ys = TString::Format("%d",yr); ys.ReplaceAll("20","");
	TFile *f = new TFile("/Users/brentstone/Dropbox/Physics/HHbbWW/plots/studyDataMCDiffs/dataMCDiffs"+ys+".root");
	doAK8YieldComp(f);
	doAK4YieldComp(f);
}
