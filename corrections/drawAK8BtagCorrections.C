#if !defined(__CINT__) || defined(__MAKECINT__)

#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include <map>

using namespace std;

class BtagPlotter {
public:
	BtagPlotter();
	void doYearComparison(TString idS);
	void doOtherComparison(TString idS);
	void doSelComparison(int year);
	void doKinComparison(int year, TString idS);
	void doComparisonWithSignal(int year, TString idS);
	bool doSlices = true;

private:
	TString path = "/Users/brentstone/Dropbox/Physics/HHbbWW/btagging/";
	TString titleS = ";jet p_{T}[GeV];jet |#eta|";
	TString denN = "incl";
	vector<TString> numNs {"loose","med"};
	vector<TString> flvs {"b","l"};
	vector<double> etabins = {0,0.2,0.4,0.6,0.8,1.2,1.6,2.4};
    vector<double> ptbins  = {20,25,30,35,40,50,60,75,100,125,150,175,200,250,300,350,400,500,600,800,1000};
	int nETA;
	int nPT;
	vector<TString> years = {"2016","2017","2018"};
    vector<TString> ids   = {"noHbbReq","reqHbb","goodMSD"};
    vector<TString> samps = {"bkg_noQCD","bkg","radion_m1000","radion_m3000"};

	TH1 * getProjectionEff(TH2 *num, TH2 *den, TString name, int firstbin, int lastbin, bool projX);
	TH2 * getTotHist(TFile *fin, TString id, TString wp);
};

TH1 * BtagPlotter::getProjectionEff(TH2 *num, TH2 *den, TString name, int firstbin, int lastbin, bool projX) {
	TH2 *hn = (TH2*)num->Clone();
	TH2 *hd = (TH2*)den->Clone();
	hd = PlotTools::rebin2D(hd,"den2D",titleS,nPT,&ptbins[0],nETA,&etabins[0]);
	hn = PlotTools::rebin2D(hn,"num2D",titleS,nPT,&ptbins[0],nETA,&etabins[0]);

	TH1 *nh=0, *nd=0;
	if(projX) {
		nh = (TH1*)hn->ProjectionX(name,firstbin,lastbin,"e");
		nd = (TH1*)hd->ProjectionX("den1D",firstbin,lastbin,"e");
	} else {
		nh = (TH1*)hn->ProjectionY(name,firstbin,lastbin,"e");
		nd = (TH1*)hd->ProjectionY("den1D",firstbin,lastbin,"e");
	}

	delete hn;
	delete hd;
	nh->Divide(nh,nd,1,1,"b");
	return nh;
}

BtagPlotter::BtagPlotter() {
	nETA = etabins.size()-1;
	nPT = ptbins.size()-1;
}

TH2 * BtagPlotter::getTotHist(TFile *fin, TString id, TString wp) {
	TH2 *h = (TH2*)fin->Get(id+"_sj1_"+wp);
	if(!h) {
		std::cout<<id+"_sj1_"+wp<<std::endl;
		throw std::invalid_argument("no hist found");
	}
	h->Add((TH2*)fin->Get(id+"_sj2_"+wp),1);
	return PlotTools::rebin2D(h,id+"_"+wp,titleS,nPT,&ptbins[0],nETA,&etabins[0]);
}

void BtagPlotter::doYearComparison(TString idS) {

    for(const auto& fl : flvs) for(const auto& num : numNs) {
    	// get slices in eta
    	if(doSlices) {
        	for(int ieta=1; ieta<=nETA; ++ieta) {
        		TString etaS = TString::Format("eta%.1fto%.1f",etabins[ieta-1],etabins[ieta]);
        		etaS.ReplaceAll(".","p");

        		Plotter *p = new Plotter();
        		for(const auto& yr : years) {
        			TFile *fin = TFile::Open(path+TString::Format("btagEffs%s.root",yr.Data()));
        	    	TH2 *hd = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,"incl");
        	    	TH2 *hn = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,num);
        	    	if(!hd || !hn) throw std::invalid_argument("missing hist");

        	    	TH1 *eff = getProjectionEff(hn,hd,fl+"_"+num,ieta,ieta,true);
        	    	p->addHist(eff,yr);

        		}
        		p->draw(false,fl+" "+num+" "+etaS);
        	}
    	}

    	// 1d full projections
		Plotter *p = new Plotter();
		Plotter *pe = new Plotter();
		for(const auto& yr : years) {
			TFile *fin = TFile::Open(path+TString::Format("btagEffs%s.root",yr.Data()));
	    	TH2 *hd = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,"incl");
	    	TH2 *hn = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,num);
	    	if(!hd || !hn) throw std::invalid_argument("missing hist");

	    	TH1 *effPT = getProjectionEff(hn,hd,fl+"_pt_"+num,1,nETA,true);
	    	TH1 *effETA = getProjectionEff(hn,hd,fl+"_eta_"+num,1,nPT,false);

	    	p->addHist(effPT,yr);
	    	pe->addHist(effETA,yr);
		}
		p->draw(false,fl+" "+num+" pt");
		pe->draw(false,fl+" "+num+" eta");
    }
}

void BtagPlotter::doSelComparison(int year) {

	TFile *fin = TFile::Open(path+TString::Format("btagEffs_%d.root",year));

    for(const auto& fl : flvs) for(const auto& num : numNs) {

    	// get slices in eta
    	if(doSlices) {
        	for(int ieta=1; ieta<=nETA; ++ieta) {
        		TString etaS = TString::Format("eta%.1fto%.1f",etabins[ieta-1],etabins[ieta]);
        		etaS.ReplaceAll(".","p");

        		Plotter *p = new Plotter();
        		for(const auto& idS : ids) {
        	    	TH2 *hd = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,"incl");
        	    	TH2 *hn = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,num);
        	    	if(!hd || !hn) throw std::invalid_argument("missing hist");

        	    	TH1 *eff = getProjectionEff(hn,hd,fl+"_"+num+"_"+etaS,ieta,ieta,true);
    //    	    	p->addHist(eff,idS,-1,1,4,20,1,true,false);
        	    	p->addHist(eff,idS);

        		}
        		p->draw(false,fl+" "+num+" "+etaS);
        	}
    	}

    	// 1d full projections
		Plotter *p = new Plotter();
		Plotter *pe = new Plotter();
		for(const auto& idS : ids) {
	    	TH2 *hd = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,"incl");
	    	TH2 *hn = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,num);
	    	if(!hd || !hn) throw std::invalid_argument("missing hist");

	    	TH1 *effPT = getProjectionEff(hn,hd,fl+"_pt_"+num,1,nETA,true);
	    	TH1 *effETA = getProjectionEff(hn,hd,fl+"_eta_"+num,1,nPT,false);

	    	p->addHist(effPT,idS);
	    	pe->addHist(effETA,idS);
		}
		p->draw(false,fl+" "+num+" pt");
		pe->draw(false,fl+" "+num+" eta");
    }
}

void BtagPlotter::doKinComparison(int year, TString idS) {

	TFile *fin = TFile::Open(path+TString::Format("btagEffs%d.root",year));

	for(const auto& fl : flvs) for(const auto& num : numNs) {

    	TH2 *hd = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,"incl");
    	TH2 *hn = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,num);
		if(!hd || !hn) throw std::invalid_argument("missing hist");

		Plotter *pp = new Plotter();
		Plotter *pe = new Plotter();

    	TH1 *effPT = getProjectionEff(hn,hd,fl+"_pt_"+num,1,nETA,true);
    	TH1 *effETA = getProjectionEff(hn,hd,fl+"_eta_"+num,1,nPT,false);

    	pp->addHist(effPT,"incl");
    	pe->addHist(effETA,"incl");

    	// slices of eta
    	for(int ieta=1; ieta<=nETA; ++ieta) {
    		TString etaS = TString::Format("%.1f < |#eta| < %.1f",etabins[ieta-1],etabins[ieta]);
    	    TH1 *eff = getProjectionEff(hn,hd,fl+"_"+num+"_"+etaS,ieta,ieta,true);
    	    pp->addHist(eff,etaS);
    	}

    	//slices of pt
    	vector<pair<pair<float,float>,TString>> ptsels = { {{20,50},"20 < p_{T} < 50"}, {{50,100},"50 < p_{T} < 100"},
    			{{100,1000},"100 < p_{T} < 1000"} };

    	for(const auto& sel : ptsels) {
    		TString str = sel.second;
    		int lowbin = hn->GetXaxis()->FindFixBin(sel.first.first + 0.1);
    		int highbin = hn->GetXaxis()->FindFixBin(sel.first.second - 0.1);

    	    TH1 *eff = getProjectionEff(hn,hd,fl+"_"+num+"_"+str,lowbin,highbin,false);
    	    pe->addHist(eff,str);
    	}

		pp->draw(false,fl+" "+num+" pt");
		pe->draw(false,fl+" "+num+" eta");
    }
}

void BtagPlotter::doComparisonWithSignal(int year, TString idS) {
	TFile *fin = TFile::Open(path+TString::Format("btagEffs%d.root",year));

    for(const auto& fl : flvs) for(const auto& num : numNs) {

    	// get slices in eta
    	if(doSlices) {
        	for(int ieta=1; ieta<=nETA; ++ieta) {
        		TString etaS = TString::Format("eta%.1fto%.1f",etabins[ieta-1],etabins[ieta]);
        		etaS.ReplaceAll(".","p");

        		Plotter *p = new Plotter();
        		for(const auto& smp : samps) {
        	    	TH2 *hd = getTotHist(fin,smp+"_"+idS+"_"+fl,"incl");
        	    	TH2 *hn = getTotHist(fin,smp+"_"+idS+"_"+fl,num);
        	    	if(!hd || !hn) throw std::invalid_argument("missing hist");

        	    	TH1 *eff = getProjectionEff(hn,hd,fl+"_"+num+"_"+etaS,ieta,ieta,true);
        	    	p->addHist(eff,idS);

        		}
        		p->draw(false,fl+" "+num+" "+etaS);
        	}
    	}

    	// 1d full projections
		Plotter *p = new Plotter();
		Plotter *pe = new Plotter();
		for(const auto& smp : samps) {
	    	TH2 *hd = getTotHist(fin,smp+"_"+idS+"_"+fl,"incl");
	    	TH2 *hn = getTotHist(fin,smp+"_"+idS+"_"+fl,num);
	    	if(!hd || !hn) throw std::invalid_argument("missing hist");

	    	TH1 *effPT = getProjectionEff(hn,hd,fl+"_pt_"+num,1,nETA,true);
	    	TH1 *effETA = getProjectionEff(hn,hd,fl+"_eta_"+num,1,nPT,false);

//	    	delete hd;
//	    	delete hn;
//	    	printf("slurm\n");

	    	p->addHist(effPT,smp);
	    	pe->addHist(effETA,smp);
		}
		p->draw(false,fl+" "+num+" pt");
		pe->draw(false,fl+" "+num+" eta");
    }
}

void BtagPlotter::doOtherComparison(TString idS) {

	std::vector<std::pair<TString,TString>> files = {{"btagEffs2018.root","new"},{"btagEffs2018_old.root","old"}};
    for(const auto& fl : flvs) for(const auto& num : numNs) {
    	// get slices in eta
    	if(doSlices) {
        	for(int ieta=1; ieta<=nETA; ++ieta) {
        		TString etaS = TString::Format("eta%.1fto%.1f",etabins[ieta-1],etabins[ieta]);
        		etaS.ReplaceAll(".","p");

        		Plotter *p = new Plotter();
        		for(const auto& f : files) {
        			TFile *fin = TFile::Open(path+f.first);
        	    	TH2 *hd = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,"incl");
        	    	TH2 *hn = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,num);
        	    	if(!hd || !hn) throw std::invalid_argument("missing hist");

        	    	TH1 *eff = getProjectionEff(hn,hd,fl+"_"+num,ieta,ieta,true);
        	    	p->addHist(eff,f.second);

        		}
        		p->draw(false,fl+" "+num+" "+etaS);
        	}
    	}

    	// 1d full projections
		Plotter *p = new Plotter();
		Plotter *pe = new Plotter();
		for(const auto& f : files) {
			TFile *fin = TFile::Open(path+f.first);
	    	TH2 *hd = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,"incl");
	    	TH2 *hn = getTotHist(fin,"bkg_noQCD_"+idS+"_"+fl,num);
	    	if(!hd || !hn) throw std::invalid_argument("missing hist");

	    	TH1 *effPT = getProjectionEff(hn,hd,fl+"_pt_"+num,1,nETA,true);
	    	TH1 *effETA = getProjectionEff(hn,hd,fl+"_eta_"+num,1,nPT,false);

	    	p->addHist(effPT,f.second);
	    	pe->addHist(effETA,f.second);
		}
		p->draw(false,fl+" "+num+" pt");
		pe->draw(false,fl+" "+num+" eta");
    }
}

#endif

void drawAK8BtagCorrections(int option, int year=0, bool doSlices = false) {

	TString idS = "ak8_reqHbb";
	BtagPlotter *plotter = new BtagPlotter();
	plotter->doSlices = doSlices;

	if(option == 0) plotter->doYearComparison(idS);
	else if(option == 1) plotter->doSelComparison(year);
	else if(option == 2) plotter->doKinComparison(year,idS);
	else if(option == 3) plotter->doComparisonWithSignal(year,idS);
	else if(option == 4) plotter->doOtherComparison(idS);

}
