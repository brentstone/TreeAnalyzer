#include "TGraph.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "HistoPlotting/include/Plotter.h"

std::vector<double> masses = {800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500};
std::map<int,TString> binIds = { {1,"mu"}, {2,"mu_unc"}, {3,"dMean"}, {4,"dMedian"}, {5,"dStdDev"} };
std::map<int,TString> yTitles = { {1,"r_{bias} #mu"}, {3,"r_{bias} mean"}, {4,"r_{bias} median"}, {5,"r_{bias} std dev"} };

std::vector<TString> rexcls = {"r1","r2","r5"};
std::map<TString,TString> rStrs = { {"r1","r_{inj} = r_{excl}"}, {"r2","r_{inj} = 2r_{excl}"}, {"r5","r_{inj} = 5r_{excl}"} };

void doBiasPlots(TString dir) {

	TFile *f = TFile::Open(dir+"/plots/biasTest.root");
	for(const auto& id : binIds) {
		if(id.first == 2) continue;

		Plotter *p = new Plotter();

		for(const auto& r : rexcls) {
			TGraphErrors *g = new TGraphErrors();

			for(unsigned im=0; im<masses.size(); ++im) {
				const auto& m = masses[im];
				TH1 *h = (TH1*)f->Get(TString::Format("biasDist_postfit_%s_%d_quants",r.Data(),int(m)));
				if(!h) continue;

				double val = h->GetBinContent(id.first);
				g->SetPoint(im,m,val);
				if(id.first == 1) g->SetPointError(im,0,h->GetBinContent(id.first+1));

//				std::cout<<"m"<<m<<" - "<<r.Data()<<std::endl<<id.second<<" = "<<val<<std::endl<<std::endl;

			}

			p->addGraph(g,rStrs[r],-1,1,2,20,1,true,true,false,"L");
		}
		if(id.first == 5) p->setMinMax(0,2);
		else p->setMinMax(-2,2);
		p->setXTitle("M_{X} [GeV]");
		p->setYTitle(yTitles[id.first]);
		p->draw(false,id.second);
	}

}
