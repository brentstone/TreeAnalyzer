#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "HistoPlotting/include/Plotter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TF1.h"

#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cctype>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include "TPRegexp.h"

struct btagLine {
	int wp;
	TString name;
	TString syst;
	int flav;
	float etamin;
	float etamax;
	float ptmin;
	float ptmax;
	float discmin;
	float discmax;
	TString form;
};

class csvData {
public:
	void read_csv(std::string filename);
	void fillBtagLine(int i, btagLine& line, std::string& colVal);

	void drawSFs(bool savePlots);
	void debugBtagLine(btagLine&);
	std::vector<btagLine> csvLines;
};

void fillHistFromTF(TH1 *h, TF1 *f, float ptmin, float ptmax) {
	int lowbin = h->FindFixBin(ptmin+0.01);
	int highbin = h->FindFixBin(ptmax-0.01);

	for(int bin = lowbin; bin <= highbin; ++bin) {
		float low = h->GetXaxis()->GetBinLowEdge(bin);
		float high = h->GetXaxis()->GetBinUpEdge(bin);
		float cont = f->Integral(low,high) / (high - low);
		h->SetBinContent(bin,cont);
	}
}

void csvData::read_csv(std::string filename){
    // Reads a CSV file into a vector of <string, vector<int>> pairs where
    // each pair represents <column name, column values>

    // Create a vector of <string, int vector> pairs to store the result
    std::vector<std::pair<std::string, std::vector<int>>> result;

    std::ifstream myFile(filename);
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line, colname;
    int val;

    if(myFile.good()) {
        std::getline(myFile, line);
        printf("%s\n",line.c_str());

        std::stringstream ss(line);

    }

    // Read data, line by line
    while(std::getline(myFile, line))
    {
//        printf("%s\n",line.c_str());

        // Create a stringstream of the current line
        std::stringstream ss(line);

        btagLine bl;

        std::vector<TString> wp, name, syst, flav, etamin, etamax, ptmin, ptmax, discmin, discmax, form;
        std::string colVal;

        // Keep track of the current column index
        int colIdx = 0;
        while (std::getline(ss,colVal,',')) {
//            printf("%s\n",colVal.c_str());
            fillBtagLine(colIdx,bl,colVal);
            colIdx++;
        }

        csvLines.push_back(bl);
//        debugBtagLine(bl);

    }

    // Close file
    myFile.close();

    return;
}

void csvData::fillBtagLine(int i, btagLine& line, std::string& colVal) {
	switch (i) {
	case 0:
		line.wp = std::stoi(colVal); break;
	case 1:
		line.name = TString(colVal); line.name.ReplaceAll(" ",""); break;
	case 2:
		line.syst = TString(colVal); line.syst.ReplaceAll(" ",""); break;
	case 3:
		line.flav = std::stoi(colVal); break;
	case 4:
		line.etamin = std::stof(colVal); break;
	case 5:
		line.etamax = std::stof(colVal); break;
	case 6:
		line.ptmin = std::stof(colVal); break;
	case 7:
		line.ptmax = std::stof(colVal); break;
	case 8:
		line.discmin = std::stof(colVal); break;
	case 9:
		line.discmax = std::stof(colVal); break;
	case 10:
		line.form = TString(colVal); line.form.ReplaceAll(" ",""); line.form.ReplaceAll("\"",""); break;
	default:
		throw std::invalid_argument("shouldn't be reaching here");
	}
}

void csvData::drawSFs(bool savePlots) {
	std::vector<int> wps = {0,1,2};
	std::vector<TString> wpS = {"loose","med","tight"};

	std::vector<int> flavs = {0,2};
	std::vector<TString> names = {"comb","incl"};
	std::vector<TString> systs = {"central","up","down"};
	float etamin = 0;
	float etamax = 2.5;
	float ptmin = 20;
	float ptmax = 1000;
	float discmin = 0;
	float discmax = 1;

	HistGetter plotter;
	std::vector<float> bins = {20,30,40,50,60,70,80,90,100,120,140,170,200,250,300,350,400,500,600,700,800,1000};

	for(const auto& sys : systs) for(const auto& fl : flavs) for(const auto& nm : names) {
		if(nm == "comb" && fl != 0) continue;
		if(nm == "incl" && fl != 2) continue;

		Plotter *p = new Plotter();

		for(const auto& wp : wps) {
			TString tit = TString::Format("h_%s_%s_%s",sys.Data(),wpS[wp].Data(),nm.Data());

			TH1F *h = new TH1F(tit,tit,bins.size()-1,&bins[0]);

			for(unsigned i=0; i<csvLines.size(); i++) {
				auto& line = csvLines[i];
				if(line.discmin < discmin) continue;
				if(line.discmax < discmax) continue;
				if(line.ptmin < ptmin) continue;
				if(line.ptmax > ptmax) continue;
				if(fabs(line.etamin) < etamin) continue;
				if(line.etamax > etamax) continue;

				if(line.wp != wp) continue;
				if(line.flav != fl) continue;
				if(line.name != nm) continue;
				if(line.syst != sys) continue;

//				std::cout<<line.ptmin<<" - "<<line.ptmax<<std::endl;
//				std::cout<<line.form<<std::endl<<std::endl;
				TString name = TString::Format("f%d",i);
				TF1 *f = new TF1(name,line.form.Data(),line.ptmin,line.ptmax);

				fillHistFromTF(h,f,line.ptmin,line.ptmax);
			}

			if(savePlots) plotter.add1D(h);
			else p->addHistLine(h,wpS[wp]);

		}
		if(savePlots) plotter.write("visBtagSFs.root");
		else p->draw(false,TString::Format("%s_%s_fl%d",nm.Data(),sys.Data(),fl));
	}


}

void csvData::debugBtagLine(btagLine& bl) {
	printf("btagLine:\nwp = %d\nname = %s\nsyst = %s\nflav = %d\nbin = (%.2f,%.2f,%.2f,%.2f,%.2f,%.2f)\nform = %s\n\n",
			bl.wp,bl.name.Data(),bl.syst.Data(),bl.flav,bl.etamin,bl.etamax,bl.ptmin,bl.ptmax,bl.discmin,bl.discmax,bl.form.Data());
}

#endif

void visBtagSFs(bool savePlots = false) {
	csvData data;
	data.read_csv("/Users/brentstone/Dropbox/Physics/HHbbWW/data/corrections/btagging/DeepJet_102XSF_WP_V1.csv");
	data.drawSFs(savePlots);

	std::cout<<data.csvLines.size()<<std::endl;

}

