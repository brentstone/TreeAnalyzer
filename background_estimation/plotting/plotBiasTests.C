#include <vector>
#include <utility>
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TArrow.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "plotTestHelper.h"

std::vector<TObject*> writeables;
class Dummy {
public:
    Dummy(const std::string& outName = "") : outName(outName) {};
    ~Dummy() {
        if(outName.size()){
            TFile * f = new TFile((outName+".root").c_str(),"recreate");
            f->cd();
            for(auto * w : writeables){
                w->Write();
                w->Print((outName +"_"+w->GetTitle() +".pdf").c_str());
            }
            f->Close();
        }
    }
    std::string outName;
};

std::vector<std::pair<int,double>> massR1s_rad = {
		{800,0.0239258},
		{900,0.0139648},
		{1000,0.00925293},
		{1200,0.00522461},
		{1400,0.00355225},
		{1600,0.00269775},
		{1800,0.00214844},
		{2000,0.00176392},
		{2500,0.0012207},
		{3000,0.000952148},
		{3500,0.000836182}
};

std::vector<std::pair<int,double>> massR1s_blk = {
		{800,0.0162598},
		{900,0.00986328},
		{1000,0.00671387},
		{1200,0.00394287},
		{1400,0.002771},
		{1600,0.00216064},
		{1800,0.00174561},
		{2000,0.00145264},
		{2500,0.00102234},
		{3000,0.000799561},
		{3500,0.000692749}
};

std::vector<std::pair<int,double>> massR1s = {};
std::vector<std::pair<int,double>> massR2s = {};
std::vector<std::pair<int,double>> massR5s = {};

void doBiasTest(std::vector<TObject*>& writeables, const std::string& limitBaseName, int sig ){

	if(sig == 0) massR1s = massR1s_rad;
	else if(sig == 1) massR1s = massR1s_blk;
	else throw std::invalid_argument("sig should be either 0 (rad) or 1 (blk), dangus");

	for(const auto& th : massR1s) {
		massR2s.push_back(std::make_pair(th.first,2*th.second));
		massR5s.push_back(std::make_pair(th.first,5*th.second));
	}

	auto getInterpQuantile = [](std::vector<double>& vec, double quant) -> double {
		// assumes quant*vec.size() is an even integer
		return 0.5*( vec[int(quant*vec.size())] + vec[int(quant*vec.size())-1]);
	};

    auto makeBiasPlot = [&](const std::string& filename, const std::string& hName, const double rV, TFitResultPtr& fitres,
    						double& dMean, double& dMedian, double& dStdDev, double& dLow, double& dUp) ->TH1*  {
        TFile * f = new TFile(filename.c_str(),"read");
        if(!f) return 0;
        TTree * tree = 0;
        f->GetObject("tree_fit_sb",tree);
        if(!tree){return 0;}
        double fr=0;
        double frerr=0;
        int numbadnll=0;
        tree->SetBranchAddress("r",&fr);
        tree->SetBranchAddress("rErr",&frerr);
        tree->SetBranchAddress("numbadnll",&numbadnll);

        std::vector<double> biases;
        for(unsigned int iE = 0; tree->GetEntry(iE); ++iE ){
            if(numbadnll!=0) continue;
        	if(frerr == 0) continue;
            if(std::fabs((fr - rV)/frerr ) > 5) continue;
            biases.push_back( (fr - rV)/frerr   );
        }
        f->Close();

        std::sort(biases.begin(),biases.end(), [](const double a, const double b){return a < b;});
        printf("biases size = %d",int(biases.size()));

        dMedian = (biases.size() % 2 == 0) ? getInterpQuantile(biases,0.5) : biases[floor(0.5*biases.size())];
        dLow = getInterpQuantile(biases,0.16);
        dUp  = getInterpQuantile(biases,0.84);

        dMean = 0;
        dStdDev = 0;
        for(const auto& b : biases) dMean += b/biases.size();
        for(const auto& b : biases) dStdDev += (b-dMean)*(b-dMean)/biases.size();
        dStdDev = std::sqrt(dStdDev);

//        std::cout<<filename<<std::endl;
//        printf("distMean = %f\ndistMedian = %f\ndistStdDev = %f\n\n",dMean,dMedian,dStdDev);

        TH1 * h = new TH1D(hName.c_str(),";signal strength pull",10, -5,5);
        h->SetDirectory(0);
        TH1::AddDirectory(kFALSE);
        for(auto b : biases) h->Fill(b);
        auto c = new TCanvas(hName.c_str(),hName.c_str());
        h->Draw();
        fitres = h->Fit("gaus","S");;
        writeables.push_back(c);

        TH1D *hStats = new TH1D((hName+"_quants").c_str(),"",9,0,9);
        hStats->SetBinContent(1,fitres->GetParams()[1]);
        hStats->SetBinContent(2,fitres->GetErrors()[1]);
        hStats->SetBinContent(3,fitres->GetParams()[2]);
        hStats->SetBinContent(4,fitres->GetErrors()[2]);
        hStats->SetBinContent(5,dMean);
        hStats->SetBinContent(6,dMedian);
        hStats->SetBinContent(7,dStdDev);
        hStats->SetBinContent(8,dLow);
        hStats->SetBinContent(9,dUp);

        writeables.push_back(hStats);

//        std::sort(biases.begin(),biases.end(), [](const double a, const double b){return a < b;});
//        double nToys = biases.size();
//        h->SetDirectory(0);
//        TH1::AddDirectory(kFALSE);
//        for(auto b : biases) h->Fill(b);
//        std::cout << hName <<" "<< biases[nToys*0.5]<<std::endl;
//        return h;

        return h;

    };

    Plotter * p = new Plotter();
    Plotter * pw = new Plotter();

    std::vector<std::vector<double>> means(massR1s.size());
    std::vector<std::vector<double>> meanErs(massR1s.size());

    std::vector<std::vector<double>> distMeans(massR1s.size());
    std::vector<std::vector<double>> distMedians(massR1s.size());
    std::vector<std::vector<double>> distStdDevs(massR1s.size());

    auto doSet = [&] (const std::string& label,const std::string& plotlabel, const std::vector<std::pair<int,double>>& massRs ){
        TGraphErrors * gr = new TGraphErrors();
        TGraphErrors * grw = new TGraphErrors();
        int iP = 0;

        for(unsigned int iM = 0; iM < massRs.size(); ++iM){
            const auto& mr = massRs[iM];
            std::cout<<mr.first<<std::endl;
            TFitResultPtr fitres ;
            double dMean=0, dSD=0, dMed=0, dLow=0, dUp=0;
            auto dist = makeBiasPlot(limitBaseName + "/biasInput_"+label+"_m"+int2Str(mr.first)+".root",
                    "biasDist_"+label+"_"+int2Str(mr.first), mr.second,fitres,dMean,dMed,dSD,dLow,dUp);
            if(dist){
                auto pars = fitres->GetParams();
                auto errs = fitres->GetErrors();
                gr->SetPoint(iP,double(mr.first),pars[1]);
                gr->SetPointError(iP,0,errs[1]);
                grw->SetPoint(iP,double(mr.first),pars[2]);
                grw->SetPointError(iP,0,errs[2]);
                means[iM].push_back(pars[1]);
                meanErs[iM].push_back(errs[1]);

                distMeans[iM].push_back(dMean);
                distMedians[iM].push_back(dMed);
                distStdDevs[iM].push_back(dSD);
            }
            ++iP;
        }
        if(iP){
            p->addGraph(gr,plotlabel);
            pw->addGraph(grw,plotlabel);
        }
    };

    doSet("postfit_r1","postfit b-model: r=excluded"  ,massR1s);
    doSet("postfit_r2","postfit b-model: r=2*excluded",massR2s);
    doSet("postfit_r5","postfit b-model: r=5*excluded",massR5s);


    p->setYTitle("signal strength bias");
    p->setXTitle((sigMCS.title).c_str());
    auto c = p->draw(false,"signalInjectTest_bias");
    c->SetTitle("signalInjectTest_bias");
    writeables.push_back(c);

    pw->setYTitle("signal strength width");
    pw->setXTitle((sigMCS.title).c_str());
    auto cw = pw->draw(false,"signalInjectTest_width");
    cw->SetTitle("signalInjectTest_width");
    writeables.push_back(cw);

    std::cout << std::endl << std::endl;
    for(unsigned int iM = 0; iM < massR1s.size(); ++iM){
        std::cout << "mx = "<<int2Str(massR1s[iM].first) << std::endl;
        for(unsigned int iR = 0; iR < means[iM].size(); ++iR ){
            std::cout << TString::Format("rexcl %d\nfitMu = %f +- %f\n",iR+1,means[iM][iR],meanErs[iM][iR]);
            std::cout << TString::Format("dMean = %f\ndMedian = %f\ndStdDev = %f\n\n",distMeans[iM][iR],distMedians[iM][iR],
            		distStdDevs[iM][iR]);
        }

    }
    std::cout << std::endl << std::endl;
}

void plotBiasTests(std::string limDir, int sig) {

	std::string outName = limDir+"/plots2/biasTest";
	doBiasTest(writeables,limDir,sig);
	Dummy d(outName);
}
