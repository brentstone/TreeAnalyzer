#include "TH2.h"
#include "TH1.h"
#include "TFile.h"
#include "TArrow.h"
#include "HistoPlotting/include/Plotter.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"

void showQuantiles(std::vector<double>& sortedToys, double dataTS) {
    auto getQuant = [&](double quant) {
        return sortedToys[int(quant*sortedToys.size())];
    };

    auto getMeanAndPval = [&]()->std::pair<double,double> {
        double cnt = 0;
	double aboveDataTS = 0;
	for(const auto& val : sortedToys) {cnt += val; if(val > dataTS) aboveDataTS++;}
	return std::make_pair(cnt/sortedToys.size(),aboveDataTS/sortedToys.size());
    };

    auto meanAndpval = getMeanAndPval();

    std::cout<<"mean     = "<<meanAndpval.first<<std::endl;
    std::cout<<"+2 sigma = "<<getQuant(0.95)<<std::endl;
    std::cout<<"+1 sigma = "<<getQuant(0.841)<<std::endl;
    std::cout<<"median   = "<<getQuant(0.5)<<std::endl;
    std::cout<<"-1 sigma = "<<getQuant(0.159)<<std::endl;
    std::cout<<"-2 sigma = "<<getQuant(0.05)<<std::endl;
    std::cout<<std::endl<<"pval = "<<meanAndpval.second<<std::endl<<std::endl;

}

double getSatTS(double nM, double nD) {
    return (nM + (nD ? (nD * std::log(nD/nM) - nD) : 0.0) );
}

std::pair<double,TH2*> calcSatTS(TH2* mod, TH2 *data, TString name) {
    if(mod->GetNbinsX() != data->GetNbinsX()) throw std::invalid_argument("unequal nbinsX!");
    if(mod->GetNbinsY() != data->GetNbinsY()) throw std::invalid_argument("unequal nbinsY!");
    
    double totTS = 0;

    TH2 *hTS = (TH2*)mod->Clone(name);

    for(unsigned iX = 1; iX<=mod->GetNbinsX(); ++iX) for(unsigned iY=1;iY<=mod->GetNbinsY();++iY) {
        double contMod = mod->GetBinContent(iX,iY);
        double contDat = data->GetBinContent(iX,iY);

        hTS->SetBinContent(iX,iY,0.0);
        hTS->SetBinError(iX,iY,0.0);

        double binLowX = mod->GetXaxis()->GetBinLowEdge(iX);
        double binLowY = mod->GetYaxis()->GetBinLowEdge(iY);

	// hack for testing a subrange
//    if(contMod < 3) continue;
//    if( binLowX >= 150 || binLowX <= 100) continue;
//    if(binLowY >= 3100 && binLowY <= 3400 ) continue;

        double ts = getSatTS(contMod,contDat);
        hTS->SetBinContent(iX,iY,2*ts);

        totTS += ts;
    }

    return std::make_pair(2*totTS,hTS);
}

void doSatGoFbyHand(TString postFitDir, bool doNonTop, std::vector<TString> exclSels = {}) {

HistGetter plotter;
std::vector<TString> selsTop = {"e_L_LP","e_T_LP","e_L_HP","e_T_HP","mu_L_LP","mu_T_LP","mu_L_HP","mu_T_HP",
				"SF_L","SF_T","OF_L","OF_T"};
std::vector<TString> selsNonTop = {"e_L_LP","e_L_HP","mu_L_LP","mu_L_HP","SF_L","OF_L"};

int nToys = 100;

TFile *f = TFile::Open(postFitDir+"/postFit.root");

double totTS = 0;
std::vector<double> toys_totTS(nToys,0.0);
std::vector<TH2*> toysHist_totTS(nToys,0);

TH2 *obsTS = 0;

for(const auto& sel : (doNonTop ? selsNonTop : selsTop)) {
    
    bool skip = false;
    for(const auto& s : exclSels) {if(sel == s) skip = true;}
    if(skip) continue;

    TH2 * data = (TH2*)f->Get("data_"+sel+"_full__MJ_MR");
    TH2 * mod  = (TH2*)f->Get("postfit_"+sel+"_full");

    auto tsinfo = calcSatTS(mod,data,sel);
    totTS += tsinfo.first;
    if(!obsTS) obsTS = (TH2*)tsinfo.second->Clone("obsTS");
    else       obsTS->Add(tsinfo.second,1);

    for(unsigned i=1; i<=nToys; ++i) {
        TString iS = TString::Format("%d",i);
        TH2 *toyData = (TH2*)f->Get("toyData_"+iS+"_"+sel+"_full__MJ_MR");
        TH2 *toyMod  = (TH2*)f->Get("toyDataFit_"+iS+"_"+sel+"_full");

        auto toyinfo = calcSatTS(toyMod,toyData,"toy"+iS+"_"+sel);
        toys_totTS[i-1] += toyinfo.first;
        if(!toysHist_totTS[i-1]) toysHist_totTS[i-1] = (TH2*)toyinfo.second->Clone(TString::Format("toyTS_%d",i));
        else                     toysHist_totTS[i-1]->Add(toyinfo.second,1);
    }
}

std::cout<<std::endl<<"data ts = "<<totTS<<std::endl;
std::sort(toys_totTS.begin(), toys_totTS.end(), [](const double a, const double b) { return a < b;});

TH1D *htoys = new TH1D("htoys","htoys",40,toys_totTS.front()-100.0,toys_totTS.back()+100.0);
for(const auto& val : toys_totTS) htoys->Fill(val);
Plotter *p = new Plotter();
p->addStackHist(htoys,"toys");

std::cout<<"Toys:"<<std::endl;
showQuantiles(toys_totTS,totTS);

p->setYTitle("N_{toys}");
p->setXTitle("GoF statistic");
TCanvas *c = p->draw(false,"toyDist");
TArrow *arr = new TArrow(totTS,10,totTS,0);
arr->SetLineWidth(2);
arr->SetLineColor(kRed);
arr->Draw();

plotter.add2D(obsTS);
for(int i=0;i<nToys;++i) plotter.add2D(toysHist_totTS[i]);
plotter.write("outSatGoF.root");

}
