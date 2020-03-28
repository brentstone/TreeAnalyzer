// systematic for di-lepton channel
{
TString pre = "/Users/brentstone/Dropbox/Physics/HHbbWW/trigger/turnons/";
vector<TString> yrs = {"2016","2017","2018"};
vector<TString> chns = {"mm_mupt27","me_mupt27","em_elpt30","ee_elpt30"};

for(const auto& y : yrs){
TFile *f = new TFile(pre+"effsAndSF_"+y+".root");

for(const auto& ch : chns){
TH2 *hn = (TH2*)f->Get("sf_"+ch+"Cross");
TH2 *hd = (TH2*)f->Get("sf_"+ch);
hn->Divide(hd);

TCanvas *c = new TCanvas(TString(ch+y),TString(ch+y));
c->cd();
hn->Draw("colztext");
}

}
}
// systematic for single-lepton channel
{
TString pre = "/Users/brentstone/Dropbox/Physics/HHbbWW/trigger/turnons/";
vector<TString> yrs = {"2016","2017","2018"};
vector<TString> eSels = {"elpt30","elpt25","elpt28","elpt32","elpt35"};
vector<TString> mSels = {"mupt27","mupt22","mupt25","mupt29","mupt32"};
vector<TString> eStrs = {"p_{T} > 30 GeV","p_{T} > 25 GeV","p_{T} > 28 GeV","p_{T} > 32 GeV","p_{T} > 35 GeV"};
vector<TString> mStrs = {"p_{T} > 27 GeV","p_{T} > 22 GeV","p_{T} > 25 GeV","p_{T} > 29 GeV","p_{T} > 32 GeV"};

for(const auto& y : yrs){
TFile *f = new TFile(pre+"effsAndSF_"+y+".root");

Plotter *pe = new Plotter(); Plotter *pm = new Plotter();

for(unsigned i=0; i<eSels.size(); ++i){
TH1 *h = (TH1*)f->Get("sf_e_"+eSels[i]);
pe->addHist(h,eStrs[i]);
}
pe->setMinMax(0.7,1.05);
pe->setBotMinMax(0.965,1.035);
pe->drawSplitRatio(0,"el"+y,0,0,"el"+y);

for(unsigned i=0; i<mSels.size(); ++i){
TH1 *h = (TH1*)f->Get("sf_m_"+mSels[i]);
pm->addHist(h,mStrs[i]);
}
pm->setMinMax(0.7,1.05);
pm->setBotMinMax(0.965,1.035);
pm->drawSplitRatio(0,"mu"+y,0,0,"mu"+y);

}
}
