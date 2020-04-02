// eff and SF for single-lepton channel
{
TString pre = "/Users/brentstone/Dropbox/Physics/HHbbWW/trigger/master/";
vector<TString> yrs = {"2016","2017","2018"};
TString elS = "e_elpt30";
TString muS = "m_mupt27";

for(const auto& y : yrs) {
TFile *f = new TFile(pre+"effsAndSF_"+y+".root");
Plotter *pe = new Plotter(); Plotter *pm = new Plotter();
TH1 *dm = (TH1*)f->Get("effDA_"+muS);
TH1 *de = (TH1*)f->Get("effDA_"+elS);
TH1 *mm = (TH1*)f->Get("effMC_"+muS);
TH1 *me = (TH1*)f->Get("effMC_"+elS);

pe->addHist(me,"MC"); pe->addHist(de,"data");
pm->addHist(mm,"MC"); pm->addHist(dm,"data");

pe->setYTitle("Efficiency"); pm->setYTitle("Efficiency");
pe->setMinMax(0.66,1.04); pm->setMinMax(0.66,1.04);
pe->setBotMinMax(0.81,1.02); pm->setBotMinMax(0.81,1.02);
pe->setYTitleBot("SF"); pm->setYTitleBot("SF");
//pe->setLegendPos(0.76087,0.951505,0.501622,0.706684); pm->setLegendPos(0.7,0.8,0,0.3);
pm->drawSplitRatio(0,"muon "+y,false,false,"muon "+y);
pe->drawSplitRatio(0,"electron "+y,false,false,"electron "+y);

}
}
// --------------------------------------------------------------------------------------

// SF for di-lepton channel
{
TString pre = "/Users/brentstone/Dropbox/Physics/HHbbWW/trigger/master/";
vector<TString> yrs = {"2016","2017","2018"};
vector<TString> chns = {"mm_mupt27","me_mupt27","em_elpt30","ee_elpt30"};

for(const auto& y : yrs){
TFile *f = new TFile(pre+"effsAndSF_"+y+".root");

for(const auto& ch : chns){
TH2 *h = (TH2*)f->Get("sf_"+ch);

TCanvas *c = new TCanvas(TString(ch+y),TString(ch+y));
c->cd();
h->GetXaxis()->SetLimits(40,1000);
h->GetYaxis()->SetLimits(10,150);
h->Draw("colztext");
}

}
}
// --------------------------------------------------------------------------------------

// systematic for single-lepton channel
{
TString pre = "/Users/brentstone/Dropbox/Physics/HHbbWW/trigger/master/";
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
pe->setYTitle("Scale factor");
pe->setMinMax(0.72,1.04);
pe->setBotMinMax(0.965,1.035);
pe->setYTitleBot("SF / SF(30)");
pe->drawSplitRatio(0,"el"+y,0,0,"el"+y);

for(unsigned i=0; i<mSels.size(); ++i){
TH1 *h = (TH1*)f->Get("sf_m_"+mSels[i]);
pm->addHist(h,mStrs[i]);
}
pm->setYTitle("Scale factor");
pm->setMinMax(0.72,1.04);
pm->setBotMinMax(0.965,1.035);
pm->setYTitleBot("SF / SF(27)");
pm->drawSplitRatio(0,"mu"+y,0,0,"mu"+y);

}
}
// --------------------------------------------------------------------------------------

// systematic for di-lepton channel
{
TString pre = "/Users/brentstone/Dropbox/Physics/HHbbWW/trigger/master/";
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
//hn->GetXaxis()->SetLimits(401,990);
//hn->GetYaxis()->SetLimits(11,139);
hn->Draw("colztext");
}

}
}
// --------------------------------------------------------------------------------------
