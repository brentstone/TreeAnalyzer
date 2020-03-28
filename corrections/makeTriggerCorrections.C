#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

using namespace std;

class TriggerEffProcessor {
public:
	TriggerEffProcessor(int era, int ept, int mpt) : era(ERA(era)), elpt(ept), mupt(mpt) {
		setFiles();
		elPreS = eraStrs[era]+"_el_tt2_passSMu";
		muPreS = eraStrs[era]+"_mu_tt2_PassSEl";
	}

	void setHTBinning(vector<double>& bins) {htbins = bins;}
	void setPTBinning(vector<double>& bins) {ptbins = bins;}
	void setPtThresh(int ept, int mpt) {
		elpt = ept;
		mupt = mpt;
	}

	vector<TH1*> getEffsAndSFs_lnuqq(bool onlyForFile);
	vector<TH2*> getEffsAndSFs_lnulnu(bool includeCross);
	vector<TH1*> getEffsForFile_2l(bool includeCross);

private:
	enum ERA {E_2016, E_2017, E_2018, E_RUN2};
	vector<TString> eraStrs = {"2016","2017","2018","Run2"};
	ERA era;
	TString elPreS;
	TString muPreS;
	int elpt;
	int mupt;
	vector<double> htbins;
	vector<double> ptbins;

	TH1 *hen_mc = 0;
	TH1 *hed_mc = 0;
	TH1 *hmn_mc = 0;
	TH1 *hmd_mc = 0;
	TH1 *hen_da = 0;
	TH1 *hed_da = 0;
	TH1 *hmn_da = 0;
	TH1 *hmd_da = 0;

	TFile *ft = 0;

	void setFiles();
	double getEffOR(const double eff1, const double eff2);
	TH2 * makeTriggerEffs2D(TString name, TH1 *h1, vector<double>& bins1, TH1 *h2, vector<double>& bins2);
	TH1 * getEff1D(TH1* num, TH1* den, vector<double> binedges, TString name, TString axS);
	TH1 * getTotHadHistMC(TFile *f, TString hS);
	TH2 * getCombPt2HistWithHT(TString name, double ht, TH2 *hist);
	void setTriggerDistributions(bool is1l, bool doHT, TString enumS, TString mnumS, bool doLumiWt=false);
	void setSumTriggerDists(TString idS, TString eNumS, TString eDenS, TString muNumS, TString muDenS, bool doLumiWt);

};

void TriggerEffProcessor::setFiles() {
	TString fPre = "/Users/brentstone/Dropbox/Physics/HHbbWW/trigger/";
	ft = TFile::Open(fPre+"master/triggerInfo_"+eraStrs[era]+".root");
//	ft = TFile::Open(fPre+"turnons/turnons"+eraStrs[era]+".root");
}

double TriggerEffProcessor::getEffOR(const double eff1, const double eff2) {
	return (eff1 + eff2 - eff1*eff2);
}

TH2 * TriggerEffProcessor::makeTriggerEffs2D(TString name, TH1 *h1, vector<double>& bins1, TH1 *h2, vector<double>& bins2) {

	static const int nbins1 = bins1.size() - 1;
	static const int nbins2 = bins2.size() - 1;

	TH2 *hh = new TH2D(name,";H_{T};p_{T}",nbins1,&bins1[0],nbins2,&bins2[0]);
	for(unsigned int i1=1; i1<=nbins1; i1++) for(unsigned int i2=1; i2<=nbins2; i2++) {
		double eff = getEffOR(h1->GetBinContent(i1),h2->GetBinContent(i2));
		hh->SetBinContent(i1,i2,eff);
	}

	return hh;

}

TH1 * TriggerEffProcessor::getEff1D(TH1* num, TH1* den, vector<double> binedges, TString name, TString axS) {

	int nbins = binedges.size() - 1;

	PlotTools::toOverflow(num);
	PlotTools::toUnderflow(num);
	PlotTools::toOverflow(den);
	PlotTools::toUnderflow(den);

	num = PlotTools::rebin(num,nbins,&binedges[0]);
	den = PlotTools::rebin(den,nbins,&binedges[0]);

	PlotTools::toOverflow(num);
	PlotTools::toUnderflow(num);
	PlotTools::toOverflow(den);
	PlotTools::toUnderflow(den);

	TH1 *rat = (TH1*)num->Clone(name);
	rat->GetXaxis()->SetTitle(axS);
	rat->Divide(num,den,1,1,"b");
	return rat;
}

TH1 * TriggerEffProcessor::getTotHadHistMC(TFile *f, TString hS) {
	vector<TString> samps = {"ttbar1_","wjets_","singlet_"};
	TH1 *h = (TH1*)f->Get(samps[0]+hS);
	h = (TH1*)h->Clone("mc_"+hS);
	for (unsigned int i=1; i<samps.size(); ++i) {
		h->Add((TH1*)f->Get(samps[i]+hS),1);
	}
	return h;
}

TH2 * TriggerEffProcessor::getCombPt2HistWithHT(TString name, double ht, TH2 *hist) {
	TH2 *hh = (TH2*)hist->Clone(name);
	for(unsigned int i1=1; i1<=hist->GetNbinsX(); ++i1) for(unsigned int i2=1; i2<=hist->GetNbinsY(); ++i2) {
		double cont = hist->GetBinContent(i1,i2);
		hh->SetBinContent(i1,i2,getEffOR(ht,cont));
	}
	return hh;
}

void TriggerEffProcessor::setTriggerDistributions(bool is1l, bool doHT, TString enumS, TString mnumS, bool doLumiWt) {
	TString idS = is1l ? "id1_" : "id2_";
	TString eVarS, mVarS;
	if(doHT) {
		eVarS = TString::Format("elpt%d_ht",elpt);
		mVarS = TString::Format("mupt%d_ht",mupt);
	} else {
		eVarS = "ht400_pt";
		mVarS = eVarS;
	}

	if(doLumiWt) {
		enumS += "_lumiwt";
		mnumS += "_lumiwt";
	}

	TString elNumS = "el_tt2_passSMuAnd"+enumS+"_"+eVarS;
	TString elDenS = "el_tt2_passSMu_"+eVarS;
	TString muNumS = "mu_tt2_passSElAnd"+mnumS+"_"+mVarS;
	TString muDenS = "mu_tt2_passSEl_"+mVarS;

//	if(era == E_RUN2) {
//		this->setSumTriggerDists(idS,elNumS,elDenS,muNumS,muDenS,doLumiWt);
//		return;
//	}

	hen_mc = (TH1*)ft->Get("ttbar2_"+idS+elNumS);
	hed_mc = (TH1*)ft->Get("ttbar2_"+idS+elDenS);
	hmn_mc = (TH1*)ft->Get("ttbar2_"+idS+muNumS);
	hmd_mc = (TH1*)ft->Get("ttbar2_"+idS+muDenS);

	if(doLumiWt) {
		elNumS.ReplaceAll("lumiwt_","");
		muNumS.ReplaceAll("lumiwt_","");
	}

	hen_da = (TH1*)ft->Get("muData_"+idS+elNumS);
	hed_da = (TH1*)ft->Get("muData_"+idS+elDenS);
	hmn_da = (TH1*)ft->Get("elData_"+idS+muNumS);
	hmd_da = (TH1*)ft->Get("elData_"+idS+muDenS);
}

void TriggerEffProcessor::setSumTriggerDists(TString idS, TString elNumS, TString elDenS, TString muNumS, TString muDenS, bool doLumiWt) {
	hen_mc = (TH1*)ft->Get("ttbar2_"+idS+"2017_"+elNumS);
	hed_mc = (TH1*)ft->Get("ttbar2_"+idS+"2017_"+elDenS);
	hmn_mc = (TH1*)ft->Get("ttbar2_"+idS+"2017_"+muNumS);
	hmd_mc = (TH1*)ft->Get("ttbar2_"+idS+"2017_"+muDenS);

	if(doLumiWt) {
		elNumS.ReplaceAll("lumiwt_","");
		muNumS.ReplaceAll("lumiwt_","");
	}

	hen_da = (TH1*)ft->Get("muData_"+idS+"2017_"+elNumS);
	hed_da = (TH1*)ft->Get("muData_"+idS+"2017_"+elDenS);
	hmn_da = (TH1*)ft->Get("elData_"+idS+"2017_"+muNumS);
	hmd_da = (TH1*)ft->Get("elData_"+idS+"2017_"+muDenS);

	vector<TString> yrs = {"2016","2018"};
	for(unsigned yr = 0; yr < yrs.size(); ++yr) {
		hen_mc->Add((TH1*)ft->Get("ttbar2_"+idS+yrs[yr]+"_"+elNumS),1);
		hed_mc->Add((TH1*)ft->Get("ttbar2_"+idS+yrs[yr]+"_"+elDenS),1);
		hmn_mc->Add((TH1*)ft->Get("ttbar2_"+idS+yrs[yr]+"_"+muNumS),1);
		hmd_mc->Add((TH1*)ft->Get("ttbar2_"+idS+yrs[yr]+"_"+muDenS),1);

		hen_da->Add((TH1*)ft->Get("muData_"+idS+yrs[yr]+"_"+elNumS),1);
		hed_da->Add((TH1*)ft->Get("muData_"+idS+yrs[yr]+"_"+elDenS),1);
		hmn_da->Add((TH1*)ft->Get("elData_"+idS+yrs[yr]+"_"+muNumS),1);
		hmd_da->Add((TH1*)ft->Get("elData_"+idS+yrs[yr]+"_"+muDenS),1);
	}

	hen_mc = (TH1*)hen_mc->Clone("ttbar2_"+idS+"Run2_"+elNumS);
	hed_mc = (TH1*)hed_mc->Clone("ttbar2_"+idS+"Run2_"+elDenS);
	hmn_mc = (TH1*)hmn_mc->Clone("ttbar2_"+idS+"Run2_"+muNumS);
	hmd_mc = (TH1*)hmd_mc->Clone("ttbar2_"+idS+"Run2_"+muDenS);

	hen_da = (TH1*)hen_da->Clone("muData_"+idS+"Run2_"+elNumS);
	hed_da = (TH1*)hed_da->Clone("muData_"+idS+"Run2_"+elDenS);
	hmn_da = (TH1*)hmn_da->Clone("elData_"+idS+"Run2_"+muNumS);
	hmd_da = (TH1*)hmd_da->Clone("elData_"+idS+"Run2_"+muDenS);

}

vector<TH1*> TriggerEffProcessor::getEffsAndSFs_lnuqq(bool onlyForFile) {
	vector<TH1*> hists;
	TString eS = TString::Format("elpt%d",elpt);
	TString mS = TString::Format("mupt%d",mupt);

	this->setTriggerDistributions(true,true,"Full","Full");
	TH1 *elEffMC = getEff1D(hen_mc,hed_mc,htbins,"effMC_e_"+eS,"H_{T}");
	TH1 *muEffMC = getEff1D(hmn_mc,hmd_mc,htbins,"effMC_m_"+mS,"H_{T}");
	TH1 *elEffDA = getEff1D(hen_da,hed_da,htbins,"effDA_e_"+eS,"H_{T}");
	TH1 *muEffDA = getEff1D(hmn_da,hmd_da,htbins,"effDA_m_"+mS,"H_{T}");

	if (!onlyForFile) {
		hists.push_back(elEffMC);
		hists.push_back(muEffMC);
		hists.push_back(elEffDA);
		hists.push_back(muEffDA);
	}

	elEffDA = (TH1*)elEffDA->Clone("sf_e_"+eS);
	elEffDA->Divide(elEffMC);
	muEffDA = (TH1*)muEffDA->Clone("sf_m_"+mS);
	muEffDA->Divide(muEffMC);

	if(onlyForFile) {
		elEffDA = (TH1*)elEffDA->Clone("electronSF_lnuqq");
		muEffDA = (TH1*)muEffDA->Clone("muonSF_lnuqq");
		hists.push_back(elEffDA);
		hists.push_back(muEffDA);
		return hists;
	}

	hists.push_back(elEffDA);
	hists.push_back(muEffDA);

	return hists;
}

vector<TH2*> TriggerEffProcessor::getEffsAndSFs_lnulnu(bool includeCross) {
	vector<TH2*> hists2D;
	TString crossS = (includeCross ? "Cross" : "");
	TString elptS = TString::Format("elpt%d",elpt);
	TString muptS = TString::Format("mupt%d",mupt);

	// leading lepton and supplemental trigger eff in di-lepton channel -----------------------------
	this->setTriggerDistributions(false,true,"Full","Full",era==E_2017?true:false);

	TH1 *elEffMC = getEff1D(hen_mc,hed_mc,htbins,"electronEffs_mcHT_lnulnu","H_{T}");
	TH1 *muEffMC = getEff1D(hmn_mc,hmd_mc,htbins,"muonEffs_mcHT_lnulnu","H_{T}");
	TH1 *elEffDA = getEff1D(hen_da,hed_da,htbins,"electronEffs_dataHT_lnulnu","H_{T}");
	TH1 *muEffDA = getEff1D(hmn_da,hmd_da,htbins,"muonEffs_dataHT_lnulnu","H_{T}");

//  get single-lepton trigger effs for second lepton -----------------------------
	this->setTriggerDistributions(false,false,"SEl"+crossS,"SMu"+crossS,(includeCross&&era==E_2017?true:false));
	TH1 *el2EffMC = getEff1D(hen_mc,hed_mc,ptbins,"electronEffs_mcPT_lnulnu","p_{T}");
	TH1 *mu2EffMC = getEff1D(hmn_mc,hmd_mc,ptbins,"muonEffs_mcPT_lnulnu","p_{T}");
	TH1 *el2EffDA = getEff1D(hen_da,hed_da,ptbins,"electronEffs_dataPT_lnulnu","p_{T}");
	TH1 *mu2EffDA = getEff1D(hmn_da,hmd_da,ptbins,"muonEffs_dataPT_lnulnu","p_{T}");

	//  get 2D effs and scale factors for di-lepton channel
	TH2 *effMC_ee = makeTriggerEffs2D("effMC_ee_"+elptS+crossS,elEffMC,htbins,el2EffMC,ptbins); hists2D.push_back(effMC_ee);
	TH2 *effMC_me = makeTriggerEffs2D("effMC_me_"+muptS+crossS,muEffMC,htbins,el2EffMC,ptbins); hists2D.push_back(effMC_me);
	TH2 *effMC_em = makeTriggerEffs2D("effMC_em_"+elptS+crossS,elEffMC,htbins,mu2EffMC,ptbins); hists2D.push_back(effMC_em);
	TH2 *effMC_mm = makeTriggerEffs2D("effMC_mm_"+muptS+crossS,muEffMC,htbins,mu2EffMC,ptbins); hists2D.push_back(effMC_mm);

	TH2 *effDA_ee = makeTriggerEffs2D("effDA_ee_"+elptS+crossS,elEffDA,htbins,el2EffDA,ptbins); hists2D.push_back(effDA_ee);
	TH2 *effDA_me = makeTriggerEffs2D("effDA_me_"+muptS+crossS,muEffDA,htbins,el2EffDA,ptbins); hists2D.push_back(effDA_me);
	TH2 *effDA_em = makeTriggerEffs2D("effDA_em_"+elptS+crossS,elEffDA,htbins,mu2EffDA,ptbins); hists2D.push_back(effDA_em);
	TH2 *effDA_mm = makeTriggerEffs2D("effDA_mm_"+muptS+crossS,muEffDA,htbins,mu2EffDA,ptbins); hists2D.push_back(effDA_mm);

	TH2 *sfee = (TH2*)effDA_ee->Clone("sf_ee_"+elptS+crossS);
	TH2 *sfem = (TH2*)effDA_em->Clone("sf_em_"+elptS+crossS);
	TH2 *sfme = (TH2*)effDA_me->Clone("sf_me_"+muptS+crossS);
	TH2 *sfmm = (TH2*)effDA_mm->Clone("sf_mm_"+muptS+crossS);

	sfee->Divide(effMC_ee); hists2D.push_back(sfee);
	sfem->Divide(effMC_em); hists2D.push_back(sfem);
	sfme->Divide(effMC_me); hists2D.push_back(sfme);
	sfmm->Divide(effMC_mm); hists2D.push_back(sfmm);

	return hists2D;
}

vector<TH1*> TriggerEffProcessor::getEffsForFile_2l(bool includeCross) {
	vector<TH1*> hists;
	TString crossS = (includeCross ? "Cross" : "");

	// leading lepton and non-lepton trigger eff in di-lepton channel -----------------------------
	this->setTriggerDistributions(false,true,"Full","Full",era==E_2017?true:false);

	TH1 *elEffMC = getEff1D(hen_mc,hed_mc,htbins,"electronEffs_mcHT_lnulnu","H_{T}");
	TH1 *muEffMC = getEff1D(hmn_mc,hmd_mc,htbins,"muonEffs_mcHT_lnulnu","H_{T}");
	TH1 *elEffDA = getEff1D(hen_da,hed_da,htbins,"electronEffs_dataHT_lnulnu","H_{T}");
	TH1 *muEffDA = getEff1D(hmn_da,hmd_da,htbins,"muonEffs_dataHT_lnulnu","H_{T}");

//  get single-lepton trigger effs for second lepton -----------------------------
	this->setTriggerDistributions(false,false,"SEl"+crossS,"SMu"+crossS,(includeCross&&era==E_2017?true:false));
	TH1 *el2EffMC = getEff1D(hen_mc,hed_mc,ptbins,"electronEffs_mcPT_lnulnu","p_{T}");
	TH1 *mu2EffMC = getEff1D(hmn_mc,hmd_mc,ptbins,"muonEffs_mcPT_lnulnu","p_{T}");
	TH1 *el2EffDA = getEff1D(hen_da,hed_da,ptbins,"electronEffs_dataPT_lnulnu","p_{T}");
	TH1 *mu2EffDA = getEff1D(hmn_da,hmd_da,ptbins,"muonEffs_dataPT_lnulnu","p_{T}");

	hists.push_back(elEffMC);
	hists.push_back(el2EffMC);
	hists.push_back(muEffMC);
	hists.push_back(mu2EffMC);
	hists.push_back(elEffDA);
	hists.push_back(el2EffDA);
	hists.push_back(muEffDA);
	hists.push_back(mu2EffDA);

	return hists;
}

#endif

void makeTriggerCorrections(int year, bool makePlotsForFileOnly = true) {

	int era;
	if(year == 2016) era = 0;
	else if(year == 2017) era = 1;
	else if(year == 2018) era = 2;
	else if(year == 0) era = 3;
	else return;

    vector<double> htbins = {100,200,250,300,350,400,450,500,550,600,650,
		700,800,900,1000,1200,2000};
    vector<double> ptbins = {5,10,15,20,25,30,35,40,50,60,70,80,90,100,150,200,300,1000};

	TriggerEffProcessor proc(era,30,27);
	proc.setPTBinning(ptbins);
	proc.setHTBinning(htbins);

	TString ys = (year == 0 ? "Run2" : TString::Format("%d",year));

	if(makePlotsForFileOnly) {
		vector<TH1*> hists1 = proc.getEffsAndSFs_lnuqq(true);
		vector<TH1*> hists2 = proc.getEffsForFile_2l(false);

		TFile *fout = new TFile("triggerSF_"+ys+".root","RECREATE");
		fout->cd();
		for(const auto& h : hists1) h->Write();
		for(const auto& h : hists2) h->Write();
		fout->Close();
		delete fout;
		return;
	}

	vector<TH2*> hists_nom       = proc.getEffsAndSFs_lnulnu(false);
	vector<TH2*> hists_nomCross  = proc.getEffsAndSFs_lnulnu(true);
	vector<TH1*> hist_nom        = proc.getEffsAndSFs_lnuqq(false);

	proc.setPtThresh(32,29);
	vector<TH2*> hists_up        = proc.getEffsAndSFs_lnulnu(false);
	vector<TH2*> hists_upCross   = proc.getEffsAndSFs_lnulnu(true);
	vector<TH1*> hist_up1        = proc.getEffsAndSFs_lnuqq(false);

	proc.setPtThresh(35,32);
	vector<TH2*> hists_upup      = proc.getEffsAndSFs_lnulnu(false);
	vector<TH2*> hists_upupCross = proc.getEffsAndSFs_lnulnu(true);
	vector<TH1*> hist_up2        = proc.getEffsAndSFs_lnuqq(false);

	proc.setPtThresh(28,25);
	vector<TH2*> hists_down      = proc.getEffsAndSFs_lnulnu(false);
	vector<TH2*> hists_downCross = proc.getEffsAndSFs_lnulnu(true);
	vector<TH1*> hist_down1      = proc.getEffsAndSFs_lnuqq(false);

	proc.setPtThresh(25,22);
	vector<TH2*> hists_downdown       = proc.getEffsAndSFs_lnulnu(false);
	vector<TH2*> hists_downdownCross  = proc.getEffsAndSFs_lnulnu(true);
	vector<TH1*> hist_down2           = proc.getEffsAndSFs_lnuqq(false);

	TFile *fout = new TFile("effsAndSF_"+ys+".root","RECREATE");
	fout->cd();
	for(const auto& h : hists_nom)       h->Write();
	for(const auto& h : hists_up)        h->Write();
	for(const auto& h : hists_down)      h->Write();
	for(const auto& h : hists_downdown)      h->Write();
	for(const auto& h : hists_upup)      h->Write();
	for(const auto& h : hists_nomCross)  h->Write();
	for(const auto& h : hists_upCross)   h->Write();
	for(const auto& h : hists_downCross) h->Write();
	for(const auto& h : hists_upupCross) h->Write();
	for(const auto& h : hists_downdownCross) h->Write();

	for(const auto& h : hist_nom) h->Write();
	for(const auto& h : hist_up1) h->Write();
	for(const auto& h : hist_up2) h->Write();
	for(const auto& h : hist_down1) h->Write();
	for(const auto& h : hist_down2) h->Write();

	fout->Close();
	delete fout;

}
