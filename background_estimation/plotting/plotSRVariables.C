
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/makePlots.C"
#include "../predTools/CutConstants.h"
#include "../../framework/Configuration/interface/FillerConstants.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/Drawing.h"
#include "HistoPlotting/include/StyleInfo.h"
#include "../predTools/makeJSON.C"

#include "TStyle.h"
#include "TGraphAsymmErrors.h"

using namespace CutConstants;
CutStr blindCut = CutStr("blindCut",std::string("(")+hbbMCS.cut+"<102||"+hbbMCS.cut+">150)");
//CutStr blindCut = CutStr("blindCut",std::string("(1.0)"));
std::vector<PlotVar> vars;
std::vector<std::string> varUnits;
std::vector<PlotSamp> samps;
float SIG_CROSS = 10; //in pb
bool preliminary = false;

std::string getTTBarSF(const std::string& filename){
    CJSON json(filename+"_ttbarSF.json");
    std::string qToW = json.getP(0).second;
    auto replace = [&](const std::string& vn, const std::string tf1n){
        std:size_t index = 0;
        while (true) {
            index = qToW.find(vn, index);
            if (index == std::string::npos) break;
            qToW.replace(index, vn.size(), tf1n);
            index += 1;
        }
    };
    replace(MOD_MS,hhMCS);
    return std::string("(")+processes[TTBAR].cut+"?("+qToW+"):1.0)";
}

void makeDataDistributions(const int channel, const std::string& name, const std::string& filename,  const std::string inputFile, const std::string& cut,bool isData){
    std::vector<PlotSel> sels;
    if (channel == 1) {
        sels.emplace_back("full"            ,preSel1.cut+"&&"+btagCats[BTAG_LMT].cut+"&&"+selCuts1[SEL1_FULL].cut);
        sels.emplace_back("fullBlind"       ,preSel1.cut+"&&"+btagCats[BTAG_LMT].cut+"&&"+selCuts1[SEL1_FULL].cut+"&&"+blindCut.cut);
        sels.emplace_back("rem_nAK4Btags"   ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+preSel1.cut+"&&"+ptomC.cut+"&&"+wjjSC.cut+"&&"+hwwLC.cut);
        sels.emplace_back("rem_wjjTau2o1"   ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+preSel1.cut+"&&"+ptomC.cut+"&&"+bV.cut+"&&"+hwwLC.cut);
        sels.emplace_back("rem_ptoM"        ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+preSel1.cut+"&&"+hwwLC.cut+"&&"+wjjSC.cut+"&&"+bV.cut);
        sels.emplace_back("rem_hwwLi"       ,blindCut.cut +"&&"+ btagCats[BTAG_LMT].cut+"&&"+preSel1.cut+"&&"+ptomC.cut+"&&"+wjjSC.cut+"&&"+bV.cut);
        sels.emplace_back("rem_hbbTag"      ,blindCut.cut +"&&"+ selCuts1[SEL1_FULL].cut);
    } else if (channel == 2) {
    	sels.emplace_back("full"            ,preSel2.cut+"&&"+btagCats[BTAG_LMT].cut+"&&"+selCuts2[SEL2_FULL].cut);
        sels.emplace_back("fullBlind"       ,preSel2.cut+"&&"+btagCats[BTAG_LMT].cut+"&&"+selCuts2[SEL2_FULL].cut+"&&"+blindCut.cut);
    	sels.emplace_back("rem_nAK4Btags"   ,blindCut.cut+"&&"+btagCats[BTAG_LMT].cut+"&&"+preSel2.cut+"&&"+dPhiC.cut+"&&"+mllV.cut+"&&"+dRC.cut+"&&"+metC.cut);
    	sels.emplace_back("rem_met"         ,blindCut.cut+"&&"+btagCats[BTAG_LMT].cut+"&&"+preSel2.cut+"&&"+dPhiC.cut+"&&"+mllV.cut+"&&"+dRC.cut+"&&"+bV.cut);
    	sels.emplace_back("rem_dilepDR"     ,blindCut.cut+"&&"+btagCats[BTAG_LMT].cut+"&&"+preSel2.cut+"&&"+dPhiC.cut+"&&"+mllV.cut+"&&"+bV.cut+"&&"+metC.cut);
    	sels.emplace_back("rem_dilepMass"   ,blindCut.cut+"&&"+btagCats[BTAG_LMT].cut+"&&"+preSel2.cut+"&&"+dPhiC.cut+"&&"+bV.cut+"&&"+dRC.cut+"&&"+metC.cut);
    	sels.emplace_back("rem_llMetDphi"   ,blindCut.cut+"&&"+btagCats[BTAG_LMT].cut+"&&"+preSel2.cut+"&&"+bV.cut+"&&"+mllV.cut+"&&"+dRC.cut+"&&"+metC.cut);
        sels.emplace_back("rem_hbbTag"      ,blindCut.cut+"&&"+selCuts2[SEL2_FULL].cut);
    }

    std::string outFileName=filename+"_"+std::to_string(channel)+"l_"+name+"_srVarDistributions.root";
    MakePlots a(inputFile,outFileName,samps,sels,vars,cut,isData ? "1.0" : nomW.cut);
}

void compilePlots(const int year, const int nLep, const std::string& prefix, const std::string& mcFile, const std::string& dataFile,  const std::vector<std::string> signalFiles,  const std::vector<std::string> signalNames){
    TFile * fd = new TFile((dataFile).c_str(),"READ");
    TFile * fm = new TFile((mcFile  ).c_str(),"READ");

    std::vector<CutStr> bkgs;
    if(nLep == 1) {
    	bkgs.push_back(processes[TTBAR]);
    	bkgs.push_back(processes[WJETS]);
    	bkgs.push_back(processes[QCD]);
    	bkgs.push_back(processes[OTHER1]);
    } else if (nLep == 2) {
    	bkgs.push_back(processes[TTBAR]);
    	bkgs.push_back(processes[ZJETS]);
    	bkgs.push_back(processes[OTHER2]);
    	SIG_CROSS /= 10;
    }

    // iV = 0-2 are search mass variables in both channels
    for(unsigned int iV = 3; iV < vars.size(); ++iV){
    	std::cout<<vars[iV].varName<<std::endl;
        Plotter * p = new Plotter;
        std::vector<Drawing::TLegendEntryDef> legEntries;

        int rebinFactor = 0;
        if(nLep == 1) {
        	if(iV == 3) rebinFactor = 4; //hbbjet
        	if(iV == 7) rebinFactor = 4; //hwwLi
        } else {
        	if(iV == 3) rebinFactor = 4; //hbbjet
        	if(iV == 5) rebinFactor = 2; //mll
        	if(iV == 7) rebinFactor = 2; //dphimetll
        	if(iV == 8) rebinFactor = 5; //met

        }
//        switch(iV){
//        case 1:
//            rebinFactor = 3;
//            break;
//        case 2:
//            rebinFactor = 4;
//            break;
//        }

        // get bkg components
        for(unsigned int iP = 0; iP < bkgs.size(); ++iP){
            TH1 * hm = 0;
            std::cout<<(bkgs[iP]+"_rem_"+vars[iV].varName+"_"+vars[iV].varName).c_str()<<std::endl;
            fm->GetObject((bkgs[iP]+"_rem_"+vars[iV].varName+"_"+vars[iV].varName).c_str(),hm);
            if(hm == 0) continue;
            if(rebinFactor) hm->Rebin(rebinFactor);
            auto * g = p->addStackHist(hm,bkgs[iP].title);
            legEntries.push_back(std::make_tuple(10+bkgs.size() +iP,g,bkgs[iP].title.c_str(),"f"));
        }

        auto * tot = (TH1*)p->getTotStack()->Clone();
        auto* errBand = new TGraphAsymmErrors();

        // get stat errors on bkg components
        for(int iB = 1; iB <= tot->GetNbinsX(); ++iB){
            double x = tot->GetBinCenter(iB);
            double y = tot->GetBinContent(iB);
            double err = tot->GetBinError(iB);
            if(iB == 1){
                y += tot->GetBinContent(0);
                err = std::sqrt(err*err + tot->GetBinError(0)*tot->GetBinError(0));
            }
            if(iB == tot->GetNbinsX()){
                y += tot->GetBinContent(iB+1);
                err = std::sqrt(err*err + tot->GetBinError(iB+1)*tot->GetBinError(iB+1));
            }
            auto setPt = [&](const int bin, float x){
                errBand->SetPoint(bin,x,y);
                errBand->SetPointError(bin,tot->GetBinWidth(iB)/2.,tot->GetBinWidth(iB)/2,err,err);
            };
            if(iB == 1) setPt(iB-1,x-tot->GetBinWidth(iB)/2.);
            setPt(iB,x);
            if(iV==1){
                const auto& hs =p->getStackHists();
                std::cout << iB <<" "<< x <<" "<< y << " "<< err;
                for(const auto& h : hs){
                    std::cout << "(" <<((const TH1*)h.obj)->GetBinContent(iB)<<","<<((const TH1*)h.obj)->GetBinError(iB)<<") ";
                }

               std::cout <<std::endl;
            }
            if(iB == tot->GetNbinsX()) setPt(iB+1,x+tot->GetBinWidth(iB)/2.);
        }

        p->clearTotStackError();

        // add error band to plot
        int fillColor = kMagenta+4;
        errBand->SetFillColor(fillColor);
        errBand->SetFillStyle(3352);
        gStyle->SetHatchesLineWidth(1);
        gStyle->SetHatchesSpacing(.5);
        auto * g_mcunc = p->addGraph(errBand,"MC stat. unc.",fillColor,1,0,20,1,false,true,false,"2");
        legEntries.push_back(std::make_tuple(2,g_mcunc,"MC stat. unc.","F"));
        legEntries.push_back(std::make_tuple(4,(TObject*)(0),"",""));

        std::vector<int> signalColors ={kSpring+10,kBlue};

        // get signals
        for(unsigned int iS = 0; iS < signalFiles.size(); ++iS){
            TFile * f = new TFile((signalFiles[iS]).c_str(),"READ");
            TH1 * hm = 0;
            f->GetObject((std::string("all_rem_")+vars[iV].varName+"_"+vars[iV].varName).c_str(),hm);
            if(hm == 0) continue;

            //Scale so that we get 1pb normalization
            hm->Scale(SIG_CROSS);
            if(rebinFactor) hm->Rebin(rebinFactor);
            hm->SetName((signalNames[iS]+vars[iV].varName).c_str());
            auto * g = p->addHistLine(hm,signalNames[iS],signalColors[iS]);
            legEntries.push_back(std::make_tuple(100+iS,g,signalNames[iS],"L"));
        }

        // get data
        TH1 * hd = 0;
        fd->GetObject((std::string("all_rem_")+vars[iV].varName+"_"+vars[iV].varName).c_str(),hd);
        if(hd == 0) continue;
        if(rebinFactor) hd->Rebin(rebinFactor);
        p->addHist(hd,"Data",kBlack,1,2,20,0.5,true,true,true);
        TGraphAsymmErrors * g_d = new TGraphAsymmErrors;
        g_d->SetLineColor  (kBlack);
        g_d->SetLineWidth  (2);
        g_d->SetLineStyle  (1);
        g_d->SetMarkerStyle(20);
        g_d->SetMarkerColor(kBlack);
        g_d->SetMarkerSize (0.5);
        legEntries.push_back(std::make_tuple(0,g_d,"Data","P E"));
        legEntries.push_back(std::make_tuple(1,(TObject*)(0),"",""));

        p->setCMSLumi(year,10);

        // customize range for each variable
        bool sup = false;
        switch(iV){
        case 1: //hbb
            p->setMinMax(0,575);
            break;
        case 2: //hh
            p->setMinMax(0.7,7500);
            p->turnOffTrailingPoissonZeros();
            break;
        case 3: // hbbTag
        	p->setMinMax(0,2500);
        	break;
        case 4://nAK4Btags
            p->setMinMax(0,20000);
            break;
        case 5://wjjTau2o1
            p->setMinMax(0,1000);
            break;
        case 6://ptoM
            p->setMinMax(0,1150);
            break;
        case 7://hwwLi
        	p->setMinMax(0,1000);
            break;
        }

//        if(nLep == 1) {
//
//        } else if(nLep == 2) {
//
//        }

        double xV =0.45;
        double yV =0.619;

//        p->addText("All categories",xV+0.0075,yV+0.2075,0.045);
        if(sup||preliminary){
            p->setCMSLumiPosition(0,1.05);
           if(sup) p->setCMSLumiExtraText("Supplementary");
           else p->setCMSLumiExtraText("Preliminary");
            p->addText("All categories",.18,0.83,0.045);
        } else {
            p->addText("All categories",.18,0.75,0.045);
        }

        //--------------------LEGEND AND TEXT------------------------------
        p->turnOffLegend();
        TLegend * legend = signalNames.size() ? new TLegend(xV,yV,xV+0.453,yV+0.25) : new TLegend(xV,yV,xV+0.445,yV+0.25);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetNColumns(2);
        std::sort(legEntries.begin(), legEntries.end(), [](const Drawing::TLegendEntryDef& a, const Drawing::TLegendEntryDef& b) {return  std::get<0>(a) < std::get<0>(b); }  );
        for(const auto& l : legEntries){
            legend->AddEntry(std::get<1>(l),std::get<2>(l),std::get<3>(l));
        }
        //--------------------LEGEND AND TEXT------------------------------


        if(signalNames.size()){
            p->addText(TString::Format("#sigma#bf{#it{#Beta}}(X #rightarrow HH) = %0.f pb",SIG_CROSS),xV+0.0075,yV-0.04,0.042);
        } else {
        }
        p->setLegendNColumns(2);

        std::string ytitle = "Events";
        if(varUnits[iV]!="-1"){
            ytitle += " / " + ASTypes::flt2Str(p->getTotStack()->GetBinWidth(1)) + " " + varUnits[iV];

        }
        p->setBotMinMax(0.05,1.95);
        p->setYTitle(ytitle);
        p->setYTitleBot("Data / MC");
//        auto * c = p->draw(false,prefix+vars[iV].varName+"_srvardists.pdf");
        auto * c = p->drawSplitRatio(-1,"stack",false,true,prefix+vars[iV].varName+"_srvardists.pdf");
//        p->yAxis()->SetTitleOffset(1.55);
        p->xAxis()->SetTitleOffset(1.05);
        c->GetPad(1)->cd();
        legend->Draw();
        c->GetPad(1)->Update();

        p->botStyle.xAxis->SetTitleOffset(1.05);


        for(unsigned int iSig = 0; iSig < signalFiles.size(); ++iSig){
            auto prim = c->GetPad(2)->GetPrimitive((signalNames[iSig]+vars[iV].varName).c_str());
            if(prim) prim->Delete();
        }

        if(iV == 1){//temp hack
            p->botStyle.xAxis->SetTitle(hbbMCS.title.c_str());
        }
        if(iV == 3){//temp hack
            p->botStyle.xAxis->SetTitle("DeepAK8-MD Z/H#rightarrowbb score");
        }
        if(nLep == 1) {
            if(iV == 5){//temp hack
                p->botStyle.xAxis->SetTitle("q#bar{q} jet #tau_{2}/#tau_{1}");
            }
            if(iV == 7){//temp hack
                p->botStyle.xAxis->SetTitle("#it{D}_{l#nuqq}");
            }
        } else if (nLep == 2) {
        	if(iV == 8) p->botStyle.xAxis->SetTitle("MET [GeV]");
        }


        if(iV == 2){
            c->GetPad(1)->SetLogy();
        }
        c->GetPad(1)->Update();
        c->Print((prefix+vars[iV].varName+"_srvardists.pdf").c_str());


    }
}

void categoryPlots(const std::string& prefix, const std::string& mcFile){
    TFile * fm = new TFile((mcFile  ).c_str(),"READ");

    for(unsigned int iV = 0; iV < vars.size(); ++iV){
        Plotter * p = new Plotter;


        //
        for(unsigned int iP = 0; iP < bkgSels.size(); ++iP){
            TH1 * hm = 0;
            fm->GetObject((bkgSels[iP]+"_full_"+vars[iV].varName).c_str(),hm);
            if(hm == 0) continue;
            p->addStackHist(hm,bkgSels[iP].title);
        }
        if(iV == 2) {
            p->rebin(4);
            p->setMinMax(0.1,10000);
        }
        if(iV== 7) p->setMinMax(20,1000000);

        if(iV == 4) p->setLegendPos(0.1588629,0.5895652,0.4397993,0.9043478);
        else p->setLegendPos(0.6555184,0.5530435,0.9364548,0.8678261);
        p->setCMSLumi();
        p->setYTitle("N. of events / bin width");
        auto * c = p->draw(false,prefix+vars[iV].varName+"_srvardists.pdf");
        p->yAxis()->SetTitleOffset(1.55);
        p->xAxis()->SetTitleOffset(1.0);
        if(iV == 7){
            p->xAxis()->SetBinLabel(1,"No loose");
            p->xAxis()->SetBinLabel(2,"One loose");
            p->xAxis()->SetBinLabel(3,"Two loose");
            p->xAxis()->SetBinLabel(4,"One medium, no loose");
            p->xAxis()->SetBinLabel(5,"One medium, one loose");
            p->xAxis()->SetBinLabel(6,"Two medium");
            p->xAxis()->SetTitle(" ");
            c->SetLogy();


        }
        if(iV == 2){
            c->SetLogy();
        }
        c->Update();
        c->Print((prefix+vars[iV].varName+"_catdists.pdf").c_str());


    }
}

#endif

void plotSRVariables(int step, int nLep, int year, int reg, std::string tree = "", std::string name = ""){
    std::string prefix = "srVarDists/";
    if(reg == REG_TOPCR){
        prefix+="TopCR_";
        hhFilename += "_TopCR";
        blindCut.cut = "(1.0)";
        selCuts1[SEL1_FULL].cut = preSel1.cut + "&&"+abV.cut+"&&"+wjjSC.cut+"&&"+hwwLC.cut;
        selCuts2[SEL2_FULL].cut  = preSel2.cut+"&&"+abV.cut+"&&"+drCrC.cut+"&&"+mllV.cut+"&&"+metC.cut+"&&"+dPhiC.cut;
    } else if(reg==REG_NONTOPCR){
        prefix+="NonTopCR_";
        hhFilename +="_NonTopCR";
        blindCut.cut = "(1.0)";
        btagCats = qgBtagCats;
    }

    bool isData = ASTypes::strFind(tree,"data");
    bool isSignal = ASTypes::strFind(tree,"spin");
    std::string cut =hhRange.cut+"&&"+hbbRange.cut ;
    std::string out = "srVarDists/"+hhFilename;
    if(ASTypes::strFind(tree,"radion")) blindCut.cut = "(1.0)";
    if(isData) cut += "&&"+blindCut.cut;

    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,
            minHbbMass,maxHbbMass,hhMCS,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins
            ,minHHMass,maxHHMass );
    vars.emplace_back(hbbMCS,std::string(";")+hbbMCS.title,hbbMCS.cut,nHbbMassBins,minHbbMass,maxHbbMass);
    vars.emplace_back(hhMCS ,std::string(";")+hhMCS.title,hhMCS.cut,nHHMassBins,minHHMass,maxHHMass );
    varUnits.emplace_back("GeV");
    varUnits.emplace_back("GeV");
    varUnits.emplace_back("GeV");

    vars.emplace_back("hbbTag","; DeepAK8-MD Z/H#rightarrowbb score","hbbTag",100,0,1);
    vars.emplace_back("nAK4Btags" ,"; N. of AK4 jet b tags","nAK4Btags",4,-0.5,3.5);
    varUnits.emplace_back("-1");
    varUnits.emplace_back("units");

    if (nLep == 1) {
        vars.emplace_back("wjjTau2o1" ,"; q#bar{q} jet #tau_{2}/#tau_{1}","wjjTau2o1",50,0,1);
        vars.emplace_back("ptoM" ,"; #it{p}_{T}/#it{m}","hwwPT/hhMass",50,0,1);
        vars.emplace_back("hwwLi","; #it{#D}_{l#nuqq}","hwwLi",100,0,20);
        varUnits.emplace_back("-1");
        varUnits.emplace_back("-1");
        varUnits.emplace_back("units");
    } else if (nLep == 2) {
    	vars.emplace_back("dilepMass",";m_{ll}","dilepMass",50,0,100);
    	vars.emplace_back("dilepDR",";#DeltaR_{ll}","dilepDR",50,0,TMath::Pi());
    	vars.emplace_back("llMetDphi",";#Delta#phi (MET,p_{ll})","llMetDphi",50,0,TMath::Pi());
    	vars.emplace_back("met",";MET)","met",100,0,500);
        varUnits.emplace_back("GeV");
        varUnits.emplace_back("units");
        varUnits.emplace_back("units");
        varUnits.emplace_back("GeV");
    }


    std::string ttSF = getTTBarSF("supportInputs/HHbb1o2l");

    samps.emplace_back("all","1.0");
    if(!isData && !isSignal){
        samps.emplace_back(bkgSels[BKG_QG],bkgSels[BKG_QG].cut);
        samps.emplace_back(bkgSels[BKG_LOSTTW],bkgSels[BKG_LOSTTW].cut+"*"+ttSF);
        samps.emplace_back(bkgSels[BKG_MW],bkgSels[BKG_MW].cut+"*"+ttSF);
        samps.emplace_back(bkgSels[BKG_MT],bkgSels[BKG_MT].cut+"*"+ttSF);

        if(nLep == 1) {
        	samps.emplace_back(processes[TTBAR],processes[TTBAR].cut+"*"+ttSF);
        	samps.emplace_back(processes[WJETS],processes[WJETS].cut);
        	samps.emplace_back(processes[QCD],processes[QCD].cut);
        	samps.emplace_back(processes[OTHER1],processes[OTHER1].cut);
        } else if(nLep == 2) {
        	samps.emplace_back(processes[TTBAR],processes[TTBAR].cut+"*"+ttSF);
        	samps.emplace_back(processes[ZJETS],processes[ZJETS].cut);
        	samps.emplace_back(processes[OTHER2],processes[OTHER2].cut);
        }
    }


    if(step==0){
        makeDataDistributions(nLep,name,out,tree,cut,isData);
    }
    if(step==1){
    	if(nLep == 1) out += "_1l";
    	else if(nLep == 2) out += "_2l";
        if(reg==REG_SR){
            compilePlots(year,nLep,prefix,out+"_mc_srVarDistributions.root",out+"_data_srVarDistributions.root",
                    {out+"_spin0_m1000_srVarDistributions.root",out+"_spin0_m3000_srVarDistributions.root"},
                    {"1 TeV X_{spin-0}","3.0 TeV X_{spin-0}"});
        } else{
            compilePlots(year,nLep,prefix,out+"_mc_srVarDistributions.root",out+"_data_srVarDistributions.root"
                    ,{},{});

        }

    }
    if(step==2) categoryPlots(prefix,out+"_mc_srVarDistributions.root");
}
