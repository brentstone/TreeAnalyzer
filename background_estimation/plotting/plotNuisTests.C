#include "../predTools/NuisPlotter.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TArrow.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include "plotTestHelper.h"
using namespace CutConstants;
using namespace ASTypes;
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

std::string getSystName (const std::string& prefix, const std::string& proc,const std::string& name, const std::string& sel){
    std::string sn = prefix;
    if(proc.size() )sn += proc;
    if(name.size() )sn += "_"+name;
    if(sel.size() )sn += "_"+sel;
    return sn;
};

//systs -> plotting name : hName
TCanvas * makeNuisPlot(const std::string& filename, const CutStr& bkg, const std::string& reg, const std::string& pName, const std::vector<std::pair<std::string,std::string>> systs, bool plotX){
	TH1::AddDirectory(kFALSE);
    TFile * f = new TFile(filename.c_str(),"read");
    std::string prefix = "shapeBkg_" + bkg + "_"+reg+"_13TeV";
    if(bkg == bkgSels[BKG_QG] || bkg == bkgSels[BKG_LOSTTW]) prefix +="_opt";
    if(f==0){
        std::cout <<"NO FILE"<<std::endl;
        return 0;
    }

    TH2 * nH = 0;
    f->GetObject((prefix+"_nom").c_str(),nH);
    if(nH==0){
        std::cout <<"NO NOM"<<std::endl;
        return 0;
    }
    nH->SetDirectory(0);

    std::vector<TH2*> sUps;
    std::vector<TH2*> sDowns;
    for(auto& s : systs){
        TH2 * uH = 0;
        TH2 * dH = 0;
        f->GetObject((prefix+"_"+s.second+"_up").c_str(),uH);
        f->GetObject((prefix+"_"+s.second+"_down").c_str(),dH);

        if(uH==0){
            std::cout <<"NO UP"<<std::endl;
            std::cout << prefix+"_"+s.second+"_up"<<std::endl;
            return 0;
        }
        if(dH==0){
            std::cout <<"NO DOWN"<<std::endl;
            return 0;
        }

        uH->SetDirectory(0);
        dH->SetDirectory(0);
        sUps.push_back(uH);
        sDowns.push_back(dH);
    }

    auto proc =[&](TH2* inH, const std::string& name) ->TH1* {
        if(plotX) return inH->ProjectionX(name.c_str());
        return inH->ProjectionY(name.c_str());
    };

    Plotter * p = new Plotter();
    p->addHistLine(proc(nH,bkg + "_"+reg+"_"+pName+"_nom"),"Nominal template"  );

    for(unsigned int iS = 0; iS < systs.size(); ++iS){
        p->addHistLine(proc(sUps[iS],bkg + "_"+reg+"_"+pName+systs[iS].second+"Up"),systs[iS].first +" unc. up",StyleInfo::getLineColor(iS+1) , 1  );
        p->addHistLine(proc(sDowns[iS],bkg + "_"+reg+"_"+pName+systs[iS].second+"Down"),systs[iS].first +" unc. down",StyleInfo::getLineColor(iS+1) , 9  );
    }
    if(plotX){
        p->setXTitle(hbbMCS.title);
        if(bkg == bkgSels[BKG_QG] || bkg == bkgSels[BKG_LOSTTW])
            p->setMinMax(0,0.18);
        else
            p->setMinMax(0,0.4);
    } else {
        p->setMinMax(0.00001,0.5);
        p->setXTitle(hhMCS.title);
    }
    p->setYTitle("Arbitrary scale");
    p->setYTitleBot("Alt. / initial temp.");
    p->setLegendPos(0.45,0.5530435,0.93,0.93);
    p->setCMSLumi();
    p->setCMSLumiPosition(0,1.05);
    p->setCMSLumiLumiText("13 TeV");
    p->setBotMinMax(0.5,1.5);

    auto c = p->drawSplitRatio(0,"stack",false,false,bkg+"_"+reg+"_"+pName);
    c->SetTitle((bkg+"_"+reg+"_"+pName).c_str());
    p->botStyle.xAxis->SetTitleOffset(1.05);
    if(!plotX){
        c->GetPad(1)->SetLogy();
        c->Update();
    }
    return c;

}

void getNuisHists(const std::string& rootFile, const std::string& outDir){

    std::string h = selCuts1[SEL1_FULL];

	for(const auto& l : lepCats) for(const auto& b : btagCats) for(const auto& p : purCats) {
		if(l == lepCats[LEP_EMU]) continue;
		if(b == btagCats[BTAG_LMT]) continue;
		if(p == purCats[PURE_I]) continue;

	    const std::string cat = l +"_"+b+"_"+p+"_"+h;
	    const std::string catNoSel = l +"_"+b+"_"+p;

	    auto systName = [&](const std::string& proc,const std::string& name, const std::string& sel = "-1")->std::string {
	        return getSystName("",proc,name, sel == "-1" ? cat : sel  );
	    };

	    auto go=[&](const CutStr& bkg, const std::vector<std::string>& nuis ){
	        std::string pdfN = "shapeBkg_"+bkg + "_"+cat+"_13TeV";
	        if(bkg == bkgSels[BKG_QG] || bkg == bkgSels[BKG_LOSTTW]) pdfN +="_opt";
	        std::cout<<pdfN<<std::endl;
	        for(const auto& s : nuis) std::cout<<s<<std::endl;
	        NuisPlotter(rootFile,pdfN,nuis,outDir+bkg+"_"+catNoSel+"_nuisHists.root");
	    };

	    go(bkgSels[BKG_QG],{
        		systName(bkgSels[BKG_QG],"PTX" ,b+"_1l"),
				systName(bkgSels[BKG_QG],"OPTX",b+"_1l"),
				systName(bkgSels[BKG_QG],"PTY",catNoSel),
				systName(bkgSels[BKG_QG],"OPTY",catNoSel)
	    });
	    go(bkgSels[BKG_LOSTTW],{
	            systName(bkgSels[BKG_LOSTTW],"PTX" ,b+"_1l"),
				systName(bkgSels[BKG_LOSTTW],"OPTX",b)+"_1l",
	            systName("top","res",catNoSel),
				systName("top","scale",catNoSel),
				systName("top","lostmw_rel_scale",b+"_1l")
	    });
	    go(bkgSels[BKG_MT],{
	            "hbb_scale_top",
				"hbb_res_top",
	            systName("top","res",catNoSel),
				systName("top","scale",catNoSel),
				systName("top","mt_rel_scale",b+"_1l")
	    });
	    go(bkgSels[BKG_MW],{
	            "hbb_scale",
				"hbb_res",
	            systName("top","res",catNoSel),
				systName("top","scale",catNoSel),
				systName("top","lostmw_rel_scale",b+"_1l")
	    });
	}

	h = selCuts2[SEL2_FULL];
	for(const auto& l : dilepCats) for(const auto& b : btagCats) {
		if(l == dilepCats[LEP_INCL]) continue;
		if(b == btagCats[BTAG_LMT]) continue;

	    const std::string cat = l +"_"+b+"_"+h;
	    const std::string catNoSel = l +"_"+b;

	    auto systName = [&](const std::string& proc,const std::string& name, const std::string& sel = "-1")->std::string {
	        return getSystName("",proc,name, sel == "-1" ? cat : sel  );
	    };

	    auto go=[&](const CutStr& bkg, const std::vector<std::string>& nuis ){
	        std::string pdfN = "shapeBkg_"+bkg + "_"+cat+"_13TeV";
	        if(bkg == bkgSels[BKG_QG] || bkg == bkgSels[BKG_LOSTTW]) pdfN +="_opt";
	        std::cout<<pdfN<<std::endl;
	        for(const auto& s : nuis) std::cout<<s<<std::endl;
	        NuisPlotter(rootFile,pdfN,nuis,outDir+bkg+"_"+catNoSel+"_nuisHists.root");
	    };

	    go(bkgSels[BKG_QG],{
	    		systName(bkgSels[BKG_QG],"PTX" ,b+"_2l"),
				systName(bkgSels[BKG_QG],"OPTX",b+"_2l"),
				systName(bkgSels[BKG_QG],"PTY",b+"_2l"),
				systName(bkgSels[BKG_QG],"OPTY",b+"_2l")
	    });
	    go(bkgSels[BKG_LOSTTW],{
	            systName(bkgSels[BKG_LOSTTW],"PTX" ,b+"_2l"),
				systName(bkgSels[BKG_LOSTTW],"OPTX",b)+"_2l",
	            systName("top","res",catNoSel),
				systName("top","scale",catNoSel),
				systName("top","lostmw_rel_scale",b+"_2l")
	    });
	    go(bkgSels[BKG_MT],{
	            "hbb_scale_top",
				"hbb_res_top",
	            systName("top","res",catNoSel),
				systName("top","scale",catNoSel),
				systName("top","mt_rel_scale",b+"_2l")
	    });
	    go(bkgSels[BKG_MW],{
	            "hbb_scale",
				"hbb_res",
	            systName("top","res",catNoSel),
				systName("top","scale",catNoSel),
				systName("top","lostmw_rel_scale",b+"_2l")
	    });
	}

}


void makeNuisPlots(const std::string& baseDir, std::vector<TObject*>& writeables){

    std::string h = selCuts1[SEL1_FULL];
	for(const auto& l : lepCats) for(const auto& b : btagCats) for(const auto& p : purCats) {
		if(l == lepCats[LEP_EMU]) continue;
		if(b == btagCats[BTAG_LMT]) continue;
		if(p == purCats[PURE_I]) continue;

	    const std::string cat = l +"_"+b+"_"+p +"_"+h;
	    const std::string catNoSel = l +"_"+b+"_"+p;

	    auto systName = [&](const std::string& proc,const std::string& name, const std::string& sel = "-1")->std::string {
	        return getSystName("",proc,name, sel == "-1" ? cat : sel  );
	    };

	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_QG]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_QG],cat,"mbb",
	    		{ {"#it{m}_{b#bar{b}} scale",systName(bkgSels[BKG_QG],"PTX",b+"_1l")},
	    		  {"#it{m}_{b#bar{b}} res." ,systName(bkgSels[BKG_QG],"OPTX",b+"_1l")} },true));

	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_QG]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_QG],cat,"mhh",
	    		{ {"#it{m}_{HH} scale" ,systName(bkgSels[BKG_QG],"PTY",catNoSel)},
	    		  {"#it{m}_{HH} res." ,systName(bkgSels[BKG_QG] ,"OPTY",catNoSel)} },false));


	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_LOSTTW]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_LOSTTW],cat,"mbb",
	    		{ {"#it{m}_{b#bar{b}} scale",systName(bkgSels[BKG_LOSTTW],"PTX",b+"_1l")},
	    		  {"#it{m}_{b#bar{b}} res."  ,systName(bkgSels[BKG_LOSTTW],"OPTX",b+"_1l")}    },true));

	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_LOSTTW]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_LOSTTW],cat,"mhh",
	    		{{"top res.",systName("top","res",catNoSel)},
	    		 {"top scale" ,systName("top","scale",catNoSel)},
				 {"rel. top scale" ,systName("top","lostmw_rel_scale",b+"_1l")}    },false));


	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_MW]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_MW],cat,"mbb",
	    		{{"b#bar{b} jet SD mass scale","hbb_scale"},
	    		 {"b#bar{b} jet SD mass res.","hbb_res"}},true));

	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_MW]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_MW],cat,"mhh",
	    		{{"top res." ,systName("top","res",catNoSel)},
	    		 {"top scale" ,systName("top","scale",catNoSel)},
				 {"rel. top scale" ,systName("top","lostmw_rel_scale",b+"_1l")}    },false));


	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_MT]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_MT],cat,"mbb",
	    		{{"b#bar{b} jet SD mass scale","hbb_scale_top"},
	    		 {"b#bar{b} jet SD mass res.","hbb_res_top"}},true));

	    //  remove rel because it is the same size as the standard scale
		//    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_MT]+"_nuisHists.root",bkgSels[BKG_MT],cat,"mhh",{{"top res." ,systName("top","res"  )}, {"top scale" ,systName("top","scale")}, {"rel. top scale" ,systName("top","mt_rel_scale",b)}    },false));
	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_MT]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_MT],cat,"mhh",
	    		{{"#it{m}_{HH} res." ,systName("top","res",catNoSel)},
	    		 {"#it{m}_{HH} scale" ,systName("top","scale",catNoSel)}   },false));
	}

    h = selCuts2[SEL2_FULL];
	for(const auto& l : dilepCats) for(const auto& b : btagCats) {
		if(l == dilepCats[LEP_INCL]) continue;
		if(b == btagCats[BTAG_LMT]) continue;

	    const std::string cat = l +"_"+b+"_"+h;
	    const std::string catNoSel = l +"_"+b;

	    auto systName = [&](const std::string& proc,const std::string& name, const std::string& sel = "-1")->std::string {
	        return getSystName("",proc,name, sel == "-1" ? cat : sel  );
	    };

	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_QG]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_QG],cat,"mbb",
	    		{ {"#it{m}_{b#bar{b}} scale" ,systName(bkgSels[BKG_QG],"PTX",b+"_2l")},
	    		  {"#it{m}_{b#bar{b}} res." ,systName(bkgSels[BKG_QG],"OPTX",b+"_2l")}    },true));

	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_QG]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_QG],cat,"mhh",
	    		{ {"#it{m}_{HH} scale" ,systName(bkgSels[BKG_QG],"PTY",b+"_2l")},
	    		  {"#it{m}_{HH} res." ,systName(bkgSels[BKG_QG] ,"OPTY",b+"_2l")}    },false));


	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_LOSTTW]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_LOSTTW],cat,"mbb",
	    		{ {"#it{m}_{b#bar{b}} scale" ,systName(bkgSels[BKG_LOSTTW],"PTX",b+"_2l")},
	    		  {"#it{m}_{b#bar{b}} res." ,systName(bkgSels[BKG_LOSTTW]  ,"OPTX",b+"_2l")}    },true));

	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_LOSTTW]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_LOSTTW],cat,"mhh",
	    		{{"top res." ,systName("top","res",catNoSel  )},
	    		 {"top scale" ,systName("top","scale",catNoSel)},
				 {"rel. top scale" ,systName("top","lostmw_rel_scale",b+"_2l")}    },false));


	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_MW]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_MW],cat,"mbb",
	    		{{"b#bar{b} jet SD mass scale","hbb_scale"},
	    		 {"b#bar{b} jet SD mass res.","hbb_res"}},true));

	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_MW]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_MW],cat,"mhh",
	    		{{"top res." ,systName("top","res",catNoSel  )},
	    		 {"top scale" ,systName("top","scale",catNoSel)},
				 {"rel. top scale" ,systName("top","lostmw_rel_scale",b+"_2l")}    },false));

	    //  remove rel because it is the same size as the standard scale
	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_MT]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_MT],cat,"mbb",
	    		{{"b#bar{b} jet SD mass scale","hbb_scale_top"},
	    		 {"b#bar{b} jet SD mass res.","hbb_res_top"}},true));

	//    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_MT]+"_nuisHists.root",bkgSels[BKG_MT],cat,"mhh",{{"top res." ,systName("top","res"  )}, {"top scale" ,systName("top","scale")}, {"rel. top scale" ,systName("top","mt_rel_scale",b)}    },false));
	    writeables.push_back(makeNuisPlot(baseDir+bkgSels[BKG_MT]+"_"+catNoSel+"_nuisHists.root",bkgSels[BKG_MT],cat,"mhh",
	    		{{"#it{m}_{HH} res." ,systName("top","res",catNoSel)},
	    		 {"#it{m}_{HH} scale" ,systName("top","scale",catNoSel)}   },false));
	}
}


void plotNuisTests(int step = 0, int inreg = REG_SR,  std::string limitBaseName = ""){
    REGION reg = REGION(inreg);

    std:: string inName =  "bkgInputs" ;
    auto srList = getSRList(reg);
    auto srListTitles = getSRListTitles(reg);
    std::string outName = limitBaseName +"/plots/";

    if(reg == REG_TOPCR){
        inName =  "bkgInputsTopCR";
        hhFilename +="_TopCR";
        outName=limitBaseName+"/plots/TopCR_";
    }
    else if(reg == REG_NONTOPCR){
        inName =  "bkgInputsNonTopCR";
        hhFilename +="_NonTopCR";
        outName=limitBaseName+"/plots/NonTopCR_";
        btagCats = qgBtagCats;
    }

    std::string filename = inName +"/"+hhFilename;
    std::string postFitFilename = limitBaseName +"/postFit.root";

    if(step==0) {//run post fit
        getNuisHists(limitBaseName+"/combined.root",outName);

    }
    if(step==1){
        makeNuisPlots(outName,writeables);
        Dummy d(outName+"nuisPlots");

    }
}







