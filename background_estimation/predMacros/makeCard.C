#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../predTools/DataCardMaker.h"
#include "../predTools/SystConstants.h"
using namespace CutConstants;
using namespace SystConstants;
using namespace ASTypes;

std::vector<double> hhMassRebins = {
        700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,
        1700,1750,1800,1850,1900,1950,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,
        3250,3400,3550,3700,3850,4000
};


std::string getSystName (const std::string& prefix, const std::string& proc,const std::string& name, const std::string& sel){
    std::string sn = prefix;

    if(proc.size() ){
        if(sn.size()) sn+="_";
        sn += proc;
    }
    if(name.size() ){
        if(sn.size()) sn+="_";
        sn += name;
    }
    if(sel.size() ){
        if(sn.size()) sn+="_";
        sn += sel;
    }
    return sn;
};

//Make search variables
void makeSearchVars(DataCardMaker& card) {
    card.addVar(MOD_MJ,100,0,1000,false);
    card.addVar(MOD_MR,1000,0,10000,false);
    card.addVar(MOD_MS,2000,true);
};

void addSDMassSysts(DataCardMaker& card) {
    card.addParamSystematic("hbb_scale",    sf_sdmass_scale_Run2-1.0,unc_sdmass_scale_Run2);
    card.addParamSystematic("hbb_scale_top",sf_sdmass_scale_Run2-1.0,0.01*sf_sdmass_scale_Run2);
    card.addParamSystematic("hbb_res",    sf_sdmass_res_Run2-1.0,unc_sdmass_res_Run2);
    card.addParamSystematic("hbb_res_top",sf_sdmass_res_Run2-1.0,unc_sdmass_res_Run2);
};

void addJECSysts(DataCardMaker& card) {
    card.addParamSystematic("unclust",0.0,0.01); // shifting MET res
    card.addParamSystematic("jes",0.0,0.01);     // shifting JES
    card.addParamSystematic("jer",0.0,0.01);     // shifting JER
    card.addParamSystematic("hem",0.0,0.01);     // HEM effect in 2018 using jets
}

void addSignalNormSysts(DataCardMaker& card, std::string sigName, bool is1l, TString fullcat) {
	bool isEorSF = is1l ? fullcat.BeginsWith("e") : fullcat.BeginsWith("SF");
	bool isLoose = fullcat.Contains("_L_");
	bool isLP = is1l ? fullcat.Contains("LP") : false;

    card.addSystematic("lumi","lnN",{{sigName,1.018}});         // from Lumi POG TWIKI
    card.addSystematic("prefire","lnN",{{sigName,1.005}});      // prefire weights (1l could be 0.4% but this is fine)

    if(is1l) {
        card.addSystematic("trigger","lnN",{{sigName,1.02}});   // trigger measurement
        card.addSystematic("pileup" ,"lnN",{{sigName,1.01}});   // pileup weights

        // lepton
    	if(isEorSF) {
    	    card.addSystematic("e_Reco" ,"lnN",{{sigName,1.005}}); // measurement -> eReco weights
    	    card.addSystematic("e_ID_1l" ,"lnN",{{sigName,1.042}});   // measurement -> eID weights
    	    card.addSystematic("e_ISO" ,"lnN",{{sigName,1.06}});   // assertion in AN
    	} else {
    	    card.addSystematic("mu_ID_1l" ,"lnN",{{sigName,1.023}});  // measurement -> muID weights
    	    card.addSystematic("mu_ISO" ,"lnN",{{sigName,1.06}});  // assertion in AN
    	}

    	// tau21
    	if(isLP) {
    		card.addSystematic("tau21_eff","lnN",{{sigName,0.727*(1-unc_tau21_LP_Run2)+(1-0.727)*(1+unc_tau21_HP_Run2)}});
    	} else {
    		card.addSystematic("tau21_eff","lnN",{{sigName,1+unc_tau21_HP_Run2}});
    	}

    } else {
        card.addSystematic("trigger","lnN",{{sigName,1.03}});   // trigger measurement
        card.addSystematic("pileup" ,"lnN",{{sigName,1.006}});   // pileup weights

        // lepton
    	if(isEorSF) {
    	    card.addSystematic("e_Reco" ,"lnN",{{sigName,1.0+(0.008*0.46)}});     // measurement -> eReco weights * ee fraction
    	    card.addSystematic("e_ID_1l" ,"lnN",{{sigName,1.0+(0.027*0.46)}});    // measurement -> eID weights * ee fraction
    	    card.addSystematic("e_ISO" ,"lnN",{{sigName,1.0+(0.06*0.46)}});       // assertion in AN * ee fraction
    	    card.addSystematic("mu_ID_2l" ,"lnN",{{sigName,1.0+(0.025*0.54)}});   // measurement -> muID weights * mumu fraction
    	    card.addSystematic("mu_ISO" ,"lnN",{{sigName,1.0+(0.04*0.54)}});      // assertion in AN * mumu fraction
    	} else {
    	    card.addSystematic("e_Reco" ,"lnN",{{sigName,1.008}});    // measurement -> eReco weights
    	    card.addSystematic("e_ID_1l" ,"lnN",{{sigName,1.026}});   // measurement -> eID weights
    	    card.addSystematic("e_ISO" ,"lnN",{{sigName,1.02}});      // assertion in AN
    	    card.addSystematic("mu_ID_2l" ,"lnN",{{sigName,1.023}});  // measurement -> muID weights
    	    card.addSystematic("mu_ISO" ,"lnN",{{sigName,1.03}});     // assertion in AN
    	}
    }

    // AK8 b jets
    if(isLoose) {
        card.addSystematic("bAk8_eff","lnN",{{sigName,1.04}});  // DeepAK8-MD Z/H->bb weights
    } else {
        card.addSystematic("bAk8_eff","lnN",{{sigName,1.055}}); // DeepAK8-MD Z/H->bb weights
    }

    // AK4 b jets
    card.addParamSystematic("bAk4_real",0.0,0.01);
    card.addParamSystematic("bAk4_fake",0.0,0.01);

}

std::string getbAk4Form(bool isReal, bool is1l, std::string l="", std::string p="") {
	std::string factor;
	std::string sysVar;
	if(is1l) {
		if(isReal) {
			sysVar = "bAk4_real";
			if(l == lepCats[LEP_MU]) {
	            // LP fit: 0.008105 -1.59657e-05*m + 1.69035e-08*m^2 - 3.1328e-12*m^3
	            // HP fit: 0.0233117 -4.64524e-05*m + 3.45413e-08*m^2 - 5.89989e-12*m^3
				if(p == purCats[PURE_LP]) factor = "(0.81-1.6*10^(-3)*"+MOD_MS+"+1.7*10^(-6)*"+MOD_MS+"^(2)-3.13*10^(-10)*"+MOD_MS+"^(3))";
				if(p == purCats[PURE_HP]) factor = "(2.33-4.65*10^(-3)*"+MOD_MS+"+3.45*10^(-6)*"+MOD_MS+"^(2)-5.9*10^(-10)*"+MOD_MS+"^(3))";
			} else {
				factor = "(1.0)"; // flat 1% for electrons
			}
		} else {
			sysVar = "bAk4_fake";
			factor = "(2.0)"; // flat 2% unc. for fake rate
		}
	} else {
		if(isReal) {
			sysVar = "bAk4_real";
			factor = "(0.5)"; // flat 0.5% unc.
		} else {
			sysVar = "bAk4_fake";
			factor = "(2.5)"; // flat 2.5% unc.
		}
	}

	std::string form = "(1.0+"+sysVar+"*"+factor+")";
	return form;
}

void go(const int insig, const std::string& filename, const std::string& mainDir, REGION reg, int channel){

std::string sigInputDir =  mainDir +  "/signalInputs/";
std::string sfPF = sigInputDir + filename;
std::string signalName = signals[insig];
std::string fPF;

std::string cmd = "combineCards.py ";

sigInputDir = mainDir + "/signalInputs/";
sfPF = sigInputDir + filename;
switch(reg){
case REG_SR:
    fPF = mainDir + "/bkgInputs/" + filename;
    break;
case REG_TOPCR:
    fPF = mainDir + "/bkgInputsTopCR/" + filename + "_TopCR";
    break;
case REG_NONTOPCR:
    fPF = mainDir + "/bkgInputsNonTopCR/" + filename + "_NonTopCR";
    break;
}

if(channel == 0 || channel == 1) {
    for(const auto& l :lepCats) for(const auto& b :btagCats) for(const auto& p :purCats) for(const auto& h :selCuts1){
        if(l == lepCats[LEP_EMU] ) continue;
        if(b == btagCats[BTAG_LMT]) continue;
        if(p == purCats[PURE_I] ) continue;
        if(h != selCuts1[SEL1_FULL] ) continue;

        const std::string cat = l +"_"+b+"_"+p +"_"+h;

        auto card = DataCardMaker(cat,"13TeV",1);
        const std::string dataCardTag = DataCardMaker::getFileNameTag(cat,"13TeV");
        cmd += std::string(" ")+ dataCardTag+"="+DataCardMaker::getOutputCardFileName(dataCardTag);

        auto fullInputName =[&](const std::string& proc, const std::string& l, const std::string& b, const std::string& p, const std::string& h, const std::string& pf) -> std::string
            {return fPF + "_"+proc +"_"+l +"_"+b+"_"+p +"_"+h +"_"+ pf; };
        auto inputName =[&](const std::string& proc, const std::string& pf) -> std::string {return fPF + "_"+proc +"_"+cat +"_"+ pf; };
        auto signalInputName =[&](const std::string& proc, const std::string& pf) -> std::string {return sfPF + "_"+proc +"_"+cat +"_"+ pf; };

        auto systName = [&](const std::string& proc,const std::string& name, const std::string& sel = "-1")->std::string {
            return getSystName("",proc,name, sel == "-1" ? l +"_"+b+"_"+p : sel  );
        };

        //Make search variables
        makeSearchVars(card);
        card.rebinY(hhMassRebins);

        //---------------------------------------------------------------------------------------------------
        //Get rates and contributions for backgrounds
        //---------------------------------------------------------------------------------------------------
        card.addFixedYieldFromFile(bkgSels[BKG_QG],1,fPF+"_"+bkgSels[BKG_QG]+"_distributions.root",bkgSels[BKG_QG]+"_"+cat+"_"+hhMCS,1.0,true);
        double rate_lostTW = card.addFixedYieldFromFile(bkgSels[BKG_LOSTTW],2,fPF+"_"+bkgSels[BKG_LOSTTW]+"_distributions.root",bkgSels[BKG_LOSTTW]+"_"+cat+"_"+hhMCS,1.0,true);
        double rate_mw =     card.addFixedYieldFromFile(bkgSels[BKG_MW],3,fPF+"_"+bkgSels[BKG_MW]+"_distributions.root",bkgSels[BKG_MW]+"_"+cat+"_"+hhMCS,1.0,true);
        double rate_mt =     card.addFixedYieldFromFile(bkgSels[BKG_MT],4,fPF+"_"+bkgSels[BKG_MT]+"_distributions.root",bkgSels[BKG_MT]+"_"+cat+"_"+hhMCS,1.0,true);

        //---------------------------------------------------------------------------------------------------
        //Add Systematics first since the param systs need to have the variables added to the workspace
        //---------------------------------------------------------------------------------------------------

        addSignalNormSysts(card,signalName,true,l+"_"+b+"_"+p); // signal normalization
        addJECSysts(card);       // JECs
        addSDMassSysts(card);    // soft-drop mass, applied to resonant mbb

        // QG bkg KDE shape systematics (X = mbb, Y = mhh -- morphing: PT = m, OPT = 1/m)
        card.addParamSystematic(systName(bkgSels[BKG_QG] ,"PTX",b+"_1l")    ,0.0,0.75);
        card.addParamSystematic(systName(bkgSels[BKG_QG] ,"OPTX",b+"_1l")   ,0.0,0.5);
        card.addParamSystematic(systName(bkgSels[BKG_QG] ,"PTY")            ,0.0,1.0);
        card.addParamSystematic(systName(bkgSels[BKG_QG] ,"OPTY")           ,0.0,1.0);

        // lost tW mbb resolution (OPT = 1/m) and scale (PT = m)
        card.addParamSystematic(systName(bkgSels[BKG_LOSTTW],"PTX" ,b+"_1l"),0.0,0.6);
        card.addParamSystematic(systName(bkgSels[BKG_LOSTTW],"OPTX",b+"_1l"),0.0,0.6);

        // Normalization (correlate top-derived bkgs, separate QG)
        card.addSystematic(systName(bkgSels[BKG_QG],"norm")  ,"lnN",{{bkgSels[BKG_QG],1.5}});
        card.addSystematic(systName("top","norm")            ,"lnN",{{bkgSels[BKG_LOSTTW],1.35},{bkgSels[BKG_MW],1.35},{bkgSels[BKG_MT],1.35}});

        //top HH resolution and scale (correlate all top-derived bkgs: lost tW, mW, mT)
        card.addParamSystematic(systName("top","res"  ) ,0.0,0.50);
        card.addParamSystematic(systName("top","scale") ,0.0,0.50);

        card.addParamSystematic(systName("top","mt_rel_scale",b+"_1l") ,0.0,0.50);
        card.addParamSystematic(systName("top","lostmw_rel_scale",b+"_1l") ,0.0,0.50);

        // avoid negative numbers for lnN systs in NonTopCR
        if(reg == REG_NONTOPCR){
            card.addSystematic(systName("top","tFrac",b+"_1l")         ,"lnN",{{bkgSels[BKG_MT],1.0+0.30},{bkgSels[BKG_MW],fabs(1.0-0.30*rate_mt/rate_mw)}});
            card.addSystematic(systName("top","lostFrac",b+"_1l")      ,"lnN",{{bkgSels[BKG_LOSTTW],1.0+0.30},{bkgSels[BKG_MW],fabs(1.0-0.30*rate_lostTW/rate_mw)}});
        } else {
            card.addSystematic(systName("top","wFrac",b+"_1l")         ,"lnN",{{bkgSels[BKG_MW],1.0+0.30},{bkgSels[BKG_MT],fabs(1.0-0.30*rate_mw/rate_mt)}});
            card.addSystematic(systName("top","lostFrac",b+"_1l")      ,"lnN",{{bkgSels[BKG_LOSTTW],1.0+0.30},{bkgSels[BKG_MT],fabs(1.0-0.30*rate_lostTW/rate_mt)}});
        }
printf("dbg0\n");
        //---------------------------------------------------------------------------------------------------
        //Signal
        //---------------------------------------------------------------------------------------------------
        StrFlts sigHHScaleSysts = {{"unclust",0.1}, {"jes",0.8}, {"jer",0.3}};
        StrFlts sigHHResSysts   = {{"unclust",1.5}, {"jes",3.0}, {"jer",4.0}, {"hem",1.0}};

        card.add2DSignalParametricShape(signalName,MOD_MJ,MOD_MR, signalInputName(signalName,"2D_fit.json"),
					{{"hbb_scale",1}},{{"hbb_res",1}},sigHHScaleSysts,sigHHResSysts, b == btagCats[BTAG_L],MOD_MS);
        printf("dbg0.1\n");

        std::string brealForm = getbAk4Form(true,true,l,p);
        printf("dbg0.2\n");

        std::string bfakeForm = getbAk4Form(false,true,l,p);
        printf("dbg0.3\n");

//        std::string tau21Form = "(1.0+tau21_PtDependence"+"*"+ (p == purCats[PURE_HP] ? "log("+MOD_MS+"/1000)" : "((0.054/0.041)*(-log("+MOD_MS+"/1000)))")+")";

        std::string jetForm = "(1.0+0.5*unclust)*(1.0+0.5*jer)*(1.0+2.0*jes)*(1.0+0.3*hem)";
        std::string uncForm = "("+brealForm+"*"+bfakeForm+"*"+jetForm+")";

        std::cout<<uncForm<<std::endl;

//            card.addParametricYieldWithUncertainty(signalName,0,signalInputName(signalName,"yield.json"),
//            		1.0,"1.0",{},MOD_MS);

        card.addParametricYieldWithUncertainty(signalName,0,signalInputName(signalName,"yield.json"),
                1.0,uncForm,{"bAk4_real","bAk4_fake","unclust","jer","jes","hem"},MOD_MS);

printf("dbg1\n");

        //---------------------------------------------------------------------------------------------------
        //QG
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts qgKDESysts;
        qgKDESysts.addSyst("PTX",{{systName(bkgSels[BKG_QG],"PTX",b+"_1l"),"1"  }});
        qgKDESysts.addSyst("OPTX",{{systName(bkgSels[BKG_QG],"OPTX",b+"_1l"),"1"  }});
        qgKDESysts.addSyst("PTY",{{systName(bkgSels[BKG_QG],"PTY"),"1"  }});
        qgKDESysts.addSyst("OPTY",{{systName(bkgSels[BKG_QG],"OPTY"),"1"  }});
        card.addHistoShapeFromFile(bkgSels[BKG_QG],{MOD_MJ,MOD_MR}, inputName(bkgSels[BKG_QG],"2D_template.root"),"histo",qgKDESysts);
        printf("dbg2\n");

        //---------------------------------------------------------------------------------------------------
        //Lost t/W
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts twKDESysts;
        twKDESysts.addSyst("PTX",{{systName(bkgSels[BKG_LOSTTW],"PTX"  ,b+"_1l"),"1"  }});
        twKDESysts.addSyst("OPTX",{{systName(bkgSels[BKG_LOSTTW],"OPTX",b+"_1l"),"1"  }});
        twKDESysts.addSyst("PTY",{{systName("top","scale"),"1"},{systName("top","lostmw_rel_scale",b+"_1l"),"1"}});
        twKDESysts.addSyst("OPTY",{{systName("top","res"  ),"1"  }});
        card.addHistoShapeFromFile(bkgSels[BKG_LOSTTW],{MOD_MJ,MOD_MR},inputName(bkgSels[BKG_LOSTTW],"2D_template.root"),"histo",twKDESysts);
        printf("dbg3\n");

        //---------------------------------------------------------------------------------------------------
        //mW
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts mwKDESysts;
        mwKDESysts.addSyst("PT",{{systName("top","scale"),"1"},{systName("top","lostmw_rel_scale",b+"_1l"),"1"}});
        mwKDESysts.addSyst("OPT",{{systName("top","res"  ),"1"  }});
        card.add1DBKGParametricShape(bkgSels[BKG_MW],MOD_MJ,inputName(bkgSels[BKG_MW],"MJJ_SFFit.json"),{{"hbb_scale",1}},{{"hbb_res",1}},MOD_MR,MOD_MJ);
        card.addHistoShapeFromFile(bkgSels[BKG_MW],{MOD_MR}, inputName(bkgSels[BKG_MW],"MVV_template.root"),"histo",mwKDESysts,false,0,MOD_MR,true);
        card.conditionalProduct(bkgSels[BKG_MW],bkgSels[BKG_MW] + "_"+MOD_MJ,MOD_MJ,bkgSels[BKG_MW] + "_"+MOD_MR);
        printf("dbg4\n");

        //---------------------------------------------------------------------------------------------------
        //mt
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts mtKDESysts;
        mtKDESysts.addSyst("PT",{{systName("top","scale"),"1"},{systName("top","mt_rel_scale",b+"_1l"),"1"}});
        mtKDESysts.addSyst("OPT",{{systName("top","res"  ),"1"  }});
        card.add1DBKGParametricShape(bkgSels[BKG_MT],MOD_MJ,inputName(bkgSels[BKG_MT],"MJJ_SFFit.json"),{{"hbb_scale_top",1}},{{"hbb_res_top",1}},MOD_MR,MOD_MJ);
        card.addHistoShapeFromFile(bkgSels[BKG_MT],{MOD_MR}, inputName(bkgSels[BKG_MT],"MVV_template.root"),"histo",mtKDESysts,false,0,MOD_MR,true);
        card.conditionalProduct(bkgSels[BKG_MT],bkgSels[BKG_MT] + "_"+MOD_MJ,MOD_MJ,bkgSels[BKG_MT]+ "_"+MOD_MR);

        //---------------------------------------------------------------------------------------------------
        //Data
        //---------------------------------------------------------------------------------------------------
        card.importBinnedData(fPF + "_data_distributions.root","data_"+cat+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});
//        card.importBinnedData(fPF + "_pd.root","data_"+cat+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});

        card.makeCard();
    }
}

if (channel == 0 || channel == 2) {
    for(const auto& l :dilepCats) for(const auto& b :btagCats) for(const auto& s :selCuts2){
        if(l == dilepCats[LEP_INCL] ) continue;
        if(b == btagCats[BTAG_LMT]) continue;
        if(s != selCuts2[SEL2_FULL] ) continue;

        const std::string cat = l +"_"+b +"_"+s;
        auto card = DataCardMaker(cat,"13TeV",1);

        const std::string dataCardTag = DataCardMaker::getFileNameTag(cat,"13TeV");
        cmd += std::string(" ")+ dataCardTag+"="+DataCardMaker::getOutputCardFileName(dataCardTag);

        auto fullInputName =[&](const std::string& proc, const std::string& l, const std::string& b, const std::string& s, const std::string& pf) -> std::string
                {return fPF + "_"+proc +"_"+l +"_"+b +"_"+s +"_"+ pf; };
        auto inputName =[&](const std::string& proc, const std::string& pf) -> std::string {return fPF + "_"+proc +"_"+cat +"_"+ pf; };
        auto signalInputName =[&](const std::string& proc, const std::string& pf) -> std::string {return sfPF + "_"+proc +"_"+cat +"_"+ pf; };

        auto systName = [&](const std::string& proc,const std::string& name, const std::string& sel = "-1")->std::string {
            return getSystName("",proc,name, sel == "-1" ? l +"_"+b : sel  );
        };

        //Make search variables
        makeSearchVars(card);
        card.rebinY(hhMassRebins);

        //---------------------------------------------------------------------------------------------------
        //Get rates and contributions for backgrounds
        //---------------------------------------------------------------------------------------------------
        card.addFixedYieldFromFile(bkgSels[BKG_QG],1,fPF+"_"+bkgSels[BKG_QG]+"_distributions.root",bkgSels[BKG_QG]+"_"+cat+"_"+hhMCS,1.0,true);
        double rate_lostTW = card.addFixedYieldFromFile(bkgSels[BKG_LOSTTW],2,fPF+"_"+bkgSels[BKG_LOSTTW]+"_distributions.root",bkgSels[BKG_LOSTTW]+"_"+cat+"_"+hhMCS,1.0,true);
        double rate_mw =     card.addFixedYieldFromFile(bkgSels[BKG_MW],3,fPF+"_"+bkgSels[BKG_MW]+"_distributions.root",bkgSels[BKG_MW]+"_"+cat+"_"+hhMCS,1.0,true);
        double rate_mt =     card.addFixedYieldFromFile(bkgSels[BKG_MT],4,fPF+"_"+bkgSels[BKG_MT]+"_distributions.root",bkgSels[BKG_MT]+"_"+cat+"_"+hhMCS,1.0,true);

        //---------------------------------------------------------------------------------------------------
        //Add Systematics first since the param systs need to have the variables added to the workspace
        //---------------------------------------------------------------------------------------------------
        addSignalNormSysts(card,signalName,false,l+"_"+b); // signal normalization
        addJECSysts(card);       // JECs
        addSDMassSysts(card);    // soft-drop mass, applied to resonant mbb

        // QG normalization and KDE shape systematics (PT = scale (m), OPT = res (1/m))
        card.addSystematic(systName(bkgSels[BKG_QG],"norm")  ,"lnN",{{bkgSels[BKG_QG],1.5}});
        card.addParamSystematic(systName(bkgSels[BKG_QG],"PTX",b+"_2l")    ,0.0,0.75);
        card.addParamSystematic(systName(bkgSels[BKG_QG],"OPTX",b+"_2l")   ,0.0,0.5);
        card.addParamSystematic(systName(bkgSels[BKG_QG],"PTY",b+"_2l")    ,0.0,1.0);
        card.addParamSystematic(systName(bkgSels[BKG_QG],"OPTY",b+"_2l")   ,0.0,1.0);

        // Lost tW mbb KDE systematics (PTX = scale (m), OPTX = res (1/m))
        card.addParamSystematic(systName(bkgSels[BKG_LOSTTW],"PTX" ,b+"_2l"),0.0,0.6);
        card.addParamSystematic(systName(bkgSels[BKG_LOSTTW],"OPTX",b+"_2l"),0.0,0.6);

        // Top normalization and HH resolution (1/m) and scale (m)
        card.addSystematic(systName("top","norm"),"lnN",{{bkgSels[BKG_LOSTTW],1.35},
        		{bkgSels[BKG_MW],1.35},{bkgSels[BKG_MT],1.35}});
        card.addParamSystematic(systName("top","res"  ) ,0.0,0.75);
        card.addParamSystematic(systName("top","scale") ,0.0,0.75);
        card.addParamSystematic(systName("top","mt_rel_scale",b+"_2l") ,0.0,0.75);
        card.addParamSystematic(systName("top","lostmw_rel_scale",b+"_2l") ,0.0,0.75);

        // Normalization

        if(reg == REG_NONTOPCR){
            card.addSystematic(systName("top","tFrac",b+"_2l")    ,"lnN",{{bkgSels[BKG_MT],1.0+0.30},{bkgSels[BKG_MW],fabs(1.0-0.30*rate_mt/rate_mw)}});
            card.addSystematic(systName("top","lostFrac",b+"_2l") ,"lnN",{{bkgSels[BKG_LOSTTW],1.0+0.30},{bkgSels[BKG_MW],fabs(1.0-0.30*rate_lostTW/rate_mw)}});
        } else {
            card.addSystematic(systName("top","wFrac",b+"_2l")    ,"lnN",{{bkgSels[BKG_MW],1.0+0.30},{bkgSels[BKG_MT],fabs(1.0-0.30*rate_mw/rate_mt)}});
            card.addSystematic(systName("top","lostFrac",b+"_2l") ,"lnN",{{bkgSels[BKG_LOSTTW],1.0+0.30},{bkgSels[BKG_MT],fabs(1.0-0.30*rate_mt/rate_lostTW)}}); // changed to mt / tw
        }

        //---------------------------------------------------------------------------------------------------
        //Signal
        //---------------------------------------------------------------------------------------------------
        StrFlts sigHHScaleSysts = {{"unclust",0.1}, {"jes",0.8}, {"jer",0.3}};
        StrFlts sigHHResSysts   = {{"unclust",1.5}, {"jes",3.0}, {"jer",4.0}, {"hem",1.0}};

        card.add2DSignalParametricShape(signalName,MOD_MJ,MOD_MR, signalInputName(signalName,"2D_fit.json"),
        		{{"hbb_scale",1}},{{"hbb_res",1}},sigHHScaleSysts,sigHHResSysts, b == btagCats[BTAG_L],MOD_MS);
//        card.add2DSignalParametricShape(signalName,MOD_MJ,MOD_MR, signalInputName(signalName,"2D_fit.json"),
//                             {{"hbb_scale",1}},{{"hbb_res",1}},{{"unclust",0.5},{"jes",1},{"jer",0.5}},{{"unclust",0.5},{"jes",2},{"jer",5}}, b == btagCats[BTAG_L],MOD_MS);

        std::string brealForm = getbAk4Form(true,false);
        std::string bfakeForm = getbAk4Form(false,false);

        std::string jetForm = "(1.0+unclust)*(1.0+jer)*(1.0+0.5*jes)";
        std::string uncForm = brealForm+"*"+jetForm;

        card.addParametricYieldWithUncertainty(signalName,0,signalInputName(signalName,"yield.json"),
                1.0,uncForm,{"bAk4_real","bAk4_fake","unclust","jer","jes","hem"},MOD_MS);

//        card.addParametricYieldWithUncertainty(signalName,0,signalInputName(signalName,"yield.json"),
//                1.0,"1.0",{},MOD_MS);
        //---------------------------------------------------------------------------------------------------
        //QG
        //---------------------------------------------------------------------------------------------------
        printf("dbg1\n");

        PDFAdder::InterpSysts qgKDESysts;
        qgKDESysts.addSyst("PTX",{{systName(bkgSels[BKG_QG],"PTX",b+"_2l"),"1"  }});
        qgKDESysts.addSyst("OPTX",{{systName(bkgSels[BKG_QG],"OPTX",b+"_2l"),"1"  }});
        qgKDESysts.addSyst("PTY",{{systName(bkgSels[BKG_QG],"PTY",b+"_2l"),"1"  }});
        qgKDESysts.addSyst("OPTY",{{systName(bkgSels[BKG_QG],"OPTY",b+"_2l"),"1"  }});
        card.addHistoShapeFromFile(bkgSels[BKG_QG],{MOD_MJ,MOD_MR}, inputName(bkgSels[BKG_QG],"2D_template.root"),"histo",qgKDESysts);
        printf("dbg2\n");

        //---------------------------------------------------------------------------------------------------
        //Lost t/W
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts twKDESysts;
        twKDESysts.addSyst("PTX",{{systName(bkgSels[BKG_LOSTTW],"PTX",b+"_2l"),"1"  }});
        twKDESysts.addSyst("OPTX",{{systName(bkgSels[BKG_LOSTTW],"OPTX",b+"_2l"),"1"  }});
        twKDESysts.addSyst("PTY",{{systName("top","scale"),"1"},{systName("top","lostmw_rel_scale",b+"_2l"),"1"}});
        twKDESysts.addSyst("OPTY",{{systName("top","res"  ),"1"  }});
        card.addHistoShapeFromFile(bkgSels[BKG_LOSTTW],{MOD_MJ,MOD_MR},inputName(bkgSels[BKG_LOSTTW],"2D_template.root"),"histo",twKDESysts);
        printf("dbg3\n");

        //---------------------------------------------------------------------------------------------------
        //mW
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts mwKDESysts;
        mwKDESysts.addSyst("PT",{{systName("top","scale"),"1"},{systName("top","lostmw_rel_scale",b+"_2l"),"1"}});
        mwKDESysts.addSyst("OPT",{{systName("top","res"  ),"1"  }});
        card.add1DBKGParametricShape(bkgSels[BKG_MW],MOD_MJ,inputName(bkgSels[BKG_MW],"MJJ_SFFit.json"),{{"hbb_scale",1}},{{"hbb_res",1}},MOD_MR,MOD_MJ);
        card.addHistoShapeFromFile(bkgSels[BKG_MW],{MOD_MR}, inputName(bkgSels[BKG_MW],"MVV_template.root"),"histo",mwKDESysts,false,0,MOD_MR,true);
        card.conditionalProduct(bkgSels[BKG_MW],bkgSels[BKG_MW] + "_"+MOD_MJ,MOD_MJ,bkgSels[BKG_MW] + "_"+MOD_MR);
        printf("dbg4\n");

        //---------------------------------------------------------------------------------------------------
        //mt
        //---------------------------------------------------------------------------------------------------
        PDFAdder::InterpSysts mtKDESysts;
        mtKDESysts.addSyst("PT",{{systName("top","scale"),"1"},{systName("top","mt_rel_scale",b+"_2l"),"1"}});
        mtKDESysts.addSyst("OPT",{{systName("top","res"  ),"1"  }});
        card.add1DBKGParametricShape(bkgSels[BKG_MT],MOD_MJ,inputName(bkgSels[BKG_MT],"MJJ_SFFit.json"),{{"hbb_scale_top",1}},{{"hbb_res_top",1}},MOD_MR,MOD_MJ);
        card.addHistoShapeFromFile(bkgSels[BKG_MT],{MOD_MR}, inputName(bkgSels[BKG_MT],"MVV_template.root"),"histo",mtKDESysts,false,0,MOD_MR,true);
        card.conditionalProduct(bkgSels[BKG_MT],bkgSels[BKG_MT] + "_"+MOD_MJ,MOD_MJ,bkgSels[BKG_MT]+ "_"+MOD_MR);
        printf("dbg5\n");

        //---------------------------------------------------------------------------------------------------
        //Data
        //---------------------------------------------------------------------------------------------------
        card.importBinnedData(fPF + "_data_distributions.root","data_"+cat+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});
//        card.importBinnedData(fPF + "_pd.root","data_"+cat+"_hbbMass_hhMass",{MOD_MJ,MOD_MR});
        printf("dbg6\n");

        card.makeCard();
    }
}

printf("dbgSLURM\n");

std::ofstream outFile("comp.sh",std::ios::out|std::ios::trunc);

outFile << cmd <<" > combinedCard.txt";
outFile.close();

std::cout << cmd <<" > combinedCard.txt"<<std::endl;

}
#endif

void makeCard(int inreg = REG_SR, int insig = RADION, int channel = 1){
	if (channel != 0 && channel != 1 && channel != 2) {
		std::cout<<"channel needs to be either 1 (single lep), 2 (dilep), or 0 (both)"<<std::endl;
		return;
	}
    REGION reg = REGION(inreg);
    if(reg == REG_NONTOPCR) btagCats = qgBtagCats;
    std::string mainDir = "../";
    go(insig,hhFilename,mainDir,reg,channel);
}
