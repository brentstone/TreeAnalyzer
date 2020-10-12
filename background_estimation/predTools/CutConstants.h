#ifndef ANALYSISTREEMAKER_BACKGROUNDESTIMATION_CUTCONSTANTS_H
#define ANALYSISTREEMAKER_BACKGROUNDESTIMATION_CUTCONSTANTS_H
#include<vector>
#include<utility>
#include <string>
#include "AnalysisSupport/Utilities/interface/Types.h"
namespace CutConstants{

class CutStr : public std::string{
public:
	CutStr() {}
    CutStr(std::string name,std::string cut) : std::string(name), cut(cut){}
    CutStr(std::string name,std::string cut, std::string title) : std::string(name),
            cut(cut), title(title){}

    std::string cut;
    std::string title;
};
std::string hhFilename    = "HHbb1o2l";

enum PROC  {TTBAR,WJETS,QCD,OTHER};
std::vector<CutStr > processes = {
        CutStr("ttbar"     ,"process==2","t#bar{t}"),
        CutStr("wjets"     ,"process==3","W+jets"),
        CutStr("qcd"       ,"process==8","Multijet"),
        CutStr("other"     ,"(process>1&&!(process==2||process==3||process==8))","Other SM")
};

enum REGION  {REG_SR, REG_TOPCR, REG_NONTOPCR};

CutStr hemW ("hemW" , "(lepChan==2&&era==2018?(!isMuon1&&lep1ETA<=-1.479&&lep1Phi>=-1.55&&lep1Phi<=-0.9?(21.08/59.74):1.0)*(!isMuon2&&lep2ETA<=-1.479&&lep2Phi>=-1.55&&lep2Phi<=-0.9?(21.08/59.74):1.0):1.0)");
CutStr nomW ("nomW"  ,  "xsec*trig_N*pu_N*lep_N*btag_N*fjbtag_N*(era==2018?1.0:prefire_N)*"+hemW.cut);

//CutStr qgWt_SR ("qgWt_SR" , "((process==3||process==8)?(era==2018?1.278:(era==2017?1.236:(era==2016?1.12:1.0))):1.0)");
//CutStr qgWt_NT ("qgWt_NT" , "((process==3||process==8)?(era==2018?0.924:(era==2017?0.955:(era==2016?0.981:1.0))):1.0)");

CutStr aQCD ("aQCD"  , "process!=8");

CutStr wjjSC("wjjBC" , "(wjjTau2o1<=0.75)");
CutStr hwwLC("hwwLC" , "(hwwLi<=11.0)");
CutStr ptomC  ("ptomC"   , "(hwwPT/hhMass>=0.3)");
CutStr bV   ("bV"    , "nAK4Btags==0");
CutStr abV  ("abV"   , "nAK4Btags!=0");

CutStr dPhiC ("dPhiC"  , "(abs(llMetDphi)<(3.14159/2))");
CutStr mllV ("mllV"  , "(dilepMass>=6&&dilepMass<=75)");
CutStr metC ("metC"  , "met>=85");
CutStr dRC   ("dRC"    , "dilepDR<=1.0");
CutStr drCrC  ("drCrC"   , "dilepDR>=0.4");

CutStr aHEM ("aHEM","((era==2018&&run>=319077)?(isMuon1==0?(lep1ETA>-1.479||lep1Phi<-1.55||lep1Phi>-0.9):1.0)&&(isMuon2==0?(lep2ETA>-1.479||lep2Phi<-1.55||lep2Phi>-0.9):1.0):1.0)");

CutStr preSel1("preSel1"  , "lepChan==1");
CutStr preSel2("preSel2"  , "(lepChan==2&&"+aHEM.cut+")");

CutStr hbbMCS("hbbMass","hbbMass","#it{m}_{b#bar{b}} [GeV]");
CutStr hbbUpMCS("hbbMassUp","hbbMassUp","#it{m}_{b#bar{b}} [GeV]");
CutStr hbbDnMCS("hbbMassDown","hbbMassDown","#it{m}_{b#bar{b}} [GeV]");
CutStr hhMCS ("hhMass" ,"hhMass","#it{m}_{HH} [GeV]");
//CutStr hhMCS ("hhMass" ,"((lepChan==1)*hhMassBasic+(lepChan==2)*hhMass)","#it{m}_{HH} [GeV]");

unsigned int nHbbMassBins   =30;
double minHbbMass = 30  ;
double maxHbbMass = 210 ;
unsigned int nHHMassBins   =86;
double minHHMass  = 700;
double maxHHMass  = 4000;
std::vector<double> hhMassBins = {
        700,725,750,775,800,825,850,875,900,925,950,975,1000,1025,1050,1075,1100,1125,1150,1175,
        1200,1225,1250,1275,1300,1325,1350,1375,1400,1425,1450,1475,1500,1525,1550,1575,1600,1625,
        1650,1675,1700,1725,1750,1775,1800,1825,1850,1875,1900,1925,1950,1975,2000,2050,2100,2150,
        2200,2250,2300,2350,2400,2450,2500,2550,2600,2650,2700,2750,2800,2850,2900,2950,3000,3050,
        3100,3175,3250,3325,3400,3475,3550,3625,3700,3775,3850,3925,4000
};


unsigned int nInclHHMassBins   =128;
double minInclHHMass  = 0   ;
double maxInclHHMass  = 5050;
std::vector<double> inclHHMassBins = {0,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,
        400,425,450,475,500,525,550,575,600,625,650,675,
        700,725,750,775,800,825,850,875,900,925,950,975,1000,1025,1050,1075,1100,1125,1150,1175,
        1200,1225,1250,1275,1300,1325,1350,1375,1400,1425,1450,1475,1500,1525,1550,1575,1600,1625,
        1650,1675,1700,1725,1750,1775,1800,1825,1850,1875,1900,1925,1950,1975,2000,2050,2100,2150,
        2200,2250,2300,2350,2400,2450,2500,2550,2600,2650,2700,2750,2800,2850,2900,2950,3000,3050,
        3100,3175,3250,3325,3400,3475,3550,3625,3700,3775,3850,3925,4000,
        4075,4150,4225,4300,4375,4450,4525,4600,4675,4750,4825,4900,4975,5050
};
unsigned int nInclHbbMassBins   =42;
double minInclHbbMass  = 0   ;
double maxInclHbbMass  = 252;


CutStr hhRange  ("hhRange", hhMCS.cut+">"+ASTypes::flt2Str(minHHMass)+
        "&&"+hhMCS.cut+"<"+ASTypes::flt2Str(maxHHMass));
CutStr hhInclRange  ("hhInclRange" , hhMCS.cut+">"+ASTypes::flt2Str(minInclHHMass)+
        "&&"+hhMCS.cut+"<"+ASTypes::flt2Str(maxInclHHMass));
CutStr hbbRange ("hbbRange", hbbMCS.cut+">"+ASTypes::flt2Str(minHbbMass)+
        "&&"+hbbMCS.cut+"<"+ASTypes::flt2Str(maxHbbMass));
CutStr hbbInclRange  ("hbbInclRange", hbbMCS.cut+">"+ASTypes::flt2Str(minInclHbbMass)+
        "&&"+hbbMCS.cut+"<"+ASTypes::flt2Str(maxInclHbbMass));

std::string getHbbBinningString(bool inclusive){
    using namespace ASTypes;
    if(inclusive)
        return int2Str(nInclHbbMassBins)+","+flt2Str(minInclHbbMass)+","+flt2Str(maxInclHbbMass);
    else
        return int2Str(nHbbMassBins)+","+flt2Str(minHbbMass)+","+flt2Str(maxHbbMass);
}
std::string getHHBinningString(bool inclusive){
    using namespace ASTypes;
    const std::vector<double> * bins = &(inclusive ? inclHHMassBins : hhMassBins);
    std::string outString="";
    for(unsigned int iB = 0; iB < bins->size(); ++iB){
        if(iB) outString += ",";
        outString += flt2Str((*bins)[iB]);
    }
    return outString;
}

enum BKGModels  {BKG_QG, BKG_LOSTTW, BKG_MW, BKG_MT};
std::vector<CutStr > bkgSels = {
        CutStr("qg"    ,"(hbbDecayType==0)","q/g bkg."),
        CutStr("losttw","((hbbDecayType>0&&hbbDecayType<=3)||hbbDecayType>5)","Lost t/W bkg."),
        CutStr("mw"     ,"(hbbDecayType==4)","#it{m}_{W} bkg."),
        CutStr("mt"     ,"(hbbDecayType==5)","#it{m}_{t} bkg.")
};

enum LEPCats  {LEP_EMU, LEP_E, LEP_MU};
std::vector<CutStr> lepCats = {
        CutStr("emu","(1.0)","e#mu"),
        CutStr("e"  ,"(isMuon1==0)","e"),
        CutStr("mu" ,"(isMuon1==1)","#mu")
};

enum DILEPCats  {LEP_INCL, LEP_SF, LEP_OF};
std::vector<CutStr> dilepCats = {
        CutStr("IF","(1.0)","incl flavor"),
        CutStr("SF","(isMuon1==isMuon2)","SF"),
        CutStr("OF","(isMuon1!=isMuon2)","OF"),
};

enum BTAGCats  {BTAG_LMT, BTAG_L, /*BTAG_M,*/ BTAG_T};
std::vector<CutStr > btagCats = {
//        CutStr("LMT","(hbbCSVCat>=4)","bLMT"),
//        CutStr("L"  ,"(hbbCSVCat==4)","bL"),
//        CutStr("M"  ,"(hbbCSVCat==5)","bM"),
//        CutStr("T"  ,"(hbbCSVCat==6)","bT")
        CutStr("LMT","hbbTag>=0.8","bLMT"),
        CutStr("L"  ,"(hbbTag>=0.8&&hbbTag<0.97)","bL"),
//        CutStr("M"  ,"(hbbTag>=0.9&&hbbTag<0.98)","bM"),
        CutStr("T"  ,"hbbTag>=0.97","bT")
};

std::vector<CutStr > qgBtagCats = {
//        CutStr("LMT","(hbbCSVCat==1)",""),
//        CutStr("L"  ,"(hbbCSVCat==1)","")
        CutStr("LMT","(hbbTag<=0.05)",""),
        CutStr("L"  ,"(hbbTag<=0.05)","")
};
CutStr inclBtagCat("I","(hbbTag>=0.0)");

enum   PURCats {PURE_I, PURE_LP, PURE_HP};
std::vector<CutStr > purCats = {
        CutStr("I","(1.0)","IP"),
        CutStr("LP" ,"(hwwLi>=2.5||(wjjTau2o1>=(era==2016?0.55:0.45)))","LP"),
        CutStr("HP" ,"(hwwLi<2.5&&(wjjTau2o1<(era==2016?0.55:0.45)))"  ,"HP")

};

enum SELCuts1  {SEL1_NONE,SEL1_LB,SEL1_LP,SEL1_LPB,SEL1_LTH,SEL1_LTMB,SEL1_FULL};
std::vector<CutStr > selCuts1 = {
        CutStr("none",preSel1.cut,"-ExB -#it{p}_{T}/#it{m} -#tau_{21} -hwwLi"),
        CutStr("lb"  ,preSel1.cut+"&&"+ptomC.cut+"&&"+wjjSC.cut+"&&"+hwwLC.cut,"-ExB"),
        CutStr("lp"  ,preSel1.cut+"&&"+wjjSC.cut+"&&"+hwwLC.cut+"&&"+bV.cut,"-#it{p}_{T}/#it{m}"),
        CutStr("lpb"  ,preSel1.cut+"&&"+wjjSC.cut+"&&"+hwwLC.cut,"-ExB -#it{p}_{T}/#it{m}"),
        CutStr("lth"  ,preSel1.cut+"&&"+ptomC.cut+"&&"+bV.cut,"-#tau_{21} -hwwLi"),
        CutStr("ltmb",preSel1.cut+"&&"+ptomC.cut,"-ExB -#tau_{21} -hwwLi"),
        CutStr("full",preSel1.cut+"&&"+ptomC.cut+"&&"+wjjSC.cut+"&&"+hwwLC.cut+"&&"+bV.cut,"")

};

enum SELCuts2  {SEL2_NONE,SEL2_MLL,SEL2_LRPB,SEL2_RPhiB,SEL2_FULL};
std::vector<CutStr > selCuts2 = {
        CutStr("none"   ,preSel2.cut,"-ExB -#it{#DeltaR}_{ll} -#it{M}_{ll} -#it{MET} -#it{#Delta#Phi}_{met,ll}"),
		CutStr("mll"    ,preSel2.cut+"&&"+mllV.cut,"-dR -ExB -Met -dphi"),
		CutStr("lrpb"   ,preSel2.cut+"&&"+mllV.cut+"&&"+metC.cut,"-dR -ExB -dphi"),
		CutStr("R_phi_b",preSel2.cut+"&&"+dRC.cut+"&&"+mllV.cut+"&&"+metC.cut,"-dphi -ExB"),
        CutStr("full"   ,preSel2.cut+"&&"+bV.cut+"&&"+dRC.cut+"&&"+dPhiC.cut+"&&"+mllV.cut+"&&"+metC.cut,"")
};

struct CatIterator{
    bool firstBin = true;
    LEPCats  l =LEP_EMU;
    BTAGCats b =BTAG_LMT;
    PURCats  p =PURE_I;
    SELCuts1  h =SEL1_NONE;
    bool is(const LEPCats  cl ) const {return l == cl;}
    bool is(const BTAGCats cb ) const {return b == cb;}
    bool is(const PURCats  cp ) const {return p == cp;}
    bool is(const SELCuts1  ch ) const {return h == ch;}
    std::string name() const {
        return lepCats[l] +"_"+btagCats[b]+"_"+purCats[p] +"_"+selCuts1[h];
    }
    std::string cut() const {
        return "("+lepCats[l].cut+"&&"+btagCats[b].cut+"&&"+purCats[p].cut +"&&"+selCuts1[h].cut+")";
    }
    bool getBin() {
        if(firstBin){
            firstBin = false;
            return true;
        }
        if(h < selCuts1.size()-1){
            h = SELCuts1(h+1);
            return true;
        }
        else if(p < purCats.size()-1){
            p= PURCats(p+1);
            h=SEL1_NONE;
            return true;
        }
        else if(b < btagCats.size()-1){
            b= BTAGCats(b+1);
            p=PURE_I;
            h=SEL1_NONE;
            return true;
        }
        else if(l < lepCats.size()-1){
            l= LEPCats(l+1);
            b=BTAG_LMT;
            p=PURE_I;
            h=SEL1_NONE;
            return true;
        }
        return false;
    }
    void reset() {
        firstBin = true;
        l =LEP_EMU;
        b =BTAG_LMT;
        p =PURE_I;
        h =SEL1_NONE;
    }
};

struct DilepCatIterator{
    bool firstBin = true;
    DILEPCats  l =LEP_INCL;
    BTAGCats b =BTAG_LMT;
    SELCuts2  s =SEL2_NONE;
    bool is(const DILEPCats  cl ) const {return l == cl;}
    bool is(const BTAGCats cb ) const {return b == cb;}
    bool is(const SELCuts2  cs ) const {return s == cs;}
    std::string name() const {
        return dilepCats[l] +"_"+btagCats[b] +"_"+selCuts2[s];
    }
    std::string cut() const {
        return "("+dilepCats[l].cut+"&&"+btagCats[b].cut +"&&"+selCuts2[s].cut+")";
    }
    bool getBin() {
        if(firstBin){
            firstBin = false;
            return true;
        }
        if(s < selCuts2.size()-1){
            s = SELCuts2(s+1);
            return true;
        }
        else if(b < btagCats.size()-1){
            b= BTAGCats(b+1);
            s=SEL2_NONE;
            return true;
        }
        else if(l < dilepCats.size()-1){
            l= DILEPCats(l+1);
            b=BTAG_LMT;
            s=SEL2_NONE;
            return true;
        }
        return false;
    }
    void reset() {
        firstBin = true;
        l =LEP_INCL;
        b =BTAG_LMT;
        s =SEL2_NONE;
    }
};

std::string getCategoryLabel(const LEPCats lep, const BTAGCats btag, const PURCats pur ){
    return lepCats[lep].title+", "+btagCats[btag].title+", "+purCats[pur].title;
}

std::string getCategoryLabel(const DILEPCats lep, const BTAGCats btag){
    return dilepCats[lep].title+", "+btagCats[btag].title;
}

std::string getCategoryLabel(const std::string& inStr, bool do1lep){


    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(inStr);
    while (std::getline(tokenStream, token, '_'))
    {tokens.push_back(token);}
    std::string title = "";
    auto getTitle = [&](const std::string in, const std::vector<CutStr >& opts ) -> std::string{
        for(const auto& c : opts) if(in == c ) return c.title;
        return "";
    };
    std::string lep,btag,pur,ex;
    if(do1lep) {
        if(tokens.size()){lep = getTitle(tokens[0],lepCats);}
        if(tokens.size()>1){btag = getTitle(tokens[1],btagCats);}
        if(tokens.size()>2){pur = getTitle(tokens[2],purCats);}
        if(tokens.size()>3){ex = getTitle(tokens[3],selCuts1);}
    } else {
        if(tokens.size()){lep = getTitle(tokens[0],dilepCats);}
        if(tokens.size()>1){btag = getTitle(tokens[1],btagCats);}
        if(tokens.size()>2){ex = getTitle(tokens[2],selCuts2);}
    }

    if(lep.size()){
        title += lep;
        if(btag.size() || pur.size()) title+=", ";
    }
    if(btag.size()){
        title += btag;
        if(pur.size())
            title += ", ";
    }
    if(pur.size()){
        title += pur;
        if(title == lepCats[LEP_EMU].title + ", "+btagCats[BTAG_LMT].title +
                ", "+purCats[PURE_I].title )
            title = "All categories";
        if(ex.size())
            title += ", ";
    }


    title += ex;
    return title;
}



std::vector<double> resPTBins = {600,700,750,800,850,900,1000,1100,1250,1500,1750,
                                 2000,2500,3000,3500,4000};

enum SIGNALS  {RADION,BLKGRAV,ALLSIGNAL};
std::vector<CutStr > signals = {
        CutStr("radHH"     ,"spin0","radion"),
        CutStr("blkHH"     ,"spin2","bulk graviton"),
        CutStr("allHH"     ,"spinE","combined radion and bulk graviton")
};
std::vector<std::vector<int> > signalMassBins = {
		{800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500},
		{800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500},
		{800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500}
};

//Constants for models when building limits
std::string MOD_MJ("MJ");
std::string MOD_MR("MR");
std::string MOD_MS("MH");

CutStr sigMCS("mx","mx","#it{m}_{#it{X}} [GeV]");
//https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG#Higgs_cross_sections_and_decay_b
const double HtoBBBF = 0.5824;
const double HtoWWBF = 0.2137;
const double HtoZZBF = 0.02619;
const double HHtobbVVBF = 2*HtoBBBF*(HtoWWBF+HtoZZBF);


}


#endif

