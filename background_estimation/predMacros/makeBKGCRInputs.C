
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "makeBKGInputs.C"
#endif

void makeBKGCRInputs(bool doTopRegion = true, int bkgToDo = BKG_QG, int channel = 0, bool doQuick = true, std::string treeDir = "../"){
	if(doTopRegion){
	    selCuts1[SEL1_NONE].cut = preSel1.cut;
//	    selCuts1[SEL1_LB].cut   = preSel1.cut+ "&&"+wjjBC.cut;
	    selCuts1[SEL1_LTH].cut   = preSel1.cut+ "&&"+abV.cut;
	    selCuts1[SEL1_LTMB].cut = preSel1.cut;
	    selCuts1[SEL1_FULL].cut = preSel1.cut + "&&"+abV.cut+"&&"+wjjSC.cut+"&&"+hwwLC.cut;

        selCuts2[SEL2_NONE].cut  = preSel2.cut;
        selCuts2[SEL2_RPhiB].cut = preSel2.cut+"&&"+abV.cut+"&&"+drCrC.cut+"&&"+mllV.cut+"&&"+metC.cut;
        selCuts2[SEL2_FULL].cut  = preSel2.cut+"&&"+abV.cut+"&&"+drCrC.cut+"&&"+mllV.cut+"&&"+metC.cut+"&&"+dPhiC.cut;

        // QG scaling for TopCR
//    	nomW.cut = nomW.cut+"*"+qgWt_SR.cut;

	    hhFilename +="_TopCR";
	    go(bkgToDo,channel,treeDir+"/bkgCompLMT/",doQuick);
	} else {
	    btagCats = qgBtagCats;
	    // QG scaling for NonTop
//		nomW.cut = nomW.cut+"*"+qgWt_NT.cut;
	    hhFilename +="_NonTopCR";
	    go(bkgToDo,channel,treeDir+"/bkgCompAB/",doQuick);
	}

}
