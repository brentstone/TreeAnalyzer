
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "makeDileptonBKGInputs.C"
#endif

void makeBKGCRInputs(bool doTopRegion = true, int bkgToDo = BKG_NOM, std::string treeDir = "../BETrees/"){
    if(doTopRegion){
        hadCuts[HAD_NONE].cut = preSel.cut;
        hadCuts[HAD_RPhiB].cut = preSel.cut+"&&"+dR.cut+"&&"+mllV.cut+"&&"+metC.cut;
        hadCuts[HAD_FULL].cut = preSel.cut+"&&"+dR.cut+"&&"+dPhi.cut+"&&"+mllV.cut+"&&"+metC.cut+"&&"+abV.cut;


//        hadCuts[HAD_NONE].cut = preSel.cut;
//        hadCuts[HAD_LB].cut   = preSel.cut+"&&"+wjjBC.cut+"&&"+exA.cut;
//        hadCuts[HAD_LT].cut   = preSel.cut+ "&&"+abV.cut+"&&"+exA.cut;
//        hadCuts[HAD_LTMB].cut = preSel.cut +"&&"+exA.cut;
//        hadCuts[HAD_FULL].cut = preSel.cut + "&&"+abV.cut+"&&"+wjjBC.cut+"&&"+exA.cut;


        hhFilename +="_TopCR";
        go(bkgToDo,treeDir);
    } else {
        btagCats = qgBtagCats;
        hhFilename +="_QGCR";
        go(bkgToDo,treeDir);
    }

}
