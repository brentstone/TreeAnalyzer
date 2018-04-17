#include "../interface/DefaultSearchRegionAnalyzer.h"

#include "TreeReaders/interface/EventReader.h"
#include "TreeReaders/interface/GenParticleReader.h"
#include "TreeReaders/interface/ElectronReader.h"
#include "TreeReaders/interface/MuonReader.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/FatJetReader.h"
#include "TreeReaders/interface/FillerConstants.h"

#include "Processors/Corrections/interface/EventWeights.h"
#include "Processors/Variables/interface/HiggsSolver.h"
#include "Processors/Variables/interface/JetKinematics.h"

#include "Processors/Variables/interface/LeptonSelection.h"
#include "Processors/Variables/interface/FatJetSelection.h"

#include "Processors/EventSelection/interface/EventSelection.h"
#include "Processors/Corrections/interface/TriggerScaleFactors.h"
#include "Processors/Corrections/interface/LeptonScaleFactors.h"
#include "Processors/Corrections/interface/BTagScaleFactors.h"
#include "Processors/Corrections/interface/FatJetScaleFactors.h"
#include "TPRegexp.h"




namespace TAna {
//--------------------------------------------------------------------------------------------------
DefaultSearchRegionAnalyzer::DefaultSearchRegionAnalyzer(std::string fileName, std::string treeName, int treeInt) : BaseTreeAnalyzer(fileName,treeName,treeInt){
    TPRegexp r1(".*m(\\d+)_[0-9]*\\..*$");
    auto match = r1.MatchS(fileName);
    const Int_t nrSubStr = match->GetLast()+1;
    if(nrSubStr>1){
        signal_mass = (((TObjString *)match->At(1))->GetString()).Atoi();
    }
    fjProc      .reset(new FatJetProcessor ()); DefaultFatJetSelections::setDefaultFatJetProcessor(*fjProc);
    leptonProc  .reset(new LeptonProcessor ()); DefaultLeptonSelections::setDefaultLeptonProcessor(*leptonProc);
    trigSFProc  .reset(new TriggerScaleFactors (dataDirectory));
    puSFProc    .reset(new PUScaleFactors (dataDirectory));
    leptonSFProc.reset(new POGLeptonScaleFactors (dataDirectory));
    ak4btagSFProc.reset(new JetBTagScaleFactors (dataDirectory));
    sjbtagSFProc.reset(new SubJetBTagScaleFactors (dataDirectory));
    hbbFJSFProc .reset(new HbbFatJetScaleFactors (dataDirectory));
    setLumi(35.922); //https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/688.html

    turnOnCorr(CORR_XSEC);
    turnOnCorr(CORR_TRIG);
    turnOnCorr(CORR_PU  );
    turnOnCorr(CORR_LEP );
    turnOnCorr(CORR_SJBTAG);
    turnOnCorr(CORR_AK4BTAG);
    turnOnCorr(CORR_SDMASS);
}
//--------------------------------------------------------------------------------------------------
DefaultSearchRegionAnalyzer::~DefaultSearchRegionAnalyzer(){}
//--------------------------------------------------------------------------------------------------
void DefaultSearchRegionAnalyzer::resetCorr() {corrections = 0;}
bool DefaultSearchRegionAnalyzer::isCorrOn(Corrections corr) const {return FillerConstants::doesPass(corrections,corr);}
void DefaultSearchRegionAnalyzer::turnOnCorr(Corrections corr) {FillerConstants::addPass(corrections,corr);}
void DefaultSearchRegionAnalyzer::turnOffCorr(Corrections corr) {FillerConstants::removePass(corrections,corr);}

//--------------------------------------------------------------------------------------------------
void DefaultSearchRegionAnalyzer::loadVariables()  {
    reader_event   =std::make_shared<EventReader>   ("event",isRealData());             load(reader_event   );
    reader_fatjet  =std::make_shared<FatJetReader>  ("ak8PuppiJet",isRealData());  load(reader_fatjet  );
    reader_fatjet_noLep=std::make_shared<FatJetReader>  ("ak8PuppiNoLepJet",isRealData());  load(reader_fatjet_noLep  );
    reader_jet_chs =std::make_shared<JetReader>     ("ak4Jet",isRealData());            load(reader_jet_chs );
    reader_jet     =std::make_shared<JetReader>     ("ak4PuppiJet",isRealData(),false);  load(reader_jet     );
    reader_electron=std::make_shared<ElectronReader>("electron");                       load(reader_electron);
    reader_muon    =std::make_shared<MuonReader>    ("muon");                           load(reader_muon    );

    if(!isRealData()){
        reader_genpart =std::make_shared<GenParticleReader>   ("genParticle");             load(reader_genpart   );
    }

    checkConfig();
}
//--------------------------------------------------------------------------------------------------
void DefaultSearchRegionAnalyzer::checkConfig()  {
    auto mkErr = [&](const std::string& reader, const std::string corr) {
        throw std::invalid_argument(std::string("DefaultSearchRegionAnalyzer::checkConfig() -> Must load ") + reader + std::string(" if you want ") + corr);};
    if(!reader_event) mkErr("event","anything");

    if(isRealData()) return;

    if(isCorrOn(CORR_XSEC) && !reader_event) mkErr("event","CORR_XSEC");
    if(isCorrOn(CORR_TRIG) && !reader_jet_chs) mkErr("ak4Jet","CORR_TRIG");
    if(isCorrOn(CORR_TRIG) && !reader_electron) mkErr("electron","CORR_TRIG");
    if(isCorrOn(CORR_TRIG) && !reader_muon) mkErr("muon","CORR_TRIG");
    if(isCorrOn(CORR_TRIG) && !reader_genpart) mkErr("genParticle","CORR_TRIG");
    if(isCorrOn(CORR_PU)   && !reader_event) mkErr("event","CORR_PU");
    if(isCorrOn(CORR_LEP) && !reader_jet_chs) mkErr("ak4Jet","CORR_LEP");
    if(isCorrOn(CORR_LEP) && !reader_electron) mkErr("electron","CORR_LEP");
    if(isCorrOn(CORR_LEP) && !reader_muon) mkErr("muon","CORR_LEP");
    if(isCorrOn(CORR_LEP) && !reader_genpart) mkErr("genParticle","CORR_LEP");
    if(isCorrOn(CORR_AK4BTAG) && !reader_jet) mkErr("ak4PuppiNoLepJet","CORR_AK4BTAG");
    if(isCorrOn(CORR_SJBTAG) && !reader_fatjet) mkErr("ak8PuppiNoLepJet","CORR_SJBTAG");
}
//--------------------------------------------------------------------------------------------------
bool DefaultSearchRegionAnalyzer::runEvent() {
    mcProc = FillerConstants::MCProcess(reader_event->process);
    if(isRealData()) smpName = "data";
    else if (mcProc == FillerConstants::SIGNAL) smpName = TString::Format("m%i",signal_mass);
    else smpName = FillerConstants::MCProcessNames[mcProc];


    //|||||||||||||||||||||||||||||| GEN PARTICLES ||||||||||||||||||||||||||||||
    if(reader_genpart){
        if(reader_event->process == FillerConstants::SIGNAL) diHiggsEvt.setDecayInfo(reader_genpart->genParticles);
        smDecayEvt.setDecayInfo(reader_genpart->genParticles);
    }

    //|||||||||||||||||||||||||||||| CHS JETS ||||||||||||||||||||||||||||||
    if(reader_jet_chs){
        jets_chs = PhysicsUtilities::selObjsMom(reader_jet_chs->jets,30);
        ht_chs = JetKinematics::ht(jets_chs);
    }


    //|||||||||||||||||||||||||||||| LEPTONS ||||||||||||||||||||||||||||||
    if(reader_electron && reader_muon){
        selectedLeptons = leptonProc->getLeptons(*reader_event,*reader_muon,*reader_electron);
        selectedLepton = selectedLeptons.size() ? selectedLeptons.front() : 0;
    }


    //|||||||||||||||||||||||||||||| FILTERS ||||||||||||||||||||||||||||||
    passEventFilters= EventSelection::passEventFilters(*reader_event);
    passTriggerPreselection= EventSelection::passTriggerPreselection(*reader_event,ht_chs,selectedLeptons);


    //|||||||||||||||||||||||||||||| FATJETS ||||||||||||||||||||||||||||||
    if(reader_fatjet && selectedLepton){
        fatjetCands = fjProc->loadFatJets(*reader_fatjet,selectedLepton);
        hbbCand     = fjProc->getHBBCand();
        wjjCand     = fjProc->getWjjCand();
        hbbCSVCat   = fjProc->getHbbCSVCat();
        wjjCSVCat   = fjProc->getWjjCSVCat();
        hbbNSJs     = fjProc->getHbbNSJs();
        wjjNSJs     = fjProc->getWjjNSJs();
        passHbbSel  = fjProc->passHbbSel();
        passWjjSel  = fjProc->passWjjSel();
    } else {
        fatjetCands.clear();
        wjjCand    =  0;
        hbbCand    =  0;
        hbbCSVCat   = BTagging::CSVSJ_INCL;
        wjjCSVCat   = BTagging::CSVSJ_INCL;
        hbbNSJs     =  0;
        wjjNSJs     =  0;
        passWjjSel =  false;
        passHbbSel =  false;
    }

    if(wjjCand){
        neutrino    = HiggsSolver::getInvisible(reader_event->met,(selectedLepton->p4() + wjjCand->p4()) );
        wlnu        =  neutrino.p4() + selectedLepton->p4();
        hWW         = wlnu.p4() + wjjCand->p4();
        wlnuDR      = PhysicsUtilities::deltaR(neutrino,*selectedLepton);
        wwDM        = PhysicsUtilities::deltaR( wlnu,*wjjCand) * hWW.pt()/125.0;
        passWlnuDR  = wlnuDR < 3.2;
        passWWDM    = wwDM < 2.0;
    } else {
        neutrino    =  MomentumF();
        wlnu        =  MomentumF();
        hWW         = MomentumF();
        wlnuDR      = 0;
        wwDM        = 0;
        passWlnuDR  = false;
        passWWDM    = false;
    }

    if(hbbCand){
        hbbMass    =   isCorrOn(CORR_SDMASS) ? hbbFJSFProc->getCorrSDMass(hbbCand) : hbbCand->sdMom().mass();
        hh =  wjjCand  ? (hWW.p4() + hbbCand->p4()) :  MomentumF();
    } else {
        hbbMass    = 0;
        hh         =  MomentumF();
    }

    //|||||||||||||||||||||||||||||| PUPPI JETS ||||||||||||||||||||||||||||||

    if(reader_jet){
        jets = PhysicsUtilities::selObjsMom(reader_jet->jets,30,2.4,[](const Jet* j){return j->passTightID();} );
        nMedBTags = PhysicsUtilities::selObjsD(jets,[](const Jet* j){return BTagging::isMediumCSVTagged(*j);} ).size();
        if(hbbCand){
            jets_HbbV = PhysicsUtilities::selObjsD(jets,[&](const Jet* j){return  PhysicsUtilities::deltaR2(*j,*hbbCand ) >= 1.2*1.2 ;} );
            nMedBTags_HbbV = PhysicsUtilities::selObjsD(jets_HbbV,[](const Jet* j){return BTagging::isMediumCSVTagged(*j);} ).size();
        } else {
            jets_HbbV = jets;
            nMedBTags_HbbV = nMedBTags;
        }
    }


    //|||||||||||||||||||||||||||||| EVENT WEIGHTS ||||||||||||||||||||||||||||||
    weight = 1;
    if(!isRealData()){
        if(isCorrOn(CORR_XSEC))
            weight *= EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi());
        if(isCorrOn(CORR_TRIG) && (smDecayEvt.promptElectrons.size() + smDecayEvt.promptMuons.size())   )
            weight *= trigSFProc->getLeptonTriggerSF(ht_chs, (selectedLepton && selectedLepton->isMuon()));
        if(isCorrOn(CORR_PU) )
            weight *= puSFProc->getCorrection(reader_event->nTruePUInts,CorrHelp::NOMINAL);
        if(isCorrOn(CORR_LEP)){
            leptonSFProc->load(smDecayEvt,selectedLeptons,&jets_chs);
            weight *= leptonSFProc->getSF();
        }
        if(isCorrOn(CORR_SJBTAG)){
            weight *= sjbtagSFProc->getSF({hbbCand});
        }
        if(isCorrOn(CORR_AK4BTAG)){
            weight *= ak4btagSFProc->getSF(jets_HbbV);
        }

    }
    return true;
}
}
