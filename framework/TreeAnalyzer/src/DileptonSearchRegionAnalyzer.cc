//#include "../interface/DileptonSearchRegionAnalyzer.h"
//
//#include "TreeReaders/interface/EventReader.h"
//#include "TreeReaders/interface/GenParticleReader.h"
//#include "TreeReaders/interface/ElectronReader.h"
//#include "TreeReaders/interface/MuonReader.h"
//#include "TreeReaders/interface/JetReader.h"
//#include "TreeReaders/interface/FatJetReader.h"
//
//#include "Processors/Corrections/interface/EventWeights.h"
//#include "Processors/Variables/interface/JetKinematics.h"
//
//#include "Processors/Variables/interface/DileptonSelection.h"
//#include "Processors/Variables/interface/FatJetSelection.h"
//
//#include "Processors/EventSelection/interface/EventSelection.h"
//#include "Processors/Corrections/interface/TriggerScaleFactors.h"
//#include "Processors/Corrections/interface/LeptonScaleFactors.h"
//#include "Processors/Corrections/interface/BTagScaleFactors.h"
//#include "Processors/Corrections/interface/FatJetScaleFactors.h"
//#include "Processors/Corrections/interface/JetAndMETCorrections.h"
//#include "TPRegexp.h"
//
//
//
//
//namespace TAna {
////--------------------------------------------------------------------------------------------------
//DileptonSearchRegionAnalyzer::DileptonSearchRegionAnalyzer(std::string fileName, std::string treeName, int treeInt, size randomSeed) : BaseTreeAnalyzer(fileName,treeName,treeInt, randomSeed){
//    TPRegexp r1(".*m(\\d+)_[0-9]*\\..*$");
//    auto match = r1.MatchS(fileName);
//    const Int_t nrSubStr = match->GetLast()+1;
//    if(nrSubStr>1){
//        signal_mass = (((TObjString *)match->At(1))->GetString()).Atoi();
//    }
//    fjProc      .reset(new FatJetProcessor ()); DefaultFatJetSelections::setDefaultFatJetProcessor(*fjProc);
//    dileptonProc.reset(new DileptonProcessor()); DefaultDileptonSelections::setDefaultDileptonProcessor(*dileptonProc);
//    trigSFProc  .reset(new TriggerScaleFactors (dataDirectory));
//    puSFProc    .reset(new PUScaleFactors (dataDirectory));
//    leptonSFProc.reset(new ActParamScaleFactors(dataDirectory));
//    ak4btagSFProc.reset(new JetBTagScaleFactors (dataDirectory));
//    sjbtagSFProc.reset(new SubJetBTagScaleFactors (dataDirectory));
//    hbbFJSFProc .reset(new HbbFatJetScaleFactors (dataDirectory));
//    topPTProc   .reset(new TopPTWeighting (dataDirectory));
//
//    JERAK4PuppiProc .reset(new JERCorrector (dataDirectory, "corrections/Summer16_25nsV1_MC_PtResolution_AK4PFPuppi.txt",randGen));;
//    JERAK4CHSProc   .reset(new JERCorrector (dataDirectory, "corrections/Summer16_25nsV1_MC_PtResolution_AK4PFCHS.txt",randGen));;
//    JERAK8PuppiProc .reset(new JERCorrector (dataDirectory, "corrections/Summer16_25nsV1_MC_PtResolution_AK8PFPuppi.txt",randGen));;
//    JESUncProc . reset(new JESUncShifter());
//    METUncProc . reset(new METUncShifter());
//
//    setLumi(35.922); //https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/688.html
//
//    turnOnCorr(CORR_XSEC);
//    turnOnCorr(CORR_TRIG);
//    turnOnCorr(CORR_PU  );
//    turnOnCorr(CORR_LEP );
//    turnOnCorr(CORR_SJBTAG);
//    turnOnCorr(CORR_AK4BTAG);
//    turnOnCorr(CORR_SDMASS);
////    turnOnCorr(CORR_TOPPT);
////    turnOnCorr(CORR_JER);
//}
////--------------------------------------------------------------------------------------------------
//DileptonSearchRegionAnalyzer::~DileptonSearchRegionAnalyzer(){}
////--------------------------------------------------------------------------------------------------
//void DileptonSearchRegionAnalyzer::resetCorr() {corrections = 0;}
//bool DileptonSearchRegionAnalyzer::isCorrOn(Corrections corr) const {return FillerConstants::doesPass(corrections,corr);}
//void DileptonSearchRegionAnalyzer::turnOnCorr(Corrections corr) {FillerConstants::addPass(corrections,corr);}
//void DileptonSearchRegionAnalyzer::turnOffCorr(Corrections corr) {FillerConstants::removePass(corrections,corr);}
//
////--------------------------------------------------------------------------------------------------
//void DileptonSearchRegionAnalyzer::loadVariables()  {
//    reader_event       =loadReader<EventReader>   ("event",isRealData());
//    reader_fatjet      =loadReader<FatJetReader>  ("ak8PuppiJet",isRealData());
//    reader_fatjet_noLep=loadReader<FatJetReader>  ("ak8PuppiNoLepJet",isRealData(),false);
//    reader_jet_chs     =loadReader<JetReader>     ("ak4Jet",isRealData());
//    reader_jet         =loadReader<JetReader>     ("ak4PuppiJet",isRealData(),false);
//    reader_electron    =loadReader<ElectronReader>("electron");
//    reader_muon        =loadReader<MuonReader>    ("muon");
//
//    if(!isRealData()){
//        reader_genpart =loadReader<GenParticleReader>   ("genParticle");
//    }
//
//    checkConfig();
//}
////--------------------------------------------------------------------------------------------------
//void DileptonSearchRegionAnalyzer::checkConfig()  {
//    auto mkErr = [&](const std::string& reader, const std::string corr) {
//        throw std::invalid_argument(std::string("DefaultSearchRegionAnalyzer::checkConfig() -> Must load ") + reader + std::string(" if you want ") + corr);};
//    if(!reader_event) mkErr("event","anything");
//
//    if(isRealData()) return;
//
//    if(isCorrOn(CORR_XSEC) && !reader_event) mkErr("event","CORR_XSEC");
//    if(isCorrOn(CORR_TRIG) && !reader_jet_chs) mkErr("ak4Jet","CORR_TRIG");
//    if(isCorrOn(CORR_TRIG) && !reader_electron) mkErr("electron","CORR_TRIG");
//    if(isCorrOn(CORR_TRIG) && !reader_muon) mkErr("muon","CORR_TRIG");
//    if(isCorrOn(CORR_TRIG) && !reader_genpart) mkErr("genParticle","CORR_TRIG");
//    if(isCorrOn(CORR_PU)   && !reader_event) mkErr("event","CORR_PU");
//    if(isCorrOn(CORR_LEP) && !reader_electron) mkErr("electron","CORR_LEP");
//    if(isCorrOn(CORR_LEP) && !reader_muon) mkErr("muon","CORR_LEP");
//    if(isCorrOn(CORR_LEP) && !reader_genpart) mkErr("genParticle","CORR_LEP");
//    if(isCorrOn(CORR_AK4BTAG) && !reader_jet) mkErr("ak4PuppiNoLepJet","CORR_AK4BTAG");
//    if(isCorrOn(CORR_SJBTAG) && !reader_fatjet) mkErr("ak8PuppiNoLepJet","CORR_SJBTAG");
//    if(isCorrOn(CORR_TOPPT) && !reader_event) mkErr("event","CORR_TOPPT");
//    if(isCorrOn(CORR_TOPPT) && !reader_genpart) mkErr("genParticle","CORR_TOPPT");
//    if(isCorrOn(CORR_JER) && !reader_fatjet) mkErr("fatjet","CORR_JER");
//    if(isCorrOn(CORR_JER) && !reader_fatjet_noLep) mkErr("fatjet_noLep","CORR_JER");
//    if(isCorrOn(CORR_JER) && !reader_jet) mkErr("jet","CORR_JER");
//    if(isCorrOn(CORR_JER) && !reader_jet_chs) mkErr("jet_chs","CORR_JER");
//
//    if(isCorrOn(CORR_JES) && !reader_fatjet) mkErr("fatjet","CORR_JES");
//    if(isCorrOn(CORR_JES) && !reader_fatjet_noLep) mkErr("fatjet_noLep","CORR_JES");
//    if(isCorrOn(CORR_JES) && !reader_jet) mkErr("jet","CORR_JES");
//    if(isCorrOn(CORR_JES) && !reader_jet_chs) mkErr("jet_chs","CORR_JES");
//}
////--------------------------------------------------------------------------------------------------
//bool DileptonSearchRegionAnalyzer::runEvent() {
//    mcProc = FillerConstants::MCProcess(*reader_event->process);
//    if(isRealData()) smpName = "data";
//    else if (mcProc == FillerConstants::SIGNAL) smpName = TString::Format("m%i",signal_mass);
//    else smpName = FillerConstants::MCProcessNames[mcProc];
//
//    //|||||||||||||||||||||||||||||| CORRECT JETS AND MET FIRST ||||||||||||||||||||||||||||||
//    if(!isRealData()){
//        if(isCorrOn(CORR_JES) ){
//            Met dummyMET =reader_event->met;
//            JESUncProc ->processJets(*reader_jet,reader_event->met);
//            JESUncProc ->processJets(*reader_jet_chs,dummyMET);
//            JESUncProc ->processFatJets(reader_fatjet_noLep->jets);
//            JESUncProc ->processFatJets(reader_fatjet->jets);
//        }
//
//        if(isCorrOn(CORR_JER) ){
//            Met dummyMET =reader_event->met;
//            JERAK4PuppiProc ->processJets(*reader_jet,reader_event->met,reader_jet_chs->genJets,*reader_event->rho);
//            JERAK4CHSProc   ->processJets(*reader_jet_chs,dummyMET,reader_jet_chs->genJets,*reader_event->rho);
//            JERAK8PuppiProc ->processFatJets(reader_fatjet_noLep->jets,std::vector<GenJet>(),*reader_event->rho);
//            JERAK8PuppiProc ->processFatJets(reader_fatjet->jets,reader_fatjet->genJets,*reader_event->rho);
//        }
//        if(isCorrOn(CORR_MET) ){
//            METUncProc->process(reader_event->met,*reader_event);
//        }
//    }
//
//    //|||||||||||||||||||||||||||||| GEN PARTICLES ||||||||||||||||||||||||||||||
//    if(reader_genpart){
//        if(*reader_event->process == FillerConstants::SIGNAL) diHiggsEvt.setDecayInfo(reader_genpart->genParticles);
//        smDecayEvt.setDecayInfo(reader_genpart->genParticles);
//    }
//
//    //|||||||||||||||||||||||||||||| CHS JETS ||||||||||||||||||||||||||||||
//    if(reader_jet_chs){
//        jets_chs = PhysicsUtilities::selObjsMom(reader_jet_chs->jets,30);
//        ht_chs = JetKinematics::ht(jets_chs);
//    }
//
//    //|||||||||||||||||||||||||||||| LEPTONS ||||||||||||||||||||||||||||||
//    if(reader_electron && reader_muon){
//        selectedDileptons = dileptonProc->getLeptons(*reader_event,*reader_muon,*reader_electron);
//    }
//
//    //|||||||||||||||||||||||||||||| FILTERS ||||||||||||||||||||||||||||||
//    passEventFilters= EventSelection::passEventFilters(*reader_event);
//    passTriggerPreselection= EventSelection::passTriggerPreselection(*reader_event,ht_puppi,selectedDileptons); // changed from ht_chs to ht_puppi due to dilepton skims
//
//    //|||||||||||||||||||||||||||||| FATJETS ||||||||||||||||||||||||||||||
//    if(reader_fatjet && selectedDileptons.size() == 2){
//
//        fjProc->loadDilepFatJet(*reader_fatjet,selectedDileptons[0],selectedDileptons[1]);
//        hbbCand     = fjProc->getDilepHbbCand();
//        hbbCSVCat   = fjProc->getDilepHbbCSVCat();
//
//		double pz = reader_event->met.pt() / TMath::Tan((selectedDileptons[0]->p4()+selectedDileptons[1]->p4()).theta());
//		pz = ((pz < 0) == ((selectedDileptons[0]->p4()+selectedDileptons[1]->p4()).pz() < 0)) ? pz : (-1)*pz;
//		ASTypes::CartLorentzVector pnunu(reader_event->met.px(),reader_event->met.py(),pz,sqrt(pow(reader_event->met.px(),2)+pow(reader_event->met.py(),2)+pz*pz+40*40));
////        HwwTestStat = hwwSolver.HwwMinimization(selectedDileptons[0]->p4(),selectedDileptons[1]->p4(),reader_event->met.p4(),&hwwInfo);
//
//        hww = selectedDileptons[0]->p4() + selectedDileptons[1]->p4() + pnunu;
//
//        if (selectedDileptons[0]->isMuon() && selectedDileptons[1]->isMuon()) dilepChan = mumu;
//        else if (selectedDileptons[0]->isElectron() && selectedDileptons[1]->isElectron()) dilepChan = ee;
//        else dilepChan = emu;
//    } else {
//        hbbCand    =  0;
//        hbbCSVCat  = BTagging::CSVSJ_INCL;
//        hww     = MomentumF();
//    }
//
//    if(hbbCand){
//        hbbMass =   isCorrOn(CORR_SDMASS) ? hbbFJSFProc->getCorrSDMass(hbbCand) : hbbCand->sdMom().mass();
//        hh = hww.p4() + hbbCand->p4();
//    } else {
//        hbbMass    = 0;
//        hh      = MomentumF();
//    }
//
//    //|||||||||||||||||||||||||||||| PUPPI JETS ||||||||||||||||||||||||||||||
//
//    if(reader_jet){
//    	jets_ht = PhysicsUtilities::selObjsMom(reader_jet->jets,30);
//        ht_puppi = JetKinematics::ht(jets_ht);
//
//        if (selectedDileptons.size() == 2) {
//        	jets_puppiNoLepOverlap = PhysicsUtilities::selObjsMom(reader_jet->jets,30,2.4,[&](const Jet* j){return (PhysicsUtilities::deltaR2(*j,*selectedDileptons[0]) > 0.4*0.4
//        			&& PhysicsUtilities::deltaR2(*j,*selectedDileptons[1]) > 0.4*0.4);} );
//        	ht_puppiNoLepOverlap = JetKinematics::ht(jets_puppiNoLepOverlap);
//        }
//
//        jets = PhysicsUtilities::selObjsMom(reader_jet->jets,30,2.4,[](const Jet* j){return j->passTightID();} );
//        nMedBTags = PhysicsUtilities::selObjsD(jets,[](const Jet* j){return BTagging::isMediumCSVTagged(*j);} ).size();
//        if(hbbCand){
//            jets_HbbV = PhysicsUtilities::selObjsD(jets,[&](const Jet* j){return  PhysicsUtilities::deltaR2(*j,*hbbCand ) >= 1.2*1.2 ;} );
//            nMedBTags_HbbV = PhysicsUtilities::selObjsD(jets_HbbV,[](const Jet* j){return BTagging::isMediumCSVTagged(*j);} ).size();
//        } else {
//            jets_HbbV = jets;
//            nMedBTags_HbbV = nMedBTags;
//        }
//    }
//
//    //|||||||||||||||||||||||||||||| EVENT WEIGHTS ||||||||||||||||||||||||||||||
//    weight = 1;
//    if(!isRealData()){
//        if(isCorrOn(CORR_XSEC))
//            weight *= EventWeights::getNormalizedEventWeight(*reader_event,xsec(),nSampEvt(),lumi());
//        if(isCorrOn(CORR_TRIG) && (selectedDileptons.size()==2) && (smDecayEvt.promptElectrons.size() + smDecayEvt.promptMuons.size())   ){
//            weight *= trigSFProc->getLeptonTriggerSF(ht_puppi, selectedDileptons[0]->isMuon());
//            }
//        if(isCorrOn(CORR_PU) )
//            weight *= puSFProc->getCorrection(*reader_event->nTruePUInts,CorrHelp::NOMINAL);
//        if(isCorrOn(CORR_LEP)){
//            leptonSFProc->load(smDecayEvt,selectedDileptons);
//            weight *= leptonSFProc->getSF();
//        }
//        if(isCorrOn(CORR_SJBTAG)){
//            weight *= sjbtagSFProc->getSF({hbbCand});
//        }
//        if(isCorrOn(CORR_AK4BTAG)){
//            weight *= ak4btagSFProc->getSF(jets_HbbV);
//        }
//        if(isCorrOn(CORR_TOPPT)){
//            weight *= topPTProc->getCorrection(mcProc,smDecayEvt);
//        }
//    }
//
//    return true;
//}
//}
