
#include "TreeReaders/interface/EventReader.h"
#include "AnalysisSupport/TreeInterface/interface/TreeReadingWrapper.h"


namespace TAna{

EventReader::EventReader(std::string branchName, bool isRealData) : BaseReader("EventReader",branchName),realData(isRealData) {};

void EventReader::setup(TreeReadingWrapper * wrapper){

    wrapper->setBranchAddressPre(branchName,"run"               , &run              , false);
    wrapper->setBranchAddressPre(branchName,"lumi"              , &lumi             , false);
    wrapper->setBranchAddressPre(branchName,"event"             , &event            , false);
    wrapper->setBranchAddressPre(branchName,"goodVtx"           , &goodVtx          , false);
    wrapper->setBranchAddressPre(branchName,"npv"               , &npv              , false);
    wrapper->setBranchAddressPre(branchName,"rho"               , &rho              , false);
    wrapper->setBranchAddressPre(branchName,"met_pt"            , &met_pt           , true );
    wrapper->setBranchAddressPre(branchName,"met_phi"           , &met_phi          , true );
    wrapper->setBranchAddressPre(branchName,"met_sig"           , &met_sig          , false);
    wrapper->setBranchAddressPre(branchName,"met_unclUp"        , &met_unclUp       , false);
    wrapper->setBranchAddressPre(branchName,"met_unclDown"      , &met_unclDown     , false);
    wrapper->setBranchAddressPre(branchName,"met_raw_pt"        , &met_raw_pt       , false);
    wrapper->setBranchAddressPre(branchName,"met_raw_phi"       , &met_raw_phi      , false);

    if(!realData){
        wrapper->setBranchAddressPre(branchName,"nTruePUInts"       , &nTruePUInts  , false);
        wrapper->setBranchAddressPre(branchName,"weight"            , &weight       , false);
        wrapper->setBranchAddressPre(branchName,"process"           , &process      , false);
    } else {
        wrapper->setBranchAddressPre(branchName,"dataset"           , &dataset      , false);
        wrapper->setBranchAddressPre(branchName,"dataRun"           , &dataRun      , false);
        weight = 1.0;
    }

    wrapper->setBranchAddressPre(branchName,"metFilterResult"   , &metFilterResult  , false);
    wrapper->setBranchAddressPre(branchName,"trigPass"          , &trigPass         , false);
    wrapper->setBranchAddressPre(branchName,"trigPrescale"      , &trigPrescale     , false);


}

void EventReader::processVars() {
    met.setP4(met_pt,float(0),met_phi,float(0));
}


}
