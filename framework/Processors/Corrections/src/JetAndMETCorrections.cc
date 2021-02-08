
#include "Processors/Corrections/interface/JetAndMETCorrections.h"
#include "DataFormats/interface/Jet.h"
#include "TreeReaders/interface/JetReader.h"
#include "TreeReaders/interface/EventReader.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"



namespace TAna {
using namespace CorrHelp;
//--------------------------------------------------------------------------------------------------
JERCorrector::JERCorrector (const std::string& dataDir,std::shared_ptr<TRandom3> rndGen, const CorrHelp::CORRTYPE cT )
: dataDir(dataDir), cT(cT),rndGen(rndGen){

}
//--------------------------------------------------------------------------------------------------
void JERCorrector::setParameters(const JetParameters& param) {

    ak8Puppi_resObj  .reset(new JMEStand::JetResolutionObject(dataDir+param.jer_AK8Puppi_resFile));
    ak8Puppi_sfObj   .reset(new JMEStand::JetResolutionObject(dataDir+param.jer_AK8Puppi_sfFile ));
    ak4CHS_resObj    .reset(new JMEStand::JetResolutionObject(dataDir+param.jer_AK4CHS_resFile  ));
    ak4CHS_sfObj     .reset(new JMEStand::JetResolutionObject(dataDir+param.jer_AK4CHS_sfFile   ));
}
//--------------------------------------------------------------------------------------------------
void JERCorrector::processJets(JetReader& jetreader,Met& met,const GenJetCollection& genjets, const float rho){
    if(cT == CORRTYPE::NONE) return;
    JMEStand::JetParameters parameters;
    parameters.setRho(std::min(40.0f,rho));
    std::vector<const GenJet*> gjptrs; gjptrs.reserve(genjets.size());
//    std::cout << "START! "<< rho <<std::endl;
    for(auto& gj: genjets) gjptrs.push_back(&gj);
    double deltaMX= 0;
    double deltaMY= 0;

    for(auto& j : jetreader.jets){
//        std::cout << j.pt() <<" "<< j.eta() <<"->";
        if(j.pt() > 7000) continue;
        if(j.pt() < 15) {
//            std::cout <<std::endl;
            continue;
        }
        parameters.setJetPt(j.pt());
        float eta = j.eta() < -4.7 ? -4.7 : (j.eta() > 4.7 ? 4.7 : j.eta());
        parameters.setJetEta(eta);
        auto oRawFact = j.toRawFactor();

        const float jres = ak4CHS_resObj->evaluateFormula(*ak4CHS_resObj
                ->getRecord(parameters),parameters);

        const float jSF =  ak4CHS_sfObj->getRecord(parameters)
                 ->getParametersValues()[getSFCount(cT)];

        auto ptSF = correctJet(&j,gjptrs,jres,jSF,0.4);

        deltaMX += ((jetreader.metUnc_rawPx)[j.index()]/oRawFact)*(1.0 - ptSF);
        deltaMY += ((jetreader.metUnc_rawPy)[j.index()]/oRawFact)*(1.0 - ptSF);
//        std::cout << std::sqrt((*jetreader.metUnc_rawPx)[j.index()]*(*jetreader.metUnc_rawPx)[j.index()] + (*jetreader.metUnc_rawPy)[j.index()]*(*jetreader.metUnc_rawPy)[j.index()] )/oRawFact <<","<<
//                (std::sqrt((*jetreader.metUnc_rawPx)[j.index()]*(*jetreader.metUnc_rawPx)[j.index()] + (*jetreader.metUnc_rawPy)[j.index()]*(*jetreader.metUnc_rawPy)[j.index()] )/oRawFact)*(1.0-ptSF) <<std::endl;
//        std::cout << j.pt() <<" "<< j.eta() <<std::endl;
    }

    ASTypes::CylLorentzVectorF deltaV(std::sqrt(deltaMX*deltaMX+deltaMY*deltaMY),0.0f, (deltaMX == 0.0 && deltaMY == 0.0) ? 0 : std::atan2(deltaMY, deltaMX),0 );
    met.p4() += deltaV;
//    std::cout << "END! "<< rho <<std::endl;
    std::sort(jetreader.jets.begin(),jetreader.jets.end(),PhysicsUtilities::greaterPT<Jet>());
}
//--------------------------------------------------------------------------------------------------
void JERCorrector::processFatJets(FatJetCollection& jets,const GenFatJetCollection& genjets, const float rho) {
    if(cT == CORRTYPE::NONE) return;
    JMEStand::JetParameters parameters;
    parameters.setRho(std::min(40.0f,rho));
    std::vector<const GenJet*> gjptrs; gjptrs.reserve(genjets.size());
    for(auto& gj: genjets) gjptrs.push_back(&gj);
    for(auto& j : jets){

        if(j.pt() < 15) continue;
        if(j.pt() > 7000) continue;
        parameters.setJetPt(j.pt());
        float eta = j.eta() < -4.7 ? -4.7 : (j.eta() > 4.7 ? 4.7 : j.eta());
        parameters.setJetEta(eta);

        const float jres = ak8Puppi_resObj->evaluateFormula(*ak8Puppi_resObj
                ->getRecord(parameters),parameters);
        const float jSF =  ak8Puppi_sfObj->getRecord(parameters)
                ->getParametersValues()[getSFCount(cT)];

        correctJet(&j,gjptrs,jres,jSF,0.8);
    }
    std::sort(jets.begin(),jets.end(),PhysicsUtilities::greaterPT<FatJet>());
}

//--------------------------------------------------------------------------------------------------
int JERCorrector::getSFCount(const CORRTYPE c) {
    switch(c){
    case NOMINAL:
        return 0;
    case DOWN:
        return 1;
    case UP:
        return 2;
    default:
        throw std::invalid_argument(std::string("JERCorrector cannot be used with NONE"));
    }
};
//--------------------------------------------------------------------------------------------------
float JERCorrector::correctJet(Jet* jet, const std::vector<const GenJet*> genjets, const float jres,
        const float resSF, const float coneR){

    double nearDR = 0;
    int idx = PhysicsUtilities::findNearestDRDeref(*jet,genjets,nearDR);
//    std::cout <<"(";
//    if(idx >= 0) std::cout <<genjets[idx]->pt()<<",";
//    else std::cout << "-1,";
//    std::cout << jres <<","<<resSF<<",";
    float ptSF = 1.0;
    // Reference: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    if(idx>= 0 && nearDR < coneR/2.0 && std::fabs(jet->pt()-genjets[idx]->pt()) <3.0*jres*jet->pt() ){
        ptSF = 1 + (resSF-1.0)*(jet->pt()-genjets[idx]->pt())/jet->pt();
//        std::cout <<(jet->pt()-genjets[idx]->pt())/jet->pt() <<","<<ptSF<<")";
    } else{
        ptSF = 1 + rndGen->Gaus(0.,jres)*std::sqrt(std::max(0.0f,resSF*resSF-1));
//        std::cout <<ptSF<<")";
    }
        if(ptSF>0){
            jet->setP4(ptSF*jet->pt(),jet->eta(),jet->phi(),ptSF*jet->mass());
            jet->setRawFactor(jet->toRawFactor()/ptSF);
        } else{
            jet->setP4(0.0f,jet->eta(),jet->phi(),0.0f);
        }

    return std::max(0.0f,ptSF);
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
JESUncShifter::JESUncShifter (const CorrHelp::CORRTYPE cT)
: cT(cT){}
//--------------------------------------------------------------------------------------------------
void JESUncShifter::processJets(JetReader& jetreader,Met& met){
    if(cT == CORRTYPE::NONE || cT == CORRTYPE::NOMINAL) return;
    double deltaMX= 0;
    double deltaMY= 0;
    for(auto& j : jetreader.jets){
        auto oRawFact = j.toRawFactor();
        auto ptSF = correctJet(&j);
        if(j.pt() >= 15) {
            deltaMX += ((jetreader.metUnc_rawPx)[j.index()]/oRawFact)*(1.0 - ptSF);
            deltaMY += ((jetreader.metUnc_rawPy)[j.index()]/oRawFact)*(1.0 - ptSF);
        }
    }
    ASTypes::CylLorentzVectorF deltaV(std::sqrt(deltaMX*deltaMX+deltaMY*deltaMY),0.0f, (deltaMX == 0.0 && deltaMY == 0.0) ? 0 : std::atan2(deltaMY, deltaMX),0 );
    met.p4() += deltaV;
    std::sort(jetreader.jets.begin(),jetreader.jets.end(),PhysicsUtilities::greaterPT<Jet>());
}
//--------------------------------------------------------------------------------------------------
void JESUncShifter::processFatJets(FatJetCollection& jets){
    if(cT == CORRTYPE::NONE || cT == CORRTYPE::NOMINAL) return;
    for(auto& j : jets){
        correctJet(&j);
        for(unsigned int iSJ= 0; iSJ < j.nSubJets(); ++iSJ)
            correctJet(&j.subJet(iSJ));

    }
    std::sort(jets.begin(),jets.end(),PhysicsUtilities::greaterPT<FatJet>());
}
//--------------------------------------------------------------------------------------------------
float JESUncShifter::correctJet(BaseRecoJet* jet) const {
    float ptSF = 1.0;
    if(cT == CORRTYPE::DOWN) ptSF = 1.0 - jet->jecUnc();
    else if(cT == CORRTYPE::UP) ptSF = 1.0 + jet->jecUnc();
    if(ptSF>0){
        jet->setP4(ptSF*jet->pt(),jet->eta(),jet->phi(),ptSF*jet->mass());
        jet->setRawFactor(jet->toRawFactor()/ptSF);
    } else{
        jet->setP4(0.0f,jet->eta(),jet->phi(),0.0f);
    }
    return std::max(0.0f,ptSF);
}
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
METUncShifter::METUncShifter (const CorrHelp::CORRTYPE cT)
: cT(cT){}
//--------------------------------------------------------------------------------------------------
void METUncShifter::process(Met& met, const EventReader& eventreader) const {
    if(cT == CORRTYPE::NONE || cT == CORRTYPE::NOMINAL) return;
    //calculate deltaMet so it can be applied pos/neg and to possibly already corrected met
    const Met met_uncUp(ASTypes::CylLorentzVectorF(eventreader.met_unclUp_pt.val(),0,eventreader.met_unclUp_phi.val(),0));
    const Met met_std(ASTypes::CylLorentzVectorF(eventreader.met_pt.val(),0,eventreader.met_phi.val(),0));
    auto deltaM = met_uncUp.p4()-met_std.p4();

    if(cT == CORRTYPE::UP) met.p4() += deltaM;
    else met.p4() -= deltaM;
}
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
HEM1516TestCorrector::HEM1516TestCorrector (const CorrHelp::CORRTYPE cT) : cT(cT){}
//--------------------------------------------------------------------------------------------------
void HEM1516TestCorrector::processJets(JetReader& jetreader, Met& met) {
    if(cT == CORRTYPE::NONE) return;
    float deltaPX = 0, deltaPY = 0;
    for(auto& j : jetreader.jets) {
    	if(j.pt() < 15) continue;
    	if(!j.passTightID()) continue;

    	// instructions given in email on Oct 20 2019
    	if(j.phi() > -1.57 && j.phi() < -0.87) {
    		float sf = 1.0;
    		if(j.eta() > -2.5 && j.eta() < -1.3) sf = 0.8;
    		else if (j.eta() > -3.0 && j.eta() < -2.5) sf = 0.65;
    		else continue;

    		const ASTypes::CartLorentzVector newMom(sf*j.px(),sf*j.py(),sf*j.pz(),sf*j.E());

    		deltaPX += (1-sf)*j.px();
    		deltaPY += (1-sf)*j.py();
    		j.setP4(newMom);
    	}
    }

    ASTypes::CartLorentzVector deltaMet(deltaPX,deltaPY,0.0f,std::sqrt(deltaPX*deltaPX + deltaPY*deltaPY));
    met.p4() += deltaMet;
    std::sort(jetreader.jets.begin(),jetreader.jets.end(),PhysicsUtilities::greaterPT<Jet>());
}
//--------------------------------------------------------------------------------------------------
void HEM1516TestCorrector::processFatJets(FatJetCollection& fatjets) {
    if(cT == CORRTYPE::NONE) return;

    for(auto& fj : fatjets) {
    	// * should be scaling down all fatjets before considering subjets
    	if(!fj.nSubJets()) continue;

    	float origSdEnergy = fj.sdMom().E();
    	for(unsigned int idx=0; idx < fj.nSubJets(); ++idx) {
    		if(fj.subJet(idx).pt() < 15) continue;
    		// should I do anything to ensure there are no leptons in this subjet??

    		if(fj.subJet(idx).phi() > -1.57 && fj.subJet(idx).phi() < -0.87) {
    			float sf = 1.0;
    			if(fj.subJet(idx).eta() > -1.3) continue;
    			if(fj.subJet(idx).eta() > -2.5 && fj.subJet(idx).eta() < -1.3) sf = 0.8;
    			else if (fj.subJet(idx).eta() > -3.0 && fj.subJet(idx).eta() < -2.5) sf = 0.65;

        		ASTypes::CartLorentzVector newMom(sf*fj.subJet(idx).px(),sf*fj.subJet(idx).py(),
        				sf*fj.subJet(idx).pz(),sf*fj.subJet(idx).E());
        		fj.subJet(idx).setP4(newMom);
    		}
    	}
    	float scaleFactor = fj.sdMom().E() / origSdEnergy;
    	if(scaleFactor < 1.0) {
    		ASTypes::CartLorentzVector newMom(scaleFactor*fj.px(),scaleFactor*fj.py(),scaleFactor*fj.pz(),scaleFactor*fj.E());
    		fj.setP4(newMom);
    	}
    }
    std::sort(fatjets.begin(),fatjets.end(),PhysicsUtilities::greaterPT<FatJet>());
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
METPhiModulationCorrector::METPhiModulationCorrector () {}
//--------------------------------------------------------------------------------------------------
void METPhiModulationCorrector::process(Met& met, const EventReader& eventreader, const bool isData) const {

	//get corrected MET from function
	auto metcorrs = METXYCorr_Met_MetPhi(met.pt(), met.phi(), eventreader, isData);

	const Met corrMet(ASTypes::CylLorentzVectorF(metcorrs.first,0,metcorrs.second,0));
	const Met initMet(ASTypes::CylLorentzVectorF(met.pt(),0,met.phi(),0));
	auto deltaM = corrMet.p4() - initMet.p4();

	met.p4() += deltaM;

}
//--------------------------------------------------------------------------------------------------
std::pair<double,double> METPhiModulationCorrector::METXYCorr_Met_MetPhi(double originalMet, double originalMet_phi, const EventReader& eventreader, const bool isData) const {

// taken from https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection_withUL17andUL18.h

  std::pair<double,double> TheXYCorr_Met_MetPhi(originalMet,originalMet_phi);

  FillerConstants::DataEra era = FillerConstants::DataEra(*eventreader.dataEra);
  int npv = *eventreader.npv;
  if(npv>100) npv=100;

  double METxcorr(0.),METycorr(0.);

  if(isData) {
//    	  std::cout<<"isrealData"<<std::endl;
//Current recommendation for 2016 and 2018 DATA
	  FillerConstants::DataRun datarun = FillerConstants::DataRun(*eventreader.dataRun);

	  if(datarun==FillerConstants::RUN2016B) {
	  	METxcorr = -(-0.0478335*npv -0.108032);
	  	METycorr = -(0.125148*npv +0.355672);
	  }
	  if(datarun==FillerConstants::RUN2016C) {
	  	METxcorr = -(-0.0916985*npv +0.393247);
	  	METycorr = -(0.151445*npv +0.114491);
	  }
	  if(datarun==FillerConstants::RUN2016D) {
	  	METxcorr = -(-0.0581169*npv +0.567316);
	  	METycorr = -(0.147549*npv +0.403088);
	  }
	  if(datarun==FillerConstants::RUN2016E) {
	  	METxcorr = -(-0.065622*npv +0.536856);
	  	METycorr = -(0.188532*npv +0.495346);
	  }
	  if(datarun==FillerConstants::RUN2016F) {
	  	METxcorr = -(-0.0313322*npv +0.39866);
	  	METycorr = -(0.16081*npv +0.960177);
	  }
	  if(datarun==FillerConstants::RUN2016G) {
	  	METxcorr = -(0.040803*npv -0.290384);
	  	METycorr = -(0.0961935*npv +0.666096);
	  }
	  if(datarun==FillerConstants::RUN2016H) {
	  	METxcorr = -(0.0330868*npv -0.209534);
	  	METycorr = -(0.141513*npv +0.816732);
	  }
// -----------------------------------------------------------------
	  if(datarun==FillerConstants::RUN2017B) {
	  	METxcorr = -(-0.19563*npv +1.51859);
	  	METycorr = -(0.306987*npv +-1.84713);
	  }
	  if(datarun==FillerConstants::RUN2017C) {
	  	METxcorr = -(-0.161661*npv +0.589933);
	  	METycorr = -(0.233569*npv +-0.995546);
	  }
	  if(datarun==FillerConstants::RUN2017D) {
	  	METxcorr = -(-0.180911*npv +1.23553);
	  	METycorr = -(0.240155*npv +-1.27449);
	  }
	  if(datarun==FillerConstants::RUN2017E) {
	  	METxcorr = -(-0.149494*npv +0.901305);
	  	METycorr = -(0.178212*npv +-0.535537);
	  }
	  if(datarun==FillerConstants::RUN2017F) {
	  	METxcorr = -(-0.165154*npv +1.02018);
	  	METycorr = -(0.253794*npv +0.75776);
	  }
// -------------------------------------------------------------------
	  if(datarun==FillerConstants::RUN2018A) {
	  	METxcorr = -(0.362865*npv -1.94505);
	  	METycorr = -(0.0709085*npv -0.307365);
	  }
	  if(datarun==FillerConstants::RUN2018B) {
	  	METxcorr = -(0.492083*npv -2.93552);
	  	METycorr = -(0.17874*npv -0.786844);
	  }
	  if(datarun==FillerConstants::RUN2018C) {
	  	METxcorr = -(0.521349*npv -1.44544);
	  	METycorr = -(0.118956*npv -1.96434);
	  }
	  if(datarun==FillerConstants::RUN2018D) {
	  	METxcorr = -(0.531151*npv -1.37568);
	  	METycorr = -(0.0884639*npv -1.57089);
	  }
  } else {
//    	  std::cout<<"isMC"<<std::endl;
// ---------------------------- MC ----------------------- //
      if(era==FillerConstants::ERA_2016) {
      	METxcorr = -(-0.195191*npv -0.170948);
      	METycorr = -(-0.0311891*npv +0.787627);
      }

      if(era==FillerConstants::ERA_2017) {
      	METxcorr = -(-0.182569*npv +0.276542);
      	METycorr = -(0.155652*npv +-0.417633);
      }

      if(era==FillerConstants::ERA_2018) {
      	METxcorr = -(0.296713*npv -0.141506);
      	METycorr = -(0.115685*npv +0.0128193);
      }
  }

  double CorrectedMET_x = originalMet *cos( originalMet_phi)+METxcorr;
  double CorrectedMET_y = originalMet *sin( originalMet_phi)+METycorr;

  double CorrectedMET = sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y);
  double CorrectedMETPhi;
  if(CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
  else if(CorrectedMET_x==0 && CorrectedMET_y<0 )CorrectedMETPhi = -TMath::Pi();
  else if(CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x);
  else if(CorrectedMET_x <0&& CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) + TMath::Pi();
  else if(CorrectedMET_x <0&& CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) - TMath::Pi();
  else CorrectedMETPhi =0;

  TheXYCorr_Met_MetPhi.first= CorrectedMET;
  TheXYCorr_Met_MetPhi.second= CorrectedMETPhi;

  return TheXYCorr_Met_MetPhi;

}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

}


