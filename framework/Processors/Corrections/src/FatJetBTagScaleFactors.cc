#include "../interface/FatJetBTagScaleFactors.h"
#include "DataFormats/interface/Jet.h"
#include "TreeReaders/interface/JetReader.h"
#include "DataFormats/interface/FatJet.h"
//#include "Processors/Corrections/interface/BTagCalibrationStandalone.h"
#include "Configuration/interface/ReaderConstants.h"


namespace TAna {
using namespace CorrHelp;
using namespace BTagging;

void FatJetBTagScaleFactors::setParameters(const FatJetParameters& params, int year) {
	sfFile = dataDir+params.fatJetBtagSFFile;
	read_csv(YR(year));
}

void FatJetBTagScaleFactors::fillSfLine(int i, sfLine& line, std::string& colVal) {
	switch (i) {
	case 0:
		line.year = yrs[colVal]; break;
	case 1:
		line.wp = colVal=="mp"?0:(colVal=="hp"?1:-1); break;
	case 2:
		line.ptmin = std::stof(colVal); break;
	case 3:
		line.ptmax = colVal=="Inf" ? float(99999) : std::stof(colVal); break;
	case 4:
		line.sf = std::stof(colVal); break;
	case 5:
		line.sfUp = line.sf + std::stof(colVal); break;
	case 6:
		line.sfDn = line.sf - std::stof(colVal); break;
	default:
		throw std::invalid_argument("shouldn't be reaching here");
	}
}

void FatJetBTagScaleFactors::read_csv(YR year) {

    std::ifstream myFile(sfFile);
    if(!myFile.is_open()) throw std::runtime_error("Could not open FatJet SF file");

    // Helper vars
    std::string line, colname;

    if(myFile.good()) {
        std::getline(myFile, line);
        std::stringstream ss(line);
    }

    // Read data, line by line
    while(std::getline(myFile, line))
    {
        std::stringstream ss(line);
        sfLine sfl;

        std::vector<TString> wp, ptmin, ptmax, sf, sfUp, sfDn;
        std::string colVal;

        int colIdx = 0;
        while (ss >> colVal) {
            fillSfLine(colIdx,sfl,colVal);
            colIdx++;
        }

        if(sfl.year != year) continue;
        sflines.push_back(sfl);
    }
}

float FatJetBTagScaleFactors::getSF(const FatJetParameters& parameters,
		const std::vector<const FatJet*>& fatJets, CorrHelp::CORRTYPE corrT) const {
	float SF = 1.0;

//	auto getsf = [&](sfLine& sfl) {
//		std::vector<float> sfs = {1.0,sfl.sfDn,sfl.sf,sfl.sfUp};
//		return sfs[corr];
//	};

	for(const auto& fj : fatJets) {
		if(!fj) continue;

		int wp = -1;
		if((fj->*parameters.getFatJetTagVal)() >= parameters.DeepAK8_TWP) wp = 1;
		else if((fj->*parameters.getFatJetTagVal)() >= parameters.DeepAK8_LWP) wp = 0;

		for(const auto& sfl : sflines) {
			if(wp != sfl.wp) continue;
			if(fj->pt() < sfl.ptmin) continue;
			if(fj->pt() > sfl.ptmax) continue;
//			SF *= getsf(sfl);
			if(corrT == CorrHelp::NOMINAL) SF *= sfl.sf;
			else if(corrT == CorrHelp::UP) SF *= sfl.sfUp;
			else if(corrT == CorrHelp::DOWN) SF *= sfl.sfDn;

//			printf("fj: pt = %.2f, dak8 = %.2f (wp=%d)\nsf = %.2f\n\n",
//					fj->pt(),fj->deep_MDZHbb(),wp,SF);
			break;
		}
	}

	return SF;
}

}

