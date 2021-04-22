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
	currYr = YR(year);
	read_csv(currYr);
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

float FatJetBTagScaleFactors::getSF(const FatJetParameters& parameters, const FatJet* fj, const bool isTT, CorrHelp::CORRTYPE corrT) const {
	float SF = 1.0;

//	auto getsf = [&](sfLine& sfl) {
//		std::vector<float> sfs = {1.0,sfl.sfDn,sfl.sf,sfl.sfUp};
//		return sfs[corr];
//	};

	if(!fj) return 1.0;
	if(isTT) return getTTBAR_dak8_SF(fj,currYr,corrT);

	int wp = -1;
	if((fj->*parameters.getFatJetTagVal)() >= parameters.DeepAK8_TWP) wp = 1;
	else if((fj->*parameters.getFatJetTagVal)() >= parameters.DeepAK8_LWP) wp = 0;

	for(const auto& sfl : sflines) {
		if(wp != sfl.wp) continue;
		if(fj->pt() < sfl.ptmin) continue;
		if(fj->pt() > sfl.ptmax) continue;
//		SF *= getsf(sfl);
		if(corrT == CorrHelp::NOMINAL) SF *= sfl.sf;
		else if(corrT == CorrHelp::UP) SF *= sfl.sfUp;
		else if(corrT == CorrHelp::DOWN) SF *= sfl.sfDn;

//			printf("fj: pt = %.2f, dak8 = %.2f (wp=%d)\nsf = %.2f\n\n",
//					fj->pt(),fj->deep_MDZHbb(),wp,SF);
		break;
	}

	return SF;
}

float FatJetBTagScaleFactors::getTTBAR_dak8_SF(const FatJet* fj, YR yr, CORRTYPE corrT) const {
	float sf = 1.0, norm = 1.0;
	float pt = fj->pt();

	// {300-600, 600-800, 800+} according to http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_279_v10.pdf (page 25)
	std::vector<float> sf16    {1.039,1.035,1.301};
	std::vector<float> sf16_up {0.061,0.105,0.325};
	std::vector<float> sf16_dn {0.058,0.098,0.266};

	std::vector<float> norm16    {0.72,0.65,0.52};
	std::vector<float> norm16_up {0.05,0.06,0.07};
	std::vector<float> norm16_dn {0.05,0.06,0.07};

	std::vector<float> sf17    {0.91,0.93,1.07};
	std::vector<float> sf17_up {0.05,0.11,0.28};
	std::vector<float> sf17_dn {0.05,0.09,0.25};

	std::vector<float> norm17    {0.85,0.87,0.74};
	std::vector<float> norm17_up {0.06,0.08,0.09};
	std::vector<float> norm17_dn {0.06,0.08,0.09};

	std::vector<float> sf18    {0.89,0.94,1.05};
	std::vector<float> sf18_up {0.04,0.08,0.21};
	std::vector<float> sf18_dn {0.05,0.08,0.19};

	std::vector<float> norm18    {0.83,0.89,0.86};
	std::vector<float> norm18_up {0.06,0.08,0.09};
	std::vector<float> norm18_dn {0.06,0.08,0.09};

	int idx = 2;
	if(pt <= 600) idx = 0;
	else if(pt <= 800) idx = 1;

	if(yr == YR_2016) {
		sf = sf16[idx];
		norm = norm16[idx];
		if(corrT == UP) {
			sf += sf16_up[idx];
			norm += norm16_up[idx];
		} else if(corrT == DOWN) {
			sf -= sf16_dn[idx];
			norm -= norm16_dn[idx];
		}
	} else if(yr == YR_2017) {
		sf = sf17[idx];
		norm = norm17[idx];
		if(corrT == UP) {
			sf += sf17_up[idx];
			norm += norm17_up[idx];
		} else if(corrT == DOWN) {
			sf -= sf17_dn[idx];
			norm -= norm17_dn[idx];
		}
	} else if(yr == YR_2018) {
		sf = sf18[idx];
		norm = norm18[idx];
		if(corrT == UP) {
			sf += sf18_up[idx];
			norm += norm18_up[idx];
		} else if(corrT == DOWN) {
			sf -= sf18_dn[idx];
			norm -= norm18_dn[idx];
		}
	}

	return sf*norm;
}

}

