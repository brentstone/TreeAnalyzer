#ifndef TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_FATJETBTAGSCALEFACTORS_H_
#define TREEANALYZER_FRAMEWORK_PROCESSORS_CORRECTIONS_INTERFACE_FATJETBTAGSCALEFACTORS_H_
#include <string>
#include "AnalysisSupport/Utilities/interface/TObjectHelper.h"
#include "AnalysisSupport/Utilities/interface/Types.h"
#include "Processors/Variables/interface/BTagging.h"
#include "Configuration/interface/ReaderConstants.h"
#include "CorrectionHelper.h"

namespace TAna {

//class Jet;
class FatJet;

class FatJetBTagScaleFactors {
protected:
	enum YR {YR_NONE, YR_2016,YR_2017,YR_2018};
	std::map<std::string,YR> yrs = { {"2016",YR_2016}, {"2017",YR_2017}, {"2018",YR_2018} };

	struct sfLine {
		YR year;
		int wp;
		float ptmin;
		float ptmax;
		float sf;
		float sfUp;
		float sfDn;
	};

	std::vector<sfLine> sflines;
	void read_csv(YR year);
	void fillSfLine(int i, sfLine& line, std::string& colVal);
	std::string dataDir;
	std::string sfFile="";

public:
	FatJetBTagScaleFactors(const std::string& dataDir) : dataDir(dataDir) {};
	~FatJetBTagScaleFactors() {};
	void setParameters(const FatJetParameters& parameters, int year);
	float getSF(const FatJetParameters& parameters, const std::vector<const FatJet*>& fatJets,
			CorrHelp::CORRTYPE corrT = CorrHelp::NOMINAL) const;

};

}



#endif
