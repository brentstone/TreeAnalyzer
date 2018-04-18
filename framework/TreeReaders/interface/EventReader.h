#ifndef TREEANALYZER_TREEREADERS_EVENTREADER_H
#define TREEANALYZER_TREEREADERS_EVENTREADER_H
#include "TreeReaders/interface/BaseReader.h"
#include "DataFormats/interface/Met.h"

using ASTypes::size8 ;
using ASTypes::size16;
using ASTypes::size64;

namespace TAna{
class EventReader: public BaseReader {
public:
    EventReader(std::string branchName, bool isRealData);
	virtual ~EventReader();
	virtual void setup(TreeReadingWrapper * wrapper);
	virtual void processVars();

	//settings
	const bool realData;

public:
	//branches from the tree
	size   run                = 0;
	size   lumi               = 0;
	size64 event              = 0;
	size8  goodVtx            = 0;
	size16 npv                = 0;
	float  rho                = 0;
	float  met_pt             = 0;
	float  met_phi            = 0;
	float  met_sig            = 0;
	float  met_unclUp_pt      = 0;
	float  met_unclUp_phi     = 0;
	float  met_raw_pt         = 0;
	float  met_raw_phi        = 0;
	float  nTruePUInts        = 0;
	float  weight             = 0;
	size8  process            = 0;
	size8  dataset            = 0;
	size8  dataRun            = 0;
    std::vector<float>* genWeights       = new std::vector<float>;

	size   metFilters         = 0;
	size64 triggerAccepts     = 0;
	size64 triggerPrescales   = 0;

	//branches from post processing
	bool   normWeightLoaded     = false;
	float  normWeight           = 0;

	//objects created in process
	Met met;
	Met rawMet;



};
}

#endif
