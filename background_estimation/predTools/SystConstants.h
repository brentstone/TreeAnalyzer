#ifndef ANALYSISTREEMAKER_BACKGROUNDESTIMATION_SYSTCONSTANTS_H
#define ANALYSISTREEMAKER_BACKGROUNDESTIMATION_SYSTCONSTANTS_H
#include<vector>
#include<utility>
#include <string>
#include "AnalysisSupport/Utilities/interface/Types.h"
namespace SystConstants{

// tau21 LP
float sf_tau21_LP_2016 = 0.95;
float sf_tau21_LP_2017 = 1.14;
float sf_tau21_LP_2018 = 1.20;
float sf_tau21_LP_Run2 = 1.116;
float unc_tau21_LP_2016 = 0.33;
float unc_tau21_LP_2017 = 0.29;
float unc_tau21_LP_2018 = 0.275;
float unc_tau21_LP_Run2 = 0.294;

// tau21 HP
float sf_tau21_HP_2016 = 1.03;
float sf_tau21_HP_2017 = 0.97;
float sf_tau21_HP_2018 = 0.98;
float sf_tau21_HP_Run2 = 0.99;
float unc_tau21_HP_2016 = 0.14;
float unc_tau21_HP_2017 = 0.06;
float unc_tau21_HP_2018 = 0.027;
float unc_tau21_HP_Run2 = 0.067;

// hbb_scale
float sf_sdmass_scale_2016 = 1.0;
float sf_sdmass_scale_2017 = 0.982;
float sf_sdmass_scale_2018 = 0.997;
float sf_sdmass_scale_Run2 = 0.993;
float unc_sdmass_scale_2016 = 0.0094;
float unc_sdmass_scale_2017 = 0.004;
float unc_sdmass_scale_2018 = 0.004;
float unc_sdmass_scale_Run2 = 0.0054;

// hbb_res
float sf_sdmass_res_2016 = 1.0;
float sf_sdmass_res_2017 = 1.09;
float sf_sdmass_res_2018 = 1.243;
float sf_sdmass_res_Run2 = 1.133;
float unc_sdmass_res_2016 = 0.2;
float unc_sdmass_res_2017 = 0.05;
float unc_sdmass_res_2018 = 0.041;
float unc_sdmass_res_Run2 = 0.086;


}


#endif

