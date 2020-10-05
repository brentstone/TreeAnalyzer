#!/bin/bash
# Usage: ./compile.sh jobs/4_30_beTreees/ compiled/
inputdir=$1
outputdir=$2

mkdir -p ${outputdir}

# BKG
hadd ${outputdir}/wjets.root ${inputdir}/out_WJetsToLNu_HT*.root
hadd ${outputdir}/zjets_m50.root ${inputdir}/out_DYJetsToLL_M-50_HT*.root
hadd ${outputdir}/ttbar_0l.root ${inputdir}/out_TTToHadronic*.root
hadd ${outputdir}/ttbar_1l.root ${inputdir}/out_TTToSemiLeptonic*.root
hadd ${outputdir}/ttbar_2l.root ${inputdir}/out_TTTo2L2Nu*.root
#hadd ${outputdir}/ttbar_m700.root ${inputdir}/out_TT_Mtt-700to1000*.root
#hadd ${outputdir}/ttbar_m1000.root ${inputdir}/out_TT_Mtt-1000toInf*.root
hadd ${outputdir}/qcd.root ${inputdir}/out_QCD_HT*.root
hadd ${outputdir}/st.root ${inputdir}/out_ST_*.root
hadd ${outputdir}/ww.root ${inputdir}/out_WW*.root
hadd ${outputdir}/wz.root ${inputdir}/out_WZ*.root
hadd ${outputdir}/zz.root ${inputdir}/out_ZZ*.root
hadd ${outputdir}/ttz.root ${inputdir}/out_TTZ*.root
hadd ${outputdir}/ttw.root ${inputdir}/out_TTW*.root
hadd ${outputdir}/tth.root ${inputdir}/out_ttH*.root

# SIGNAL
for mx in 800 900 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000 4500
do
  hadd ${outputdir}/radion_bbVV_m${mx}.root       ${inputdir}/out_Radion_hh_hVV*M-${mx}*[0-9].root
  hadd ${outputdir}/radion_bbtautau_m${mx}.root   ${inputdir}/out_Radion_hh_htata*M${mx}*[0-9].root
  hadd ${outputdir}/bulkgrav_bbVV_m${mx}.root     ${inputdir}/out_BulkGravTohhTohVV*M-${mx}*[0-9].root
  hadd ${outputdir}/bulkgrav_bbtautau_m${mx}.root ${inputdir}/out_BulkGrav_hh_htata*M${mx}*[0-9].root

  for sys in HEM JESUp JESDOWN JERUp JERDown METUp METDOWN
  do
    hadd ${outputdir}/radion_bbVV_m${mx}_${sys}.root       ${inputdir}/out_Radion_hh_hVV*M-${mx}*_${sys}.root
    hadd ${outputdir}/radion_bbtautau_m${mx}_${sys}.root   ${inputdir}/out_Radion_hh_htata*M${mx}*_${sys}.root
    hadd ${outputdir}/bulkgrav_bbVV_m${mx}_${sys}.root     ${inputdir}/out_BulkGravTohhTohVV*M-${mx}*_${sys}.root
    hadd ${outputdir}/bulkgrav_bbtautau_m${mx}_${sys}.root ${inputdir}/out_BulkGrav_hh_htata*M${mx}*_${sys}.root
  done
done

yr=$(echo ${inputdir} | tr -dc '1,6-8')

# DATA AND LOW MASS DYJETSTOLL
if [ ${yr} -eq 16 ]
then
	hadd ${outputdir}/zjets_m5to50.root ${inputdir}/out_DYJetsToLL_M-5to50_HT*.root

	hadd ${outputdir}/data_mu-2016b.root      ${inputdir}/out_SingleMuon_Run2016B*.root
	hadd ${outputdir}/data_mu-2016c.root      ${inputdir}/out_SingleMuon_Run2016C*.root
	hadd ${outputdir}/data_mu-2016d.root      ${inputdir}/out_SingleMuon_Run2016D*.root
	hadd ${outputdir}/data_mu-2016e.root      ${inputdir}/out_SingleMuon_Run2016E*.root
	hadd ${outputdir}/data_mu-2016f.root      ${inputdir}/out_SingleMuon_Run2016F*.root
	hadd ${outputdir}/data_mu-2016g.root      ${inputdir}/out_SingleMuon_Run2016G*.root
	hadd ${outputdir}/data_mu-2016h.root      ${inputdir}/out_SingleMuon_Run2016H*.root
	hadd ${outputdir}/data_el-2016b.root      ${inputdir}/out_SingleElectron_Run2016B*.root
	hadd ${outputdir}/data_el-2016c.root      ${inputdir}/out_SingleElectron_Run2016C*.root
	hadd ${outputdir}/data_el-2016d.root      ${inputdir}/out_SingleElectron_Run2016D*.root
	hadd ${outputdir}/data_el-2016e.root      ${inputdir}/out_SingleElectron_Run2016E*.root
	hadd ${outputdir}/data_el-2016f.root      ${inputdir}/out_SingleElectron_Run2016F*.root
	hadd ${outputdir}/data_el-2016g.root      ${inputdir}/out_SingleElectron_Run2016G*.root
	hadd ${outputdir}/data_el-2016h.root      ${inputdir}/out_SingleElectron_Run2016H*.root
	hadd ${outputdir}/data_phot-2016b.root      ${inputdir}/out_SinglePhoton_Run2016B*.root
	hadd ${outputdir}/data_phot-2016c.root      ${inputdir}/out_SinglePhoton_Run2016C*.root
	hadd ${outputdir}/data_phot-2016d.root      ${inputdir}/out_SinglePhoton_Run2016D*.root
	hadd ${outputdir}/data_phot-2016e.root      ${inputdir}/out_SinglePhoton_Run2016E*.root
	hadd ${outputdir}/data_phot-2016f.root      ${inputdir}/out_SinglePhoton_Run2016F*.root
	hadd ${outputdir}/data_phot-2016g.root      ${inputdir}/out_SinglePhoton_Run2016G*.root
	hadd ${outputdir}/data_phot-2016h.root      ${inputdir}/out_SinglePhoton_Run2016H*.root
	hadd ${outputdir}/data_met-2016b.root      ${inputdir}/out_MET_Run2016B*.root
	hadd ${outputdir}/data_met-2016c.root      ${inputdir}/out_MET_Run2016C*.root
	hadd ${outputdir}/data_met-2016d.root      ${inputdir}/out_MET_Run2016D*.root
	hadd ${outputdir}/data_met-2016e.root      ${inputdir}/out_MET_Run2016E*.root
	hadd ${outputdir}/data_met-2016f.root      ${inputdir}/out_MET_Run2016F*.root
	hadd ${outputdir}/data_met-2016g.root      ${inputdir}/out_MET_Run2016G*.root
	hadd ${outputdir}/data_met-2016h.root      ${inputdir}/out_MET_Run2016H*.root
	hadd ${outputdir}/data_jetht-2016b.root      ${inputdir}/out_JetHT_Run2016B*.root
	hadd ${outputdir}/data_jetht-2016c.root      ${inputdir}/out_JetHT_Run2016C*.root
	hadd ${outputdir}/data_jetht-2016d.root      ${inputdir}/out_JetHT_Run2016D*.root
	hadd ${outputdir}/data_jetht-2016e.root      ${inputdir}/out_JetHT_Run2016E*.root
	hadd ${outputdir}/data_jetht-2016f.root      ${inputdir}/out_JetHT_Run2016F*.root
	hadd ${outputdir}/data_jetht-2016g.root      ${inputdir}/out_JetHT_Run2016G*.root
	hadd ${outputdir}/data_jetht-2016h.root      ${inputdir}/out_JetHT_Run2016H*.root
elif [ ${yr} -eq 17 ]
then
	hadd ${outputdir}/zjets_m4to50.root ${inputdir}/out_DYJetsToLL_M-4to50_HT*.root

	hadd ${outputdir}/data_mu-2017b.root      ${inputdir}/out_SingleMuon_Run2017B*.root
	hadd ${outputdir}/data_mu-2017c.root      ${inputdir}/out_SingleMuon_Run2017C*.root
	hadd ${outputdir}/data_mu-2017d.root      ${inputdir}/out_SingleMuon_Run2017D*.root
	hadd ${outputdir}/data_mu-2017e.root      ${inputdir}/out_SingleMuon_Run2017E*.root
	hadd ${outputdir}/data_mu-2017f.root      ${inputdir}/out_SingleMuon_Run2017F*.root
	hadd ${outputdir}/data_el-2017b.root      ${inputdir}/out_SingleElectron_Run2017B*.root
	hadd ${outputdir}/data_el-2017c.root      ${inputdir}/out_SingleElectron_Run2017C*.root
	hadd ${outputdir}/data_el-2017d.root      ${inputdir}/out_SingleElectron_Run2017D*.root
	hadd ${outputdir}/data_el-2017e.root      ${inputdir}/out_SingleElectron_Run2017E*.root
	hadd ${outputdir}/data_el-2017f.root      ${inputdir}/out_SingleElectron_Run2017F*.root
	hadd ${outputdir}/data_phot-2017b.root      ${inputdir}/out_SinglePhoton_Run2017B*.root
	hadd ${outputdir}/data_phot-2017c.root      ${inputdir}/out_SinglePhoton_Run2017C*.root
	hadd ${outputdir}/data_phot-2017d.root      ${inputdir}/out_SinglePhoton_Run2017D*.root
	hadd ${outputdir}/data_phot-2017e.root      ${inputdir}/out_SinglePhoton_Run2017E*.root
	hadd ${outputdir}/data_phot-2017f.root      ${inputdir}/out_SinglePhoton_Run2017F*.root
	hadd ${outputdir}/data_met-2017b.root      ${inputdir}/out_MET_Run2017B*.root
	hadd ${outputdir}/data_met-2017c.root      ${inputdir}/out_MET_Run2017C*.root
	hadd ${outputdir}/data_met-2017d.root      ${inputdir}/out_MET_Run2017D*.root
	hadd ${outputdir}/data_met-2017e.root      ${inputdir}/out_MET_Run2017E*.root
	hadd ${outputdir}/data_met-2017f.root      ${inputdir}/out_MET_Run2017F*.root
	hadd ${outputdir}/data_jetht-2017b.root      ${inputdir}/out_JetHT_Run2017B*.root
	hadd ${outputdir}/data_jetht-2017c.root      ${inputdir}/out_JetHT_Run2017C*.root
	hadd ${outputdir}/data_jetht-2017d.root      ${inputdir}/out_JetHT_Run2017D*.root
	hadd ${outputdir}/data_jetht-2017e.root      ${inputdir}/out_JetHT_Run2017E*.root
	hadd ${outputdir}/data_jetht-2017f.root      ${inputdir}/out_JetHT_Run2017F*.root
elif [ ${yr} -eq 18 ]
then
	hadd ${outputdir}/zjets_m4to50.root ${inputdir}/out_DYJetsToLL_M-4to50_HT*.root

	hadd ${outputdir}/data_mu-2018a.root      ${inputdir}/out_SingleMuon_Run2018A*.root
	hadd ${outputdir}/data_mu-2018b.root      ${inputdir}/out_SingleMuon_Run2018B*.root
	hadd ${outputdir}/data_mu-2018c.root      ${inputdir}/out_SingleMuon_Run2018C*.root
	hadd ${outputdir}/data_mu-2018d.root      ${inputdir}/out_SingleMuon_Run2018D*.root
	hadd ${outputdir}/data_el-2018a.root      ${inputdir}/out_EGamma_Run2018A*.root
	hadd ${outputdir}/data_el-2018b.root      ${inputdir}/out_EGamma_Run2018B*.root
	hadd ${outputdir}/data_el-2018c.root      ${inputdir}/out_EGamma_Run2018C*.root
	hadd ${outputdir}/data_el-2018d.root      ${inputdir}/out_EGamma_Run2018D*.root
	hadd ${outputdir}/data_met-2018a.root      ${inputdir}/out_MET_Run2018A*.root
	hadd ${outputdir}/data_met-2018b.root      ${inputdir}/out_MET_Run2018B*.root
	hadd ${outputdir}/data_met-2018c.root      ${inputdir}/out_MET_Run2018C*.root
	hadd ${outputdir}/data_met-2018d.root      ${inputdir}/out_MET_Run2018D*.root
	hadd ${outputdir}/data_jetht-2018a.root      ${inputdir}/out_JetHT_Run2018A*.root
	hadd ${outputdir}/data_jetht-2018b.root      ${inputdir}/out_JetHT_Run2018B*.root
	hadd ${outputdir}/data_jetht-2018c.root      ${inputdir}/out_JetHT_Run2018C*.root
	hadd ${outputdir}/data_jetht-2018d.root      ${inputdir}/out_JetHT_Run2018D*.root
else
	echo "COULD NOT FIND NUMBER IN 16, 17, OR 18"
fi

