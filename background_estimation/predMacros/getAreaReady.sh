#!/bin/bash
# ./getAreaReady.sh /uscms/home/bwstone/nobackup/TreeAnalyzer/ betrees 0
inputdir=$1
remoteFld=$2
copyFromRemote=$3

skimLoc="/Users/brentstone/Dropbox/Physics/GitRepos/TreeAnalyzer/background_estimation/predMacros/skimTree.C"
OUTSRV="bwstone@cmslpc-sl7.fnal.gov:"

for fld in 2016 2017 2018 Run2
do
  mkdir -p ${fld}
  mkdir -p ${fld}/bkgInputs
  mkdir -p ${fld}/bkgInputsTopCR
  mkdir -p ${fld}/bkgInputsNonTopCR
  mkdir -p ${fld}/signalInputs
  mkdir -p ${fld}/supportInputs
  mkdir -p ${fld}/bkgCompLMT
  mkdir -p ${fld}/bkgCompAB
  mkdir -p ${fld}/signalPieces
  rm ${fld}/betrees_mc.root
  rm ${fld}/betrees_data.root
done

if [ ${copyFromRemote} -eq 1 ]
then
  for fld in 2016 2017 2018
  do
    yr=$(echo ${fld} | tr -dc '1,6-8')
    rsync -azP ${OUTSRV}${inputdir}/${remoteFld}${yr}/*.root ${fld}/
  done
fi

wait

for fld in 2016 2017 2018
do
  mkdir -p ${fld}/mcPieces
  mkdir -p ${fld}/dataPieces
  mv ${fld}/data*.root ${fld}/dataPieces/
  mv ${fld}/radion*.root ${fld}/signalPieces/
  mv ${fld}/bulkgrav*.root ${fld}/signalPieces/
  mv ${fld}/*.root ${fld}/mcPieces/
  hadd -f ${fld}/betrees_mc.root ${fld}/mcPieces/*.root
  hadd -f ${fld}/betrees_data.root ${fld}/dataPieces/*.root

  rm ${fld}/signalPieces/spin*.root
  for mx in 800 900 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000 4500
  do
  	hadd -f ${fld}/signalPieces/spin0_m${mx}.root ${fld}/signalPieces/radion_bbVV_m${mx}.root ${fld}/signalPieces/radion_bbtautau_m${mx}.root
  	hadd -f ${fld}/signalPieces/spin2_m${mx}.root ${fld}/signalPieces/bulkgrav_bbVV_m${mx}.root ${fld}/signalPieces/bulkgrav_bbtautau_m${mx}.root
  	hadd -f ${fld}/signalPieces/spinE_m${mx}.root ${fld}/signalPieces/spin*_m${mx}.root

  	for sys in HEM JESUp JESDOWN JERUp JERDown METUp METDOWN
  	do
  	  hadd -f ${fld}/signalPieces/spin0_m${mx}_${sys}.root ${fld}/signalPieces/radion_bbVV_m${mx}_${sys}.root ${fld}/signalPieces/radion_bbtautau_m${mx}_${sys}.root
  	  hadd -f ${fld}/signalPieces/spin2_m${mx}_${sys}.root ${fld}/signalPieces/bulkgrav_bbVV_m${mx}_${sys}.root ${fld}/signalPieces/bulkgrav_bbtautau_m${mx}_${sys}.root
  	  hadd -f ${fld}/signalPieces/spinE_m${mx}_${sys}.root ${fld}/signalPieces/spin*_m${mx}_${sys}.root
  	done

  done

  RCMD="root -l -b -q '${skimLoc}(\"${fld}/betrees_mc.root\",\"${fld}/bkgCompLMT/betrees\",\"hbbTag>=0.8\",true)' &"
  eval ${RCMD}
  RCMD="root -l -b -q '${skimLoc}(\"${fld}/betrees_mc.root\",\"${fld}/bkgCompAB/betrees\",\"hbbTag<=0.05\",true)' &"
  eval ${RCMD}
  
  wait
  
done

rm Run2/signalPieces/spin*.root
rm Run2/betrees_*.root

hadd Run2/betrees_mc.root 201*/betrees_mc.root
hadd Run2/betrees_data.root 201*/betrees_data.root

rm Run2/signalPieces/spin*.root
for mx in 800 900 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000 4500
do
  hadd -f Run2/signalPieces/spin0_m${mx}.root 201?/signalPieces/spin0_m${mx}.root
  hadd -f Run2/signalPieces/spin2_m${mx}.root 201?/signalPieces/spin2_m${mx}.root
  hadd -f Run2/signalPieces/spinE_m${mx}.root Run2/signalPieces/spin*_m${mx}.root
  
  for sys in HEM JESUp JESDOWN JERUp JERDown METUp METDOWN
  do
  	hadd -f Run2/signalPieces/spin0_m${mx}_${sys}.root 201?/signalPieces/spin0_m${mx}_${sys}.root
  	hadd -f Run2/signalPieces/spin2_m${mx}_${sys}.root 201?/signalPieces/spin2_m${mx}_${sys}.root
  	hadd -f Run2/signalPieces/spinE_m${mx}_${sys}.root Run2/signalPieces/spin*_m${mx}_${sys}.root
  done

done

for bk in qg mw mt losttw
do
  hadd -f Run2/bkgCompLMT/betrees_${bk}.root 201?/bkgCompLMT/betrees_${bk}.root
  hadd -f Run2/bkgCompAB/betrees_${bk}.root 201?/bkgCompAB/betrees_${bk}.root
done

