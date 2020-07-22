#!/bin/bash
# ./getAreaReady.sh /uscms/home/bwstone/nobackup/TreeAnalyzer/ betrees 
inputdir=$1
remoteFld=$2

skimLoc="/Users/brentstone/Dropbox/Physics/GitRepos/TreeAnalyzer/background_estimation/predMacros/skimTree.C"
OUTSRV="bwstone@cmslpc111.fnal.gov:"

for fld in 2016 2017 2018 Run2
do
  mkdir ${fld}
  mkdir ${fld}/bkgInputs
  mkdir ${fld}/bkgInputsTopCR
  mkdir ${fld}/bkgInputsNonTopCR
  mkdir ${fld}/signalInputs
  mkdir ${fld}/supportInputs
  mkdir ${fld}/bkgCompLMT
  mkdir ${fld}/bkgCompAB
done

for fld in 2016 2017 2018
do
  yr=$(echo ${fld} | tr -dc '1,6-8')
  rsync -azP ${OUTSRV}${inputdir}/${remoteFld}${yr}/*.root ${fld}/
done

wait

for fld in 2016 2017 2018
do
  mkdir ${fld}/mcPieces
  mkdir ${fld}/dataPieces
  mkdir ${fld}/signalPieces
  mv ${fld}/data*.root ${fld}/dataPieces/
  mv ${fld}/radion*.root ${fld}/signalPieces/
  mv ${fld}/bulkgrav*.root ${fld}/signalPieces/
  mv ${fld}/*.root ${fld}/mcPieces/
  hadd ${fld}/betrees_mc.root ${fld}/mcPieces/*.root
  hadd ${fld}/betrees_data.root ${fld}/dataPieces/*.root

  for mx in 800 900 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000 4500
  do
  	hadd ${fld}/signalPieces/spin0_m${mx}.root ${fld}/signalPieces/radion_bbVV_m${mx}.root ${fld}/signalPieces/radion_bbtautau_m${mx}.root
  	hadd ${fld}/signalPieces/spin2_m${mx}.root ${fld}/signalPieces/bulkgrav_bbVV_m${mx}.root ${fld}/signalPieces/bulkgrav_bbtautau_m${mx}.root
  	hadd ${fld}/signalPieces/spinE_m${mx}.root ${fld}/signalPieces/spin*_m${mx}.root
  done

done

for fld in 2016 2017 2018
do
  RCMD="root -l -b -q '${skimLoc}(\"${fld}/betrees_mc.root\",\"${fld}/bkgCompLMT/betrees\",\"hbbTag>=0.8\",true)'"
  eval ${RCMD}
  RCMD="root -l -b -q '${skimLoc}(\"${fld}/betrees_mc.root\",\"${fld}/bkgCompAB/betrees\",\"hbbTag<=0.05\",true)'"
  eval ${RCMD}
done

hadd Run2/betrees_mc.root 201?/betrees_mc.root
hadd Run2/betrees_data.root 201?/betrees_data.root
mkdir Run2/signalPieces
for mx in 800 900 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000 4500
do
  hadd Run2/signalPieces/spin0_m${mx}.root 201?/signalPieces/spin0_m${mx}.root
  hadd Run2/signalPieces/spin2_m${mx}.root 201?/signalPieces/spin2_m${mx}.root
  hadd Run2/signalPieces/spinE_m${mx}.root Run2/signalPieces/spin*_m${mx}.root
done

for bk in qg mw mt losttw
do
  hadd Run2/bkgCompLMT/betrees_${bk}.root 201?/bkgCompLMT/betrees_${bk}.root
  hadd Run2/bkgCompAB/betrees_${bk}.root 201?/bkgCompAB/betrees_${bk}.root
done


