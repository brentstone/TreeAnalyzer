#!/bin/bash
# ./getAreaReady.sh /uscms/home/nmccoll/nobackup/2011-04-15-susyra2/rel_HbbWW/work/analyzer_running/compiled/ trees/ /Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/background_estimation/predMacros/skimTree.C

startDir=$1
macroLoc=$2

cd ${startDir}/signalInputs
RCMD="root -b -q '${macroLoc}/makeSignalInputs.C+(3)'"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeSignalInputs.C+(0)' &"
eval $RCMD
cd ../signalInputsNoCond
RCMD="root -b -q '${macroLoc}/makeSignalInputs.C+(0)' &"
eval $RCMD
wait
RCMD="root -b -q '${macroLoc}/makeSignalInputs.C+(2)' &"
eval $RCMD
cd ../signalInputs
RCMD="root -b -q '${macroLoc}/makeSignalInputs.C+(1)' &"
eval $RCMD
cd ../

cd ${startDir}/bkgInputs
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(5)'"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(0)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(1)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(2)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(3)' &"
eval $RCMD
cd .. 

cd ${startDir}/bkgInputsCR
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(5)'"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(0)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(1)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(2)' &"
eval $RCMD
RCMD="root -b -q '${macroLoc}/makeBKGCRInputs.C+(3)' &"
eval $RCMD
cd .. 

wait
cd ${startDir}/bkgInputs
RCMD="root -b -q '${macroLoc}/makeBKGInputs.C+(4)' &"
eval $RCMD
cd .. 

