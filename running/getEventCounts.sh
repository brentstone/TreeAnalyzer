#!/bin/bash
python running/getEventCounts.py -i /uscms/home/bwstone/nobackup/AnalysisTreeMaker/CMSSW_10_2_16_patch2/src/AnalysisTreeMaker/TreeMaker/run/datasets_2016.conf -o dataset_confs/procdatasets2016_brent.conf -d /eos/uscms/store/user/bstone/HHbbWW/ntuples/2016/ &
sleep 5

python running/getEventCounts.py -i /uscms/home/bwstone/nobackup/AnalysisTreeMaker/CMSSW_10_2_16_patch2/src/AnalysisTreeMaker/TreeMaker/run/datasets_2016.conf -o dataset_confs/procdatasets2016_nick.conf -d /eos/uscms/store/user/nmccoll/HHWWbb_trees/10_21_19_ntuples_2016/ &
sleep 5

python running/getEventCounts.py -i /uscms/home/bwstone/nobackup/AnalysisTreeMaker/CMSSW_10_2_16_patch2/src/AnalysisTreeMaker/TreeMaker/run/datasets_2016.conf -o dataset_confs/procdatasets2016_lpchh.conf -d /eos/uscms/store/user/lpchh/HHWWbb_trees/10_21_19_ntuples_2016/ &
sleep 5

python running/getEventCounts.py -i /uscms/home/bwstone/nobackup/AnalysisTreeMaker/CMSSW_10_2_16_patch2/src/AnalysisTreeMaker/TreeMaker/run/datasets_2017.conf -o dataset_confs/procdatasets2017_brent.conf -d /eos/uscms/store/user/bstone/HHbbWW/ntuples/2017/ &
sleep 5

python running/getEventCounts.py -i /uscms/home/bwstone/nobackup/AnalysisTreeMaker/CMSSW_10_2_16_patch2/src/AnalysisTreeMaker/TreeMaker/run/datasets_2017.conf -o dataset_confs/procdatasets2017_nick.conf -d /eos/uscms/store/user/nmccoll/HHWWbb_trees/10_5_19_ntuples_2017/ &
sleep 5

python running/getEventCounts.py -i /uscms/home/bwstone/nobackup/AnalysisTreeMaker/CMSSW_10_2_16_patch2/src/AnalysisTreeMaker/TreeMaker/run/datasets_2017.conf -o dataset_confs/procdatasets2017_lpchh.conf -d /eos/uscms/store/user/lpchh/HHWWbb_trees/10_5_19_ntuples_2017/ &
sleep 5

python running/getEventCounts.py -i /uscms/home/bwstone/nobackup/AnalysisTreeMaker/CMSSW_10_2_16_patch2/src/AnalysisTreeMaker/TreeMaker/run/datasets_2018.conf -o dataset_confs/procdatasets2018_brent.conf -d /eos/uscms/store/user/bstone/HHbbWW/ntuples/2018/ &
sleep 5

python running/getEventCounts.py -i /uscms/home/bwstone/nobackup/AnalysisTreeMaker/CMSSW_10_2_16_patch2/src/AnalysisTreeMaker/TreeMaker/run/datasets_2018.conf -o dataset_confs/procdatasets2018_nick.conf -d /eos/uscms/store/user/nmccoll/HHWWbb_trees/10_8_19_ntuples_2018/ &
sleep 5

python running/getEventCounts.py -i /uscms/home/bwstone/nobackup/AnalysisTreeMaker/CMSSW_10_2_16_patch2/src/AnalysisTreeMaker/TreeMaker/run/datasets_2018.conf -o dataset_confs/procdatasets2018_lpchh.conf -d /eos/uscms/store/user/lpchh/HHWWbb_trees/10_8_19_ntuples_2018/ &
