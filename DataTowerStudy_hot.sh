#!/bin/bash
# 'DataTowerStudy_hot.sh'
#
# Use this to run 'DataTowerStudy_hot.C'
# in batch mode.

iFile=\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/merged/pp200r9.epCalc2.root\"
oFile=\"pp200r9.trgEtaVsClustAndPtCut.root\"
root -b -q DataTowerStudy_hot.C\("true","$iFile","$oFile"\)
