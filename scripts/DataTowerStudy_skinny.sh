#!/bin/bash
# 'DataTowerStudy_skinny.sh'
#
# Use this to run 'DataTowerStudy_skinny.C'
# in batch mode.

iFile=\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/merged/pp200r9.merge.root\"
oFile=\"pp200r9.merge.root\"

root -b -q DataTowerStudy_skinny.C\("true","$iFile","$oFile"\)
