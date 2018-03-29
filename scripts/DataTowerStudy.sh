#!/bin/bash
# 'DataTowerStudy.sh'
#
# Use this to run 'DataTowerStudy.C'
# in batch mode.

iFile=\"TowerResponseStudy/Run9/pp200r9.epCalc2.root\"
oFile=\"pp200r9.highTwr.d8m9y2017.root\"

root -b -q DataTowerStudy.C\("true","$iFile","$oFile"\)
