#!/bin/bash
# 'EmbeddingTowerStudy_geant.sh'
#
# Use this to run 'EmbeddingTowerStudy_geant.C'
# in batch mode.

iFile=\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/embedding/pt5_-1.ePcalc2.root\"
oFile=\"pp200r12pt5g.gntPi0.root\"
root -b -q EmbeddingTowerStudy_geant.C\("true","$iFile","$oFile"\)
