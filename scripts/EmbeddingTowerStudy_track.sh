#!/bin/bash
# 'EmbeddingTowerStudy_track.sh'
#
# Use this to run 'EmbeddingTowerStudy_skinny.C'
# in batch mode.

iFile=\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt5.match.root\"
oFile=\"pp200r9pt5u.et9vz55track.root\"
root -b -q EmbeddingTowerStudy_track.C\("true","$iFile","$oFile"\)
