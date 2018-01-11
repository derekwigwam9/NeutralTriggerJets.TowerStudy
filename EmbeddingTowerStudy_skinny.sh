#!/bin/bash
# 'EmbeddingTowerStudy_skinny.sh'
#
# Use this to run 'EmbeddingTowerStudy_skinny.C'
# in batch mode.

iFile=\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt11.match.root\"
oFile=\"pp200r9pt11u.et9vz55hightwr.root\"
root -b -q EmbeddingTowerStudy_skinny.C\("true","$iFile","$oFile"\)
