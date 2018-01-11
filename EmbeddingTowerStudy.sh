#!/bin/bash
# 'EmbeddingTowerStudy.sh'
#
# Use this to run 'EmbeddingTowerStudy.C'
# in batch mode.

iFile=\"../../Embedding/Run12pp/MuDstMatching/output/merged/pt5_-1.ePcalc2.root\"
oFile=\"test2.root\"

root -b -q EmbeddingTowerStudy.C\("true","$iFile","$oFile"\)
