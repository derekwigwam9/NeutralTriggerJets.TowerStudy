#!/bin/bash
# 'SimpleSpectrumHad.sh'
# Derek Anderson
#
# Use this to run 'SimpleSpectrumHad.C'
# in batch mode.

declare -a input
declare -a output

# input files
input[0]="\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt5.match.root\""
input[1]="\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt7.match.root\""
input[2]="\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt9.match.root\""
input[3]="\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt11.match.root\""
input[4]="\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt15.match.root\""
input[5]="\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt25.match.root\""
input[6]="\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/run9embedding/pt35.match.root\""

# output files
output[0]="\"pp200r9pt5.simpleAllTrks.et9vz55nHot41had.root\""
output[1]="\"pp200r9pt7.simpleAllTrks.et9vz55nHot41had.root\""
output[2]="\"pp200r9pt9.simpleAllTrks.et9vz55nHot41had.root\""
output[3]="\"pp200r9pt11.simpleAllTrks.et9vz55nHot41had.root\""
output[4]="\"pp200r9pt15.simpleAllTrks.et9vz55nHot41had.root\""
output[5]="\"pp200r9pt25.simpleAllTrks.et9vz55nHot41had.root\""
output[6]="\"pp200r9pt35.simpleAllTrks.et9vz55nHot41had.root\""

# run macro
for iFile in `seq 0 6`; do
  root -b -q "SimpleSpectrumHad.C(true, ${input[iFile]}, ${output[iFile]})"
done

unset input
unset output
