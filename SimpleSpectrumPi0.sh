#!/bin/bash
# 'SimpleSpectrumPi0.sh'
# Derek Anderson
#
# Use this to run 'SimpleSpectrumPi0.sh'
# in batch mode.

# input, output files
input="\"/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/merged/pp200r9.epCalc2.root\""
output="\"pp200r9.simpleAllTrks.et9vz55nHot41pi0.root\""

# run script
root -b -q "SimpleSpectrumPi0.C(true, $input, $output)"
