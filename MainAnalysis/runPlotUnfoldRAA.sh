#!/bin/bash

files=(output/plotUnfoldedSpectra_NSuperBayes0_akCs10PU3PFFlowJetAnalyzer_ak10PFJetAnalyzer_20180913.root output/plotUnfoldedSpectra_NSuperBayes0_akCs3PU3PFFlowJetAnalyzer_ak3PFJetAnalyzer_20180913.root output/plotUnfoldedSpectra_NSuperBayes0_akCs4PU3PFFlowJetAnalyzer_ak4PFJetAnalyzer_20180913.root output/plotUnfoldedSpectra_NSuperBayes0_akCs6PU3PFFlowJetAnalyzer_ak6PFJetAnalyzer_20180913.root output/plotUnfoldedSpectra_NSuperBayes0_akCs8PU3PFFlowJetAnalyzer_ak8PFJetAnalyzer_20180913.root)

fileStr=""

for i in "${files[@]}"
do
    echo $i
    fileStr="$i,$fileStr"
done

echo ""
./bin/plotUnfoldedRAA.exe "$fileStr"
echo ""

wait

echo "runPlotUnfoldRAA.sh complete!"