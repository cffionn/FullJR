#!/bin/bash

runDir=$FULLJRDIR/MainAnalysis
curDir=$PWD

if [[ "$runDir" == "$curDir" ]]
then
    dummyVal=0
else
    echo "Please run from $PWD, as bash bash/runResponse.sh. exit 1"
    exit 1
fi

DATE=`date +%Y%m%d`

mkdir -p logs
mkdir -p logs/$DATE

jetsPP=(ak3PFJetAnalyzer ak4PFJetAnalyzer ak6PFJetAnalyzer ak8PFJetAnalyzer ak10PFJetAnalyzer)
jetsPbPb=(akCs3PU3PFFlowJetAnalyzer akCs4PU3PFFlowJetAnalyzer akCs6PU3PFFlowJetAnalyzer akCs8PU3PFFlowJetAnalyzer akCs10PU3PFFlowJetAnalyzer)

#jetsPP=(ak3PFJetAnalyzer)
#jetsPbPb=(akCs3PU3PFFlowJetAnalyzer)

for i in "${jetsPP[@]}"
do
    ./bin/makeJetResponseTree.exe paths/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_$i.txt 1 1.0 >& logs/$DATE/responsePP_$i.log &
done

for i in "${jetsPbPb[@]}"
do
    ./bin/makeJetResponseTree.exe paths/Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbPb_MCDijet_20180521_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180608_SVM_$i.txt 0 1.0 >& logs/$DATE/responsePbPb_$i.log &
done
