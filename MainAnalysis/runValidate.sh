#!/bin/bash

jetsPP=(ak3PFJetAnalyzer ak4PFJetAnalyzer ak6PFJetAnalyzer ak8PFJetAnalyzer ak10PFJetAnalyzer)
jetsPbPb=(akCs3PU3PFFlowJetAnalyzer akCs4PU3PFFlowJetAnalyzer akCs6PU3PFFlowJetAnalyzer akCs8PU3PFFlowJetAnalyzer akCs10PU3PFFlowJetAnalyzer)

PbPbMCFilePre=output/Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbPb_MCDijet_20180521_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180608_SVM_
PbPbMCFilePost=_FracNEntries1p00_JetResponse_20180911.root

PPMCFilePre=output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_
PPMCFilePost=_FracNEntries1p00_JetResponse_20180911.root

pos=0
for i in "${jetsPP[@]}"
do
    echo "Processing $i, ${jetsPbPb[$pos]}..."

    ./bin/validateJetResponse.exe $PbPbMCFilePre${jetsPbPb[$pos]}$PbPbMCFilePost >& logs/validatePbPb_${jetsPbPb[$pos]}.log &

    ./bin/validateJetResponse.exe  $PPMCFilePre$i$PPMCFilePost >& logs/validatePP_$i.log &

    pos=$((pos + 1))

    wait
done

echo "runValidate.sh Complete"