#!/bin/bash

ppVals=(ak3PF ak4PF ak6PF ak8PF ak10PF)
pbpbVals=(akCs3PU3PFFlow akCs4PU3PFFlow akCs6PU3PFFlow akCs8PU3PFFlow akCs10PU3PFFlow)

pbpbFilePre=output/Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbPb_MCDijet_20180521_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180608_SVM_
pbpbFilePost=JetAnalyzer_FracNEntries1p00_JetResponse_20180830.root

ppFilePre=output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_
ppFilePost=JetAnalyzer_FracNEntries1p00_JetResponse_20180830.root

for i in "${pbpbVals[@]}"
do
    ./bin/plotJetResponse.exe $pbpbFilePre$i$pbpbFilePost >& logs/plotJetResPbPb_$i.log &
done

for i in "${ppVals[@]}"
do
    ./bin/plotJetResponse.exe $ppFilePre$i$ppFilePost >& logs/plotJetResPP_$i.log &
done

wait

echo "runPlotJetRes.sh complete"