#!/bin/bash

ppVals=(ak3PF ak4PF ak6PF ak8PF ak10PF)
pbpbVals=(akCs3PU3PFFlow akCs4PU3PFFlow akCs6PU3PFFlow akCs8PU3PFFlow akCs10PU3PFFlow)

PPFilePre=output/HiForestAOD_HighPtJet80_HLTJet80_LargeRO_Pythia6_Dijet_pp502_MCDijet_20180712_Exc_UnfoldRawData_NSuperBayes0_
PPFilePost=_20180828.root

PbPbFilePre=output/HiForestAOD_HIHardProbes_HLTJet100_AllR_Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbP_UnfoldRawData_NSuperBayes0_
PbPbFilePost=_20180828.root

PPFileOut=$PPFilePre\CombinedAlgos$PPFilePost
PbPbFileOut=$PbPbFilePre\CombinedAlgos$PbPbFilePost

#editing here 20180828
filesPbPbStr=""
for i in "${pbpbVals[@]}"
do
    filesPbPbStr="$filesPbPbStr $PbPbFilePre$i$PbPbFilePost"
done

filesPpStr=""
for i in "${ppVals[@]}"
do
    filesPPStr="$filesPPStr $PPFilePre$i$PPFilePost"
done


./bin/combineFiles.exe $PbPbFileOut $filesPbPbStr >& logs/comboPbPb.log &
./bin/combineFiles.exe $PPFileOut $filesPPStr >& logs/comboPP.log &

wait

echo "combined.sh complete1"