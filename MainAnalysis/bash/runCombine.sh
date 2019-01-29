#!/bin/bash

runDir=$FULLJRDIR/MainAnalysis
curDir=$PWD

if [[ "$runDir" == "$curDir" ]]
then
    dummyVal=0
else
    echo "Please run from $PWD, as bash bash/runUnfold.sh. exit 1"
    exit 1
fi

dateStrPP=20190129
dateStrPbPb=20190129

DATE=`date +%Y%m%d`

mkdir -p logs
mkdir -p logs/$DATE

ppVals=(ak3PF ak4PF ak6PF ak8PF ak10PF)
pbpbVals=(akCs3PU3PFFlow akCs4PU3PFFlow akCs6PU3PFFlow akCs8PU3PFFlow akCs10PU3PFFlow)

PPFilePre=output/"$dateStrPP"/HiForestAOD_HighPtJet80_HLTJet80_LargeRO_Pythia6_Dijet_pp502_MCDijet_20180712_Exc_UnfoldRawData_NSuperBayes0_
PPFilePost=_"$dateStrPP".root

PbPbFilePre=output/"$dateStrPbPb"/HiForestAOD_HIHardProbes_HLTJet100_AllR_Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbP_UnfoldRawData_NSuperBayes0_
PbPbFilePost=_"$dateStrPbPb".root

PPFileOut="$PPFilePre"CombinedAlgos"$PPFilePost"
PbPbFileOut="$PbPbFilePre"CombinedAlgos"$PbPbFilePost"

filesPbPbStr=""
for i in "${pbpbVals[@]}"
do
    echo "Adding algo $i..."
    newFile="$PbPbFilePre$i$PbPbFilePost"
    if [[ -f $newFile ]]
    then
	filesPbPbStr="$filesPbPbStr $newFile"
    else
	echo " File $newFile for algo $i not found. continue"
    fi
done

filesPpStr=""
for i in "${ppVals[@]}"
do
    echo "Adding algo $i..."
    newFile="$PPFilePre$i$PPFilePost"
    if [[ -f $newFile ]]
    then
	filesPPStr="$filesPPStr $PPFilePre$i$PPFilePost"
    else
        echo " File $newFile for algo $i not found. continue"
    fi

done

if [[ $filesPbPbStr == "" ]]
then
    echo "No files found for PbPb! no combo run"
else
    ./bin/combineFiles.exe $PbPbFileOut $filesPbPbStr >& logs/$DATE/comboPbPb.log &
fi

if [[ $filesPPStr == "" ]]
then
    echo "No files found for PP! no combo run"
else
    ./bin/combineFiles.exe $PPFileOut $filesPPStr >& logs/$DATE/comboPP.log &
fi

wait

echo "runCombine.sh complete!"
