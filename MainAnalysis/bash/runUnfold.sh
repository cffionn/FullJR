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

if [[ -f $runDir/bin/unfoldRawData.exe ]]
then
    dummyVal=0
else
    echo "Necessary executable ./bin/unfoldRawData.exe is missing. run make. exit 1"
    exit 1
fi


dateStrPP=20190129
dateStrPbPb=20190129

DATE=`date +%Y%m%d`

mkdir -p logs
mkdir -p logs/$DATE

ppVals=(ak3PF ak4PF ak6PF ak8PF ak10PF)
pbpbVals=(akCs3PU3PFFlow akCs4PU3PFFlow akCs6PU3PFFlow akCs8PU3PFFlow akCs10PU3PFFlow)

PPDataFilePre=output/"$dateStrPP"/HiForestAOD_HighPtJet80_HLTJet80_LargeRO_Pythia6_Dijet_pp502_MCDijet_20180712_Exc_ProcessRawData_
PPDataFilePost=JetAnalyzer_"$dateStrPP".root
PPResFilePre=output/"$dateStrPP"/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_
PPResFilePost=JetAnalyzer_FracNEntries1p00_JetResponse_"$dateStrPP".root

PbPbDataFilePre=output/"$dateStrPbPb"/HiForestAOD_HIHardProbes_HLTJet100_AllR__Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbP_ProcessRawData_
PbPbDataFilePost=JetAnalyzer_"$dateStrPbPb".root
PbPbResFilePre=output/"$dateStrPbPb"/Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbPb_MCDijet_20180521_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180608_SVM_
PbPbResFilePost=JetAnalyzer_FracNEntries1p00_JetResponse_"$dateStrPbPb".root

for i in "${ppVals[@]}"
do
    echo "Processing $i..."
    fileNameData=$PPDataFilePre$i$PPDataFilePost
    fileNameMC=$PPResFilePre$i$PPResFilePost

    if [[ -f $fileNameData ]]
    then
	if [[ -f $fileNameMC ]]
	then
	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i >& logs/$DATE/unfoldPP_$i.log &
	else
	    echo " Missing $fileNameMC, continue on $i"
	fi
    else
	echo " Missing $fileNameData, continue on $i"
    fi
done

for i in "${pbpbVals[@]}"
do
    echo "Processing $i..."
    fileNameData=$PbPbDataFilePre$i$PbPbDataFilePost
    fileNameMC=$PbPbResFilePre$i$PbPbResFilePost

    if [[ -f $fileNameData ]]
    then
        if [[ -f $fileNameMC ]]
        then
	    ./bin/unfoldRawData.exe $PbPbDataFilePre$i$PbPbDataFilePost $PbPbResFilePre$i$PbPbResFilePost $i >& logs/$DATE/unfoldPbPb_$i.log &
	else
            echo " Missing $fileNameMC, continue on $i"
        fi
    else
        echo " Missing $fileNameData, continue on $i"
    fi
done

wait

echo "runUnfold.sh Complete!"
