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


dateStrPPMC=20190308
dateStrPbPbMC=20190308
dateStrPPData=20190308
dateStrPbPbData=20190308

DATE=`date +%Y%m%d`

mkdir -p logs
mkdir -p logs/$DATE

ppVals=(ak3PF ak4PF ak6PF ak8PF ak10PF)
pbpbVals=(akCs3PU3PFFlow akCs4PU3PFFlow akCs6PU3PFFlow akCs8PU3PFFlow akCs10PU3PFFlow)

PPDataFilePre=output/"$dateStrPPData"/HiForestAOD_HighPtJet80_HLTJet80_LargeROR_PtCut110_AbsEta5_20190220_11Lumi_190220_221659_561_OutOf561_MERGED_ProcessRawData_
PPDataFilePost=JetAnalyzer_"$dateStrPPData".root
PPResFilePre=output/"$dateStrPPMC"/combinedResponse_
PPResFilePost=_"$dateStrPPMC".root

PbPbDataFilePre=output/"$dateStrPbPbData"/HiForestAOD_HIHardProbes_HLTJet100_AllR_PtCut140_AbsEta3_20180626_21LumiPer_180626_152510_1050_OutOf1050_MERGED_ProcessRawData_
PbPbDataFilePost=JetAnalyzer_"$dateStrPbPbData".root
PbPbResFilePre=output/"$dateStrPbPbMC"/combinedResponse_
PbPbResFilePost=_"$dateStrPbPbMC".root

for i in "${ppVals[@]}"
do
    echo "Processing $i..."
    jetTrunc=${i%PF*}

    fileNameData=$PPDataFilePre$i$PPDataFilePost
    fileNameMC=$PPResFilePre$jetTrunc$PPResFilePost

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
    jetTrunc=${i%PU*}

    echo $jetTrunc

    fileNameData=$PbPbDataFilePre$i$PbPbDataFilePost
    fileNameMC=$PbPbResFilePre$jetTrunc$PbPbResFilePost

    if [[ -f $fileNameData ]]
    then
        if [[ -f $fileNameMC ]]
        then
	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i >& logs/$DATE/unfoldPbPb_$i.log &
	else
            echo " Missing $fileNameMC, continue on $i"
        fi
    else
        echo " Missing $fileNameData, continue on $i"
    fi
done

wait

echo "runUnfold.sh Complete!"
