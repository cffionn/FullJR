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

dateStrPPMC=20190604
dateStrPbPbMC=20190604
dateStrPPData=20190604
dateStrPbPbData=20190604

DATE=`date +%Y%m%d`

mkdir -p /data/cmcginn/logs
mkdir -p /data/cmcginn/logs/$DATE

ppVals=(ak2PF ak3PF ak4PF ak6PF ak8PF ak10PF)
pbpbVals=(akCs2PU3PFFlow akCs3PU3PFFlow akCs4PU3PFFlow akCs6PU3PFFlow akCs8PU3PFFlow akCs10PU3PFFlow)

#ppVals=(ak10PF)
#pbpbVals=(akCs10PU3PFFlow)

PPDataFilePre=output/"$dateStrPPData"/dataPP_ProcessRawData_
PPDataFilePost=JetAnalyzer_"$dateStrPPData".root
PPResFilePre=output/"$dateStrPPMC"/combinedResponse_
PPResFilePost=_"$dateStrPPMC".root

PbPbDataFilePre=output/"$dateStrPbPbData"/dataPbPb_ProcessRawData_
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
	    dummy=0
#	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBinsPP.txt 0 100 >& /data/cmcginn/logs/$DATE/unfoldPP_NoClean_100_$i.log &
	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBinsPP.txt 0 1000 >& /data/cmcginn/logs/$DATE/unfoldPP_NoClean_1000_$i.log &
#	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBinsPP.txt 1 1000 >& /data/cmcginn/logs/$DATE/unfoldPP_NoClean_1000_$i.log &
#	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBinsPP.txt 0 1000 >& /data/cmcginn/logs/$DATE/unfoldPP_NoClean_1000_$i.log &
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
#	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBins.txt 0 100 >& /data/cmcginn/logs/$DATE/unfoldPbPb_NoClean_100_$i.log &
#	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBins_NOSPLIT.txt 0 100 >& /data/cmcginn/logs/$DATE/unfoldPbPb_NoClean_100_NOSPLIT_$i.log &
#	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBins_SPLIT.txt 0 100 >& /data/cmcginn/logs/$DATE/unfoldPbPb_NoClean_100_SPLIT_$i.log &
	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBins.txt 0 1000 >& /data/cmcginn/logs/$DATE/unfoldPbPb_NoClean_1000_$i.log &
#	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBins.txt 1 100 >& /data/cmcginn/logs/$DATE/unfoldPbPb_NoClean_100_$i.log &
#	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBins_SPLIT.txt 0 100 >& /data/cmcginn/logs/$DATE/unfoldPbPb_NoClean_100_SPLIT_$i.log &
	else
            echo " Missing $fileNameMC, continue on $i"
        fi
    else
        echo " Missing $fileNameData, continue on $i"
    fi
done

wait

#for i in "${ppVals[@]}"
#do
#    echo "Processing $i..."
#    jetTrunc=${i%PF*}
#
#    fileNameData=$PPDataFilePre$i$PPDataFilePost
#    fileNameMC=$PPResFilePre$jetTrunc$PPResFilePost
#
#    if [[ -f $fileNameData ]]
#    then
#	if [[ -f $fileNameMC ]]
#	then
#	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBinsPP.txt 1 >& logs/$DATE/unfoldPP_Clean_$i.log &
#	else
#	    echo " Missing $fileNameMC, continue on $i"
#	fi
#    else
#	echo " Missing $fileNameData, continue on $i"
#    fi
#done
#
#for i in "${pbpbVals[@]}"
#do
#    echo "Processing $i..."
#    jetTrunc=${i%PU*}
#
#    echo $jetTrunc
#
#    fileNameData=$PbPbDataFilePre$i$PbPbDataFilePost
#    fileNameMC=$PbPbResFilePre$jetTrunc$PbPbResFilePost
#
#    if [[ -f $fileNameData ]]
#    then
#        if [[ -f $fileNameMC ]]
#        then
#	    ./bin/unfoldRawData.exe $fileNameData $fileNameMC $i tables/overrideBins.txt 1 >& logs/$DATE/unfoldPbPb_Clean_$i.log &
#	else
#            echo " Missing $fileNameMC, continue on $i"
#        fi
#    else
#        echo " Missing $fileNameData, continue on $i"
#    fi
#done

wait

echo "runUnfold.sh Complete!"
