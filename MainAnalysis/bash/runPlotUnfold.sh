#!/bin/bash

runDir=$FULLJRDIR/MainAnalysis
curDir=$PWD

if [[ "$runDir" == "$curDir" ]]
then
    dummyVal=0
else
    echo "Please run from $PWD, as bash bash/runPlotUnfold.sh. exit 1"
    exit 1
fi

if [[ -f $runDir/bin/plotUnfoldedSpectra.exe ]]
then
    dummyVal=0
else
    echo "Necessary executable ./bin/plotUnfoldedSpectra.exe is missing. run make. exit 1"
    exit 1
fi

dateStrPP=20190129
dateStrPbPb=20190129

DATE=`date +%Y%m%d`

mkdir -p logs
mkdir -p logs/$DATE

pbpbVals=(akCs3PU3PFFlow akCs4PU3PFFlow akCs6PU3PFFlow akCs8PU3PFFlow akCs10PU3PFFlow CombinedAlgos)
ppVals=(ak3PF ak4PF ak6PF ak8PF ak10PF CombinedAlgos)

fileOutPbPbPre=output/"$dateStrPbPb"/HiForestAOD_HIHardProbes_HLTJet100_AllR_Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbP_UnfoldRawData_NSuperBayes0_
fileOutPbPbPost=_"$dateStrPbPb".root

fileOutPPPre=output/"$dateStrPP"/HiForestAOD_HighPtJet80_HLTJet80_LargeRO_Pythia6_Dijet_pp502_MCDijet_20180712_Exc_UnfoldRawData_NSuperBayes0_
fileOutPPPost=_"$dateStrPP".root

pos=0
for i in "${pbpbVals[@]}"
do
    echo "Processing $i, ${ppVals[$pos]}"
    
    fileNamePP=$fileOutPPPre${ppVals[$pos]}$fileOutPPPost
    fileNamePbPb=$fileOutPbPbPre$i$fileOutPbPbPost

    if [[ -f $fileNamePP ]]
    then
	if [[ -f $fileNamePbPb ]]
	then
	    ./bin/plotUnfoldedSpectra.exe $fileNamePP $fileNamePbPb >& logs/$DATE/plotUnfold_${ppVals[$pos]}_$i.log &
	else
	    echo " $fileNamePbPb not found, continuing on $i, ${ppVals[$pos]}"
	fi
    else
	echo " $fileNamePP not found, continuing on $i, ${ppVals[$pos]}"
    fi

    pos=$((pos + 1))    
done

wait

echo "runPlotUnfold.sh complete!"
