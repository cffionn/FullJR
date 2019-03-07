#!/bin/bash

runDir=$FULLJRDIR/MainAnalysis
curDir=$PWD

if [[ "$runDir" == "$curDir" ]]
then
    dummyVal=0
else
    echo "Please run from $PWD, as bash bash/runPlotUnfoldAll.sh. exit 1"
    exit 1
fi

if [[ -f $runDir/bin/plotUnfoldedAll.exe ]]
then
    dummyVal=0
else
    echo "Necessary executable ./bin/plotUnfoldedAll.exe is missing. run make. exit 1"
    exit 1
fi

dateStrPP=20190306
dateStrPbPb=20190306

DATE=`date +%Y%m%d`

mkdir -p logs
mkdir -p logs/$DATE

#pbpbVals=(akCs3PU3PFFlow akCs4PU3PFFlow akCs6PU3PFFlow akCs8PU3PFFlow akCs10PU3PFFlow CombinedAlgos)
#ppVals=(ak3PF ak4PF ak6PF ak8PF ak10PF CombinedAlgos)
pbpbVals=(akCs3PU3PFFlow akCs4PU3PFFlow akCs6PU3PFFlow akCs8PU3PFFlow akCs10PU3PFFlow)
ppVals=(ak3PF ak4PF ak6PF ak8PF ak10PF)

fileOutPbPbPre=output/"$dateStrPbPb"/HiForestAOD_HIHardProbes_HLTJet100_AllR_PtCut140_AbsEta3_20180626_21LumiPer_180626_152510_1050_OutOf1050_MERGED_UnfoldRawData_NSuperBayes0_
fileOutPbPbPost=_"$dateStrPbPb".root

fileOutPPPre=output/"$dateStrPP"/HiForestAOD_HighPtJet80_HLTJet80_LargeROR_PtCut110_AbsEta5_20190220_11Lumi_190220_221659_561_OutOf561_MERGED_UnfoldRawData_NSuperBayes0_
fileOutPPPost=_"$dateStrPP".root

atlasFile=output/HEPData-ins1673184-v1-root.root

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
	    ./bin/plotUnfoldedAll.exe $fileNamePP $fileNamePbPb $atlasFile $i >& logs/$DATE/plotUnfoldAll_${ppVals[$pos]}_$i.log &
	else
	    echo " $fileNamePbPb not found, continuing on $i, ${ppVals[$pos]}"
	fi
    else
	echo " $fileNamePP not found, continuing on $i, ${ppVals[$pos]}"
    fi

    pos=$((pos + 1))    
done

wait

echo "runPlotUnfoldAll.sh complete!"
