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

dateStrPP=20190606
dateStrPbPb=20190606

DATE=`date +%Y%m%d`

mkdir -p logs
mkdir -p logs/$DATE

pbpbVals=(CombinedAlgos)
ppVals=(CombinedAlgos)

fileOutPbPbPre=output/"$dateStrPbPb"/dataPbPb_UnfoldRawData_
fileOutPbPbPost=_NoClean_NToy1000_"$dateStrPbPb".root

fileOutPPPre=output/"$dateStrPP"/dataPP_UnfoldRawData_
fileOutPPPost=_NoClean_NToy1000_"$dateStrPP".root

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
	    ./bin/plotUnfoldedAll.exe $fileNamePP $fileNamePbPb $atlasFile $i tables/overrideBinsPlot.txt >& logs/$DATE/plotUnfoldAll_${ppVals[$pos]}_$i.log &
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
