#!/bin/bash

runDir=$FULLJRDIR/MainAnalysis
curDir=$PWD

if [[ "$runDir" == "$curDir" ]]
then
    dummyVal=0
else
    echo "Please run from $PWD, as bash bash/runCombineResponse.sh. exit 1"
    exit 1
fi

if [[ -f $runDir/bin/combineResponse.exe ]]
then
    dummyVal=0
else
    echo "Necessary executable ./bin/combineResponse.exe is missing. run make. exit 1"
    exit 1
fi


DATE=`date +%Y%m%d`

mkdir -p logs
mkdir -p logs/$DATE

FILEDATE=20190301

jets=(ak3PFJetAnalyzer ak4PFJetAnalyzer ak6PFJetAnalyzer ak8PFJetAnalyzer ak10PFJetAnalyzer akCs3PU3PFFlowJetAnalyzer akCs4PU3PFFlowJetAnalyzer akCs6PU3PFFlowJetAnalyzer akCs8PU3PFFlowJetAnalyzer akCs10PU3PFFlowJetAnalyzer)

for i in "${jets[@]}"
do
    files=""
    for j in output/$FILEDATE/*$i*Response*.root
    do
	files="$files$j,"
    done
#    echo $files
    echo "$runDir/bin/combineResponse.exe $files $i"
done


wait
echo "bash/runCombineResponse.sh Complete!"
