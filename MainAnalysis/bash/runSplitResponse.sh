#!/bin/bash

runDir=$FULLJRDIR/MainAnalysis
curDir=$PWD

if [[ "$runDir" == "$curDir" ]]
then
    dummyVal=0
else
    echo "Please run from $PWD, as bash bash/runResponse.sh. exit 1"
    exit 1
fi

if [[ -f $runDir/bin/makeJetResponseTree.exe ]]
then
    dummyVal=0
else
    echo "Necessary executable ./bin/makeJetResponseTree.exe is missing. run make. exit 1"
    exit 1
fi

jobNumber=20
if [[ -d /home/cfmcginn ]]
then
    jobNumber=12
fi

DATE=`date +%Y%m%d`
TIMESTART=`date +%H%M%S`

mkdir -p logs
mkdir -p logs/$DATE


for i in paths/*/*SPLIT*.txt
do
    pbpbBool=1
    if [[ $i == *"akCs"* ]]
    then
	pbpbBool=0
    fi

    logName=${i%.txt}
    while [[ $logName == *"/"* ]]
    do
	logName=${logName#*/}
    done

    ./bin/makeJetResponseTree.exe $i $pbpbBool 1.0 >& logs/$DATE/response_$logName.log &

    counts=$(ps | grep makeJet | wc -l)
    while [[ $counts -ge $jobNumber ]]
    do
	sleep 10
	counts=$(ps | grep makeJet | wc -l)
    done
done

wait

rVals=(ak3 ak4 ak6 ak8 ak10 akCs3 akCs4 akCs6 akCs8 akCs10)
for r in "${rVals[@]}"
do
    files=""

    count=$(ls output/$DATE/*$r*SPLIT*.root | wc -l)

    if [[ $count -eq 0 ]]
    then
	continue
    fi
    
    for i in output/$DATE/*$r*SPLIT*.root
    do
	files=$files$i,
    done
    
    ./bin/combineResponse.exe $files $r >& logs/$DATE/combineResponse_$r.log &
done

wait
	     
TIMEEND=`date +%H%M%S`

echo "Time start: $TIMESTART"
echo "Time end: $TIMEEND"

echo "bash/runSplitResponse.sh Complete!"
