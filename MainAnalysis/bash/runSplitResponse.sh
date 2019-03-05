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

jobNumber=16
if [[ -d /home/cfmcginn ]]
then
    jobNumber=10
fi

DATE=`date +%Y%m%d`
TIMESTART=`date +%H%M%S`

mkdir -p logs
mkdir -p logs/$DATE

filesToProcess=()
counts=$(ls paths/*SPLIT*.txt | wc -l)
if [[ $counts -ge 1 ]]
then
    for i in paths/*SPLIT*.txt
    do
	filesToProcess+=($i)
    done
fi    

counts=$(ls paths/*/*SPLIT*.txt | wc -l)
if [[ $counts -ge 1 ]]
then
    for i in paths/*/*SPLIT*.txt
    do
	filesToProcess+=($i)
    done
fi


for i in "${filesToProcess[@]}"
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

    echo "./bin/makeJetResponseTree.exe $i $pbpbBool 1.0 >& logs/$DATE/response_$logName.log &"
    exit 1

    counts=$(ps | grep makeJet | wc -l)
    while [[ $counts -ge $jobNumber ]]
    do
	sleep 10
	counts=$(ps | grep makeJet | wc -l)
    done
done

wait

DATEEND=`date +%Y%m%d`

rVals=(ak3 ak4 ak6 ak8 ak10 akCs3 akCs4 akCs6 akCs8 akCs10)
for r in "${rVals[@]}"
do
    files=""

    count=$(ls output/$DATE/*$r*SPLIT*.root | wc -l)
    if [[ -d output/$DATEEND ]]
    then
	counts2=$(ls output/$DATEEND/*$r*SPLIT*.root | wc -l)
	counts=$((counts + $counts2))
    fi


    if [[ $count -eq 0 ]]
    then
	continue
    fi
    
    for i in output/$DATE/*$r*SPLIT*.root
    do
	files=$files$i,
    done
    
    if [[ -d output/$DATEEND ]]
    then
	for i in output/$DATEEND/*$r*SPLIT*.root
	do
	    files=$files$i,
	done	
    fi

    ./bin/combineResponse.exe $files $r >& logs/$DATE/combineResponse_$r.log &
    jobCounts=$(ps | grep combineR | wc -l)
    if [[ $jobCounts -ge 5 ]]
    then
	sleep 10
	jobCounts=$(ps | grep combineR | wc -l)
    fi
done

wait
	     
TIMEEND=`date +%H%M%S`

echo "Time start: $TIMESTART"
echo "Time end: $TIMEEND"

echo "WARNING: DATE START \'$DATE\' not same as DATE END \'$DATEEND\'"

echo "bash/runSplitResponse.sh Complete!"
