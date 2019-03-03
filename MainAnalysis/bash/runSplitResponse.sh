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


DATE=`date +%Y%m%d`

mkdir -p logs
mkdir -p logs/$DATE

counts=0
for i in paths/*/*.txt
do
    counts=$((counts + 1))
done

if [[ $counts -ne 60 ]]
then
    echo "$counts should equal 60. exit 1"
    exit 1
fi

for i in paths/*/*.txt
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

    ./bin/makeJetResponseTree.exe $i $pbpbBool 1.0 >& logs/$DATE/$logName.log &

    counts=$(ps | grep makeJet | wc -l)
    while [[ $counts -ge 16 ]]
    do
	sleep 10
	counts=$(ps | grep makeJet | wc -l)
    done
done

wait
echo "bash/runSplitResponse.sh Complete!"
