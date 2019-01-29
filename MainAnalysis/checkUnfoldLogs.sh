#!/bin/bash


for i in logs/unfold*.log
do
    cat $i
    echo "Printed $i"
    sleep 3
done
