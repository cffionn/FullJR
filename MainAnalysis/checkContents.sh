#!/bin/bash

files=(output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_AllAlgos_FracNEntries1p00_JetResponse_20180823.root output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_FracNEntries1p00_JetResponse_20180823.root)

fileStr=""
for i in "${files[@]}"
do
    echo $i
#    ./bin/checkFileNClassContents.exe $i
    fileStr="$fileStr $i"
done

./bin/checkFileNClassContents.exe $fileStr