#!/bin/bash

mkdir -p logs

./bin/makeDeriveFlatResponse.exe paths/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM.txt >& logs/flatPP.log &

./bin/makeDeriveFlatResponse.exe paths/Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbPb_MCDijet_20180521_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180608_SVM.txt >& logs/flat.log &