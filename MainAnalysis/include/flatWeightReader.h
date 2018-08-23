#ifndef FLATWEIGHTREADER_H
#define FLATWEIGHTREADER_H

//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/stringUtil.h"

class flatWeightReader{
 public:
  flatWeightReader(std::string inWeightFileName, cutPropagator inCutProp);
  bool Init(std::string inWeightFileName, cutPropagator inCutProp);
  void Clean();
  void HistSetNull();
  double getJtWeight(std::string algo, int cent, double jtpt, double jteta);


  std::string weightFileName = "";
  cutPropagator cutPropIn;
  cutPropagator cutPropWeight;
  TFile* weightFile_p=NULL;
  
  static const Int_t nMaxCentBins = 10;
  static const Int_t nMaxJtAbsEtaBins = 6;
  static const Int_t nMaxJtAlgos = 20;

  TH1D* genJtPt_h[nMaxJtAlgos][nMaxCentBins][nMaxJtAbsEtaBins];

  Bool_t isPP = true;

  Int_t nCentBins = -1;
  std::vector<Int_t> centBinsLow;
  std::vector<Int_t> centBinsHi;

  Int_t nJtAlgos = -1;
  std::vector<std::string> jtAlgos;

  Int_t nJtAbsEtaBins = -1;
  std::vector<Double_t> jtAbsEtaBinsLow;
  std::vector<Double_t> jtAbsEtaBinsHi;
};

flatWeightReader::flatWeightReader(std::string inWeightFileName, cutPropagator inCutProp)
{
  HistSetNull();
  Clean();
  Init(inWeightFileName, inCutProp);
  return;
}

bool flatWeightReader::Init(std::string inWeightFileName, cutPropagator inCutProp)
{
  HistSetNull();

  if(inWeightFileName.find(".root") == std::string::npos || !checkFile(inWeightFileName)){
    std::cout << "flatWeightReader: WARNING, weightFileName \'" << inWeightFileName << "\' is not a valid root file. return false" << std::endl;
    return false;
  }
  
  weightFileName = inWeightFileName;
  weightFile_p = new TFile(weightFileName.c_str(), "READ");
  cutPropWeight.Clean();
  cutPropWeight.GetAllVarFromFile(weightFile_p);
  cutPropIn = inCutProp;

  bool cutsMatch = true;
  if(!cutPropIn.CheckJtAbsEtaMax(cutPropWeight)) cutsMatch = false;
  else if(!cutPropIn.CheckNJtAbsEtaBins(cutPropWeight)) cutsMatch = false;
  else if(!cutPropIn.CheckJtAbsEtaBinsLow(cutPropWeight)) cutsMatch = false;
  else if(!cutPropIn.CheckJtAbsEtaBinsHi(cutPropWeight)) cutsMatch = false;
  else if(!cutPropIn.CheckNPthats(cutPropWeight)) cutsMatch = false;
  else if(!cutPropIn.CheckPthats(cutPropWeight)) cutsMatch = false; 
  else if(!cutPropIn.CheckPthatWeights(cutPropWeight)) cutsMatch = false;
  //  else if(!cutPropIn.CheckNJtAlgos(cutPropWeight)) cutsMatch = false;
  //  else if(!cutPropIn.CheckJtAlgos(cutPropWeight)) cutsMatch = false;
  else if(!cutPropIn.CheckIsPP(cutPropWeight)) cutsMatch = false;
  else if(!cutPropIn.CheckNCentBins(cutPropWeight)) cutsMatch = false;
  else if(!cutPropIn.CheckCentBinsLow(cutPropWeight)) cutsMatch = false;
  else if(!cutPropIn.CheckCentBinsHi(cutPropWeight)) cutsMatch = false;
 
  if(!cutsMatch){
    std::cout << "flatWeightReader: cutPropagator for weight file \'" << weightFileName << "\' does not match that of the input cutPropagator. return false" << std::endl;
    Clean();
    return false;
  }

  isPP = cutPropIn.GetIsPP();

  nCentBins = cutPropIn.GetNCentBins();
  centBinsLow = cutPropIn.GetCentBinsLow();
  centBinsHi = cutPropIn.GetCentBinsHi();

  nJtAlgos = cutPropIn.GetNJtAlgos();
  jtAlgos = cutPropIn.GetJtAlgos();

  nJtAbsEtaBins = cutPropIn.GetNJtAbsEtaBins();
  jtAbsEtaBinsLow = cutPropIn.GetJtAbsEtaBinsLow();
  jtAbsEtaBinsHi = cutPropIn.GetJtAbsEtaBinsHi();

  if(nCentBins > nMaxCentBins){
    std::cout << "nCentBins, " << nCentBins << ", from weightFile \'" << weightFileName << "\' is greater than the max allowed value, " << nMaxCentBins << ". return false" << std::endl;
    Clean();
    return false;
  }
 
  if(nJtAbsEtaBins > nMaxJtAbsEtaBins){
    std::cout << "nJtAbsEtaBins, " << nJtAbsEtaBins << ", from weightFile \'" << weightFileName << "\' is greater than the max allowed value, " << nMaxJtAbsEtaBins << ". return false" << std::endl;
    Clean();
    return false;
  }

  if(nJtAlgos > nMaxJtAlgos){
    std::cout << "nJtAlgos, " << nJtAlgos << ", from weightFile \'" << weightFileName << "\' is greater than the max allowed value, " << nMaxJtAlgos << ". return false" << std::endl;
    Clean();
    return false;
  }


  for(Int_t jI = 0; jI < nJtAlgos; ++jI){
    std::string jtName = jtAlgos.at(jI);
    jtName = jtName.substr(0, jtName.find("/"));

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
      if(isPP) centStr = "PP";

      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

	genJtPt_h[jI][cI][aI] = (TH1D*)weightFile_p->Get((jtName + "/genJtPt_" + jtName + "_" + centStr + "_" + jtAbsEtaStr + "_h").c_str());
      }
    }
  }


  return true;
}

void flatWeightReader::Clean()
{
  HistSetNull();

  weightFileName = "";
  cutPropIn.Clean();
  cutPropWeight.Clean();
  if(weightFile_p!=NULL){
    weightFile_p->Close();
    delete weightFile_p;
    weightFile_p = NULL;
  }

  isPP = true;

  nCentBins = -1;
  centBinsLow.clear();
  centBinsHi.clear();

  nJtAlgos = -1;
  jtAlgos.clear();

  nJtAbsEtaBins = -1;
  jtAbsEtaBinsLow.clear();
  jtAbsEtaBinsHi.clear();

  return;
}

void flatWeightReader::HistSetNull()
{
  for(Int_t jI = 0; jI < nMaxJtAlgos; ++jI){
    for(Int_t cI = 0; cI < nMaxCentBins; ++cI){
      for(Int_t aI = 0; aI < nMaxJtAbsEtaBins; ++aI){
	genJtPt_h[jI][cI][aI] = NULL;
      }
    }
  }

  return;
}


double flatWeightReader::getJtWeight(std::string algo, int cent, double jtpt, double jteta)
{
  int algPos = -1;
  int centPos = -1;
  int etaPos = -1;

  jteta = TMath::Abs(jteta);

  for(int jI = 0; jI < nJtAlgos; ++jI){
    if(isStrSame(algo, jtAlgos.at(jI))){
      algPos = jI;
      break;
    }    
  }

  if(algPos == -1){
    std::cout << "flatWeightReader: Warning algPos == -1 for \'" << algo << "\', return weight of 1" << std::endl;
    std::cout << " Options were: ";
    for(int jI = 0; jI < nJtAlgos; ++jI){
      std::cout << jtAlgos.at(jI) << ", ";
    }
    std::cout << std::endl;
    return 1.;
  }

  if(isPP) centPos = 0;
  else{
    for(int cI = 0; cI < nCentBins; ++cI){
      if(cent >= centBinsLow.at(cI) && cent < centBinsHi.at(cI)){
	centPos = cI;
	break;
      }    
    }
  }
    
  if(centPos == -1){
    std::cout << "flatWeightReader: Warning centPos == -1, return weight of 1" << std::endl;
    return 1.;
  }


  for(int aI = 0; aI < nJtAbsEtaBins; ++aI){
    if(jteta >= jtAbsEtaBinsLow.at(aI) && jteta < jtAbsEtaBinsHi.at(aI)){
      etaPos = aI;
      break;
    }    
  }

  if(etaPos == -1){
    if(jteta == jtAbsEtaBinsHi.at(jtAbsEtaBinsHi.size()-1)) etaPos = jtAbsEtaBinsHi.size()-2;
    else{
      std::cout << "flatWeightReader: Warning etaPos == -1 for eta, " << jteta << ", return weight of 1" << std::endl;
      return 1.;
    }
  }

  int binPos = -1;
  for(Int_t bI = 0; bI < genJtPt_h[algPos][centPos][etaPos]->GetNbinsX(); ++bI){
    if(jtpt >= genJtPt_h[algPos][centPos][etaPos]->GetBinLowEdge(bI+1) && jtpt < genJtPt_h[algPos][centPos][etaPos]->GetBinLowEdge(bI+2)){
      binPos = bI;
      break;
    }
  }

  if(binPos == -1){
    if(jtpt < genJtPt_h[algPos][centPos][etaPos]->GetBinLowEdge(1)) binPos = 1;
    else if(jtpt >= genJtPt_h[algPos][centPos][etaPos]->GetBinLowEdge(genJtPt_h[algPos][centPos][etaPos]->GetNbinsX()+1)) binPos = genJtPt_h[algPos][centPos][etaPos]->GetNbinsX();
    else{
      std::cout << "flatWeightReader: Warning binPos == -1, return weight of 1" << std::endl;
      return 1.; 
    }
  }

  double weight = 1./(Double_t)genJtPt_h[algPos][centPos][etaPos]->GetBinContent(binPos+1);
  return weight;
}

#endif
