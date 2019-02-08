//cpp dependencies
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT dependencies
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/smallOrLargeR.h"
#include "MainAnalysis/include/systFunctions.h"

//Non-Local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/vanGoghPalette.h"

template <class T>
bool compWithWarning(const std::string typeStr, T a, T b)
{
  bool val = a == b;
  if(!val) std::cout << typeStr << " \'" << a << "\' from one file does not match \'" << b << "\' from other. fail";
  return val;
}

ULong64_t getKey(ULong64_t jI, ULong64_t cI, ULong64_t idI, ULong64_t rI, ULong64_t aI, ULong64_t sI){return jI + 100*cI + 10000*idI + 1000000*rI + 100000000*aI + 10000000000*sI;}

ULong64_t getKey(ULong64_t jI, ULong64_t cI, ULong64_t idI, ULong64_t rI, ULong64_t aI, ULong64_t sI, ULong64_t bI){return jI + 100*cI + 10000*idI + 1000000*rI + 100000000*aI + 10000000000*sI + 10000000000000*bI;}

int posInStrVect(std::string strToCheck, std::string front, std::vector<std::string> vect, std::string back)
{
  int pos = -1;
  for(unsigned i = 0; i < vect.size(); ++i){
    if(vect[i].size() == 0) continue;
    if(strToCheck.find(front + vect[i] + back) != std::string::npos){
      pos = i;
      break;
    }
  }

  return pos;
}

int posInStrVectExact(std::string inStr, std::vector<std::string> inVect)
{
  Int_t pos = -1;
  for(unsigned int i = 0; i < inVect.size(); ++i){
    if(isStrSame(inStr, inVect[i])){pos = i; break;}
  }
  return pos;
}



void defineCanv(TCanvas* canv_p)
{
  canv_p->SetTopMargin(0.08);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.18);
  canv_p->SetBottomMargin(0.1);

  return;
}

template <class T>
bool setGeneralPtBins(const Int_t nGeneralPtBins, Double_t generalPtBins[], std::vector<T> binVect)
{
  if(nGeneralPtBins+1 < (Int_t)binVect.size()){
    std::cout << "testGeneralPtBins input binVect size \'" << binVect.size() << "\' is greater than the generalBins max+1, \'" << nGeneralPtBins << "+1\'. return false"  << std::endl;
    return false;
  }
  for(unsigned int vI = 0; vI < binVect.size(); ++vI){generalPtBins[vI] = binVect[vI];}
  return true;
}


bool macroHistToSubsetHist(TH1D* macroHist_p, TH1D* subsetHist_p)
{
  if(macroHist_p->GetBinLowEdge(1) - subsetHist_p->GetBinLowEdge(1) > 1){
    std::cout << "macroHistToSubsetHist Warning: macroHist low edge \'" << macroHist_p->GetBinLowEdge(1) << "\' is greater than low edge of subset hist, \'" << subsetHist_p->GetBinLowEdge(1) << "\', return false" << std::endl;
    return false;
  }
  if(macroHist_p->GetBinLowEdge(macroHist_p->GetNbinsX()+1) - subsetHist_p->GetBinLowEdge(subsetHist_p->GetNbinsX()+1) < -1){
    std::cout << "macroHistToSubsetHist Warning: macroHist high edge \'" << macroHist_p->GetBinLowEdge(1) << "\' is less than high edge of subset hist, \'" << subsetHist_p->GetBinLowEdge(1) << "\', return false" << std::endl;
    return false;
  }

  for(Int_t bIX = 0; bIX < subsetHist_p->GetNbinsX(); ++bIX){
    Double_t center = subsetHist_p->GetBinCenter(bIX+1);
    
    Int_t pos2 = -1;
    for(Int_t bIX2 = 0; bIX2 < macroHist_p->GetNbinsX(); ++bIX2){
      if(center >= macroHist_p->GetBinLowEdge(bIX2+1) && center < macroHist_p->GetBinLowEdge(bIX2+2)){
	pos2 = bIX2;
	break;
      }
    }
    
    subsetHist_p->SetBinContent(bIX+1, macroHist_p->GetBinContent(pos2+1));
    subsetHist_p->SetBinError(bIX+1, macroHist_p->GetBinError(pos2+1));
  }  
  return true;
}

void divHistByWidth(TH1D* hist_p)
{
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    Double_t binWidth = hist_p->GetBinWidth(bIX+1);
    hist_p->SetBinContent(bIX+1, hist_p->GetBinContent(bIX+1)/binWidth);
    hist_p->SetBinError(bIX+1, hist_p->GetBinError(bIX+1)/binWidth);
  }
  return;
}

void scaleHist(TH1D* hist_p, Double_t scaleFactor)
{
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    hist_p->SetBinContent(bIX+1, hist_p->GetBinContent(bIX+1)*scaleFactor);
    hist_p->SetBinError(bIX+1, hist_p->GetBinError(bIX+1)*scaleFactor);
  }
  return;
}

void scaleVect(std::vector<Double_t>* vect, Double_t scaleFactor)
{
  for(ULong64_t vI = 0; vI < (ULong64_t)vect->size(); ++vI){    
    (*vect)[vI] *= scaleFactor;
  }
  return;
}

Double_t getHistMax(TH1D* hist_p)
{
  Double_t max = hist_p->GetMinimum();
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    if(max < hist_p->GetBinContent(bIX+1)) max = hist_p->GetBinContent(bIX+1);
  }
  return max;
}

Double_t getHistMin(TH1D* hist_p)
{
  Double_t min = hist_p->GetMaximum();
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    if(min > hist_p->GetBinContent(bIX+1)) min = hist_p->GetBinContent(bIX+1);
  }
  return min;
}

Double_t getHistMinGTZero(TH1D* hist_p)
{
  Double_t min = hist_p->GetMaximum();
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    if(hist_p->GetBinContent(bIX+1) <= 0) continue;
    if(min > hist_p->GetBinContent(bIX+1)) min = hist_p->GetBinContent(bIX+1);
  }
  return min;
}


void createRAA(TH1D* pbpbHist_p, TH1D* ppHist_p)
{
  for(Int_t bIX = 0; bIX < pbpbHist_p->GetNbinsX(); ++bIX){
    Double_t val = pbpbHist_p->GetBinContent(bIX+1)/ppHist_p->GetBinContent(bIX+1);
    Double_t relErr1 = pbpbHist_p->GetBinError(bIX+1)/pbpbHist_p->GetBinContent(bIX+1);
    Double_t relErr2 = ppHist_p->GetBinError(bIX+1)/ppHist_p->GetBinContent(bIX+1);
    Double_t relErr = TMath::Sqrt(relErr1*relErr1 + relErr2*relErr2);
    Double_t err = val*relErr;

    pbpbHist_p->SetBinContent(bIX+1, val);
    pbpbHist_p->SetBinError(bIX+1, err);
  }
  
  return;
}


int plotUnfoldedAll(const std::string inFileNamePP, const std::string inFileNamePbPb)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);

  gStyle->SetOptStat(0);

  kirchnerPalette kPalette;

  const Double_t lumiFactor = getLumiFactor();
  const Double_t nMBEvents = getNMBEvents();

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  std::string outFileName = "output/" + dateStr + "/plotUnfoldedAll_";
  if(checkFile(outFileName + dateStr + ".root")) outFileName = outFileName + "UPDATED_" + dateStr + ".root";
  else outFileName = outFileName + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE"); 
  
  const Int_t nFiles = 2;
  TFile* inFile_p[nFiles];
  std::vector<bool> isPbPb = {true, false};
  std::vector<std::string> fileNames = {inFileNamePbPb, inFileNamePP};
  std::vector<cutPropagator*> cutProps;

  std::vector<std::string> jtAlgosPbPb;
  std::vector<std::string> jtAlgosPP;
  std::vector<std::string> dirListPbPb;
  std::vector<std::string> dirListPP;

  std::vector<std::string> histTagsPbPb;
  std::vector<std::string> histTagsPP;
  std::vector<int> histBestBayesPbPb;
  std::vector<int> histBestBayesPP;
  std::map<std::string, int> histTagsToBestBayesPbPb;
  std::map<std::string, int> histTagsToBestBayesPP;
  std::map<ULong64_t, ULong64_t> histKeyToBestBayesKeyPbPb;
  std::map<ULong64_t, ULong64_t> histKeyToBestBayesKeyPP;

  Int_t nCentBins = -1;
  std::vector<Int_t> centBinsLow;
  std::vector<Int_t> centBinsHi;
  std::vector<std::string> centBinsStr;
  std::vector<Double_t> centBinsWidth;

  Int_t nID = -1;
  std::vector<std::string> idStr;
  std::vector<bool> goodID;

  Int_t nResponseMod = -1;
  std::vector<double> responseMod;
  std::vector<std::string> responseModStr;
  std::vector<bool> goodResponseMod;

  Int_t nJtAbsEtaBins = -1;  
  std::vector<double> jtAbsEtaBinsLow;
  std::vector<double> jtAbsEtaBinsHi;
  std::vector<std::string> jtAbsEtaBinsStr;
  std::vector<bool> goodJtAbsEtaBins;
  std::vector<double> jtAbsEtaBinsWidth;

  Int_t nSystInFile = -1;
  std::vector<std::string> systStrInFile;
  //  Int_t nSystOutFile = 4;
  std::vector<std::string> systStrOutFile = {"LumiUp", "LumiDown", "TAAUp", "TAADown"};

  std::map<std::string, std::string> systToCombo = {{"JECUpMC", "JECMC"}, {"JECDownMC", "JECMC"}, {"JECUpData", "JECData"}, {"JECDownData", "JECData"}, {"JECUpUE", "JECUE"}, {"JECDownUE", "JECUE"}, {"PriorUp1PowerPthat", "Prior1PowerPthat"}, {"PriorDown1PowerPthat", "Prior1PowerPthat"}, {"LumiUp", "Lumi"}, {"LumiDown", "Lumi"}, {"TAAUp", "TAA"}, {"TAADown", "TAA"}};
  std::vector<std::string> reducedSystStr;

  std::vector<std::string> corrSystStr = {"JECUpMC", "JECDownMC", "JECUpData", "JECDownData", "JECUpUE", "JECDownUE", "PriorUp1PowerPthat", "PriorDown1PowerPthat"};
  
  Int_t nBayes = -1;

  std::vector<int> bayesVal;
  std::vector<std::string> bayesValStr;

  //Use this to do any or all pt binnings
  const Int_t nGeneralPtBins = 100;
  Double_t generalPtBins[nGeneralPtBins+1];

  Int_t nGenJtPtBinsSmallR;
  Int_t nGenJtPtBinsLargeR;
  std::vector<Double_t> genJtPtBinsSmallRTemp;
  std::vector<Double_t> genJtPtBinsLargeRTemp;
  std::vector<Double_t> genJtPtBinsSmallR;
  std::vector<Double_t> genJtPtBinsLargeR;
  Int_t nRecoJtPtBinsSmallR;
  Int_t nRecoJtPtBinsLargeR;
  std::vector<Double_t> recoJtPtBinsSmallRTemp;
  std::vector<Double_t> recoJtPtBinsLargeRTemp;

  std::cout << "Processing inputs..." << std::endl;
  unsigned int pos = 0;
  for(auto const & fName : fileNames){
    inFile_p[pos] = new TFile(fName.c_str(), "READ");

    cutProps.push_back(new cutPropagator());
    cutProps[pos]->Clean();
    cutProps[pos]->GetAllVarFromFile(inFile_p[pos]);        
    
    if(isPbPb[pos]){
      dirListPbPb = returnRootFileContentsList(inFile_p[pos], "TDirectoryFile", "JetAnalyzer", 1);
      histTagsPbPb = cutProps[pos]->GetHistTag();
      histBestBayesPbPb = cutProps[pos]->GetHistBestBayes();

      nCentBins = cutProps[pos]->GetNCentBins();
      centBinsLow = cutProps[pos]->GetCentBinsLow();
      centBinsHi = cutProps[pos]->GetCentBinsHi();
      jtAlgosPbPb = cutProps[pos]->GetJtAlgos();
    }      
    else{
      dirListPP = returnRootFileContentsList(inFile_p[pos], "TDirectoryFile", "JetAnalyzer", 1);
      histTagsPP = cutProps[pos]->GetHistTag();
      histBestBayesPP = cutProps[pos]->GetHistBestBayes();

      jtAlgosPP = cutProps[pos]->GetJtAlgos();      
    }
    
    if(pos == 0){
      nID = cutProps[pos]->GetNID();
      nResponseMod = cutProps[pos]->GetNResponseMod();
      nSystInFile = cutProps[pos]->GetNSyst();
      nJtAbsEtaBins = cutProps[pos]->GetNJtAbsEtaBins();
      nBayes = cutProps[pos]->GetNBayes();

      idStr = cutProps[pos]->GetIdStr();
      responseMod = cutProps[pos]->GetResponseMod();
      jtAbsEtaBinsLow = cutProps[pos]->GetJtAbsEtaBinsLow();
      jtAbsEtaBinsHi = cutProps[pos]->GetJtAbsEtaBinsHi();
      systStrInFile = cutProps[pos]->GetSystStr();
      bayesVal = cutProps[pos]->GetBayesVal();

      nGenJtPtBinsSmallR = cutProps[pos]->GetNGenJtPtBinsSmallR();
      nGenJtPtBinsLargeR = cutProps[pos]->GetNGenJtPtBinsLargeR();
      genJtPtBinsSmallRTemp = cutProps[pos]->GetGenJtPtBinsSmallR();
      genJtPtBinsLargeRTemp = cutProps[pos]->GetGenJtPtBinsLargeR();
      nRecoJtPtBinsSmallR = cutProps[pos]->GetNRecoJtPtBinsSmallR();
      nRecoJtPtBinsLargeR = cutProps[pos]->GetNRecoJtPtBinsLargeR();
      recoJtPtBinsSmallRTemp = cutProps[pos]->GetRecoJtPtBinsSmallR();
      recoJtPtBinsLargeRTemp = cutProps[pos]->GetRecoJtPtBinsLargeR();
    }
    else{
      if(!compWithWarning("nID", nID, cutProps[pos]->GetNID())) return 1;
      if(!compWithWarning("nResponseMod", nResponseMod, cutProps[pos]->GetNResponseMod())) return 1;
      if(!compWithWarning("nSystInFile", nSystInFile, cutProps[pos]->GetNSyst())) return 1;
      if(!compWithWarning("nJtAbsEtaBins", nJtAbsEtaBins, cutProps[pos]->GetNJtAbsEtaBins())) return 1;

      bool propMatch = cutProps[0]->CheckPropagatorsMatch(*(cutProps[pos]), true, false, false);

      if(!propMatch) std::cout << " Mismatch in propagator 0 to " << pos << std::endl;
      else std::cout << " Good propagator " << pos << std::endl;
    }

    inFile_p[pos]->Close();
    delete inFile_p[pos];
    ++pos;
  }

  smallOrLargeR rReader;
  std::cout << "Checking binning against rReader.." << std::endl;
  if(!rReader.CheckNGenJtPtBinsSmallR(nGenJtPtBinsSmallR)) return 1;
  else if(!rReader.CheckNGenJtPtBinsLargeR(nGenJtPtBinsLargeR)) return 1;
  else if(!rReader.CheckGenJtPtBinsSmallR(genJtPtBinsSmallRTemp)) return 1;
  else if(!rReader.CheckGenJtPtBinsLargeR(genJtPtBinsLargeRTemp)) return 1;
  else if(!rReader.CheckNRecoJtPtBinsSmallR(nRecoJtPtBinsSmallR)) return 1;
  else if(!rReader.CheckNRecoJtPtBinsLargeR(nRecoJtPtBinsLargeR)) return 1;
  else if(!rReader.CheckRecoJtPtBinsSmallR(recoJtPtBinsSmallRTemp)) return 1;
  else if(!rReader.CheckRecoJtPtBinsLargeR(recoJtPtBinsLargeRTemp)) return 1;
  else std::cout << " All bins good!" << std::endl;

  const Int_t nXVals = 5;
  const Int_t xVals[nXVals] = {200, 400, 600, 800, 1000};
  
  const Double_t minValSmallR = 200;
  const Double_t minValLargeR = 300;
  const Double_t maxValSmallR = 1000;
  const Double_t maxValLargeR = 1000;
  bool minSmallRFound = false;
  bool minLargeRFound = false;
  bool maxSmallRFound = false;
  bool maxLargeRFound = false;
  
  for(Int_t gI = 0; gI < nGenJtPtBinsSmallR+1; ++gI){
    if(TMath::Abs(genJtPtBinsSmallRTemp[gI] - minValSmallR) < 1.) minSmallRFound = true;
    if(TMath::Abs(genJtPtBinsSmallRTemp[gI] - maxValSmallR) < 1.) maxSmallRFound = true;

    if(minSmallRFound) genJtPtBinsSmallR.push_back(genJtPtBinsSmallRTemp[gI]);
    if(minSmallRFound && maxSmallRFound) break;
  }

  for(Int_t gI = 0; gI < nGenJtPtBinsLargeR+1; ++gI){
    if(TMath::Abs(genJtPtBinsLargeRTemp[gI] - minValLargeR) < 1.) minLargeRFound = true;
    if(TMath::Abs(genJtPtBinsLargeRTemp[gI] - maxValLargeR) < 1.) maxLargeRFound = true;

    if(minLargeRFound) genJtPtBinsLargeR.push_back(genJtPtBinsLargeRTemp[gI]);
    if(minLargeRFound && maxLargeRFound) break;
  }

  if(!minSmallRFound){
    std::cout << "Min smallR val \'" << minValSmallR << "\' is not found in bins: ";
    for(Int_t gI = 0; gI < nGenJtPtBinsSmallR; ++gI){
      std::cout << genJtPtBinsSmallRTemp[gI] << ", ";
    }
    std::cout << genJtPtBinsSmallRTemp[nGenJtPtBinsSmallR] << ". return 1" << std::endl;
    return 1;
  }
  if(!maxSmallRFound){
    std::cout << "Max smallR val \'" << maxValSmallR << "\' is not found in bins: ";
    for(Int_t gI = 0; gI < nGenJtPtBinsSmallR; ++gI){
      std::cout << genJtPtBinsSmallRTemp[gI] << ", ";
    }
    std::cout << genJtPtBinsSmallRTemp[nGenJtPtBinsSmallR] << ". return 1" << std::endl;
    return 1;
  }

  if(!minLargeRFound){
    std::cout << "Min largeR val \'" << minValLargeR << "\' is not found in bins: ";
    for(Int_t gI = 0; gI < nGenJtPtBinsLargeR; ++gI){
      std::cout << genJtPtBinsLargeRTemp[gI] << ", ";
    }
    std::cout << genJtPtBinsLargeRTemp[nGenJtPtBinsLargeR] << ". return 1" << std::endl;
    return 1;
  }
  if(!maxLargeRFound){
    std::cout << "Max largeR val \'" << maxValLargeR << "\' is not found in bins: ";
    for(Int_t gI = 0; gI < nGenJtPtBinsLargeR; ++gI){
      std::cout << genJtPtBinsLargeRTemp[gI] << ", ";
    }
    std::cout << genJtPtBinsLargeRTemp[nGenJtPtBinsLargeR] << ". return 1" << std::endl;
    return 1;
  }

  for(unsigned int sI = 0; sI < systStrInFile.size(); ++sI){
    if(systToCombo.count(systStrInFile[sI]) > 0){
      std::string systTemp = systToCombo[systStrInFile[sI]];
      bool stringIsFound = false;
      for(auto const & sI2 : reducedSystStr){
	if(isStrSame(sI2, systTemp)) stringIsFound = true;
      }
      if(!stringIsFound) reducedSystStr.push_back(systToCombo[systStrInFile[sI]]);
    }
    else reducedSystStr.push_back(systStrInFile[sI]);
  }
  Int_t nReducedSyst = reducedSystStr.size();

  std::cout << "Full systematics: " << systStrInFile.size() << std::endl;
  for(unsigned int sI = 1; sI < systStrInFile.size(); ++sI){
    std::cout << " " << systStrInFile[sI] << std::endl;
  }

  std::cout << "Reduced systematics: " << reducedSystStr.size() << std::endl;
  for(unsigned int sI = 1; sI < reducedSystStr.size(); ++sI){
    std::cout << " " << reducedSystStr[sI] << std::endl;
  }

  for(unsigned int hI = 0; hI < histTagsPbPb.size(); ++hI){
    if(histTagsPbPb[hI].find("nHistDim") != std::string::npos) continue;
    histTagsToBestBayesPbPb[histTagsPbPb[hI]] = histBestBayesPbPb[hI];
  }
  for(unsigned int hI = 0; hI < histTagsPP.size(); ++hI){
    if(histTagsPP[hI].find("nHistDim") != std::string::npos) continue;
    histTagsToBestBayesPP[histTagsPP[hI]] = histBestBayesPP[hI];
  }

  for(auto const & res : responseMod){responseModStr.push_back("ResponseMod" + prettyString(res, 2, true));}

  for(unsigned int aI = 0; aI < jtAbsEtaBinsLow.size(); ++aI){
    jtAbsEtaBinsStr.push_back("AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true));
  }

  for(auto const & b : bayesVal){bayesValStr.push_back("Bayes" + std::to_string(b));}

  for(unsigned int cI = 0; cI < centBinsLow.size(); ++cI){
    centBinsStr.push_back("Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]));
    centBinsWidth.push_back(centBinsHi[cI] - centBinsLow[cI]);
  }

  for(unsigned int aI = 0; aI < jtAbsEtaBinsLow.size(); ++aI){
    jtAbsEtaBinsStr.push_back("JtAbsEta" + std::to_string(jtAbsEtaBinsLow[aI]) + "to" + std::to_string(jtAbsEtaBinsHi[aI]));
    jtAbsEtaBinsWidth.push_back(jtAbsEtaBinsHi[aI] - jtAbsEtaBinsLow[aI]);
  }

  std::cout << "nCentBins: " << nCentBins << std::endl;
  std::cout << "nBayes: " << nBayes << std::endl;

  std::cout << "PbPb Algos..." << std::endl;
  for(auto const & algo : jtAlgosPbPb){
    std::cout << " " << algo << std::endl;
  }
  std::cout << "PbPb Dir..." << std::endl;
  for(auto const & dir : dirListPbPb){
    std::cout << " " << dir << std::endl;
  }

  std::cout << "PP Algos..." << std::endl;
  for(auto const & algo : jtAlgosPP){
    std::cout << " " << algo << std::endl;
  }
  std::cout << "PP Dir..." << std::endl;
  for(auto const & dir : dirListPP){
    std::cout << " " << dir << std::endl;
  }

  const UInt_t nJtAlgos = jtAlgosPbPb.size();
  //nJtAlgos + nCentBins*100 + nMaxID*10000 + nMaxResponseMod*1000000 + nMaxJtAbsEtaBins*100000000 + nMaxSyst*10000000000

  std::map<ULong64_t, ULong64_t> keyToVectPosPbPb;
  std::map<ULong64_t, ULong64_t> keyPbPbToKeyPP;
  UInt_t nKeyPbPb = 0;
  std::map<ULong64_t, ULong64_t> keyToVectPosPP;
  UInt_t nKeyPP = 0;
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
      goodID.push_back(false);

      for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	goodResponseMod.push_back(false);

	for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	  goodJtAbsEtaBins.push_back(false);

	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    for(ULong64_t bI = 0; bI < (ULong64_t)nBayes; ++bI){

	      ULong64_t keyPP = getKey(jI, 0, idI, rI, aI, sI, bI);
	      keyToVectPosPP[keyPP] = nKeyPP;
	      ++nKeyPP;

	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		ULong64_t keyPbPb = getKey(jI, cI, idI, rI, aI, sI, bI);
		keyToVectPosPbPb[keyPbPb] = nKeyPbPb;
		keyPbPbToKeyPP[keyPbPb] = keyPP;
		++nKeyPbPb;
	      }	      
	    }
	  }
	}
      }
    }
  }
  std::cout << "nKeys: " << nKeyPbPb << std::endl;
  std::vector<TH1D*> histVectPbPb;
  std::vector<TH1D*> histVectPP;

  histVectPbPb.reserve(nKeyPbPb);
  for(ULong64_t kI = 0; kI < nKeyPbPb; ++kI){histVectPbPb.push_back(NULL);}

  histVectPP.reserve(nKeyPP);
  for(ULong64_t kI = 0; kI < nKeyPP; ++kI){histVectPP.push_back(NULL);}

  std::cout << "Getting all hists...." << std::endl;
  for(Int_t fI = 0; fI < nFiles; ++fI){
    std::vector<std::string> dirs;
    if(isPbPb[fI]) dirs = dirListPbPb;
    else dirs = dirListPP;

    inFile_p[fI] = new TFile(fileNames[fI].c_str(), "READ");

    for(auto const & dir : dirs){
      TDirectory* dir_p = (TDirectory*)inFile_p[fI]->Get(dir.c_str());
      TIter next(dir_p->GetListOfKeys());
      TKey* key = NULL;
    
      int algoPos = -1;
      std::vector<std::string> algos;
      if(isPbPb[fI]) algos = jtAlgosPbPb;
      else algos = jtAlgosPP;

      for(UInt_t algoI = 0; algoI < algos.size(); ++algoI){
	if(algos[algoI].find(dir) != std::string::npos){
	  algoPos = algoI;
	  break;
	}
      }
    
      while((key = (TKey*)next())){
	const std::string name = key->GetName();
	const std::string className = key->GetClassName();
	if(className.find("TH1") == std::string::npos) continue;
	
	int idPos = posInStrVect(name, "_", idStr, "_");
	int modPos = posInStrVect(name, "_", responseModStr, "_");
	int absEtaPos = posInStrVect(name, "_", jtAbsEtaBinsStr, "_");
	int systPos = posInStrVect(name, "_", systStrInFile, "_");
	if(systPos < 0 && name.find("FlatPrior") == std::string::npos) systPos = 0;
	int bayesPos = posInStrVect(name, "_", bayesValStr, "_");
	int centPos = 0;
	if(isPbPb[fI]) centPos = posInStrVect(name, "_", centBinsStr, "_");

	if(algoPos < 0) std::cout << "Missing algoPos in name \'" << name << "\'" << std::endl;
	if(idPos < 0) std::cout << "Missing idPos in name \'" << name << "\'" << std::endl;
	if(modPos < 0) std::cout << "Missing modPos in name \'" << name << "\'" << std::endl;
	if(absEtaPos < 0) std::cout << "Missing absEtaPos in name \'" << name << "\'" << std::endl;
	if(systPos < 0) std::cout << "Missing systPos in name \'" << name << "\'" << std::endl;
	if(bayesPos < 0) std::cout << "Missing bayesPos in name \'" << name << "\'" << std::endl;
	if(centPos < 0) std::cout << "Missing centPos in name \'" << name << "\'" << std::endl;	

	goodID[idPos] = true;
	goodResponseMod[modPos] = true;
	goodJtAbsEtaBins[absEtaPos] = true;

	ULong64_t idKey = getKey(algoPos, centPos, idPos, modPos, absEtaPos, systPos, bayesPos);
	if(isPbPb[fI]){
	  ULong64_t vectPos = keyToVectPosPbPb[idKey];	  
	  histVectPbPb[vectPos] = (TH1D*)key->ReadObj();
	}
	else{
	  ULong64_t vectPos = keyToVectPosPP[idKey];
	  histVectPP[vectPos] = (TH1D*)key->ReadObj();
	}
      }
    }
  }

  std::cout << "Processing best bayes PbPb..." << std::endl;
  for(auto const & tagToBB : histTagsToBestBayesPbPb){
    int algoPos = posInStrVect(tagToBB.first, "", dirListPbPb, "_");
    int idPos = posInStrVect(tagToBB.first, "_", idStr, "_");
    int modPos = posInStrVect(tagToBB.first, "_", responseModStr, "_");
    int absEtaPos = posInStrVect(tagToBB.first, "_", jtAbsEtaBinsStr, "_");
    int systPos = posInStrVect(tagToBB.first, "_", systStrInFile, "");
    if(systPos < 0 && tagToBB.first.find("FlatPrior") == std::string::npos) systPos = 0;
    int centPos = posInStrVect(tagToBB.first, "_", centBinsStr, "_");
    int bayesPos = tagToBB.second;
    
    if(algoPos < 0) std::cout << "Missing algoPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(idPos < 0) std::cout << "Missing idPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(modPos < 0) std::cout << "Missing modPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(absEtaPos < 0) std::cout << "Missing absEtaPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(systPos < 0) std::cout << "Missing systPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(centPos < 0) std::cout << "Missing centPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;	
    if(bayesPos < 0){
      std::cout << "Missing bayesPos in \'" << tagToBB.first << "\'. set to default 10" << std::endl;
      bayesPos = 10;
    }

    ULong64_t idKey = getKey(algoPos, centPos, idPos, modPos, absEtaPos, systPos);
    ULong64_t idKeyBB = getKey(algoPos, centPos, idPos, modPos, absEtaPos, systPos, (ULong64_t)bayesPos);
    histKeyToBestBayesKeyPbPb[idKey] = idKeyBB;
  }

  std::cout << "Processing best bayes PP..." << std::endl;
  for(auto const & tagToBB : histTagsToBestBayesPP){
    int algoPos = posInStrVect(tagToBB.first, "", dirListPP, "_");
    int idPos = posInStrVect(tagToBB.first, "_", idStr, "_");
    int modPos = posInStrVect(tagToBB.first, "_", responseModStr, "_");
    int absEtaPos = posInStrVect(tagToBB.first, "_", jtAbsEtaBinsStr, "_");
    int systPos = posInStrVect(tagToBB.first, "_", systStrInFile, "");
    if(systPos < 0 && tagToBB.first.find("FlatPrior") == std::string::npos) systPos = 0;

    int centPos = 0;
    int bayesPos = tagToBB.second;
    
    if(algoPos < 0) std::cout << "Missing algoPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(idPos < 0) std::cout << "Missing idPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(modPos < 0) std::cout << "Missing modPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(absEtaPos < 0) std::cout << "Missing absEtaPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(systPos < 0) std::cout << "Missing systPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(centPos < 0) std::cout << "Missing centPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;	
    if(bayesPos < 0){
      std::cout << "Missing bayesPos in \'" << tagToBB.first << "\'. set to default 10" << std::endl;
      bayesPos = 10;
    }

    ULong64_t idKey = getKey(algoPos, centPos, idPos, modPos, absEtaPos, systPos);
    ULong64_t idKeyBB = getKey(algoPos, centPos, idPos, modPos, absEtaPos, systPos, (ULong64_t)bayesPos);
    histKeyToBestBayesKeyPP[idKey] = idKeyBB;
  }
  
  //  const UInt_t nDummyStr = 256;
  //  char dummyStr[nDummyStr];
  //  std::cin.get(dummyStr, nDummyStr);

  std::cout << "Traditional Spectra" << std::endl;
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    std::cout << "Jet Algo: " << jtAlgosPbPb[jI] << "/" << jtAlgosPP[jI] << std::endl;
    const Int_t rVal = getRVal(dirListPbPb.at(jI));
    const std::string rValStr = getRValStr(dirListPbPb.at(jI));
    const bool isSmallR = rReader.GetIsSmallR(rVal);

    Int_t nBinsTemp = -1;
    if(isSmallR){
      nBinsTemp = genJtPtBinsSmallR.size()-1;
      setGeneralPtBins(nGeneralPtBins, generalPtBins, genJtPtBinsSmallR);
    }
    else{
      nBinsTemp = genJtPtBinsLargeR.size()-1;
      setGeneralPtBins(nGeneralPtBins, generalPtBins, genJtPtBinsLargeR);
    }
    const Int_t nBins = nBinsTemp;

    for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
      if(!goodID[idI]) continue;

      for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	if(!goodResponseMod[rI]) continue;

	for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	  if(!goodJtAbsEtaBins[aI]) continue;	 

	  TLatex* label_p = new TLatex();
	  //	  label_p->SetNDC();
	  label_p->SetTextFont(43);
	  label_p->SetTextSize(16);
	  
	  TLegend* leg_p = new TLegend(0.25, 0.12, 0.4, 0.35);
	  leg_p->SetBorderSize(0);
	  leg_p->SetFillColor(0);
	  leg_p->SetFillStyle(0);
	  leg_p->SetTextFont(43);
	  leg_p->SetTextSize(16);

	  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);	  
	  defineCanv(canv_p);
	  
	  std::vector<TH1D*> histVectPbPbClones;
	  histVectPbPbClones.reserve(nSystInFile*nCentBins);
	  std::vector<TH1D*> histVectPPClones;	  
	  histVectPPClones.reserve(nSystInFile);

	  //Grab the subset histograms
	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    ULong64_t idKey = getKey(jI, 0, idI, rI, aI, sI);
	    ULong64_t idKeyBB = histKeyToBestBayesKeyPP[idKey];
	    ULong64_t vectPos = keyToVectPosPP[idKeyBB];
	    std::string tempName = std::string(histVectPP[vectPos]->GetName()) + "_Clone";
	    histVectPPClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];#frac{1}{#LTTAA#GT} #frac{1}{N_{evt}} #frac{d^{2}N_{Jet}}{dp_{T}d#eta}  and  #frac{d^{2}#sigma_{Jet}}{dp_{T}d#eta} [nb/GeV]", nBins, generalPtBins));
	    if(!macroHistToSubsetHist(histVectPP[vectPos], histVectPPClones[sI])) return 1;

	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      idKey = getKey(jI, cI, idI, rI, aI, sI);
	      idKeyBB = histKeyToBestBayesKeyPbPb[idKey];
	      vectPos = keyToVectPosPbPb[idKeyBB];
	      tempName = std::string(histVectPbPb[vectPos]->GetName()) + "_Clone";
	      histVectPbPbClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];#frac{1}{#LTTAA#GT} #frac{1}{N_{evt}} #frac{d^{2}N_{Jet}}{dp_{T}d#eta}  and  #frac{d^{2}#sigma_{Jet}}{dp_{T}d#eta} [nb/GeV]", nBins, generalPtBins));
	      if(!macroHistToSubsetHist(histVectPbPb[vectPos], histVectPbPbClones[cI + sI*nCentBins])) return 1;
	    }
	  }

	  //Div by bin widths;
          for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    divHistByWidth(histVectPPClones[sI]);
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      divHistByWidth(histVectPbPbClones[cI + sI*nCentBins]);
	    }
	  }

	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){	  
	    scaleHist(histVectPPClones[sI], 1./(2.*jtAbsEtaBinsWidth[aI]));
	    scaleHist(histVectPPClones[sI], 1./(lumiFactor));
	  }

	  Double_t scaleFactor = 1.;
	  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	    Double_t tempTAAFactor = getTAAScaleFactor(centBinsStr[cI]);
	    scaleFactor *= 10.;
	    
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	      scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./(nMBEvents*tempTAAFactor));
	      scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./(2.*jtAbsEtaBinsWidth[aI]));
	      scaleHist(histVectPbPbClones[cI + sI*nCentBins], 100./centBinsWidth[cI]);
	      scaleHist(histVectPbPbClones[cI + sI*nCentBins], scaleFactor);
	    }
	  }

	  //Setup systematics
	  std::vector<std::vector<Double_t> > reducedSystPP;
	  std::vector<std::vector<Double_t> > reducedSystPbPb;
	  std::vector<Double_t> reducedSystSumPP;
	  std::vector<std::vector<Double_t> > reducedSystSumPbPb;
	  for(ULong64_t sI = 0; sI < (ULong64_t)nReducedSyst; ++sI){
	    reducedSystPP.push_back({});
	    for(Int_t bI = 0; bI < nBins; ++bI){
	      if(sI == 0) reducedSystSumPP.push_back(0.0);
	      reducedSystPP[sI].push_back(1000000000000.);
	    }

            for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      if(sI == 0) reducedSystSumPbPb.push_back({});
	      reducedSystPbPb.push_back({});

	      for(Int_t bI = 0; bI < nBins; ++bI){
		if(sI == 0) reducedSystSumPbPb[cI].push_back(0.0);
		reducedSystPbPb[cI + sI*nCentBins].push_back(1000000000000.);
	      }
	    }
	  }

	  //Extract systematics
	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    ULong64_t pos = 0;
	    if(systToCombo.count(systStrInFile[sI]) > 0) pos = posInStrVectExact(systToCombo[systStrInFile[sI]], reducedSystStr);
	    else pos = posInStrVectExact(systStrInFile[sI], reducedSystStr);

	    for(Int_t bI = 0; bI < nBins; ++bI){
	      Double_t delta = TMath::Abs(histVectPPClones[sI]->GetBinContent(bI+1) - histVectPPClones[0]->GetBinContent(bI+1));
	      (reducedSystPP[pos])[bI] = TMath::Min((reducedSystPP[pos])[bI], delta);
	    }

            for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      for(Int_t bI = 0; bI < nBins; ++bI){
		Double_t delta = TMath::Abs(histVectPbPbClones[cI + sI*nCentBins]->GetBinContent(bI+1) - histVectPbPbClones[cI]->GetBinContent(bI+1));

		(reducedSystPbPb[cI + pos*nCentBins])[bI] = TMath::Min((reducedSystPbPb[cI + pos*nCentBins])[bI], delta);
	      }
	    }
	  }

	  for(Int_t bI = 0; bI < nBins; ++bI){
	    for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
	      if(reducedSystStr[sI].find("PriorFlat") != std::string::npos) continue;
	      reducedSystSumPP[bI] = TMath::Sqrt(reducedSystSumPP[bI]*reducedSystSumPP[bI] + (reducedSystPP[sI])[bI]*(reducedSystPP[sI])[bI]);
	    }
	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
		if(reducedSystStr[sI].find("PriorFlat") != std::string::npos) continue;
		(reducedSystSumPbPb[cI])[bI] = TMath::Sqrt((reducedSystSumPbPb[cI])[bI]*(reducedSystSumPbPb[cI])[bI] + (reducedSystPbPb[cI + sI*nCentBins])[bI]*(reducedSystPbPb[cI + sI*nCentBins])[bI]);
	      }
	    }
	  }
	  
	  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	    std::cout << "CentBins: " << centBinsStr[cI] << std::endl;
	    for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
	      std::cout << " SYST CHECK: " << reducedSystStr[sI] << std::endl;
	      for(Int_t bI = 0; bI < nBins; ++bI){
		std::cout <<"  " << bI << ": " << (reducedSystPbPb[cI + nCentBins*sI])[bI] << std::endl;
	      }	    
	    }
	  }
	  

	  //Grab Max and Min vals...	  
	  Double_t min = getHistMinGTZero(histVectPPClones[0]);
	  Double_t max = getHistMax(histVectPPClones[0]);
	  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	    Double_t tempMin = getHistMinGTZero(histVectPbPbClones[cI]);
	    Double_t tempMax = getHistMax(histVectPbPbClones[cI]);

	    if(tempMin < min) min = tempMin;
	    if(tempMax > max) max = tempMax;
	  }

	  max = getNearestFactor10Up(max, 1);
	  min = getNearestFactor10Down(min, 2);

	  histVectPPClones[0]->SetMaximum(max);
	  histVectPPClones[0]->SetMinimum(min);
	  
	  canv_p->cd();
	  histVectPPClones[0]->SetMarkerColor(kPalette.getColor(getColorPosFromCent("", true)));
	  histVectPPClones[0]->SetLineColor(kPalette.getColor(getColorPosFromCent("", true)));
	  histVectPPClones[0]->SetMarkerStyle(getStyleFromCent("", true));
	  histVectPPClones[0]->SetMarkerSize(1);
	  centerTitles(histVectPPClones[0]);

	  histVectPPClones[0]->GetXaxis()->SetTitleOffset(histVectPPClones[0]->GetXaxis()->GetTitleOffset()*1.3);
	  histVectPPClones[0]->GetYaxis()->SetTitleOffset(2.2);
	  histVectPPClones[0]->GetXaxis()->SetLabelColor(0);
	  
	  histVectPPClones[0]->DrawCopy("HIST E1 P");	  
	  gPad->SetLogx();
	  
	  drawSyst(canv_p, histVectPPClones[0], reducedSystSumPP, histVectPPClones[0]->GetBinLowEdge(1), histVectPPClones[0]->GetBinLowEdge(histVectPPClones[0]->GetNbinsX()+1));
	  histVectPPClones[0]->DrawCopy("HIST E1 P SAME");	  
	  

	  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	    histVectPbPbClones[cI]->SetMarkerColor(kPalette.getColor(getColorPosFromCent(centBinsStr[cI], false)));
	    histVectPbPbClones[cI]->SetLineColor(kPalette.getColor(getColorPosFromCent(centBinsStr[cI], false)));
	    histVectPbPbClones[cI]->SetMarkerStyle(getStyleFromCent(centBinsStr[cI], false));
	    histVectPbPbClones[cI]->SetMarkerSize(1);
	    histVectPbPbClones[cI]->DrawCopy("HIST E1 P SAME");
	    drawSyst(canv_p, histVectPbPbClones[cI], reducedSystSumPbPb[cI], histVectPbPbClones[cI]->GetBinLowEdge(1), histVectPbPbClones[cI]->GetBinLowEdge(histVectPbPbClones[cI]->GetNbinsX()+1));
	    histVectPbPbClones[cI]->DrawCopy("HIST E1 P SAME");
	  }
	  
	  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	    ULong64_t pos = nCentBins - 1 - cI;
	    std::string centStr = std::to_string(centBinsLow[pos]) + "-" + std::to_string(centBinsHi[pos]) + "%";
	    leg_p->AddEntry(histVectPbPbClones[pos], (centStr + " #times 10^{" + std::to_string((Int_t)pos+1) + "}").c_str(), "F P L");
	  }

	  leg_p->AddEntry(histVectPPClones[0], "pp #times 10^{0}", "F P L");
		
	  canv_p->cd();
	  gPad->SetLogy();

	  Double_t firstXVal = -1;
	  Double_t secondXVal = -1;
	  for(Int_t xI = 0; xI < nXVals; ++xI){
	    if(xVals[xI] < histVectPPClones[0]->GetBinLowEdge(1)) continue;
	    if(xVals[xI] >= histVectPPClones[0]->GetBinLowEdge(histVectPPClones[0]->GetNbinsX()+1)) continue;

	    if(firstXVal < 0) firstXVal = xVals[xI];
	    else if(secondXVal < 0) secondXVal = xVals[xI];
	    label_p->DrawLatex(xVals[xI], min/3., std::to_string(xVals[xI]).c_str());
	  }
	  label_p->DrawLatex(firstXVal + 10, max/5., "#bf{CMS}");
	  label_p->DrawLatex(firstXVal + 5, max*2., "#sqrt{s_{NN}} = 5.02 TeV, PbPb 404 #mub^{-1}, pp 27.3 pb^{-1}");
	  label_p->DrawLatex(secondXVal, max/5., ("anti-k_{t} R=" + rValStr + " jets").c_str());
	  label_p->DrawLatex(secondXVal, max/25., "|#eta_{jets}| #LT 2");
	  
	  leg_p->Draw("SAME");
	  gPad->RedrawAxis();
	  const std::string saveName = "pdfDir/" + dateStr + "/spectra_" + dirListPbPb[jI] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + dateStr+ ".pdf";
	  canv_p->SaveAs(saveName.c_str());
	  delete canv_p;
	  delete leg_p;

	  delete label_p;

	  outFile_p->cd();

	  scaleFactor = 1.;
	  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	    scaleFactor *= 10.;
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	      scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./scaleFactor);
	    }
	  }
	  
	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    histVectPPClones[sI]->Write("", TObject::kOverwrite);
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      histVectPbPbClones[cI + sI*nCentBins]->Write("", TObject::kOverwrite);
	    }
	  }
	
	  
	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    delete histVectPPClones[sI];
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      delete histVectPbPbClones[cI + sI*nCentBins];
	    }
	  }

	  histVectPPClones.clear();
	  histVectPbPbClones.clear();
	}
      }
    }
  }  
  std::cout << "End Traditional Spectra" << std::endl;

  std::cout << "Traditional RAA" << std::endl;
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    std::cout << "Jet Algo: " << jtAlgosPbPb[jI] << "/" << jtAlgosPP[jI] << std::endl;
    const Int_t rVal = getRVal(dirListPbPb.at(jI));
    const std::string rValStr = getRValStr(dirListPbPb.at(jI));
    const bool isSmallR = rReader.GetIsSmallR(rVal);
    
    Int_t nBinsTemp = -1;
    if(isSmallR){
      nBinsTemp = genJtPtBinsSmallR.size()-1;
      setGeneralPtBins(nGeneralPtBins, generalPtBins, genJtPtBinsSmallR);
    }
    else{
      nBinsTemp = genJtPtBinsLargeR.size()-1;
      setGeneralPtBins(nGeneralPtBins, generalPtBins, genJtPtBinsLargeR);
    }
    const Int_t nBins = nBinsTemp;

    for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
      if(!goodID[idI]) continue;

      for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	if(!goodResponseMod[rI]) continue;

	for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	  if(!goodJtAbsEtaBins[aI]) continue;	 

	  TLatex* label_p = new TLatex();
	  //	  label_p->SetNDC();
	  label_p->SetTextFont(43);
	  label_p->SetTextSize(16);
	  
	  TLegend* leg_p = new TLegend(0.25, 0.12, 0.4, 0.35);
	  leg_p->SetBorderSize(0);
	  leg_p->SetFillColor(0);
	  leg_p->SetFillStyle(0);
	  leg_p->SetTextFont(43);
	  leg_p->SetTextSize(16);

	  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);	  
	  defineCanv(canv_p);
	  
	  std::vector<TH1D*> histVectPbPbClones;
	  histVectPbPbClones.reserve(nSystInFile*nCentBins);
	  std::vector<TH1D*> histVectPPClones;	  
	  histVectPPClones.reserve(nSystInFile);

	  //Grab the subset histograms
	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    ULong64_t idKey = getKey(jI, 0, idI, rI, aI, sI);
	    ULong64_t idKeyBB = histKeyToBestBayesKeyPP[idKey];
	    ULong64_t vectPos = keyToVectPosPP[idKeyBB];
	    std::string tempName = "raa_" + std::string(histVectPP[vectPos]->GetName()) + "_Clone";
	    histVectPPClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];R_{AA}", nBins, generalPtBins));
	    if(!macroHistToSubsetHist(histVectPP[vectPos], histVectPPClones[sI])) return 1;

	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      idKey = getKey(jI, cI, idI, rI, aI, sI);
	      idKeyBB = histKeyToBestBayesKeyPbPb[idKey];
	      vectPos = keyToVectPosPbPb[idKeyBB];
	      tempName = "raa_" + std::string(histVectPbPb[vectPos]->GetName()) + "_Clone";
	      histVectPbPbClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];R_{AA}", nBins, generalPtBins));
	      if(!macroHistToSubsetHist(histVectPbPb[vectPos], histVectPbPbClones[cI + sI*nCentBins])) return 1;
	    }
	  }

	  //Div by bin widths;
          for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    divHistByWidth(histVectPPClones[sI]);
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      divHistByWidth(histVectPbPbClones[cI + sI*nCentBins]);
	    }
	  }

	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){	  
	    scaleHist(histVectPPClones[sI], 1./(2.*jtAbsEtaBinsWidth[aI]));
	    scaleHist(histVectPPClones[sI], 1./(lumiFactor));
	  }

	  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	    Double_t tempTAAFactor = getTAAScaleFactorNB(centBinsStr[cI]);
	    
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	      scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./(tempTAAFactor*nMBEvents));
	      //	      scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./tempTAAFactor);
	      scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./(2.*jtAbsEtaBinsWidth[aI]));
	      scaleHist(histVectPbPbClones[cI + sI*nCentBins], 100./centBinsWidth[cI]);

	      bool isSystCorr = false;
	      for(ULong64_t sI2 = 0; sI2 < (ULong64_t)corrSystStr.size(); ++sI2){
		if(isStrSame(corrSystStr[sI2], systStrInFile[sI])){
		  isSystCorr = true;
		  break;
		}
	      }

	      if(isSystCorr) createRAA(histVectPbPbClones[cI + sI*nCentBins], histVectPPClones[sI]);
	      else createRAA(histVectPbPbClones[cI + sI*nCentBins], histVectPPClones[0]);
	    }
	  }
	  
	  

	  //Setup systematics
	  std::vector<std::vector<Double_t> > reducedSystPbPb;
	  std::vector<std::vector<Double_t> > reducedSystSumPbPb;
	  for(ULong64_t sI = 0; sI < (ULong64_t)nReducedSyst; ++sI){
            for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      if(sI == 0) reducedSystSumPbPb.push_back({});
	      reducedSystPbPb.push_back({});

	      for(Int_t bI = 0; bI < nBins; ++bI){
		if(sI == 0) reducedSystSumPbPb[cI].push_back(0.0);
		reducedSystPbPb[cI + sI*nCentBins].push_back(1000000000000.);
	      }
	    }
	  }

	  //Extract systematics
	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    ULong64_t pos = 0;
	    if(systToCombo.count(systStrInFile[sI]) > 0) pos = posInStrVectExact(systToCombo[systStrInFile[sI]], reducedSystStr);
	    else pos = posInStrVectExact(systStrInFile[sI], reducedSystStr);

            for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      for(Int_t bI = 0; bI < nBins; ++bI){
		Double_t delta = TMath::Abs(histVectPbPbClones[cI + sI*nCentBins]->GetBinContent(bI+1) - histVectPbPbClones[cI]->GetBinContent(bI+1));

		(reducedSystPbPb[cI + pos*nCentBins])[bI] = TMath::Min((reducedSystPbPb[cI + pos*nCentBins])[bI], delta);
	      }
	    }
	  }

	  for(Int_t bI = 0; bI < nBins; ++bI){	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
		if(reducedSystStr[sI].find("PriorFlat") != std::string::npos) continue;
		(reducedSystSumPbPb[cI])[bI] = TMath::Sqrt((reducedSystSumPbPb[cI])[bI]*(reducedSystSumPbPb[cI])[bI] + (reducedSystPbPb[cI + sI*nCentBins])[bI]*(reducedSystPbPb[cI + sI*nCentBins])[bI]);
	      }
	    }
	  }
	  
	  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	    std::cout << "CentBins: " << centBinsStr[cI] << std::endl;
	    for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
	      std::cout << " SYST CHECK: " << reducedSystStr[sI] << std::endl;
	      for(Int_t bI = 0; bI < nBins; ++bI){
		std::cout <<"  " << bI << ": " << (reducedSystPbPb[cI + nCentBins*sI])[bI] << std::endl;
	      }	    
	    }
	  }
	  

	  //Grab Max and Min vals...	  
	  Double_t min = getHistMinGTZero(histVectPbPbClones[0]);
	  Double_t max = getHistMax(histVectPbPbClones[0]);
	  for(ULong64_t cI = 1; cI < (ULong64_t)nCentBins; ++cI){
	    Double_t tempMin = getHistMinGTZero(histVectPbPbClones[cI]);
	    Double_t tempMax = getHistMax(histVectPbPbClones[cI]);

	    if(tempMin < min) min = tempMin;
	    if(tempMax > max) max = tempMax;
	  }

	  max *= 1.5;
	  min = 0.0;

	  histVectPbPbClones[0]->SetMaximum(max);
	  histVectPbPbClones[0]->SetMinimum(min);
	  
	  canv_p->cd();
	  histVectPbPbClones[0]->SetMarkerColor(kPalette.getColor(getColorPosFromCent(centBinsStr[0], false)));
	  histVectPbPbClones[0]->SetLineColor(kPalette.getColor(getColorPosFromCent(centBinsStr[0], false)));
	  histVectPbPbClones[0]->SetMarkerStyle(getStyleFromCent(centBinsStr[0], false));
	  histVectPbPbClones[0]->SetMarkerSize(1);
	  centerTitles(histVectPbPbClones[0]);

	  histVectPbPbClones[0]->GetXaxis()->SetTitleOffset(histVectPbPbClones[0]->GetXaxis()->GetTitleOffset()*1.3);
	  histVectPbPbClones[0]->GetYaxis()->SetTitleOffset(2.2);
	  histVectPbPbClones[0]->GetXaxis()->SetLabelColor(0);

	  histVectPbPbClones[0]->Print("ALL");
	  histVectPbPbClones[0]->DrawCopy("HIST E1 P");	  
	  gPad->SetLogx();
	  
	  drawSyst(canv_p, histVectPbPbClones[0], reducedSystSumPbPb[0], histVectPbPbClones[0]->GetBinLowEdge(1), histVectPbPbClones[0]->GetBinLowEdge(histVectPbPbClones[0]->GetNbinsX()+1));
	  histVectPbPbClones[0]->DrawCopy("HIST E1 P SAME");	  
	  

	  for(ULong64_t cI = 1; cI < (ULong64_t)nCentBins; ++cI){
	    histVectPbPbClones[cI]->SetMarkerColor(kPalette.getColor(getColorPosFromCent(centBinsStr[cI], false)));
	    histVectPbPbClones[cI]->SetLineColor(kPalette.getColor(getColorPosFromCent(centBinsStr[cI], false)));
	    histVectPbPbClones[cI]->SetMarkerStyle(getStyleFromCent(centBinsStr[cI], false));
	    histVectPbPbClones[cI]->SetMarkerSize(1);
	    histVectPbPbClones[cI]->DrawCopy("HIST E1 P SAME");
	    drawSyst(canv_p, histVectPbPbClones[cI], reducedSystSumPbPb[cI], histVectPbPbClones[cI]->GetBinLowEdge(1), histVectPbPbClones[cI]->GetBinLowEdge(histVectPbPbClones[cI]->GetNbinsX()+1));
	    histVectPbPbClones[cI]->DrawCopy("HIST E1 P SAME");
	  }
	  
	  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	    ULong64_t pos = nCentBins - 1 - cI;
	    std::string centStr = std::to_string(centBinsLow[pos]) + "-" + std::to_string(centBinsHi[pos]) + "%";
	    leg_p->AddEntry(histVectPbPbClones[pos], centStr.c_str(), "F P L");
	  }
		
	  canv_p->cd();

	  Double_t firstXVal = -1;
	  Double_t secondXVal = -1;
	  for(Int_t xI = 0; xI < nXVals; ++xI){
	    if(xVals[xI] < histVectPPClones[0]->GetBinLowEdge(1)) continue;
	    if(xVals[xI] >= histVectPPClones[0]->GetBinLowEdge(histVectPPClones[0]->GetNbinsX()+1)) continue;

	    if(firstXVal < 0) firstXVal = xVals[xI];
	    else if(secondXVal < 0) secondXVal = xVals[xI];
	    label_p->DrawLatex(xVals[xI], min - (max - min)/25., std::to_string(xVals[xI]).c_str());
	  }
	  label_p->DrawLatex(firstXVal + 10, max*9.5/10., "#bf{CMS}");
	  label_p->DrawLatex(firstXVal + 5, max*10.25/10., "#sqrt{s_{NN}} = 5.02 TeV, PbPb 404 #mub^{-1}, pp 27.3 pb^{-1}");
	  label_p->DrawLatex(secondXVal, max*9.5/10., ("anti-k_{t} R=" + rValStr + " jets").c_str());
	  label_p->DrawLatex(secondXVal, max*8.75/10., "|#eta_{jets}| #LT 2");
	  
	  leg_p->Draw("SAME");
	  gPad->RedrawAxis();
	  const std::string saveName = "pdfDir/" + dateStr + "/raa_" + dirListPbPb[jI] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + dateStr+ ".pdf";
	  canv_p->SaveAs(saveName.c_str());
	  delete canv_p;
	  delete leg_p;

	  delete label_p;

	  outFile_p->cd();

	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      histVectPbPbClones[cI + sI*nCentBins]->Write("", TObject::kOverwrite);
	    }
	  }
	
	  
	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    delete histVectPPClones[sI];
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      delete histVectPbPbClones[cI + sI*nCentBins];
	    }
	  }

	  histVectPPClones.clear();
	  histVectPbPbClones.clear();
	}
      }
    }
  }
    
  std::cout << "End Traditional RAA" << std::endl;
  
  for(Int_t fI = 0; fI < nFiles; ++fI){
    inFile_p[fI]->Close();
    delete inFile_p[fI];
  }

  for(unsigned int i = 0; i < cutProps.size(); ++i){
    delete (cutProps[i]);
  }

  outFile_p->Close();
  delete outFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/plotUnfoldedAll.exe <inFileNamePP> <inFileNamePbPb>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotUnfoldedAll(argv[1], argv[2]);
  return retVal;
}
