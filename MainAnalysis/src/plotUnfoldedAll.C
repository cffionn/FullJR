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
#include "MainAnalysis/include/canvNDCToXY.h"
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/macroHistToSubsetHist.h"
#include "MainAnalysis/include/smallOrLargeR.h"
#include "MainAnalysis/include/systFunctions.h"

//Non-Local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/fileUtilities.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vanGoghPalette.h"

template <class T>
bool compWithWarning(const std::string typeStr, T a, T b)
{
  bool val = a == b;
  if(!val) std::cout << typeStr << " \'" << a << "\' from one file does not match \'" << b << "\' from other. fail";
  return val;
}

ULong64_t getKey(ULong64_t binsI, ULong64_t jI, ULong64_t cI, ULong64_t idI, ULong64_t rI, ULong64_t aI, ULong64_t sI){return binsI + jI*10 + 1000*cI + 100000*idI + 10000000*rI + 1000000000*aI + 100000000000*sI;}

ULong64_t getKey(ULong64_t binsI, ULong64_t jI, ULong64_t cI, ULong64_t idI, ULong64_t rI, ULong64_t aI, ULong64_t sI, ULong64_t bI){return binsI + 10*jI + 1000*cI + 100000*idI + 10000000*rI + 1000000000*aI + 100000000000*sI + 100000000000000*bI;}

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

void divHistByWidth(TH1D* hist_p)
{
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    Double_t binWidth = hist_p->GetBinWidth(bIX+1);
    hist_p->SetBinContent(bIX+1, hist_p->GetBinContent(bIX+1)/binWidth);
    hist_p->SetBinError(bIX+1, hist_p->GetBinError(bIX+1)/binWidth);
  }
  return;
}

void divHistByWidth(std::vector<TH1D*> hist_p)
{
  for(unsigned int sI = 0; sI < hist_p.size(); ++sI){
    for(Int_t bIX = 0; bIX < hist_p[sI]->GetNbinsX(); ++bIX){
      Double_t binWidth = hist_p[sI]->GetBinWidth(bIX+1);
      hist_p[sI]->SetBinContent(bIX+1, hist_p[sI]->GetBinContent(bIX+1)/binWidth);
      hist_p[sI]->SetBinError(bIX+1, hist_p[sI]->GetBinError(bIX+1)/binWidth);
    }
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

template <class T>
bool binRCheck(bool maxOrMin, std::string maxOrMinStr, Double_t maxOrMinVal, Int_t nBins, std::vector<T> bins)
{
  if(!maxOrMin){
    std::cout << maxOrMinStr << " largeR val \'" << maxOrMinVal << "\' is not found in bins: ";
    for(Int_t gI = 0; gI < nBins; ++gI){
      std::cout << bins[gI] << ", ";
    }
    std::cout << bins[nBins] << ". return 1" << std::endl;
    return false;
  }
  return true;
}

void labelStandard(TLatex* label_p)
{
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);
  
  return;
}

void legStandard(TLegend* leg_p)
{
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(16);
  
  return;
}

void setupSyst(const ULong64_t nReducedSyst, const Int_t nBins, const ULong64_t nCentBins, std::vector<std::vector<Double_t> >* systPbPb, std::vector<std::vector<Double_t> >* systSumPbPb, std::vector<std::vector<Double_t> >* systPP=NULL, std::vector<std::vector<Double_t> >* systSumPP=NULL)
{
  for(ULong64_t sI = 0; sI < nReducedSyst; ++sI){
    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
      if(sI == 0) systSumPbPb->push_back({});
      systPbPb->push_back({});
      
      for(Int_t bI = 0; bI < nBins; ++bI){
	if(sI == 0) (*systSumPbPb)[cI].push_back(0.0);
	(*systPbPb)[cI + sI*nCentBins].push_back(1000000000000.);
      }
    }
  }

  if(systPP != NULL){
    for(ULong64_t sI = 0; sI < nReducedSyst; ++sI){
      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	if(sI == 0) systSumPP->push_back({});
	systPP->push_back({});
	
	for(Int_t bI = 0; bI < nBins; ++bI){
	  if(sI == 0) (*systSumPP)[cI].push_back(0.0);
	  (*systPP)[cI + sI*nCentBins].push_back(1000000000000.);
	}
      }
    }
  }

  return;
}

void defaultPlotSet(TH1D* hist_p, std::string centStr = "")
{
  bool isPP = true;
  if(centStr.size() != 0) isPP = false;
  
  kirchnerPalette kPalette;

  hist_p->SetMarkerColor(kPalette.getColor(getColorPosFromCent(centStr, isPP)));
  hist_p->SetLineColor(kPalette.getColor(getColorPosFromCent(centStr, isPP)));
  hist_p->SetMarkerStyle(getStyleFromCent(centStr, isPP));
  hist_p->SetMarkerSize(1);
  centerTitles(hist_p);

  hist_p->GetXaxis()->SetTitleOffset(hist_p->GetXaxis()->GetTitleOffset()*1.3);
  hist_p->GetYaxis()->SetTitleOffset(2.2);
  hist_p->GetXaxis()->SetLabelColor(0);

  return;
}

void writeAll(TFile* outFile_p, std::vector<TH1D*> hist_p)
{
  outFile_p->cd();
  for(auto const & hist : hist_p){
    hist->Write("", TObject::kOverwrite);
  }
  return;
}

void deleteAll(TFile* outFile_p, std::vector<TH1D*>* hist_p)
{
  outFile_p->cd();
  for(unsigned int hI = 0; hI < hist_p->size(); ++hI){
    delete (*hist_p)[hI];
  }
  return;
}

std::string labelAbsEta(std::string absEtaString)
{
  absEtaString.replace(absEtaString.find("AbsEta"), 6, "");
  while(absEtaString.find("p") != std::string::npos){absEtaString.replace(absEtaString.find("p"), 1, ".");}
  absEtaString.replace(absEtaString.find("to"), 2, " < |#eta_{Jet}| < ");
  while(absEtaString.substr(0,1).find(" ") != std::string::npos){absEtaString.replace(0,1,"");}
  if(absEtaString.substr(0,3).find("0.0") != std::string::npos){
    unsigned int endPos = absEtaString.find("<");
    bool goodToRemove = true;
    for(unsigned int iter = 3; iter < endPos; ++iter){
      std::string tempStr = absEtaString.substr(iter,1);
      if(tempStr.find(" ") != std::string::npos) continue;
      if(tempStr.find("0") != std::string::npos) continue;
      if(tempStr.find("<") != std::string::npos) continue;

      goodToRemove = false;
      break;
    }

    if(goodToRemove) absEtaString.replace(0, absEtaString.find("<")+1, "");
  }
  while(absEtaString.substr(0,1).find(" ") != std::string::npos){absEtaString.replace(0,1,"");}

  unsigned int startPos = absEtaString.find("<");
  unsigned int endPos = absEtaString.rfind("<");

  if(startPos != endPos){
    std::string tempStr = absEtaString.substr(0, startPos);
    absEtaString.replace(0, startPos, "");
    while(isStrSame(tempStr.substr(tempStr.size()-1, 1), " ") || isStrSame(tempStr.substr(tempStr.size()-1, 1), "0")){
      tempStr = tempStr.substr(0,tempStr.size()-1);
    }
    tempStr = tempStr.substr(0,tempStr.size()-1);
    absEtaString = tempStr + " " + absEtaString;
  }

  std::string tempStr = absEtaString.substr(endPos, absEtaString.size());
  absEtaString.replace(endPos, absEtaString.size(), "");
  while(isStrSame(tempStr.substr(tempStr.size()-1, 1), " ") || isStrSame(tempStr.substr(tempStr.size()-1, 1), "0")){
    tempStr = tempStr.substr(0,tempStr.size()-1);
  }
  tempStr = tempStr.substr(0,tempStr.size()-1);
  absEtaString = absEtaString + tempStr;
  
  return absEtaString;
}

void drawAllLabels(TCanvas* canv_p, TH1D* hist_p, TLatex* label_p, const Int_t nXVals, const Int_t xVals[], const std::string rValStr, const std::string absEtaStr)
{
  canvNDCToXY labelAid(canv_p, hist_p);

  canv_p->cd();	  
  for(Int_t xI = 0; xI < nXVals; ++xI){
    if(xVals[xI] < hist_p->GetBinLowEdge(1)) continue;
    if(xVals[xI] >= hist_p->GetBinLowEdge(hist_p->GetNbinsX()+1)) continue;
    
    label_p->DrawLatex(labelAid.getXRelFromAbs(xVals[xI], false), canv_p->GetBottomMargin()-0.035, std::to_string(xVals[xI]).c_str());
  }
  label_p->DrawLatex(canv_p->GetLeftMargin()+0.04, 1.0-canv_p->GetTopMargin()-0.05, "#bf{CMS}");
  label_p->DrawLatex(canv_p->GetLeftMargin(), 1.0-canv_p->GetTopMargin()+0.02, "#sqrt{s_{NN}} = 5.02 TeV, PbPb 404 #mub^{-1}, pp 27.4 pb^{-1}");
  label_p->DrawLatex(0.5, 1.0-canv_p->GetTopMargin()-0.05, ("anti-k_{t} R=" + rValStr + " jets").c_str());
  label_p->DrawLatex(0.5, 1.0-canv_p->GetTopMargin()-0.10, labelAbsEta(absEtaStr).c_str());

  return;
}

int plotUnfoldedAll(const std::string inFileNamePP, const std::string inFileNamePbPb, const std::string inATLASFileName="")
{
  if(!fileIsGood(inFileNamePP, ".root")) return 1;
  if(!fileIsGood(inFileNamePbPb, ".root")) return 1;
  if(inATLASFileName.size() != 0 && !fileIsGood(inATLASFileName, ".root")) return 1;
    
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

  const Int_t nMaxCentBins = 4;
  Int_t nCentBins = -1;
  std::vector<Int_t> centBinsLow;
  std::vector<Int_t> centBinsHi;
  std::vector<std::string> centBinsStr;
  std::vector<Double_t> centBinsWidth;

  Int_t nID = -1;
  std::vector<std::string> idStr;
  std::vector<bool> goodID;

  const Int_t nSmallLargeBins = 2;
  std::vector<std::string> smallLargeBinsStr = {"SmallBins", "LargeBins"};

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

  Int_t nGenJtPtSmallBinsSmallR[nMaxCentBins];
  Int_t nGenJtPtSmallBinsLargeR[nMaxCentBins];
  Int_t nGenJtPtLargeBinsSmallR[nMaxCentBins];
  Int_t nGenJtPtLargeBinsLargeR[nMaxCentBins];
  Int_t nRecoJtPtBinsSmallR[nMaxCentBins];
  Int_t nRecoJtPtBinsLargeR[nMaxCentBins];

  std::vector<Double_t> genJtPtSmallBinsSmallRTemp;
  std::vector<Double_t> genJtPtSmallBinsLargeRTemp;
  std::vector<Double_t> genJtPtLargeBinsSmallRTemp;
  std::vector<Double_t> genJtPtLargeBinsLargeRTemp;
  std::vector<Double_t> genJtPtSmallBinsSmallR;
  std::vector<Double_t> genJtPtSmallBinsLargeR;
  std::vector<Double_t> genJtPtLargeBinsSmallR;
  std::vector<Double_t> genJtPtLargeBinsLargeR;
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

      nGenJtPtSmallBinsSmallR[0] = cutProps[pos]->GetNGenJtPtSmallBinsSmallRCent0to10();
      nGenJtPtSmallBinsLargeR[0] = cutProps[pos]->GetNGenJtPtSmallBinsLargeRCent0to10();
      nGenJtPtLargeBinsSmallR[0] = cutProps[pos]->GetNGenJtPtLargeBinsSmallRCent0to10();
      nGenJtPtLargeBinsLargeR[0] = cutProps[pos]->GetNGenJtPtLargeBinsLargeRCent0to10();
      nRecoJtPtBinsSmallR[0] = cutProps[pos]->GetNRecoJtPtBinsSmallRCent0to10();
      nRecoJtPtBinsLargeR[0] = cutProps[pos]->GetNRecoJtPtBinsLargeRCent0to10();

      nGenJtPtSmallBinsSmallR[1] = cutProps[pos]->GetNGenJtPtSmallBinsSmallRCent10to30();
      nGenJtPtSmallBinsLargeR[1] = cutProps[pos]->GetNGenJtPtSmallBinsLargeRCent10to30();
      nGenJtPtLargeBinsSmallR[1] = cutProps[pos]->GetNGenJtPtLargeBinsSmallRCent10to30();
      nGenJtPtLargeBinsLargeR[1] = cutProps[pos]->GetNGenJtPtLargeBinsLargeRCent10to30();
      nRecoJtPtBinsSmallR[1] = cutProps[pos]->GetNRecoJtPtBinsSmallRCent10to30();
      nRecoJtPtBinsLargeR[1] = cutProps[pos]->GetNRecoJtPtBinsLargeRCent10to30();

      nGenJtPtSmallBinsSmallR[2] = cutProps[pos]->GetNGenJtPtSmallBinsSmallRCent30to50();
      nGenJtPtSmallBinsLargeR[2] = cutProps[pos]->GetNGenJtPtSmallBinsLargeRCent30to50();
      nGenJtPtLargeBinsSmallR[2] = cutProps[pos]->GetNGenJtPtLargeBinsSmallRCent30to50();
      nGenJtPtLargeBinsLargeR[2] = cutProps[pos]->GetNGenJtPtLargeBinsLargeRCent30to50();
      nRecoJtPtBinsSmallR[2] = cutProps[pos]->GetNRecoJtPtBinsSmallRCent30to50();
      nRecoJtPtBinsLargeR[2] = cutProps[pos]->GetNRecoJtPtBinsLargeRCent30to50();

      nGenJtPtSmallBinsSmallR[3] = cutProps[pos]->GetNGenJtPtSmallBinsSmallRCent50to90();
      nGenJtPtSmallBinsLargeR[3] = cutProps[pos]->GetNGenJtPtSmallBinsLargeRCent50to90();
      nGenJtPtLargeBinsSmallR[3] = cutProps[pos]->GetNGenJtPtLargeBinsSmallRCent50to90();
      nGenJtPtLargeBinsLargeR[3] = cutProps[pos]->GetNGenJtPtLargeBinsLargeRCent50to90();
      nRecoJtPtBinsSmallR[3] = cutProps[pos]->GetNRecoJtPtBinsSmallRCent50to90();
      nRecoJtPtBinsLargeR[3] = cutProps[pos]->GetNRecoJtPtBinsLargeRCent50to90();

      genJtPtSmallBinsSmallRTemp = cutProps[pos]->GetGenJtPtSmallBinsSmallR();
      genJtPtSmallBinsLargeRTemp = cutProps[pos]->GetGenJtPtSmallBinsLargeR();
      genJtPtLargeBinsSmallRTemp = cutProps[pos]->GetGenJtPtLargeBinsSmallR();
      genJtPtLargeBinsLargeRTemp = cutProps[pos]->GetGenJtPtLargeBinsLargeR();
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

  if(nCentBins > nMaxCentBins){
    std::cout << "nCentBins \'" << nCentBins << "\' is less than nMaxCentBins \'" << nMaxCentBins << "\'. return 1" << std::endl;
    return 1;
  }

  smallOrLargeR rReader;
  std::cout << "Checking binning against rReader.." << std::endl;

  
  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent0to10(nGenJtPtSmallBinsSmallR[0])) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent0to10(nGenJtPtSmallBinsLargeR[0])) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent0to10(nGenJtPtLargeBinsSmallR[0])) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent0to10(nGenJtPtLargeBinsLargeR[0])) return 1;
  if(!rReader.CheckNRecoJtPtBinsSmallRCent0to10(nRecoJtPtBinsSmallR[0])) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent0to10(nRecoJtPtBinsLargeR[0])) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent10to30(nGenJtPtSmallBinsSmallR[1])) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent10to30(nGenJtPtSmallBinsLargeR[1])) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent10to30(nGenJtPtLargeBinsSmallR[1])) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent10to30(nGenJtPtLargeBinsLargeR[1])) return 1;
  if(!rReader.CheckNRecoJtPtBinsSmallRCent10to30(nRecoJtPtBinsSmallR[1])) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent10to30(nRecoJtPtBinsLargeR[1])) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent30to50(nGenJtPtSmallBinsSmallR[2])) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent30to50(nGenJtPtSmallBinsLargeR[2])) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent30to50(nGenJtPtLargeBinsSmallR[2])) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent30to50(nGenJtPtLargeBinsLargeR[2])) return 1;
  if(!rReader.CheckNRecoJtPtBinsSmallRCent30to50(nRecoJtPtBinsSmallR[2])) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent30to50(nRecoJtPtBinsLargeR[2])) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent50to90(nGenJtPtSmallBinsSmallR[3])) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent50to90(nGenJtPtSmallBinsLargeR[3])) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent50to90(nGenJtPtLargeBinsSmallR[3])) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent50to90(nGenJtPtLargeBinsLargeR[3])) return 1;
  if(!rReader.CheckNRecoJtPtBinsSmallRCent50to90(nRecoJtPtBinsSmallR[3])) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent50to90(nRecoJtPtBinsLargeR[3])) return 1;
  if(!rReader.CheckGenJtPtSmallBinsSmallR(genJtPtSmallBinsSmallRTemp)) return 1;
  if(!rReader.CheckGenJtPtSmallBinsLargeR(genJtPtSmallBinsLargeRTemp)) return 1;
  if(!rReader.CheckGenJtPtLargeBinsSmallR(genJtPtLargeBinsSmallRTemp)) return 1;
  if(!rReader.CheckGenJtPtLargeBinsLargeR(genJtPtLargeBinsLargeRTemp)) return 1;
  if(!rReader.CheckRecoJtPtBinsSmallR(recoJtPtBinsSmallRTemp)) return 1;
  if(!rReader.CheckRecoJtPtBinsLargeR(recoJtPtBinsLargeRTemp)) return 1;

  std::cout << " All bins good!" << std::endl;

  const Int_t nXVals = 5;
  const Int_t xVals[nXVals] = {200, 400, 600, 800, 1000};
  
  const Int_t nXValsReduced = 3;
  const Int_t xValsReduced[nXValsReduced] = {300, 400, 500};
  
  const Double_t minValSmallR = 200;
  const Double_t minValLargeR = 300;
  const Double_t maxValSmallR = 1000;
  const Double_t maxValLargeR = 1000;
  bool minSmallBinsSmallRFound = false;
  bool minSmallBinsLargeRFound = false;
  bool maxSmallBinsSmallRFound = false;
  bool maxSmallBinsLargeRFound = false;
  bool minLargeBinsSmallRFound = false;
  bool minLargeBinsLargeRFound = false;
  bool maxLargeBinsSmallRFound = false;
  bool maxLargeBinsLargeRFound = false;
  
  for(unsigned int gI = 0; gI < genJtPtSmallBinsSmallRTemp.size(); ++gI){
    if(TMath::Abs(genJtPtSmallBinsSmallRTemp[gI] - minValSmallR) < 1.) minSmallBinsSmallRFound = true;
    if(TMath::Abs(genJtPtSmallBinsSmallRTemp[gI] - maxValSmallR) < 1.) maxSmallBinsSmallRFound = true;

    if(minSmallBinsSmallRFound) genJtPtSmallBinsSmallR.push_back(genJtPtSmallBinsSmallRTemp[gI]);
    if(minSmallBinsSmallRFound && maxSmallBinsSmallRFound) break;
  }

  for(unsigned int gI = 0; gI < genJtPtSmallBinsLargeRTemp.size(); ++gI){
    if(TMath::Abs(genJtPtSmallBinsLargeRTemp[gI] - minValLargeR) < 1.) minSmallBinsLargeRFound = true;
    if(TMath::Abs(genJtPtSmallBinsLargeRTemp[gI] - maxValLargeR) < 1.) maxSmallBinsLargeRFound = true;

    if(minSmallBinsLargeRFound) genJtPtSmallBinsLargeR.push_back(genJtPtSmallBinsLargeRTemp[gI]);
    if(minSmallBinsLargeRFound && maxSmallBinsLargeRFound) break;
  }

  for(unsigned int gI = 0; gI < genJtPtLargeBinsSmallRTemp.size(); ++gI){
    if(TMath::Abs(genJtPtLargeBinsSmallRTemp[gI] - minValSmallR) < 1.) minLargeBinsSmallRFound = true;
    if(TMath::Abs(genJtPtLargeBinsSmallRTemp[gI] - maxValSmallR) < 1.) maxLargeBinsSmallRFound = true;

    if(minLargeBinsSmallRFound) genJtPtLargeBinsSmallR.push_back(genJtPtLargeBinsSmallRTemp[gI]);
    if(minLargeBinsSmallRFound && maxLargeBinsSmallRFound) break;
  }

  for(unsigned int gI = 0; gI < genJtPtLargeBinsLargeRTemp.size(); ++gI){
    if(TMath::Abs(genJtPtLargeBinsLargeRTemp[gI] - minValLargeR) < 1.) minLargeBinsLargeRFound = true;
    if(TMath::Abs(genJtPtLargeBinsLargeRTemp[gI] - maxValLargeR) < 1.) maxLargeBinsLargeRFound = true;

    if(minLargeBinsLargeRFound) genJtPtLargeBinsLargeR.push_back(genJtPtLargeBinsLargeRTemp[gI]);
    if(minLargeBinsLargeRFound && maxLargeBinsLargeRFound) break;
  }

  if(!binRCheck(minSmallBinsSmallRFound, "Min", minValSmallR, genJtPtSmallBinsSmallRTemp.size() - 1, genJtPtSmallBinsSmallRTemp)) return 1;
  if(!binRCheck(maxSmallBinsSmallRFound, "Max", maxValSmallR, genJtPtSmallBinsSmallRTemp.size() - 1, genJtPtSmallBinsSmallRTemp)) return 1;
  if(!binRCheck(minSmallBinsLargeRFound, "Min", minValLargeR, genJtPtSmallBinsLargeRTemp.size() - 1, genJtPtSmallBinsLargeRTemp)) return 1;
  if(!binRCheck(maxSmallBinsLargeRFound, "Max", maxValLargeR, genJtPtSmallBinsLargeRTemp.size() - 1, genJtPtSmallBinsLargeRTemp)) return 1;
  if(!binRCheck(minLargeBinsSmallRFound, "Min", minValSmallR, genJtPtLargeBinsSmallRTemp.size() - 1, genJtPtLargeBinsSmallRTemp)) return 1;
  if(!binRCheck(maxLargeBinsSmallRFound, "Max", maxValSmallR, genJtPtLargeBinsSmallRTemp.size() - 1, genJtPtLargeBinsSmallRTemp)) return 1;
  if(!binRCheck(minLargeBinsLargeRFound, "Min", minValLargeR, genJtPtLargeBinsLargeRTemp.size() - 1, genJtPtLargeBinsLargeRTemp)) return 1;
  if(!binRCheck(maxLargeBinsLargeRFound, "Max", maxValLargeR, genJtPtLargeBinsLargeRTemp.size() - 1, genJtPtLargeBinsLargeRTemp)) return 1;

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

  const UInt_t nJtAlgos = jtAlgosPbPb.size();
  std::vector<Int_t> rValI;
  std::vector<std::string> rValStr;
  std::vector<Bool_t> isSmallR;

  std::map<ULong64_t, ULong64_t> keyToVectPos;
  UInt_t nKey = 0;
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    rValI.push_back(getRVal(dirListPbPb.at(jI)));
    rValStr.push_back(getRValStr(dirListPbPb.at(jI)));
    isSmallR.push_back(rReader.GetIsSmallR(rValI[jI]));

    for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
      goodID.push_back(false);
      for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	goodResponseMod.push_back(false);
	for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	  goodJtAbsEtaBins.push_back(false);
	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	    for(ULong64_t bI = 0; bI < (ULong64_t)nBayes; ++bI){
	      for(ULong64_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  ULong64_t key = getKey(binsI, jI, cI, idI, rI, aI, sI, bI);
		  keyToVectPos[key] = nKey;
		  ++nKey;
		}	      
	      }
	    }
	  }
	}
      }
    }
  }
  std::cout << "nKeys: " << nKey << std::endl;
  std::vector<TH1D*> histVectPbPb;
  std::vector<TH1D*> histVectPP;

  histVectPbPb.reserve(nKey);
  for(ULong64_t kI = 0; kI < nKey; ++kI){histVectPbPb.push_back(NULL);}

  histVectPP.reserve(nKey);
  for(ULong64_t kI = 0; kI < nKey; ++kI){histVectPP.push_back(NULL);}

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
	
	int binsPos = posInStrVect(name, "_", smallLargeBinsStr, "_");
	int idPos = posInStrVect(name, "_", idStr, "_");
	int modPos = posInStrVect(name, "_", responseModStr, "_");
	int absEtaPos = posInStrVect(name, "_", jtAbsEtaBinsStr, "_");	int systPos = posInStrVect(name, "_", systStrInFile, "_");
	if(systPos < 0 && name.find("PriorFlat") == std::string::npos) systPos = 0;
	int bayesPos = posInStrVect(name, "_", bayesValStr, "_");
	int centPos = posInStrVect(name, "_", centBinsStr, "_");

	if(algoPos < 0) std::cout << __LINE__ << ": Missing algoPos in name \'" << name << "\'" << std::endl;
	if(idPos < 0) std::cout << __LINE__ << ": Missing idPos in name \'" << name << "\'" << std::endl;
	if(binsPos < 0) std::cout << __LINE__ << ": Missing binsPos in name \'" << name << "\'" << std::endl;
	if(modPos < 0) std::cout << __LINE__ << ": Missing modPos in name \'" << name << "\'" << std::endl;
	if(absEtaPos < 0) std::cout << __LINE__ << ": Missing absEtaPos in name \'" << name << "\'" << std::endl;
	if(systPos < 0) std::cout << __LINE__ << ": Missing systPos in name \'" << name << "\'" << std::endl;
	if(bayesPos < 0) std::cout << __LINE__ << ": Missing bayesPos in name \'" << name << "\'" << std::endl;
	if(centPos < 0) std::cout << __LINE__ << ": Missing centPos in name \'" << name << "\'" << std::endl;	

	goodID[idPos] = true;
	goodResponseMod[modPos] = true;
	goodJtAbsEtaBins[absEtaPos] = true;

	ULong64_t idKey = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos, bayesPos);
	ULong64_t vectPos = keyToVectPos[idKey];	  

	if(isPbPb[fI]) histVectPbPb[vectPos] = (TH1D*)key->ReadObj();
	else histVectPP[vectPos] = (TH1D*)key->ReadObj();
      }
    }
  }

  std::cout << "Processing best bayes PbPb..." << std::endl;
  for(auto const & tagToBB : histTagsToBestBayesPbPb){
    int binsPos = posInStrVect(tagToBB.first, "_", smallLargeBinsStr, "_");
    int algoPos = posInStrVect(tagToBB.first, "", dirListPbPb, "_");
    int idPos = posInStrVect(tagToBB.first, "_", idStr, "_");
    int modPos = posInStrVect(tagToBB.first, "_", responseModStr, "_");
    int absEtaPos = posInStrVect(tagToBB.first, "_", jtAbsEtaBinsStr, "_");
    int systPos = posInStrVect(tagToBB.first, "_", systStrInFile, "");
    if(systPos < 0 && tagToBB.first.find("PriorFlat") == std::string::npos) systPos = 0;
    int centPos = posInStrVect(tagToBB.first, "_", centBinsStr, "_");
    int bayesPos = tagToBB.second;
    
    if(binsPos < 0) std::cout << __LINE__ << ": Missing binsPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(algoPos < 0) std::cout << __LINE__ << ": Missing algoPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(idPos < 0) std::cout << __LINE__ << ": Missing idPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(modPos < 0) std::cout << __LINE__ << ": Missing modPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(absEtaPos < 0) std::cout << __LINE__ << ": Missing absEtaPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(systPos < 0) std::cout << __LINE__ << ": Missing systPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(centPos < 0) std::cout << __LINE__ << ": Missing centPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;	
    if(bayesPos < 0){
      std::cout << __LINE__ << ": Missing bayesPos in \'" << tagToBB.first << "\'. set to default 10" << std::endl;
      bayesPos = 10;
    }

    ULong64_t idKey = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos);
    ULong64_t idKeyBB = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos, (ULong64_t)bayesPos);
    histKeyToBestBayesKeyPbPb[idKey] = idKeyBB;
  }

  std::cout << "Processing best bayes PP..." << std::endl;
  for(auto const & tagToBB : histTagsToBestBayesPP){
    int binsPos = posInStrVect(tagToBB.first, "", smallLargeBinsStr, "_");
    int algoPos = posInStrVect(tagToBB.first, "", dirListPP, "_");
    int idPos = posInStrVect(tagToBB.first, "_", idStr, "_");
    int modPos = posInStrVect(tagToBB.first, "_", responseModStr, "_");
    int absEtaPos = posInStrVect(tagToBB.first, "_", jtAbsEtaBinsStr, "_");
    int systPos = posInStrVect(tagToBB.first, "_", systStrInFile, "");
    if(systPos < 0 && tagToBB.first.find("PriorFlat") == std::string::npos) systPos = 0;
    int centPos = posInStrVect(tagToBB.first, "_", centBinsStr, "_");

    int bayesPos = tagToBB.second;
    
    if(binsPos < 0) std::cout << __LINE__ << ": Missing binsPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(algoPos < 0) std::cout << __LINE__ << ": Missing algoPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(idPos < 0) std::cout << __LINE__ << ": Missing idPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(modPos < 0) std::cout << __LINE__ << ": Missing modPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(absEtaPos < 0) std::cout << __LINE__ << ": Missing absEtaPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(systPos < 0) std::cout << __LINE__ << ": Missing systPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(centPos < 0) std::cout << __LINE__ << ": Missing centPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;	
    if(bayesPos < 0){
      std::cout << __LINE__ << ": Missing bayesPos in \'" << tagToBB.first << "\'. set to default 10" << std::endl;
      bayesPos = 10;
    }

    ULong64_t idKey = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos);
    ULong64_t idKeyBB = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos, (ULong64_t)bayesPos);
    histKeyToBestBayesKeyPP[idKey] = idKeyBB;
  }

  std::cout << "Traditional Spectra" << std::endl;
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    std::cout << "Jet Algo: " << jtAlgosPbPb[jI] << "," << jtAlgosPP[jI] << std::endl;

    for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
      Int_t nBins[nMaxCentBins];
      Int_t nBinsTemp = -1;
      if(isSmallR[jI]){
	nBinsTemp = genJtPtLargeBinsSmallR.size()-1;
	setGeneralPtBins(nGeneralPtBins, generalPtBins, genJtPtLargeBinsSmallR);
      }
      else{
	nBinsTemp = genJtPtLargeBinsLargeR.size()-1;
	setGeneralPtBins(nGeneralPtBins, generalPtBins, genJtPtLargeBinsLargeR);
      }

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
	Int_t nBinsDelta = (genJtPtLargeBinsSmallRTemp.size()-1) - rReader.GetSmallOrLargeRNBins(isSmallR[jI], true, false, centStr);
	nBins[cI] = nBinsTemp - nBinsDelta;
      }

      for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
	if(!goodID[idI]) continue;
	
	for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	  if(!goodResponseMod[rI]) continue;
	  
	  for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	    if(!goodJtAbsEtaBins[aI]) continue;	 

	    TLatex* label_p = new TLatex();
	    labelStandard(label_p);
	    
	    TLegend* leg_p = new TLegend(0.25, 0.12, 0.4, 0.35);
	    legStandard(leg_p);
	    
	    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);	  
	    defineCanv(canv_p);
	    
	    std::vector<TH1D*> histVectPbPbClones;
	    histVectPbPbClones.reserve(nSystInFile*nCentBins);
	    std::vector<TH1D*> histVectPPClones;	  
	    histVectPPClones.reserve(nSystInFile*nCentBins);
	    
	    //Grab the subset histograms
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){

		ULong64_t idKey = getKey(binsI, jI, cI, idI, rI, aI, sI);
		ULong64_t idKeyBBPP = histKeyToBestBayesKeyPP[idKey];
		ULong64_t idKeyBBPbPb = histKeyToBestBayesKeyPbPb[idKey];
		ULong64_t vectPos = keyToVectPos[idKeyBBPP];
		std::string tempName = std::string(histVectPP[vectPos]->GetName()) + "_Clone";
		histVectPPClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];#frac{1}{#LTTAA#GT} #frac{1}{N_{evt}} #frac{d^{2}N_{Jet}}{dp_{T}d#eta}  and  #frac{d^{2}#sigma_{Jet}}{dp_{T}d#eta} [nb/GeV]", nBins[cI], generalPtBins));
	      
		if(!macroHistToSubsetHist(histVectPP[vectPos], histVectPPClones[cI + sI*nCentBins])) return 1;

		vectPos = keyToVectPos[idKeyBBPbPb];
		tempName = std::string(histVectPbPb[vectPos]->GetName()) + "_Clone";
		histVectPbPbClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];#frac{1}{#LTTAA#GT} #frac{1}{N_{evt}} #frac{d^{2}N_{Jet}}{dp_{T}d#eta}  and  #frac{d^{2}#sigma_{Jet}}{dp_{T}d#eta} [nb/GeV]", nBins[cI], generalPtBins));

		if(!macroHistToSubsetHist(histVectPbPb[vectPos], histVectPbPbClones[cI + sI*nCentBins])) return 1;
	      }
	    }
	    
	    //Div by bin widths;
	    divHistByWidth(histVectPPClones);
	    divHistByWidth(histVectPbPbClones);
      	    
	    Double_t scaleFactor = 1.;
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){

	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){	  
		scaleHist(histVectPPClones[cI + sI*nCentBins], 1./(2.*jtAbsEtaBinsWidth[aI]*lumiFactor));
	      }

	      scaleFactor *= 10.;
	      Double_t totalFactor = getTAAScaleFactorNB(centBinsStr[cI])*nMBEvents*2.*jtAbsEtaBinsWidth[aI]*centBinsWidth[cI]/(100.*scaleFactor);
	      
	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
		scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./totalFactor);
	      }
	    }

	    //Setup systematics
	    std::vector<std::vector<Double_t> > reducedSystPP, reducedSystPbPb, reducedSystSumPP, reducedSystSumPbPb;
	    setupSyst(nReducedSyst, nBins[0], nCentBins, &reducedSystPbPb, &reducedSystSumPbPb, &reducedSystPP, &reducedSystSumPP);

	    //Extract systematics
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	      ULong64_t pos = 0;
	      if(systToCombo.count(systStrInFile[sI]) > 0) pos = posInStrVectExact(systToCombo[systStrInFile[sI]], reducedSystStr);
	      else pos = posInStrVectExact(systStrInFile[sI], reducedSystStr);
	      
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		for(Int_t bI = 0; bI < nBins[cI]; ++bI){
		  Double_t delta = TMath::Abs(histVectPPClones[cI + sI*nCentBins]->GetBinContent(bI+1) - histVectPPClones[cI]->GetBinContent(bI+1));
		  (reducedSystPP[cI + pos*nCentBins])[bI] = TMath::Min((reducedSystPP[cI + pos*nCentBins])[bI], delta);

		  delta = TMath::Abs(histVectPbPbClones[cI + sI*nCentBins]->GetBinContent(bI+1) - histVectPbPbClones[cI]->GetBinContent(bI+1));
		  
		  (reducedSystPbPb[cI + pos*nCentBins])[bI] = TMath::Min((reducedSystPbPb[cI + pos*nCentBins])[bI], delta);
		}
	      }
	    }

	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      for(Int_t bI = 0; bI < nBins[cI]; ++bI){
		for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
		  if(reducedSystStr[sI].find("PriorFlat") != std::string::npos) continue;

		  (reducedSystSumPP[cI])[bI] = TMath::Sqrt((reducedSystSumPP[cI])[bI]*(reducedSystSumPP[cI])[bI] + (reducedSystPP[cI + sI*nCentBins])[bI]*(reducedSystPP[cI + sI*nCentBins])[bI]);

		  (reducedSystSumPbPb[cI])[bI] = TMath::Sqrt((reducedSystSumPbPb[cI])[bI]*(reducedSystSumPbPb[cI])[bI] + (reducedSystPbPb[cI + sI*nCentBins])[bI]*(reducedSystPbPb[cI + sI*nCentBins])[bI]);
		}
	      }
	    }

	    //Grab Max and Min vals...	  
	    Double_t max = 10000;
	    Double_t min = 0.00000001;

	    histVectPPClones[0]->SetMaximum(max);
	    histVectPPClones[0]->SetMinimum(min);

	    canv_p->cd();
	    defaultPlotSet(histVectPPClones[0], "");	    
	    histVectPPClones[0]->DrawCopy("HIST E1 P");	  
	    canvNDCToXY labelAid(canv_p, histVectPPClones[0]);
	    //	    gPad->SetLogx();	    
	    //	    drawSyst(canv_p, histVectPPClones[0], reducedSystSumPP[0], histVectPPClones[0]->GetBinLowEdge(1), histVectPPClones[0]->GetBinLowEdge(histVectPPClones[0]->GetNbinsX()+1));
	    histVectPPClones[0]->DrawCopy("HIST E1 P SAME");	  

	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      defaultPlotSet(histVectPbPbClones[cI], centBinsStr[cI]);
	      histVectPbPbClones[cI]->DrawCopy("HIST E1 P SAME");
	      //	      drawSyst(canv_p, histVectPbPbClones[cI], reducedSystSumPbPb[cI], histVectPbPbClones[cI]->GetBinLowEdge(1), histVectPbPbClones[cI]->GetBinLowEdge(histVectPbPbClones[cI]->GetNbinsX()+1));
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
	  
	    drawAllLabels(canv_p, histVectPPClones[0], label_p, nXVals, xVals, rValStr[jI], jtAbsEtaBinsStr[aI]);
	    
	    leg_p->Draw("SAME");
	    TFile* inATLASFile_p = NULL;
	    if(inATLASFileName.size() != 0 && rValI[jI] == 4){
	      inATLASFile_p = new TFile(inATLASFileName.c_str(), "READ");
	      TH1D* histPP_p = (TH1D*)inATLASFile_p->Get("Table 5/Hist1D_y1");
	      TH1D* histPbPb5090_p = (TH1D*)inATLASFile_p->Get("Table 11/Hist1D_y1");
	      
	      histPP_p->SetMarkerSize(1);
	      histPP_p->SetMarkerStyle(24);
	      histPP_p->SetMarkerColor(1);
	      histPP_p->SetFillColor(0);
	      histPP_p->SetLineColor(1);
	      histPP_p->SetLineWidth(1);
	      for(Int_t bIX = 0; bIX < histPP_p->GetNbinsX(); ++bIX){
		histPP_p->SetBinError(bIX+1, 0);
	      }
	      scaleHist(histPP_p, 2.8/2.0);
	      histPP_p->DrawCopy("HIST E1 SAME");
	      histPP_p->DrawCopy("HIST E1 P SAME");

	      histPbPb5090_p->SetMarkerSize(1);
	      histPbPb5090_p->SetMarkerStyle(25);
	      histPbPb5090_p->SetMarkerColor(1);
	      histPbPb5090_p->SetFillColor(0);
	      histPbPb5090_p->SetLineColor(1);
	      histPbPb5090_p->SetLineWidth(1);
	      for(Int_t bIX = 0; bIX < histPbPb5090_p->GetNbinsX(); ++bIX){
		histPbPb5090_p->SetBinError(bIX+1, 0);
	      }

	      scaleFactor = 1;
	      for(Int_t scaleI = 0; scaleI < nCentBins; ++scaleI){
		scaleFactor *= 10;
	      }
	      scaleHist(histPbPb5090_p, scaleFactor*2.8/2.0);
	      histPbPb5090_p->DrawCopy("HIST E1  SAME");	     
	      histPbPb5090_p->DrawCopy("HIST E1 P SAME");	     

	      leg_p->AddEntry(histPP_p, "Matched ATLAS PP, 50-90%", "L P");
	    }

	    gPad->RedrawAxis();
	    const std::string saveName = "pdfDir/" + dateStr + "/spectra_" + dirListPbPb[jI] + "_" + smallLargeBinsStr[binsI] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + dateStr+ ".pdf";
	    canv_p->SaveAs(saveName.c_str());
	    delete canv_p;
	    delete leg_p;

	    if(inATLASFile_p != NULL){
	      inATLASFile_p->Close();
	      delete inATLASFile_p;
	    }
	    
	    delete label_p;
	    
	    outFile_p->cd();
	    
	    scaleFactor = 1.;
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      scaleFactor *= 10.;
	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
		scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./scaleFactor);
	      }
	    }

	    writeAll(outFile_p, histVectPPClones);
	    writeAll(outFile_p, histVectPbPbClones);
	    deleteAll(outFile_p, &histVectPPClones);
	    deleteAll(outFile_p, &histVectPbPbClones);	    
	    histVectPPClones.clear();
	    histVectPbPbClones.clear();
	  }
	}
      }
    }
  }  
  std::cout << "End Traditional Spectra" << std::endl;
  
  std::cout << "Traditional RAA" << std::endl;
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    std::cout << "Jet Algo: " << jtAlgosPbPb[jI] << "/" << jtAlgosPP[jI] << std::endl;

    for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
      Int_t nBins[nMaxCentBins];
      Int_t nBinsTemp = -1;
      if(isSmallR[jI]){
	nBinsTemp = genJtPtLargeBinsSmallR.size()-1;
	setGeneralPtBins(nGeneralPtBins, generalPtBins, genJtPtLargeBinsSmallR);
      }
      else{
	nBinsTemp = genJtPtLargeBinsLargeR.size()-1;
	setGeneralPtBins(nGeneralPtBins, generalPtBins, genJtPtLargeBinsLargeR);
      }

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
        Int_t nBinsDelta = (genJtPtLargeBinsSmallRTemp.size()-1) - rReader.GetSmallOrLargeRNBins(isSmallR[jI], true, false, centStr);
        nBins[cI] = nBinsTemp - nBinsDelta;
      }

      for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
	if(!goodID[idI]) continue;
	
	for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	  if(!goodResponseMod[rI]) continue;
	  
	  for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	    if(!goodJtAbsEtaBins[aI]) continue;	 

	    TLatex* label_p = new TLatex();
	    labelStandard(label_p);
	    
	    TLegend* leg_p = new TLegend(0.25, 0.12, 0.4, 0.35);
	    legStandard(leg_p);

	    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);	  
	    defineCanv(canv_p);
	    
	    std::vector<TH1D*> histVectPbPbClones;
	    histVectPbPbClones.reserve(nSystInFile*nCentBins);
	    std::vector<TH1D*> histVectPPClones;	  
	    histVectPPClones.reserve(nSystInFile*nCentBins);
	    
	    //Grab the subset histograms
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		ULong64_t idKey = getKey(binsI, jI, cI, idI, rI, aI, sI);
		ULong64_t idKeyBBPP = histKeyToBestBayesKeyPP[idKey];
		ULong64_t idKeyBBPbPb = histKeyToBestBayesKeyPbPb[idKey];
		ULong64_t vectPos = keyToVectPos[idKeyBBPP];
		std::string tempName = "raa_" + std::string(histVectPP[vectPos]->GetName()) + "_Clone";
		histVectPPClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];R_{AA}", nBins[cI], generalPtBins));

		if(!macroHistToSubsetHist(histVectPP[vectPos], histVectPPClones[cI + sI*nCentBins])) return 1;	      

		idKey = getKey(binsI, jI, cI, idI, rI, aI, sI);
		vectPos = keyToVectPos[idKeyBBPbPb];
		tempName = "raa_" + std::string(histVectPbPb[vectPos]->GetName()) + "_Clone";

		histVectPbPbClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];R_{AA}", nBins[cI], generalPtBins));

		if(!macroHistToSubsetHist(histVectPbPb[vectPos], histVectPbPbClones[cI + sI*nCentBins])) return 1;
	      }
	    }

	    //Div by bin widths;
	    divHistByWidth(histVectPPClones);
	    divHistByWidth(histVectPbPbClones);
	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      Double_t tempTAAFactor = getTAAScaleFactorNB(centBinsStr[cI]);	      

	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
		Double_t totalFactor = tempTAAFactor*nMBEvents*2.*jtAbsEtaBinsWidth[aI]*centBinsWidth[cI]/100.;
		scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./totalFactor);
		scaleHist(histVectPPClones[cI + sI*nCentBins], 1./(2.*jtAbsEtaBinsWidth[aI]*lumiFactor));

		bool isSystCorr = false;
		for(ULong64_t sI2 = 0; sI2 < (ULong64_t)corrSystStr.size(); ++sI2){
		  if(isStrSame(corrSystStr[sI2], systStrInFile[sI])){
		    isSystCorr = true;
		    break;
		  }
		}
		
		if(isSystCorr) createRAA(histVectPbPbClones[cI + sI*nCentBins], histVectPPClones[cI + sI*nCentBins]);
		else createRAA(histVectPbPbClones[cI + sI*nCentBins], histVectPPClones[cI]);
	      }
	    }
	    	
	    //Setup systematics
	    std::vector<std::vector<Double_t> > reducedSystPbPb, reducedSystSumPbPb;
	    setupSyst(nReducedSyst, nBins[0], nCentBins, &reducedSystPbPb, &reducedSystSumPbPb);
	    
	    //Extract systematics
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
	      ULong64_t pos = 0;
	      if(systToCombo.count(systStrInFile[sI]) > 0) pos = posInStrVectExact(systToCombo[systStrInFile[sI]], reducedSystStr);
	      else pos = posInStrVectExact(systStrInFile[sI], reducedSystStr);
	    
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		for(Int_t bI = 0; bI < nBins[cI]; ++bI){
		  Double_t delta = TMath::Abs(histVectPbPbClones[cI + sI*nCentBins]->GetBinContent(bI+1) - histVectPbPbClones[cI]->GetBinContent(bI+1));
		  
		  (reducedSystPbPb[cI + pos*nCentBins])[bI] = TMath::Min((reducedSystPbPb[cI + pos*nCentBins])[bI], delta);
		}
	      }
	    }
	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      for(Int_t bI = 0; bI < nBins[cI]; ++bI){	    
		for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
		  if(reducedSystStr[sI].find("PriorFlat") != std::string::npos) continue;
		  (reducedSystSumPbPb[cI])[bI] = TMath::Sqrt((reducedSystSumPbPb[cI])[bI]*(reducedSystSumPbPb[cI])[bI] + (reducedSystPbPb[cI + sI*nCentBins])[bI]*(reducedSystPbPb[cI + sI*nCentBins])[bI]);
		}
	      }
	    }
	    
	    //Grab Max and Min vals...	  
	    Double_t min = 0.0;
	    Double_t max = 1.2;
	    histVectPbPbClones[0]->SetMaximum(max);
	    histVectPbPbClones[0]->SetMinimum(min);
	    
	    canv_p->cd();
	    defaultPlotSet(histVectPbPbClones[0], centBinsStr[0]);
	    histVectPbPbClones[0]->DrawCopy("HIST E1 P");	  
	    canvNDCToXY labelAid(canv_p, histVectPbPbClones[0]);
	    gPad->SetLogx();
	    
	    drawSyst(canv_p, histVectPbPbClones[0], reducedSystSumPbPb[0], histVectPbPbClones[0]->GetBinLowEdge(1), histVectPbPbClones[0]->GetBinLowEdge(histVectPbPbClones[0]->GetNbinsX()+1));
	    histVectPbPbClones[0]->DrawCopy("HIST E1 P SAME");

	    for(ULong64_t cI = 1; cI < (ULong64_t)nCentBins; ++cI){
	      defaultPlotSet(histVectPbPbClones[cI], centBinsStr[cI]);
	      histVectPbPbClones[cI]->DrawCopy("HIST E1 P SAME");
	      drawSyst(canv_p, histVectPbPbClones[cI], reducedSystSumPbPb[cI], histVectPbPbClones[cI]->GetBinLowEdge(1), histVectPbPbClones[cI]->GetBinLowEdge(histVectPbPbClones[cI]->GetNbinsX()+1));
	      histVectPbPbClones[cI]->DrawCopy("HIST E1 P SAME");
	    }
	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      ULong64_t pos = nCentBins - 1 - cI;
	      std::string centStr = std::to_string(centBinsLow[pos]) + "-" + std::to_string(centBinsHi[pos]) + "%";
	      leg_p->AddEntry(histVectPbPbClones[pos], centStr.c_str(), "F P L");
	    }
	    
	    //	    drawAllLabels(canv_p, histVectPPClones[0], label_p, nXVals, xVals, rValStr[jI], jtAbsEtaBinsStr[aI]);
	  	    
	    leg_p->Draw("SAME");
	    gPad->RedrawAxis();
	    const std::string saveName = "pdfDir/" + dateStr + "/raa_" + dirListPbPb[jI] + "_" + smallLargeBinsStr[binsI] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + dateStr+ ".pdf";
	    canv_p->SaveAs(saveName.c_str());
	    delete canv_p;
	    delete leg_p;
	    
	    delete label_p;
	    
	    outFile_p->cd();
	    
	    writeAll(outFile_p, histVectPbPbClones);
	    deleteAll(outFile_p, &histVectPPClones);
	    deleteAll(outFile_p, &histVectPbPbClones);	    
	    histVectPPClones.clear();
	    histVectPbPbClones.clear();
	  }
	}
      }
    }
  }
    
  std::cout << "End Traditional RAA" << std::endl;
  
  std::cout << "Reduced RAA" << std::endl;
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    std::cout << "Jet Algo: " << jtAlgosPbPb[jI] << "/" << jtAlgosPP[jI] << std::endl;
    
    for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
      const Int_t nBins = 1;
      generalPtBins[0] = 300;
      generalPtBins[1] = 600;

      for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
	if(!goodID[idI]) continue;
	
	for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	  if(!goodResponseMod[rI]) continue;
	  
	  for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	    if(!goodJtAbsEtaBins[aI]) continue;	 
	    
	    TLatex* label_p = new TLatex();
	    labelStandard(label_p);
	    
	    TLegend* leg_p = new TLegend(0.25, 0.12, 0.4, 0.35);
	    legStandard(leg_p);
	    
	    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);	  
	    defineCanv(canv_p);
	    
	    std::vector<TH1D*> histVectPbPbClones;
	    histVectPbPbClones.reserve(nSystInFile*nCentBins);
	    std::vector<TH1D*> histVectPPClones;	  
	    histVectPPClones.reserve(nSystInFile*nCentBins);
	    
	    //Grab the subset histograms
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){	      
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		ULong64_t idKey = getKey(binsI, jI, cI, idI, rI, aI, sI);
		ULong64_t idKeyBBPP = histKeyToBestBayesKeyPP[idKey];
		ULong64_t idKeyBBPbPb = histKeyToBestBayesKeyPbPb[idKey];
		ULong64_t vectPos = keyToVectPos[idKeyBBPP];
		std::string tempName = "raaReduced_" + std::string(histVectPP[vectPos]->GetName()) + "_Clone";
		histVectPPClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];R_{AA}", nBins, generalPtBins));
		if(!macroHistToSubsetHist(histVectPP[vectPos], histVectPPClones[cI + sI*nCentBins])) return 1;

		vectPos = keyToVectPos[idKeyBBPbPb];
		tempName = "raaReduced_" + std::string(histVectPbPb[vectPos]->GetName()) + "_Clone";
		histVectPbPbClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];R_{AA}", nBins, generalPtBins));
		if(!macroHistToSubsetHist(histVectPbPb[vectPos], histVectPbPbClones[cI + sI*nCentBins])) return 1;
	      }
	    }
	    
	    //Div by bin widths;
	    divHistByWidth(histVectPPClones);
	    divHistByWidth(histVectPbPbClones);
	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      Double_t tempTAAFactor = getTAAScaleFactorNB(centBinsStr[cI]);
	      
	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){
		scaleHist(histVectPPClones[cI + sI*nCentBins], 1./(2.*jtAbsEtaBinsWidth[aI]*lumiFactor));

		Double_t totalFactor = tempTAAFactor*nMBEvents*2.*jtAbsEtaBinsWidth[aI]*centBinsWidth[cI]/100.;
		scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./totalFactor);
		
		bool isSystCorr = false;
		for(ULong64_t sI2 = 0; sI2 < (ULong64_t)corrSystStr.size(); ++sI2){
		  if(isStrSame(corrSystStr[sI2], systStrInFile[sI])){
		    isSystCorr = true;
		    break;
		  }
		}
		
		if(isSystCorr) createRAA(histVectPbPbClones[cI + sI*nCentBins], histVectPPClones[cI + sI*nCentBins]);
		else createRAA(histVectPbPbClones[cI + sI*nCentBins], histVectPPClones[cI]);
	      }
	    }
	    	    
	    //Setup systematics
	    std::vector<std::vector<Double_t> > reducedSystPbPb, reducedSystSumPbPb;
	    setupSyst(nReducedSyst, nBins, nCentBins, &reducedSystPbPb, &reducedSystSumPbPb);	    
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
	    	    
	    //Grab Max and Min vals...	  
	    Double_t max = 1.2;
	    Double_t min = 0.0;	    
	    histVectPbPbClones[0]->SetMaximum(max);
	    histVectPbPbClones[0]->SetMinimum(min);
	    
	    canv_p->cd();
	    defaultPlotSet(histVectPbPbClones[0], centBinsStr[0]);
	    histVectPbPbClones[0]->DrawCopy("HIST E1 P");	  
	    canvNDCToXY labelAid(canv_p, histVectPbPbClones[0]);
	    gPad->SetLogx();	    
	    drawSyst(canv_p, histVectPbPbClones[0], reducedSystSumPbPb[0], histVectPbPbClones[0]->GetBinLowEdge(1), histVectPbPbClones[0]->GetBinLowEdge(histVectPbPbClones[0]->GetNbinsX()+1));
	    histVectPbPbClones[0]->DrawCopy("HIST E1 P SAME");	  
  	    	    
	    for(ULong64_t cI = 1; cI < (ULong64_t)nCentBins; ++cI){
	      defaultPlotSet(histVectPbPbClones[cI], centBinsStr[cI]);
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
	    drawAllLabels(canv_p, histVectPPClones[0], label_p, nXValsReduced, xValsReduced, rValStr[jI], jtAbsEtaBinsStr[aI]);
	    
	    leg_p->Draw("SAME");
	    gPad->RedrawAxis();
	    const std::string saveName = "pdfDir/" + dateStr + "/raaReduced_" + dirListPbPb[jI] + "_" + smallLargeBinsStr[binsI] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + dateStr+ ".pdf";
	    canv_p->SaveAs(saveName.c_str());
	    delete canv_p;
	    delete leg_p;
	    
	    delete label_p;
	    
	    outFile_p->cd();

	    writeAll(outFile_p, histVectPbPbClones);
	    deleteAll(outFile_p, &histVectPPClones);
	    deleteAll(outFile_p, &histVectPbPbClones);	    
	    histVectPPClones.clear();
	    histVectPbPbClones.clear();
	  }
	}
      }
    }
  }
  std::cout << "End Reduced RAA" << std::endl;
  
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
  if(argc < 3 || argc > 4){
    std::cout << "Usage: ./bin/plotUnfoldedAll.exe <inFileNamePP> <inFileNamePbPb> <inATLASFileName-opt>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 3) retVal += plotUnfoldedAll(argv[1], argv[2]);
  else if(argc == 4) retVal += plotUnfoldedAll(argv[1], argv[2], argv[3]);
  return retVal;
}
