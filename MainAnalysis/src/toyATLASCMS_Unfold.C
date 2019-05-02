//cpp
#include <iostream>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TStyle.h"

//Local dependencies
#include "MainAnalysis/include/macroHistToSubsetHist.h"
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/plotUtilities.h"

//RooUnfoldResponse
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldResponse.h"

void matrixPrint(TH2D* hist_p, bool doYNorm)
{
  std::cout << "PRINTING, " << hist_p->GetName() << std::endl;
  for(Int_t bIY = 0; bIY < hist_p->GetYaxis()->GetNbins(); ++bIY){
    Double_t total = 0.0;
    for(Int_t bIX = 0; bIX < hist_p->GetXaxis()->GetNbins(); ++bIX){
      total += hist_p->GetBinContent(bIX+1, bIY+1);
    }

    for(Int_t bIX = 0; bIX < hist_p->GetXaxis()->GetNbins(); ++bIX){
      if(doYNorm) std::cout << hist_p->GetBinContent(bIX+1, bIY+1)/total << ",";
      else std::cout << hist_p->GetBinContent(bIX+1, bIY+1)<< ",";
    }
    std::cout << std::endl;
  }

  return;
}

void setPrior(TH2D* res_p, TH1D* prior_p)
{
  //  prior_p->Scale(res_p->Integral()/prior_p->Integral());

  Double_t origScaleFact = res_p->Integral()/prior_p->Integral();
  
  for(Int_t bIY = 0; bIY < res_p->GetYaxis()->GetNbins(); ++bIY){
    Double_t total = 0.0;
    for(Int_t bIX = 0; bIX < res_p->GetXaxis()->GetNbins(); ++bIX){
      total += res_p->GetBinContent(bIX+1, bIY+1);
    }

    if(total <= TMath::Power(10, -20)) continue;

    total = prior_p->GetBinContent(bIY+1)/total;
    for(Int_t bIX = 0; bIX < res_p->GetXaxis()->GetNbins(); ++bIX){
      Double_t val = res_p->GetBinContent(bIX+1, bIY+1)*total*origScaleFact;
      Double_t err = res_p->GetBinError(bIX+1, bIY+1)*total*origScaleFact;

      res_p->SetBinContent(bIX+1, bIY+1, val);
      res_p->SetBinError(bIX+1, bIY+1, err);
    }
  }
  
  return;
}

std::vector<double> binsStrToVect(std::string binsStr)
{
  binsStr = binsStr + ",";
  while(binsStr.find(",,") != std::string::npos){
    binsStr.replace(binsStr.find(",,"), 2, ",");
  }
  if(binsStr.size() == 1) return {};

  std::vector<double> binVect;
  while(binsStr.find(",") != std::string::npos){
    binVect.push_back(std::stod(binsStr.substr(0, binsStr.find(","))));
    binsStr.replace(0, binsStr.find(",")+1, "");
  }

  
  return binVect;
}

bool replaceBins(std::vector<double> inBinsVect, Int_t* nBins, Double_t bins[])
{  
  for(unsigned int bI = 0; bI < inBinsVect.size(); ++bI){
    bool isBinGood = false;

    for(Int_t bI2 = 0; bI2 < (*nBins)+1; ++bI2){
      if(TMath::Abs(inBinsVect[bI] - bins[bI2]) < 1.){
	isBinGood = true;
	break;
      }
    }

    if(!isBinGood) return false;
  }
  
  (*nBins) = inBinsVect.size()-1;
  for(unsigned int bI = 0; bI < inBinsVect.size(); ++bI){
    bins[bI] = inBinsVect[bI];
  }

  return true;
}

void histPrint(TH1D* hist_p)
{
  std::cout << hist_p->GetName() << ": ";
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    std::cout << hist_p->GetBinContent(bIX+1) << ", ";
  }
  std::cout << std::endl;
  
  return;
}

void doUnfold(TH2D* matrix_p, TH1D* prior_p, TH1D* data_p, Int_t nBayes, Int_t bayesVals[], std::string tagStr, TH1D* unfolded_p[], Double_t chi2[], TH2D* newMatrix_p, cppWatch* inWatch)
{
  TH2D* matrixClone_p = (TH2D*)matrix_p->Clone("matrixClone_h");

  const Int_t nMaxBins = 1000;
  Double_t recoBins[nMaxBins+1];
  Double_t genBins[nMaxBins+1];
  Int_t nRecoBins = -1;
  Int_t nGenBins = -1;
  
  for(Int_t bIY = 0; bIY < matrixClone_p->GetYaxis()->GetNbins()+1; ++bIY){
    genBins[bIY] = matrixClone_p->GetYaxis()->GetBinLowEdge(bIY+1);
    ++nGenBins;
  }

  for(Int_t bIX = 0; bIX < matrixClone_p->GetXaxis()->GetNbins()+1; ++bIX){
    recoBins[bIX] = matrixClone_p->GetXaxis()->GetBinLowEdge(bIX+1);
    ++nRecoBins;
  }

  TH1D* gen_p = new TH1D("gen_h", "", nGenBins, genBins);
  TH1D* reco_p = new TH1D("reco_h", "", nRecoBins, recoBins);
  //  if(newMatrix_p != NULL){


  setPrior(matrixClone_p, prior_p);
  if(tagStr.find("Random0") != std::string::npos && tagStr.find("Target0") != std::string::npos){
    std::cout << "PRIOR CHECK (response, prior): " << tagStr << std::endl;
    
    //    matrixPrint(matrixClone_p, false);
    //    histPrint(prior_p);
  }  
  
  if(newMatrix_p != NULL){
    for(Int_t bIX = 0; bIX < matrixClone_p->GetXaxis()->GetNbins(); ++bIX){
      for(Int_t bIY = 0; bIY < matrixClone_p->GetYaxis()->GetNbins(); ++bIY){
	newMatrix_p->SetBinContent(bIX+1, bIY+1, matrixClone_p->GetBinContent(bIX+1, bIY+1));
	newMatrix_p->SetBinError(bIX+1, bIY+1, 0.0);
      }
    }
  }
  
  
  RooUnfoldResponse* rooRes_p = NULL;
  RooUnfoldBayes* rooBayes_p = NULL;
  
  macroHistToSubsetHistX(matrixClone_p, reco_p, true);
  macroHistToSubsetHistY(matrixClone_p, gen_p, true);

  TH1D* tempUnfold_h = NULL;

  /*
  rooRes_p = new RooUnfoldResponse(reco_p, gen_p, matrixClone_p);
  rooBayes_p = new RooUnfoldBayes(rooRes_p, data_p, 3, false, "temp");
  rooBayes_p->SetVerbose(-1);
  rooBayes_p->SetNToys(1000);
  tempUnfold_h = (TH1D*)rooBayes_p->Hreco(RooUnfold::kCovToy);
  setPrior(matrixClone_p, tempUnfold_h);
  
  delete rooRes_p;
  delete rooBayes_p;
  delete tempUnfold_h;
  */

  //  std::cout << "RUNNING" << std::endl;
  for(Int_t bI = 0; bI < nBayes; ++bI){   
    TH2D* matrixClone2_p = (TH2D*)matrixClone_p->Clone("matrixClone2_p");    
    macroHistToSubsetHistX(matrixClone2_p, reco_p, true);
    macroHistToSubsetHistY(matrixClone2_p, gen_p, true);
        
    rooRes_p = new RooUnfoldResponse(reco_p, gen_p, matrixClone2_p);
    rooBayes_p = new RooUnfoldBayes(rooRes_p, data_p, bayesVals[bI], false, "temp");
    rooBayes_p->SetVerbose(-1);
    //    rooBayes_p->SetNToys(1000);

    inWatch->start();
    tempUnfold_h = (TH1D*)rooBayes_p->Hreco(RooUnfold::kCovToy);
    inWatch->stop();
    
    //    chi2[bI] = rooBayes_p->getChi2();
    if(true/*bI == 0*/) chi2[bI] = 100.;
    else{
      Double_t chi2Val;
      Int_t npdf, igood;
      chi2[bI] = tempUnfold_h->Chi2TestX(unfolded_p[bI-1], chi2Val, npdf, igood, "WW");
    }
    
    unfolded_p[bI] = new TH1D((tagStr + "_unfoldedBayes" + std::to_string(bayesVals[bI]) + "_h").c_str(), ";pT;Counts", nGenBins, genBins);

    for(Int_t tI = 0; tI < tempUnfold_h->GetNbinsX(); ++tI){
      unfolded_p[bI]->SetBinContent(tI+1, tempUnfold_h->GetBinContent(tI+1));
      unfolded_p[bI]->SetBinError(tI+1, tempUnfold_h->GetBinError(tI+1));
    }

    delete matrixClone2_p;
    delete tempUnfold_h;
    delete rooBayes_p;
    delete rooRes_p;
  }   

  delete gen_p;
  delete reco_p;
  delete matrixClone_p;

  return;
}

void ratioPrint(TH1D* histNum_p, TH1D* histDenom_p, bool doNorm=false)
{
  std::cout << "RATIO OF " << histNum_p->GetName() << "/" << histDenom_p->GetName() << ": ";
  for(Int_t bIX = 0; bIX < histNum_p->GetXaxis()->GetNbins(); ++bIX){
    if(doNorm) std::cout << (histNum_p->GetBinContent(bIX+1)*histDenom_p->Integral())/(histDenom_p->GetBinContent(bIX+1)*histNum_p->Integral()) << ",";
    else std::cout << (histNum_p->GetBinContent(bIX+1))/(histDenom_p->GetBinContent(bIX+1)) << ",";
  }
  std::cout << std::endl;

  return;
}

void th2ToTH1(TH2D* th2_p, TH1D** th1_p, bool doYAxis)
{
  const Int_t nMaxBins = 100;
  Int_t nBins = th2_p->GetXaxis()->GetNbins();
  if(doYAxis) nBins = th2_p->GetYaxis()->GetNbins();

  Double_t bins[nMaxBins+1];
  for(Int_t bI = 0; bI < nBins+1; ++bI){
    if(doYAxis) bins[bI] = th2_p->GetYaxis()->GetBinLowEdge(bI+1);
    else bins[bI] = th2_p->GetXaxis()->GetBinLowEdge(bI+1);
  }

  std::string name = th2_p->GetName();
  if(doYAxis) name = name + "_YAxis";
  else name = name + "_XAxis";
  (*th1_p) = new TH1D(name.c_str(), "", nBins, bins);

  if(doYAxis) macroHistToSubsetHistY(th2_p, *th1_p, true);
  else macroHistToSubsetHistX(th2_p, *th1_p, true);

  return;
}

void ratioPrint(TH2D* histNum_p, TH1D* histDenom_p, bool doYAxisNum, bool doNorm=false)
{
  TH1D* histNum2_p = NULL;
  th2ToTH1(histNum_p, &histNum2_p, doYAxisNum);
  ratioPrint(histNum2_p, histDenom_p, doNorm);
  delete histNum2_p;

  return;
}

void relErrPrint(TH1D* hist_p)
{
  std::cout << "REL ERR OF " << hist_p->GetName() << ": ";
  for(Int_t bIX = 0; bIX < hist_p->GetXaxis()->GetNbins(); ++bIX){
    if(hist_p->GetBinContent(bIX+1) <= TMath::Power(10, -20)) std::cout << "0,";
    else std::cout << hist_p->GetBinError(bIX+1)/hist_p->GetBinContent(bIX+1) << ",";
  }
  std::cout << std::endl;

  return;
}


std::string lowerToUpper(std::string inStr)
{
  const std::string lower = "abcdefghijklmnopqrstuvwxyz";
  const std::string upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

  for(unsigned int i = 0; i < inStr.size(); ++i){
    if(lower.find(inStr.substr(i, 1)) != std::string::npos){
      inStr.replace(i, 1, upper.substr(lower.find(inStr.substr(i, 1)), 1));
    }
  }
  
  return inStr;
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


int toyATLASCMS_Unfold(const std::string inFileName, const std::string startTag, bool doClean, const std::string inFileName2 = "", std::string newGenBins = "", std::string newRecoBins = "", std::string newReducedBins = "", std::string inDataFile = "") 
{
  cppWatch total;
  total.start();

  cppWatch unfoldTime;
  cppWatch unfoldTimeRand;
  cppWatch justUnfoldTimeRand;
  cppWatch unfoldTimeNonRand;
  cppWatch justUnfoldTimeNonRand;

  kirchnerPalette kPal;

  Int_t nCol = 4;
  const Int_t nStyle = 5;
  Int_t styles[nStyle] = {24, 25, 28, 27, 46};


  std::cout << "LINE: " << __LINE__ << std::endl;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  Long64_t nFill = -1;
  std::string nFillStr = inFileName.substr(0, inFileName.find("_Seed"));
  nFillStr.replace(0, nFillStr.find("NFill")+5, "");
  nFill = std::stol(nFillStr);

  std::cout << "LINE: " << __LINE__ << std::endl;

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  bool doFile2 = inFileName2.size() != 0 && checkFile(inFileName2);

  std::string file2Str = "NoFile2";
  if(doFile2) file2Str = "YesFile2";

  std::string cleanStr = "Clean";
  if(!doClean) cleanStr = "NoClean";

  TFile* inFile2_p = NULL;
  if(doFile2) inFile2_p = new TFile(inFileName2.c_str(), "READ");

  const Int_t nPowers = 2;//3;
  Double_t powers[nPowers] = {7., 8.};

  const Int_t nRandom = 1000;
  const Int_t nTarget = 1;
  const Double_t targetPoint = 500.;
  const Double_t errTarget = 200.;
  //  const Double_t errTarget = 200./20.;

  //  if(startTag.find("R10") != std::string::npos) ++powers[0];
  
  TH2D* response_h[nPowers];
  TH2D* response2_h[nPowers];
  TH2D* response_NewPrior_h[nPowers][nPowers];
  TH1D* genSpectraPower_p[nPowers];
  //  TH1D* genSpectraTruncPower_p[nPowers];
  //  TH1D* recoSpectraPower_p[nPowers];
  TH2D* recoGenSpectraTruncPower_p[nPowers];
  //  TH2D* recoGenSpectraTruncPower_Random_p[nPowers][nRandom][nTarget];

  
  TH2D* response_Rebin_h[nPowers];
  TH1D* genSpectraPower_Rebin_p[nPowers];
  TH1D* genSpectraTruncPower_Rebin_p[nPowers];
  TH1D* recoSpectraPower_Rebin_p[nPowers];
  TH2D* recoGenSpectraTruncPower_Rebin_p[nPowers];
  TH2D* recoGenSpectraTruncPower_Random_Rebin_p[nPowers][nRandom][nTarget];
  TH1D* genSpectraTruncPower_Random_Rebin_p[nPowers][nRandom][nTarget];
  TH1D* recoSpectraPower_Random_Rebin_p[nPowers][nRandom][nTarget];
    
  std::string startTag2 = lowerToUpper(startTag.substr(0, startTag.find("R"))) + startTag.substr(startTag.find("R"), startTag.size());;
  
  for(Int_t pI = 0; pI < nPowers; ++pI){
    std::string powStr = "Power" + prettyString(powers[pI], 1, true);
    response_h[pI] = (TH2D*)inFile_p->Get((startTag + "Response" + powStr + "_h").c_str());      
    response2_h[pI] = NULL;
    if(doFile2) response2_h[pI] = (TH2D*)inFile2_p->Get((startTag + "Response" + powStr + "_h").c_str());      

    genSpectraPower_p[pI] = (TH1D*)inFile_p->Get(("genSpectra" + powStr + "_h").c_str());
    //    genSpectraTruncPower_p[pI] = (TH1D*)inFile_p->Get(("genSpectraTrunc" + powStr + "_h").c_str());
    //    recoSpectraPower_p[pI] = (TH1D*)inFile_p->Get(("recoSpectra" + powStr + "_h").c_str());
    recoGenSpectraTruncPower_p[pI] = (TH2D*)inFile_p->Get(("recoGenSpectra" + startTag2 + "Trunc" + powStr + "_h").c_str());

    for(Int_t rI = 0; rI < nRandom; ++rI){
      for(Int_t tI = 0; tI < nTarget; ++tI){
	//	recoGenSpectraTruncPower_Random_p[pI][rI][tI] = (TH2D*)inFile_p->Get(("recoGenSpectra" + startTag2 + "Trunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_h").c_str());
      }
    }
  }

  const Int_t nMaxBins = 1000;
  Int_t nGenBins = genSpectraPower_p[0]->GetXaxis()->GetNbins();
  Double_t genBins[nMaxBins];
  for(Int_t bI = 0; bI < genSpectraPower_p[0]->GetXaxis()->GetNbins()+1; ++bI){
    genBins[bI] = genSpectraPower_p[0]->GetXaxis()->GetBinLowEdge(bI+1);
  }

  Int_t nRecoBins = response_h[0]->GetXaxis()->GetNbins();
  Double_t recoBins[nMaxBins];
  for(Int_t bI = 0; bI < response_h[0]->GetXaxis()->GetNbins()+1; ++bI){
    recoBins[bI] = response_h[0]->GetXaxis()->GetBinLowEdge(bI+1);
  }

  std::vector<double> binsVect = binsStrToVect(newGenBins);
  if(binsVect.size() != 0){
    std::cout << "BINSVECT: ";
    for(unsigned int bI = 0; bI < binsVect.size(); ++bI){
      std::cout << binsVect[bI] << ", ";
    }
    std::cout << std::endl;

    std::cout << "NGENBINS: ";
    for(int bI = 0; bI < nGenBins+1; ++bI){
      std::cout << genBins[bI] << ", ";
    }
    std::cout << std::endl;

    if(!replaceBins(binsVect, &nGenBins, genBins)){
      std::cout << "GIVEN nGenBins doesnt represent a subset of macrobins. return 1" << std::endl;
      return 1;
    }
  }

  binsVect = binsStrToVect(newRecoBins);
  if(binsVect.size() != 0){
    if(!replaceBins(binsVect, &nRecoBins, recoBins)){
      std::cout << "GIVEN nRecoBins doesnt represent a subset of macrobins. return 1" << std::endl;
      return 1;
    }
  }

  /*
  newGenBins = newGenBins + ",";
  while(newGenBins.find(",,") != std::string::npos){
    newGenBins.replace(newGenBins.find(",,"), 2, ",");
  }
  
  if(newGenBins.size() != 1){
    std::vector<int> newGenBinsVect;

    while(newGenBins.find(",") != std::string::npos){
      newGenBinsVect.push_back(std::stod(newGenBins.substr(0, newGenBins.find(","))));
      newGenBins.replace(0, newGenBins.find(",")+1, "");
    }

    bool allBinsGood = true;

    for(unsigned int gI = 0; gI < newGenBinsVect.size(); ++gI){
      bool binIsGood = false;
      
      for(Int_t gI2 = 0; gI2 < nGenBins+1; ++gI2){
	if(TMath::Abs(newGenBinsVect[gI] - genBins[gI2]) < 1.){
	  binIsGood = true;
	  break;
	}
      }

      if(!binIsGood){
	allBinsGood = false;
	break;
      }
    }

    if(!allBinsGood){
      std::cout << "INPUT NEW BINS ARE INVALID. return 1" << std::endl;
      return 1;
    }
    else{
      nGenBins = newGenBinsVect.size()-1;
      for(unsigned int gI = 0; gI < newGenBinsVect.size(); ++gI){
	genBins[gI] = newGenBinsVect[gI];
      }
    }
  }
  
  newRecoBins = newRecoBins + ",";
  while(newRecoBins.find(",,") != std::string::npos){
    newRecoBins.replace(newRecoBins.find(",,"), 2, ",");
  }
  
  if(newRecoBins.size() != 1){
    std::vector<int> newRecoBinsVect;

    while(newRecoBins.find(",") != std::string::npos){
      newRecoBinsVect.push_back(std::stod(newRecoBins.substr(0, newRecoBins.find(","))));
      newRecoBins.replace(0, newRecoBins.find(",")+1, "");
    }

    bool allBinsGood = true;

    for(unsigned int gI = 0; gI < newRecoBinsVect.size(); ++gI){
      bool binIsGood = false;
      
      for(Int_t gI2 = 0; gI2 < nRecoBins+1; ++gI2){
	if(TMath::Abs(newRecoBinsVect[gI] - recoBins[gI2]) < 1.){
	  binIsGood = true;
	  break;
	}
      }

      if(!binIsGood){
	allBinsGood = false;
	break;
      }
    }

    if(!allBinsGood){
      std::cout << "INPUT NEW BINS ARE INVALID. return 1" << std::endl;
      return 1;
    }
    else{
      nRecoBins = newRecoBinsVect.size()-1;
      for(unsigned int gI = 0; gI < newRecoBinsVect.size(); ++gI){
	recoBins[gI] = newRecoBinsVect[gI];
      }
    }
  }
  */

  Int_t nGeneralBins = 0;
  Double_t generalBins[nMaxBins];
  for(Int_t bI = 0; bI < nGenBins+1; ++bI){    
    if(genBins[bI] < recoBins[0]) continue;
    if(genBins[bI] > recoBins[nRecoBins]) continue;
    
    generalBins[nGeneralBins] = genBins[bI];
    ++nGeneralBins;
  }
  --nGeneralBins;

  binsVect = binsStrToVect(newReducedBins);
  if(binsVect.size() != 0){
    if(!replaceBins(binsVect, &nGeneralBins, generalBins)){
      std::cout << "GIVEN nGeneralBins doesnt represent a subset of macrobins. return 1" << std::endl;
      return 1;
    }
  }

  Int_t nRecoBins2 = TMath::Min(10, nRecoBins);
  Double_t recoBins2[nMaxBins+1];
  for(Int_t rI = 0; rI < nRecoBins2+1; ++rI){
    recoBins2[rI] = recoBins[rI];
  }

  Int_t nGenBins2 = TMath::Min(20, nGenBins);
  Double_t genBins2[nMaxBins+1];
  for(Int_t rI = 0; rI < nGenBins2+1; ++rI){
    genBins2[rI] = genBins[rI];
  }
  
  for(Int_t pI = 0; pI < nPowers; ++pI){
    std::string powStr = "Power" + prettyString(powers[pI], 1, true);
    response_Rebin_h[pI] = new TH2D((startTag + "Response" + powStr + "_Rebin_h").c_str(), "", nRecoBins, recoBins, nGenBins, genBins);      
    genSpectraPower_Rebin_p[pI] = new TH1D(("genSpectra" + powStr + "_Rebin_h").c_str(), "", nGenBins, genBins);
    genSpectraTruncPower_Rebin_p[pI] = new TH1D(("genSpectra" + startTag2 + "Trunc" + powStr + "_Rebin_h").c_str(), "", nGenBins, genBins);
    recoGenSpectraTruncPower_Rebin_p[pI] = new TH2D(("recoGenSpectra" + startTag2 + "Trunc" + powStr + "_Rebin_h").c_str(), "", nRecoBins, recoBins, nGenBins, genBins);
    recoSpectraPower_Rebin_p[pI] = new TH1D(("reco" + startTag2 + "Spectra" + powStr + "_Rebin_h").c_str(), "", nRecoBins, recoBins);

    if(!doFile2) macroHistToSubsetHist(response_h[pI], response_Rebin_h[pI], true);
    else macroHistToSubsetHist(response2_h[pI], response_Rebin_h[pI], true);
    
    macroHistToSubsetHist(genSpectraPower_p[pI], genSpectraPower_Rebin_p[pI], true);    
    macroHistToSubsetHist(recoGenSpectraTruncPower_p[pI], recoGenSpectraTruncPower_Rebin_p[pI], true);

    for(Int_t rI = 0; rI < nRandom; ++rI){
      for(Int_t tI = 0; tI < nTarget; ++tI){
	recoGenSpectraTruncPower_Random_Rebin_p[pI][rI][tI] = new TH2D(("recoGenSpectra" + startTag2 + "Trunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_Rebin_h").c_str(), "", nRecoBins, recoBins, nGenBins, genBins);
	
	genSpectraTruncPower_Random_Rebin_p[pI][rI][tI] = new TH1D(("genSpectra" + startTag2 + "Trunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_Rebin_h").c_str(), "", nGenBins, genBins);
	recoSpectraPower_Random_Rebin_p[pI][rI][tI] = new TH1D(("reco" + startTag2 + "Spectra" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_Rebin_h").c_str(), "", nRecoBins, recoBins);

	Double_t totalValTarget = 0.0;
	Double_t totalErrTarget = 0.0;
	Double_t rat = 1.0;
	while(rat > 1./TMath::Sqrt(errTarget*TMath::Power(10, tI))){
	  Double_t valX, valY;
	  response_Rebin_h[pI]->GetRandom2(valX, valY);
	  
	  if(valX > targetPoint){
	    ++totalValTarget;
	    totalErrTarget = TMath::Sqrt(totalErrTarget*totalErrTarget + 1*1);
	  }
	  
	  recoGenSpectraTruncPower_Random_Rebin_p[pI][rI][tI]->Fill(valX, valY);
	  if(totalValTarget <= TMath::Power(10, -20)) rat = 1.0;
          else rat = totalErrTarget/totalValTarget;	  
	}
      }
    }
    
    for(Int_t rI = 0; rI < nRandom; ++rI){
      for(Int_t tI = 0; tI < nTarget; ++tI){
	//	macroHistToSubsetHist(recoGenSpectraTruncPower_Random_p[pI][rI][tI], recoGenSpectraTruncPower_Random_Rebin_p[pI][rI][tI], true);
      }
    }

    if(doClean){
      for(Int_t bIX = 0; bIX < nRecoBins; ++bIX){
	for(Int_t bIY = 0; bIY < nGenBins; ++bIY){
	  Double_t val = recoGenSpectraTruncPower_Rebin_p[pI]->GetBinContent(bIX+1, bIY+1);
	  Double_t err = recoGenSpectraTruncPower_Rebin_p[pI]->GetBinError(bIX+1, bIY+1);
	  if(val <= TMath::Power(10, -20)) continue;

	  if(err/val > 1./TMath::Sqrt(10)){
	    recoGenSpectraTruncPower_Rebin_p[pI]->SetBinContent(bIX+1, bIY+1, 0.0);
	    recoGenSpectraTruncPower_Rebin_p[pI]->SetBinError(bIX+1, bIY+1, 0.0);
	  }

	  for(Int_t rI = 0; rI < nRandom; ++rI){
	    for(Int_t tI = 0; tI < nTarget; ++tI){
	      val = recoGenSpectraTruncPower_Random_Rebin_p[pI][rI][tI]->GetBinContent(bIX+1, bIY+1);
	      err = recoGenSpectraTruncPower_Random_Rebin_p[pI][rI][tI]->GetBinError(bIX+1, bIY+1);
	      if(val <= TMath::Power(10, -20)) continue;
	      
	      if(err/val > 1./TMath::Sqrt(10)){
		recoGenSpectraTruncPower_Random_Rebin_p[pI][rI][tI]->SetBinContent(bIX+1, bIY+1, 0.0);
		recoGenSpectraTruncPower_Random_Rebin_p[pI][rI][tI]->SetBinError(bIX+1, bIY+1, 0.0);
	      }
	    }	    
	  }
	}
      }
    }

    macroHistToSubsetHistY(recoGenSpectraTruncPower_Rebin_p[pI], genSpectraTruncPower_Rebin_p[pI], true);
    macroHistToSubsetHistX(recoGenSpectraTruncPower_Rebin_p[pI], recoSpectraPower_Rebin_p[pI], true);

    for(Int_t rI = 0; rI < nRandom; ++rI){
      for(Int_t tI = 0; tI < nTarget; ++tI){
	macroHistToSubsetHistY(recoGenSpectraTruncPower_Random_Rebin_p[pI][rI][tI], genSpectraTruncPower_Random_Rebin_p[pI][rI][tI], true);
	macroHistToSubsetHistX(recoGenSpectraTruncPower_Random_Rebin_p[pI][rI][tI], recoSpectraPower_Random_Rebin_p[pI][rI][tI], true);
      }
    }

    //    macroHistToSubsetHist(genSpectraTruncPower_p[pI], genSpectraTruncPower_Rebin_p[pI], true);
    //    macroHistToSubsetHist(recoSpectraPower_p[pI], recoSpectraPower_Rebin_p[pI], true);

    /*
    macroHistToSubsetHistY(response_Rebin_h[pI], genSpectraTruncPower_Rebin_p[pI], true);       
    genSpectraTruncPower_Rebin_p[pI]->Print("ALL");

    std::cout << "LINE: " << __LINE__ << std::endl;

    macroHistToSubsetHist(genSpectraTruncPower_p[pI], genSpectraTruncPower_Rebin_p[pI], true);
    

    std::cout << "LINE: " << __LINE__ << std::endl;
    */
    /*
    macroHistToSubsetHistX(response_Rebin_h[pI], recoSpectraPower_Rebin_p[pI], true);       
    macroHistToSubsetHistY(response_Rebin_h[pI], genSpectraTruncPower_Rebin_p[pI], true);       
    */
    std::string powOrigStr = "OrigPower" + prettyString(powers[pI], 1, true);
      
    for(Int_t pI2 = TMath::Max(0, pI-1); pI2 < TMath::Min(nPowers, pI+1); ++pI2){
      std::string powPriorStr = "PriorPower" + prettyString(powers[pI2], 1, true);
      response_NewPrior_h[pI][pI2] = new TH2D((startTag + "Response_NewPrior_" + powOrigStr + "_" + powPriorStr + "_h").c_str(), "", nRecoBins, recoBins, nGenBins, genBins);
    }
  }

  
  TH2D* subsetRes_p = new TH2D("subsetRes_h", ";Reco Normalized to Unity;Gen", nRecoBins2, recoBins2, nGenBins2, genBins2);
  TH2D* subsetResErr_p = new TH2D("subsetResErr_h", ";Reco;Gen", nRecoBins2, recoBins2, nGenBins2, genBins2);
  macroHistToSubsetHist(response_Rebin_h[0], subsetRes_p, true, true, true); 
  macroHistToSubsetHist(response_Rebin_h[0], subsetResErr_p, true, true, true); 

  centerTitles(subsetRes_p);
  centerTitles(subsetResErr_p);
  
  std::string saveStr = startTag + "_Gen" + prettyString(genBins[0],1,true) + "to" + prettyString(genBins[nGenBins], 1, true) + "_Reco" + prettyString(recoBins[0],1,true) + "to" + prettyString(recoBins[nRecoBins], 1, true) + "_nFill" + std::to_string(nFill) + "_" + file2Str + "_" + cleanStr + "_" + dateStr + ".pdf";
  
  for(Int_t bIX = 0; bIX < nRecoBins2; ++bIX){
    Double_t total = 0.0;

    for(Int_t bIY = 0; bIY < nGenBins2; ++bIY){
      total += subsetRes_p->GetBinContent(bIX+1, bIY+1);      
    }

    if(total <= TMath::Power(10, -20)) continue;
    
    for(Int_t bIY = 0; bIY < nGenBins2; ++bIY){
      Double_t val = subsetRes_p->GetBinContent(bIX+1, bIY+1)/total;
      Double_t err = subsetRes_p->GetBinError(bIX+1, bIY+1)/total;

      if(doClean && err/val > 1./TMath::Sqrt(10)){
	val = 0.0;
      }

      subsetRes_p->SetBinContent(bIX+1, bIY+1, val);
      subsetRes_p->SetBinError(bIX+1, bIY+1, err);

      if(val <= TMath::Power(10, -20)){
	subsetResErr_p->SetBinContent(bIX+1, bIY+1, 0.0);
	subsetResErr_p->SetBinError(bIX+1, bIY+1, 0.0);
      }
      else{
	subsetResErr_p->SetBinContent(bIX+1, bIY+1, err/val);
	subsetResErr_p->SetBinError(bIX+1, bIY+1, 0.0);
      }
    }
  }

  if(doClean){
    for(Int_t pI = 0; pI < nPowers; ++pI){
      for(Int_t bIX = 0; bIX < nRecoBins; ++bIX){
	for(Int_t bIY = 0; bIY < nGenBins; ++bIY){
	  
	  Double_t val = response_Rebin_h[pI]->GetBinContent(bIX+1, bIY+1);
	  Double_t err = response_Rebin_h[pI]->GetBinError(bIX+1, bIY+1);

	  if(val <= TMath::Power(10, -20)) continue;
	  
	  if(err/val > 1./TMath::Sqrt(10)){
	    response_Rebin_h[pI]->SetBinContent(bIX+1, bIY+1, 0.0);
	    response_Rebin_h[pI]->SetBinError(bIX+1, bIY+1, 0.0);
	  }
	}
      }
    }
  }
  
  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  checkMakeDir("pdfDir/" + dateStr + "/nFill" + std::to_string(nFill));

  TCanvas* canv_p = new TCanvas("canv_p", "canv_p", 450/**2*/, 450);
  //  canv_p->Divide(2, 1);
  //  canv_p->cd(1);
  
  subsetRes_p->SetMaximum(1.0);
  subsetRes_p->SetMinimum(0.0);
  subsetRes_p->DrawCopy("COL TEXT");
  gStyle->SetPaintTextFormat("1.3f");
  gPad->SetLogx();
  gStyle->SetOptStat(0);
  /*
  canv_p->cd(2);
  subsetResErr_p->SetMaximum(1.0);
  subsetResErr_p->SetMinimum(0.0);
  subsetResErr_p->DrawCopy("COL TEXT");
  gStyle->SetPaintTextFormat("1.3f");
  gPad->SetLogx();
  gStyle->SetOptStat(0);
  */
  quietSaveAs(canv_p, "pdfDir/" + dateStr + "/nFill" + std::to_string(nFill) + "/subsetRes_" + saveStr);
  delete canv_p;
  
  delete subsetRes_p;
  delete subsetResErr_p;

  /*  
  for(Int_t pI = 0; pI < nPowers; ++pI){
    std::string powOrigStr = "OrigPower" + prettyString(powers[pI], 1, true);

    for(Int_t pI2 = 0; pI2 < nPowers; ++pI2){
      std::string powPriorStr = "PriorPower" + prettyString(powers[pI2], 1, true);
      response_NewPrior_h[pI][pI2] = new TH2D((startTag + "Response_NewPrior_" + powOrigStr + "_" + powPriorStr + "_h").c_str(), "", nRecoBins, recoBins, nGenBins, genBins);
    }
  }
  */

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  checkMakeDir("pdfDir/" + dateStr + "/nFill" + std::to_string(nFill));

  std::cout << "STATS COMPARE" << std::endl;

  for(Int_t pI = 0; pI < nPowers; ++pI){
    TCanvas* canv_p = new TCanvas("canv_p", "canv_p", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.16);
    canv_p->SetBottomMargin(0.16);
    
    TH1D* tempHist_p = new TH1D("tempHist_h", (";" + std::string(recoSpectraPower_Rebin_p[pI]->GetXaxis()->GetTitle()) + ";" + std::string(recoSpectraPower_Rebin_p[pI]->GetYaxis()->GetTitle())).c_str(), nRecoBins, recoBins);
    tempHist_p->SetMarkerStyle(styles[0%nStyle]);
    tempHist_p->SetMarkerColor(kPal.getColor(0%nCol));
    tempHist_p->SetLineColor(kPal.getColor(0%nCol));

    for(Int_t bIX = 0; bIX < tempHist_p->GetNbinsX(); ++ bIX){
      double val = recoSpectraPower_Rebin_p[pI]->GetBinContent(bIX+1);
      double err = recoSpectraPower_Rebin_p[pI]->GetBinError(bIX+1);

      std::cout << "pI, bIX: " << pI << ", " << bIX << ": " << val << ", " << err << std::endl;
      tempHist_p->SetBinContent(bIX+1, err/val);
      tempHist_p->SetBinError(bIX+1, 0.0);
    }
    tempHist_p->SetMaximum(1.);
    tempHist_p->SetMinimum(0.00001);
    tempHist_p->GetYaxis()->SetTitle("Rel. Stat. Error");
    tempHist_p->GetXaxis()->SetTitle("Jet Reco. p_{T}");
    centerTitles(tempHist_p);
    tempHist_p->DrawCopy("HIST E1");
    
    for(Int_t rI = 0; rI < nRandom; ++rI){
      for(Int_t tI = 0; tI < nTarget; ++tI){
	for(Int_t bIX = 0; bIX < tempHist_p->GetNbinsX(); ++ bIX){
	  double val = recoSpectraPower_Random_Rebin_p[pI][rI][tI]->GetBinContent(bIX+1);
	  double err = recoSpectraPower_Random_Rebin_p[pI][rI][tI]->GetBinError(bIX+1);
	  tempHist_p->SetBinContent(bIX+1, err/val);
	  tempHist_p->SetBinError(bIX+1, 0.0);
	}
	
	tempHist_p->SetMarkerStyle(styles[(rI+1)%nStyle]);
	tempHist_p->SetMarkerColor(kPal.getColor((rI+1)%nCol));
	tempHist_p->SetLineColor(kPal.getColor((rI+1)%nCol));
	
	tempHist_p->DrawCopy("HIST E1 SAME");
      }
    }

    if(inDataFile.size() != 0){
      std::string histName = "akCs10PU3PFFlowJetAnalyzer/jtPtRaw_General_akCs10PU3PFFlowJetAnalyzer_PbPb_Cent0to10_LightMUAndCHID_AbsEta0p0to2p0_h";

      if(inDataFile.find("Cs4") != std::string::npos){
	while(histName.find("Cs10") != std::string::npos){
	  histName.replace(histName.find("Cs10"), 4, "Cs4");
	}
      }
      
      TFile* dataFile_p = new TFile(inDataFile.c_str(), "READ");
      TH1D* tempHist2_p = (TH1D*)dataFile_p->Get(histName.c_str());
      macroHistToSubsetHist(tempHist2_p, tempHist_p);

      tempHist_p->SetMarkerColor(1);
      tempHist_p->SetMarkerStyle(20);
      tempHist_p->SetLineColor(1);

      for(Int_t bIX = 0; bIX < tempHist_p->GetNbinsX(); ++ bIX){
        double val = tempHist_p->GetBinContent(bIX+1);
        double err = tempHist_p->GetBinError(bIX+1);
	tempHist_p->SetBinContent(bIX+1, err/val);
        tempHist_p->SetBinError(bIX+1, 0.0);
      }

      tempHist_p->DrawCopy("*HIST E1 SAME");
      
      dataFile_p->Close();
      delete dataFile_p;
    }
  
    gPad->SetLogy();
    quietSaveAs(canv_p, "pdfDir/" + dateStr + "/nFill" + std::to_string(nFill) + "/statsCompare_Power" + prettyString(powers[pI], 1, true) + "_" + dateStr + ".pdf");
    delete canv_p;
    delete tempHist_p;
  }
  
  std::cout << "STARTED UNFOLDING" << std::endl;
  //Unfolding testing
  const Int_t nBayes = 5;
  Int_t bayesVals[nBayes] = {1,2,4,10,100};

  //ONE FOR UNFOLDING MATRIX YADA YADA
  TH1D* powerUnfolds[nPowers][nPowers][nBayes];
  TH1D* powerUnfoldsGenRats[nPowers][nPowers][nBayes];
  
  Double_t chi2[nPowers][nPowers][nBayes];
  
  TH1D* powerUnfoldsRandom[nPowers][nPowers][nRandom][nTarget][nBayes];
  TH1D* powerUnfoldsGenRatsRandom[nPowers][nPowers][nRandom][nTarget][nBayes];
  TH1D* powerUnfoldsBinDeltasRandom[nPowers][nPowers][nTarget][nBayes][nMaxBins];
  Double_t chi2Random[nPowers][nPowers][nRandom][nTarget][nBayes];

  for(Int_t pI = 0; pI < nPowers; ++pI){
    std::cout << " Unfold power " << powers[pI] << std::endl;
    for(Int_t pI2 = TMath::Max(pI-1, 0); pI2 < TMath::Min(pI+1, nPowers); ++pI2){
      
      if(powers[pI] == 8 && powers[pI2] == 7 && false){
	genSpectraTruncPower_Rebin_p[pI]->SetBinContent(1, 0.0);
	genSpectraTruncPower_Rebin_p[pI]->SetBinContent(2, 0.000784);
	genSpectraTruncPower_Rebin_p[pI]->SetBinContent(3, 0.003745);
	genSpectraTruncPower_Rebin_p[pI]->SetBinContent(4, 0.009974);
	genSpectraTruncPower_Rebin_p[pI]->SetBinContent(5, 0.006371);
	genSpectraTruncPower_Rebin_p[pI]->SetBinContent(6, 0.001784);
	genSpectraTruncPower_Rebin_p[pI]->SetBinContent(7, 0.000517);
	genSpectraTruncPower_Rebin_p[pI]->SetBinContent(8, 0.000249);
	genSpectraTruncPower_Rebin_p[pI]->SetBinContent(9, 0.000002);
      }
      

      unfoldTime.start();
      unfoldTimeNonRand.start();
      doUnfold(response_Rebin_h[pI2], genSpectraTruncPower_Rebin_p[pI], recoSpectraPower_Rebin_p[pI2], nBayes, bayesVals, startTag + "UnfoldPower" + prettyString(powers[pI], 1, true) + "_SpectPower" + prettyString(powers[pI2], 1, true), powerUnfolds[pI][pI2], chi2[pI][pI2], response_NewPrior_h[pI][pI2], &justUnfoldTimeNonRand); 
      unfoldTime.stop();
      unfoldTimeNonRand.stop();

      for(Int_t rI = 0; rI < nRandom; ++rI){
	for(Int_t tI = 0; tI < nTarget; ++tI){
	  if(rI == 0){
	    TH1D* genHighStat_p = new TH1D("genHighStat_h", "", nGenBins, genBins);	
	    macroHistToSubsetHist(genSpectraTruncPower_Rebin_p[pI], genHighStat_p, true);
	    
	    TH1D* recoHighStat_p = new TH1D("recoHighStat_h", "", nRecoBins, recoBins);	
	    macroHistToSubsetHist(recoSpectraPower_Rebin_p[pI], recoHighStat_p, true);
	    
	    TH1D* genLowerStat_p = new TH1D("genLowerStat_h", "", nGenBins, genBins);	
	    macroHistToSubsetHist(genSpectraTruncPower_Random_Rebin_p[pI2][rI][tI], genLowerStat_p, true);
	    
	    TH1D* recoLowerStat_p = new TH1D("recoLowerStat_h", "", nRecoBins, recoBins);	
	    macroHistToSubsetHist(recoSpectraPower_Random_Rebin_p[pI2][rI][tI], recoLowerStat_p, true);
	    
	    
	    Double_t scaleFactor = recoHighStat_p->Integral()/recoLowerStat_p->Integral();
	    genLowerStat_p->Scale(scaleFactor);
	    recoLowerStat_p->Scale(scaleFactor);
	    
	    divHistByWidth(recoLowerStat_p);
	    divHistByWidth(recoHighStat_p);
	    
	    divHistByWidth(genLowerStat_p);
	    divHistByWidth(genHighStat_p);
	    
	    Double_t par1Guess = 7;
	    Double_t par0Guess = recoLowerStat_p->GetBinContent(1)/TMath::Power(1./recoLowerStat_p->GetBinCenter(1), par1Guess);
	    
	    TF1* fitRecoLowerStat_p = new TF1("fitRecoLowerStat_p", "[0]*TMath::Power(1./x, [1])", recoBins[0], recoBins[nRecoBins]);
	    fitRecoLowerStat_p->SetParameter(0, par0Guess);
	    fitRecoLowerStat_p->SetParameter(1, par1Guess);
	    
	    
	    TF1* fitRecoHighStat_p = new TF1("fitRecoHighStat_p", "[0]*TMath::Power(1./x, [1])", recoBins[0], recoBins[nRecoBins]);
	    fitRecoHighStat_p->SetParameter(0, par0Guess);
	    fitRecoHighStat_p->SetParameter(1, par1Guess);
	    recoHighStat_p->Fit("fitRecoHighStat_p", "M", "", recoBins[0], recoBins[nRecoBins]);
	    
	    Int_t startBinPos = -1;
	    Int_t endBinPos = -1;
	    
	    for(Int_t bI = 0; bI < nGenBins+1; ++bI){
	      if(TMath::Abs(genBins[bI] - recoBins[0]) < 1.) startBinPos = bI;	 
	      if(TMath::Abs(genBins[bI] - recoBins[nRecoBins]) < 1.) endBinPos = bI;	 
	    }
	    
	    TF1* fitGenHighStat_p = new TF1("fitGenHighStat_p", "[0]*TMath::Power(1./x, [1])", genBins[startBinPos+1], genBins[endBinPos]);
	    fitGenHighStat_p->SetParameter(0, par0Guess);
	    fitGenHighStat_p->SetParameter(1, par1Guess);      
	  
	    genHighStat_p->Fit("fitGenHighStat_p", "M", "", genBins[startBinPos+1], genBins[endBinPos]);
	    
	    TF1* fitGenLowerStat_p = new TF1("fitGenLowerStat_p", "[0]*TMath::Power(1./x, [1])", genBins[startBinPos+1], genBins[endBinPos]);
	    fitGenLowerStat_p->SetParameter(0, par0Guess);
	    fitGenLowerStat_p->SetParameter(1, par1Guess);      
	    
	    genLowerStat_p->Fit("fitGenLowerStat_p", "M", "", genBins[startBinPos+1], genBins[endBinPos]);
	  	    
	    TCanvas* canv_p = new TCanvas("canv_p", "canv_p", 450, 450);
	    canv_p->SetTopMargin(0.01);
	    canv_p->SetRightMargin(0.01);
	    
	    Double_t min = TMath::Min(getMinGTZero(recoLowerStat_p), getMinGTZero(recoHighStat_p));
	    min = TMath::Min(min, getMinGTZero(genLowerStat_p));
	    min = TMath::Min(min, getMinGTZero(genHighStat_p));
	    
	    Double_t max = TMath::Max(recoLowerStat_p->GetMaximum(), recoHighStat_p->GetMaximum());
	    max = TMath::Max(max, genLowerStat_p->GetMaximum());
	    max = TMath::Max(max, genHighStat_p->GetMaximum());
	    
	    std::cout << "MINMAX: " << min << ", " << max << std::endl;
	    genLowerStat_p->SetMinimum(min/5.);
	    genLowerStat_p->SetMaximum(max*5.);
	    
	    genLowerStat_p->SetMarkerColor(4);
	    genLowerStat_p->SetLineColor(4);
	    genLowerStat_p->SetLineStyle(2);
	    genLowerStat_p->DrawCopy("HIST E1");
	    
	    recoLowerStat_p->SetMarkerColor(2);
	    recoLowerStat_p->SetLineColor(2);
	    recoLowerStat_p->SetLineStyle(2);
	    recoLowerStat_p->DrawCopy("HIST E1 SAME");
	    
	    genHighStat_p->SetMarkerColor(4);
	    genHighStat_p->SetLineColor(4);
	    genHighStat_p->DrawCopy("HIST E1 SAME");
	    
	    recoHighStat_p->SetMarkerColor(2);
	    recoHighStat_p->SetLineColor(2);
	    recoHighStat_p->DrawCopy("HIST E1 SAME");
	    
	    fitRecoLowerStat_p->SetMarkerColor(2);
	    fitRecoLowerStat_p->SetLineColor(2);
	    fitRecoLowerStat_p->SetLineStyle(2);
	    
	    fitGenLowerStat_p->SetMarkerColor(4);
	    fitGenLowerStat_p->SetLineColor(4);
	    fitGenLowerStat_p->SetLineStyle(2);
	    
	    fitRecoHighStat_p->SetMarkerColor(2);
	    fitRecoHighStat_p->SetLineColor(2);
	    
	    fitGenHighStat_p->SetMarkerColor(4);
	    fitGenHighStat_p->SetLineColor(4);
	    
	    fitRecoLowerStat_p->Draw("SAME");
	    fitRecoHighStat_p->Draw("SAME");
	    fitGenLowerStat_p->Draw("SAME");
	    fitGenHighStat_p->Draw("SAME");
	    
	    
	    TLatex* label_p = new TLatex();
	    label_p->SetNDC();
	    label_p->SetTextFont(43);
	    label_p->SetTextSize(10);
	    
	    label_p->DrawLatex(0.7, 0.95, ("GenLower=" + prettyString(fitGenLowerStat_p->GetParameter(1), 2, false)).c_str());
	    label_p->DrawLatex(0.7, 0.89, ("GenHigh=" + prettyString(fitGenHighStat_p->GetParameter(1), 2, false)).c_str());
	    
	    label_p->DrawLatex(0.7, 0.83, ("RecoLower=" + prettyString(fitRecoLowerStat_p->GetParameter(1), 2, false)).c_str());
	    label_p->DrawLatex(0.7, 0.77, ("RecoHigh=" + prettyString(fitRecoHighStat_p->GetParameter(1), 2, false)).c_str());
	    
	    delete label_p;
      
	    gPad->SetLogy();
	    quietSaveAs(canv_p, "pdfDir/" + dateStr + "/" + "nFill" + std::to_string(nFill) + "/fits_UnfoldPower" + prettyString(powers[pI], 1, true) + "_SpectPower" + prettyString(powers[pI2], 1, true) + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_" + dateStr + ".pdf");
	    delete canv_p;

	    
	    delete genHighStat_p;
	    delete recoHighStat_p;
	    delete genLowerStat_p;
	    delete recoLowerStat_p;
	    delete fitRecoLowerStat_p;	  
	    delete fitGenHighStat_p;	  
	    delete fitRecoHighStat_p;	  
	    delete fitGenLowerStat_p;	  
	  }
	  
	  unfoldTime.start();
	  unfoldTimeRand.start();
	  doUnfold(response_Rebin_h[pI2], genSpectraTruncPower_Rebin_p[pI], recoSpectraPower_Random_Rebin_p[pI2][rI][tI], nBayes, bayesVals, startTag + "UnfoldPower" + prettyString(powers[pI], 1, true) + "_SpectPower" + prettyString(powers[pI2], 1, true) + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI), powerUnfoldsRandom[pI][pI2][rI][tI], chi2Random[pI][pI2][rI][tI], NULL, &justUnfoldTimeRand);
	  unfoldTime.stop();
	  unfoldTimeRand.stop();

	}
      }
      
      for(Int_t bI = 0; bI < nBayes; ++bI){
	powerUnfoldsGenRats[pI][pI2][bI] = new TH1D((startTag + "PowerUnfoldsGenRats_UnfoldPow" + prettyString(powers[pI], 1, true) + "_SpectPow" + prettyString(powers[pI2], 1, true) + "_Bayes" + std::to_string(bayesVals[bI]) + "_h").c_str(), ";Jet pT;Unfolded/Gen", nGeneralBins, generalBins);
	centerTitles(powerUnfoldsGenRats[pI][pI2][bI]);
	
	TH1D* tempHist_p = new TH1D("tempHist_p", "", nGeneralBins, generalBins);
	TH1D* tempHist2_p = new TH1D("tempHist2_p", "", nGeneralBins, generalBins);
	macroHistToSubsetHist(powerUnfolds[pI][pI2][bI], tempHist_p);
	macroHistToSubsetHist(genSpectraTruncPower_Rebin_p[pI2], tempHist2_p);
	
	for(Int_t gI = 0; gI < nGeneralBins; ++gI){
	  Double_t val = 1.0;
	  if(tempHist_p->GetBinContent(gI+1) <= TMath::Power(10,-20)){
	    if(tempHist_p->GetBinContent(gI+1) >= TMath::Power(10,-20)){
	      val = 0.0;
	    }
	  }
	  else{
	    val = tempHist_p->GetBinContent(gI+1)/tempHist2_p->GetBinContent(gI+1);
	  }
	  
	  powerUnfoldsGenRats[pI][pI2][bI]->SetBinContent(gI+1, val);
	  powerUnfoldsGenRats[pI][pI2][bI]->SetBinError(gI+1, 0.0);
	}

	delete tempHist_p;
	delete tempHist2_p;
      }

      for(Int_t rI = 0; rI < nRandom; ++rI){
	for(Int_t tI = 0; tI < nTarget; ++tI){
	  for(Int_t bI = 0; bI < nBayes; ++bI){
	    powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI] = new TH1D((startTag + "PowerUnfoldsGenRats_UnfoldPow" + prettyString(powers[pI], 1, true) + "_SpectPow" + prettyString(powers[pI2], 1, true) + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_Bayes" + std::to_string(bayesVals[bI]) + "_h").c_str(), ";Jet pT;Unfolded/Gen", nGeneralBins, generalBins);
	    centerTitles(powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]);

	    if(rI == 0){
	      for(Int_t bIX = 0; bIX < nGeneralBins; ++bIX){
		powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX] = new TH1D((startTag + "PowerUnfoldsBinDeltas_UnfoldPow" + prettyString(powers[pI], 1, true) + "_SpectPow" + prettyString(powers[pI2], 1, true) + "_Target" + std::to_string(tI) + "_Bayes" + std::to_string(bayesVals[bI]) + "_GeneralBin" + prettyString(generalBins[bIX], 1, true) + "to" + prettyString(generalBins[bIX+1], 1, true) + "_h").c_str(), ";Reco/Gen;Counts", 25, 0., 2.);
	      }
	    }

	    TH1D* tempHist_p = new TH1D("tempHist_p", "", nGeneralBins, generalBins);
	    TH1D* tempHist2_p = new TH1D("tempHist2_p", "", nGeneralBins, generalBins);
	    macroHistToSubsetHist(powerUnfoldsRandom[pI][pI2][rI][tI][bI], tempHist_p);
	    macroHistToSubsetHist(genSpectraTruncPower_Random_Rebin_p[pI2][rI][tI], tempHist2_p);

	    for(Int_t gI = 0; gI < nGeneralBins; ++gI){
	      Double_t val = 1.0;
	      if(tempHist_p->GetBinContent(gI+1) <= TMath::Power(10,-20)){
		if(tempHist_p->GetBinContent(gI+1) >= TMath::Power(10,-20)){
		  val = 0.0;
		}
	      }
	      else{
		val = tempHist_p->GetBinContent(gI+1)/tempHist2_p->GetBinContent(gI+1);
	      }
	      
	      powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->SetBinContent(gI+1, val);
	      powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->SetBinError(gI+1, 0.0);
	    }

	    delete tempHist_p;
	    delete tempHist2_p;
	  }
	}      
      }
    }
  }

  for(Int_t pI = 0; pI < nPowers; ++pI){
    for(Int_t pI2 = TMath::Max(pI-1, 0); pI2 < TMath::Min(pI+1, nPowers); ++pI2){

      for(Int_t tI = 0; tI < nTarget; ++tI){
	for(Int_t bI = 0; bI < nBayes; ++bI){
	  for(Int_t bIX = 0; bIX < nGeneralBins; ++bIX){

	    for(Int_t rI = 0; rI < nRandom; ++rI){
	      Double_t val = powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->GetBinContent(bIX+1);
	      
	      if(val > generalBins[nGeneralBins]) val = (generalBins[nGeneralBins-1] + generalBins[nGeneralBins])/2.;
	      if(val < 0) val = (generalBins[1] + generalBins[0])/2.;

	      powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->Fill(val);
	    }
	  }
	}
      }
    }
  }

    
  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);

  std::string outFileName = inFileName;
  outFileName.replace(outFileName.find(".root"), 5, "");
  outFileName = outFileName + "_" + startTag;
  outFileName = outFileName + "_Unfold.root";

  while(outFileName.find("/") != std::string::npos){
    outFileName.replace(0, outFileName.find("/")+1, "");
  }
  
  TFile* outFile_p = new TFile(("output/" + dateStr + "/" + outFileName).c_str(), "RECREATE");
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);


  Double_t max = -1;
  Double_t min = 10;

  for(Int_t pI = 0; pI < nPowers; ++pI){
    for(Int_t pI2 = TMath::Max(0, pI-1); pI2 < TMath::Min(pI+1, nPowers); ++pI2){

      for(Int_t bI = 0; bI < nBayes; ++bI){
	if(bayesVals[bI] != 100) continue;

	if(max < getMax(powerUnfoldsGenRats[pI][pI2][bI])) max = getMax(powerUnfoldsGenRats[pI][pI2][bI]);
	if(min >= getMin(powerUnfoldsGenRats[pI][pI2][bI])) min = getMin(powerUnfoldsGenRats[pI][pI2][bI]);
      }
    }
  }

  //  Double_t interval = max-min;
  //  max += interval/10.;
  //  min -= interval/10.;  
  max = 1.5;
  min = 0.5;
    
  for(Int_t pI = 0; pI < nPowers; ++pI){
    response_Rebin_h[pI]->Write("", TObject::kOverwrite);
    
    genSpectraPower_Rebin_p[pI]->Write("", TObject::kOverwrite);
    genSpectraTruncPower_Rebin_p[pI]->Write("", TObject::kOverwrite);
    recoSpectraPower_Rebin_p[pI]->Write("", TObject::kOverwrite);

    for(Int_t rI = 0; rI < nRandom; ++rI){
      for(Int_t tI = 0; tI < nTarget; ++tI){
	genSpectraTruncPower_Random_Rebin_p[pI][rI][tI]->Write("", TObject::kOverwrite);
	recoSpectraPower_Random_Rebin_p[pI][rI][tI]->Write("", TObject::kOverwrite);
      }
    }
      
    if(!doFile2){
      canv_p = new TCanvas("canv_p", "canv_p", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);

      genSpectraPower_Rebin_p[pI]->SetMarkerStyle(24);
      genSpectraPower_Rebin_p[pI]->SetMarkerSize(1);
      genSpectraPower_Rebin_p[pI]->SetMarkerColor(1);
      genSpectraPower_Rebin_p[pI]->SetLineColor(1);
      
      genSpectraTruncPower_Rebin_p[pI]->SetMarkerStyle(25);
      genSpectraTruncPower_Rebin_p[pI]->SetMarkerSize(1);
      genSpectraTruncPower_Rebin_p[pI]->SetMarkerColor(4);
      genSpectraTruncPower_Rebin_p[pI]->SetLineColor(4);
      
      recoSpectraPower_Rebin_p[pI]->SetMarkerStyle(28);
      recoSpectraPower_Rebin_p[pI]->SetMarkerSize(1);
      recoSpectraPower_Rebin_p[pI]->SetMarkerColor(2);
      recoSpectraPower_Rebin_p[pI]->SetLineColor(2);

      TLegend* leg_p = new TLegend(0.45, 0.75, 0.9, 0.95);
      leg_p->SetBorderSize(0);
      leg_p->SetFillStyle(0);
      leg_p->SetTextFont(43);
      leg_p->SetTextSize(16);

      TLatex* label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(43);
      label_p->SetTextSize(10);

      Double_t total = 0.0;
      for(Int_t bIY = 0; bIY < recoGenSpectraTruncPower_Rebin_p[pI]->GetNbinsY(); ++bIY){
	total += recoGenSpectraTruncPower_Rebin_p[pI]->GetBinContent(1, bIY+1);
      }
      
      Int_t pos = -1;
      Double_t total2 = 0.0;
      for(Int_t bIY = 0; bIY < recoGenSpectraTruncPower_Rebin_p[pI]->GetNbinsY(); ++bIY){
	total2 += recoGenSpectraTruncPower_Rebin_p[pI]->GetBinContent(1, bIY+1);
	if(total2/total > 0.01){
	  pos = bIY;
	  break;
	}
      }
      
      Double_t center = recoGenSpectraTruncPower_Rebin_p[pI]->GetYaxis()->GetBinCenter(pos+1);
      Double_t lineMin = genSpectraPower_Rebin_p[pI]->GetMinimum();
      Double_t lineMax = genSpectraPower_Rebin_p[pI]->GetMaximum();
      
      std::cout << "POS: " << pos << ", " << center << std::endl;
      
      const Int_t nRecoSubsetBins = 1;
      Double_t recoSubsetBins[nRecoSubsetBins+1];
      recoSubsetBins[0] = recoBins[0];
      recoSubsetBins[1] = recoBins[1];
      TH2D* subset2D_p = new TH2D("subset2D_p", "", nRecoSubsetBins, recoSubsetBins, nGenBins, genBins);
      macroHistToSubsetHist(recoGenSpectraTruncPower_Rebin_p[pI], subset2D_p, true);
      TH1D* subset1D_p = new TH1D("subset1D_p", "", nGenBins, genBins);
      macroHistToSubsetHistY(subset2D_p, subset1D_p, true);
      
      subset1D_p->SetMarkerColor(4);
      subset1D_p->SetLineColor(4);
          
      genSpectraPower_Rebin_p[pI]->DrawCopy("HIST E1");
      genSpectraTruncPower_Rebin_p[pI]->DrawCopy("HIST E1 SAME");
      recoSpectraPower_Rebin_p[pI]->DrawCopy("HIST E1 SAME");
      subset1D_p->DrawCopy("*HIST E1 SAME");
      
      gPad->SetLogy();
      gPad->SetLogx();
      gStyle->SetOptStat(0);
      
      TLine* line_p = new TLine();    
      line_p->SetLineColor(4);
      line_p->SetLineStyle(2);
      line_p->DrawLine(center, lineMin, center, lineMax);    
      delete line_p;

      leg_p->AddEntry(genSpectraPower_Rebin_p[pI], "Gen.", "P L");
      leg_p->AddEntry(genSpectraTruncPower_Rebin_p[pI], "Gen. w/ Reco.", "P L");
      leg_p->AddEntry(recoSpectraPower_Rebin_p[pI], "Reco.", "P L");
      leg_p->AddEntry(subset1D_p, "Gen. From 1st Reco. Bin", "P L");
      
      leg_p->Draw("SAME");
      label_p->DrawLatex(0.5, 0.65, ("Reco pT > " + prettyString(recoBins[0],1,false)).c_str());

      quietSaveAs(canv_p, "pdfDir/" + dateStr + "/nFill" + std::to_string(nFill) + "/" + genSpectraTruncPower_Rebin_p[pI]->GetName() + "_" + saveStr);
      delete canv_p;

      delete leg_p;
      delete label_p;
      delete subset2D_p;
      delete subset1D_p;
    }

    for(Int_t pI2 = TMath::Max(pI-1,0); pI2 < TMath::Min(pI+1, nPowers); ++pI2){
      response_NewPrior_h[pI][pI2]->Write("", TObject::kOverwrite);

      canv_p = new TCanvas("canv_p", "canv_p", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);

      TLegend* leg_p = new TLegend(0.7, 0.15, 0.9, 0.35);
      leg_p->SetBorderSize(0);
      leg_p->SetFillStyle(0);
      
      for(Int_t bI = 0; bI < nBayes; ++bI){
	powerUnfolds[pI][pI2][bI]->Write("", TObject::kOverwrite);

	if(bayesVals[bI] < 6){
	  powerUnfoldsGenRats[pI][pI2][bI]->SetMaximum(max);
	  powerUnfoldsGenRats[pI][pI2][bI]->SetMinimum(min);

	  powerUnfoldsGenRats[pI][pI2][bI]->SetMarkerStyle(styles[bI%nStyle]);
	  powerUnfoldsGenRats[pI][pI2][bI]->SetMarkerColor(kPal.getColor(bI%nCol));
	  powerUnfoldsGenRats[pI][pI2][bI]->SetLineColor(kPal.getColor(bI%nCol));

	  if(bI != 0) powerUnfoldsGenRats[pI][pI2][bI]->DrawCopy("HIST E1 P SAME");
	  else powerUnfoldsGenRats[pI][pI2][bI]->DrawCopy("HIST E1 P");
	  
	  if(bayesVals[bI] == 3){
	    Double_t mean = 0.0;
	    
	    for(Int_t bIX = 0; bIX < powerUnfoldsGenRats[pI][pI2][bI]->GetNbinsX(); ++bIX){
	      mean += powerUnfoldsGenRats[pI][pI2][bI]->GetBinContent(bIX+1);
	    }
	    mean /= (Double_t)powerUnfoldsGenRats[pI][pI2][bI]->GetNbinsX();
	    
	    Double_t sigma = 0.0;
	    
	    for(Int_t bIX = 0; bIX < powerUnfoldsGenRats[pI][pI2][bI]->GetNbinsX(); ++bIX){
	      Double_t val = powerUnfoldsGenRats[pI][pI2][bI]->GetBinContent(bIX+1) - mean;
	      sigma += val*val;
	    }
	    
	    sigma /= (Double_t)(powerUnfoldsGenRats[pI][pI2][bI]->GetNbinsX() - 1);
	    sigma = TMath::Sqrt(sigma);
	    
	    TLatex* label_p = new TLatex();
	    label_p->SetNDC();
	    label_p->SetTextFont(43);
	    label_p->SetTextSize(10);
	    label_p->DrawLatex(0.6, 0.95, ("#mu=" + prettyString(mean, 3, false)).c_str());
	    label_p->DrawLatex(0.6, 0.89, ("#sigma=" + prettyString(sigma, 3, false)).c_str());
	    Double_t statsRes = response_Rebin_h[pI2]->Integral()*100;
	    Double_t statsReco = recoSpectraPower_Rebin_p[pI2]->Integral()*100;
	    
	    label_p->DrawLatex(0.6, 0.83, ("IntRes=" + prettyStringE(statsRes, 3, false)).c_str());
	    label_p->DrawLatex(0.6, 0.77, ("IntReco=" + prettyStringE(statsReco, 3, false)).c_str());
	    label_p->DrawLatex(0.6, 0.71, ("Reco pT > " + prettyString(recoBins[0], 1, false)).c_str());


	    Double_t startVal = 0.95;
	    Double_t sum1 = 0.0;
	    Double_t sum2 = 0.0;
	    
	    for(Int_t bIX = 0; bIX < powerUnfolds[pI][pI2][bI]->GetNbinsX(); ++bIX){
	      Double_t val1 = powerUnfolds[pI][pI2][bI]->GetBinContent(bIX+1);
	      Double_t val2 = genSpectraTruncPower_Rebin_p[pI2]->GetBinContent(bIX+1);
	      sum1 += val1;
	      sum2 += val2;
	      
	      label_p->DrawLatex(0.16, startVal, (std::to_string(bIX) + ": " + prettyString(val1, 6, false) +  "/" + prettyString(val2, 6, false) + "=" + prettyString(val1/val2, 6, false)).c_str());		  
	      startVal -= 0.03;
	    }
	    
	    label_p->DrawLatex(0.16, startVal, ("SUM: " + prettyString(sum1, 6, false) +  "/" + prettyString(sum2, 6, false) + "=" + prettyString(sum1/sum2, 6, false)).c_str());		  	    
	    delete label_p;
	  }
	  leg_p->AddEntry(powerUnfoldsGenRats[pI][pI2][bI], ("Bayes=" + std::to_string(bayesVals[bI]) + "," + prettyString(chi2[pI][pI2][bI], 3, false)).c_str(), "P L");
	}
      
	powerUnfoldsGenRats[pI][pI2][bI]->Write("", TObject::kOverwrite);
      }

      gStyle->SetOptStat(0);
      leg_p->Draw("SAME");
      
      quietSaveAs(canv_p, "pdfDir/" + dateStr + "/nFill" + std::to_string(nFill) + "/" + powerUnfoldsGenRats[pI][pI2][nBayes-1]->GetName() + "_" + saveStr);	  
      delete canv_p;
      delete leg_p;
    }

    for(Int_t pI2 = TMath::Max(pI-1,0); pI2 < TMath::Min(pI+1, nPowers); ++pI2){
      for(Int_t rI = 0; rI < nRandom; ++rI){
	for(Int_t tI = 0; tI < nTarget; ++tI){
	  for(Int_t bI = 0; bI < nBayes; ++bI){
	    powerUnfoldsRandom[pI][pI2][rI][tI][bI]->Write("", TObject::kOverwrite);
	    powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->Write("", TObject::kOverwrite);

	    if(rI == 0){	  
	      for(Int_t bIX = 0; bIX < nGeneralBins; ++bIX){
		powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->Write("", TObject::kOverwrite);
	      }
	    }
	  }
	}
      }
    }
	    
    for(Int_t pI2 = TMath::Max(pI-1,0); pI2 < TMath::Min(pI+1, nPowers); ++pI2){
      for(Int_t rI = 0; rI < 1; ++rI){
	for(Int_t tI = 0; tI < nTarget; ++tI){
	  canv_p = new TCanvas("canv_p", "canv_p", 450, 450);
	  canv_p->SetTopMargin(0.01);
	  canv_p->SetRightMargin(0.01);
	
	  TLegend* leg_p = new TLegend(0.7, 0.15, 0.9, 0.35);
	  leg_p->SetBorderSize(0);
	  leg_p->SetFillStyle(0);
	  
	  for(Int_t bI = 0; bI < nBayes; ++bI){	    
	    if(bayesVals[bI] < 6){
	      
	      powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->SetMaximum(max);
	      powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->SetMinimum(min);
	      
	      powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->SetMarkerStyle(styles[bI%nStyle]);
	      powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->SetMarkerColor(kPal.getColor(bI%nCol));
	      powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->SetLineColor(kPal.getColor(bI%nCol));
	      
	      if(bI != 0) powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->DrawCopy("HIST E1 P SAME");
	      else powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->DrawCopy("HIST E1 P");
	      
	      if(bayesVals[bI] == 3){
		Double_t mean = 0.0;
		
		for(Int_t bIX = 0; bIX < powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->GetNbinsX(); ++bIX){
		  mean += powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->GetBinContent(bIX+1);
		}
		mean /= (Double_t)powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->GetNbinsX();
		
		Double_t sigma = 0.0;
		
		for(Int_t bIX = 0; bIX < powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->GetNbinsX(); ++bIX){
		  Double_t val = powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->GetBinContent(bIX+1) - mean;
		  sigma += val*val;
		}
		
		sigma /= (Double_t)(powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI]->GetNbinsX() - 1);
		sigma = TMath::Sqrt(sigma);
		
		TLatex* label_p = new TLatex();
		label_p->SetNDC();
		label_p->SetTextFont(43);
		label_p->SetTextSize(10);
		label_p->DrawLatex(0.6, 0.95, ("#mu=" + prettyString(mean, 3, false)).c_str());
		label_p->DrawLatex(0.6, 0.89, ("#sigma=" + prettyString(sigma, 3, false)).c_str());
		Double_t statsRes = response_Rebin_h[pI2]->Integral()*100;
		Double_t statsReco = recoSpectraPower_Rebin_p[pI2]->Integral()*100;
		
		label_p->DrawLatex(0.6, 0.83, ("IntRes=" + prettyStringE(statsRes, 3, false)).c_str());
		label_p->DrawLatex(0.6, 0.77, ("IntReco=" + prettyStringE(statsReco, 3, false)).c_str());
		label_p->DrawLatex(0.6, 0.71, ("Reco pT > " + prettyString(recoBins[0], 1, false)).c_str());
	  
		Double_t startVal = 0.95;
		Double_t sum1 = 0.0;
		Double_t sum2 = 0.0;

		for(Int_t bIX = 0; bIX < powerUnfoldsRandom[pI][pI2][rI][tI][bI]->GetNbinsX(); ++bIX){
		  Double_t val1 = powerUnfoldsRandom[pI][pI2][rI][tI][bI]->GetBinContent(bIX+1);
		  Double_t val2 = genSpectraTruncPower_Random_Rebin_p[pI2][rI][tI]->GetBinContent(bIX+1);
		  sum1 += val1;
		  sum2 += val2;
		  
		  label_p->DrawLatex(0.16, startVal, (std::to_string(bIX) + ": " + prettyString(val1, 6, false) +  "/" + prettyString(val2, 6, false) + "=" + prettyString(val1/val2, 6, false)).c_str());		  
		  startVal -= 0.03;
		}
		
		label_p->DrawLatex(0.16, startVal, ("SUM: " + prettyString(sum1, 6, false) +  "/" + prettyString(sum2, 6, false) + "=" + prettyString(sum1/sum2, 6, false)).c_str());		  

		delete label_p;
	      }
	      leg_p->AddEntry(powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI], ("Bayes=" + std::to_string(bayesVals[bI])  + "," + prettyString(chi2Random[pI][pI2][rI][tI][bI], 3, false)).c_str(), "P L");
	    }
	  }

	  gStyle->SetOptStat(0);
	  leg_p->Draw("SAME");
	
	  //	  canv_p->SaveAs(("pdfDir/" + dateStr + "/nFill" + std::to_string(nFill) + "/" + powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][nBayes-1]->GetName() + "_" + saveStr).c_str());	  
	  delete canv_p;
	  delete leg_p;
	}
      }
    }
  }

  for(Int_t pI = 0; pI < nPowers; ++pI){
    for(Int_t pI2 = TMath::Max(0, pI-1); pI2 < TMath::Min(nPowers, pI+1); ++pI2){
      for(Int_t tI = 0; tI < nTarget; ++tI){
	for(Int_t bIX = 0; bIX < nGeneralBins; ++bIX){
	  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
	  canv_p->SetTopMargin(0.01);
	  canv_p->SetRightMargin(0.01);

	  TLegend* leg_p = new TLegend(0.60, 0.60, 0.95, 0.95);
	  leg_p->SetTextFont(43);
	  leg_p->SetTextSize(10);
	  leg_p->SetBorderSize(0);
	  leg_p->SetFillStyle(0);
	  leg_p->SetFillColor(0);
	  
	  for(Int_t bI = 0; bI < nBayes; ++bI){	    
	    powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->Sumw2();
	    powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->Scale(1./powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->Integral());
	    powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->SetMaximum(1.1);
	    powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->SetMinimum(0.0);

	    powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->SetMarkerColor(kPal.getColor(bI%nCol));
	    powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->SetLineColor(kPal.getColor(bI%nCol));
	    powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->SetMarkerStyle(styles[bI%nStyle]);
	    powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->SetMarkerSize(1.);
	    
	    if(bI == 0) powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->DrawCopy("HIST E1 P");
	    else powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->DrawCopy("HIST E1 P SAME");

	    std::string legStr = std::to_string(bayesVals[bI]) + ", #mu=" + prettyString(powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->GetMean(), 3, false) + ", #sigma=" + prettyString(powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX]->GetStdDev(), 3, false);
	    leg_p->AddEntry(powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX], legStr.c_str(), "P L");
	  }

	  leg_p->Draw("SAME");
	  
	  std::string saveName = "pdfDir/" + dateStr + "/nFill" + std::to_string(nFill) + "/" + powerUnfoldsBinDeltasRandom[pI][pI2][tI][0][bIX]->GetName() + "_" + saveStr;
	  quietSaveAs(canv_p, saveName);
	    
	  delete canv_p;	  
	}
      }
    }
  }
  
  
  for(Int_t pI = 0; pI < nPowers; ++pI){
    delete response_Rebin_h[pI];
    delete genSpectraPower_Rebin_p[pI];
    delete genSpectraTruncPower_Rebin_p[pI];
    delete recoSpectraPower_Rebin_p[pI];

    for(Int_t rI = 0; rI < nRandom; ++rI){
      for(Int_t tI = 0; tI < nTarget; ++tI){
	delete genSpectraTruncPower_Random_Rebin_p[pI][rI][tI];
	delete recoSpectraPower_Random_Rebin_p[pI][rI][tI];
      }
    }
    
    for(Int_t pI2 = TMath::Max(pI-1,0); pI2 < TMath::Min(pI+1, nPowers); ++pI2){
      delete response_NewPrior_h[pI][pI2];

      for(Int_t bI = 0; bI < nBayes; ++bI){
	delete powerUnfolds[pI][pI2][bI];
	delete powerUnfoldsGenRats[pI][pI2][bI];
      }	

      for(Int_t rI = 0; rI < nRandom; ++rI){
	for(Int_t tI = 0; tI < nTarget; ++tI){
	  for(Int_t bI = 0; bI < nBayes; ++bI){
	    delete powerUnfoldsRandom[pI][pI2][rI][tI][bI];
	    delete powerUnfoldsGenRatsRandom[pI][pI2][rI][tI][bI];

	    if(rI == 0){
	      for(Int_t bIX = 0; bIX < nGeneralBins; ++bIX){
		delete powerUnfoldsBinDeltasRandom[pI][pI2][tI][bI][bIX];
	      }
	    }
	  }
	}
      }
    }
  }
  
  outFile_p->Close();
  delete outFile_p;

  if(doFile2){
    inFile2_p->Close();
    delete inFile2_p;
  }
  
  inFile_p->Close();
  delete inFile_p;

  total.stop();
  std::cout << "Total timing: " << total.totalCPU() << ", " << total.totalWall() << std::endl;
  std::cout << "Unfold timing: " << unfoldTime.totalCPU() << ", " << unfoldTime.totalWall() << std::endl;
  std::cout << "Unfold NonRand timing: " << unfoldTimeNonRand.totalCPU() << ", " << unfoldTimeNonRand.totalWall() << std::endl;
  std::cout << "JustUnfold NonRand timing: " << justUnfoldTimeNonRand.totalCPU() << ", " << justUnfoldTimeNonRand.totalWall() << std::endl;
  std::cout << "Unfold Rand timing: " << unfoldTimeRand.totalCPU() << ", " << unfoldTimeRand.totalWall() << std::endl;
  std::cout << "JustUnfold Rand timing: " << justUnfoldTimeRand.totalCPU() << ", " << justUnfoldTimeRand.totalWall() << std::endl;
  
  return 0;
}


int main(int argc, char* argv[])
{
  if(argc < 4 || argc > 9){
    std::cout << "Usage: ./bin/toyATLASCMS_Unfold.exe <inFileName> <startTag> <doClean> <inFileName2> <newGenBins> <newRecoBins> <newReducedBins> <inDataFile>" << std::endl;
    return 1;
  }

  for(int aI = 0; aI < argc; ++aI){
    std::cout << " " << aI << ": " << argv[aI] << std::endl;
  }
  
  int retVal = 0;
  if(argc == 4) retVal += toyATLASCMS_Unfold(argv[1], argv[2], std::stoi(argv[3]));
  else if(argc == 5) retVal += toyATLASCMS_Unfold(argv[1], argv[2], std::stoi(argv[3]), argv[4]);
  else if(argc == 6) retVal += toyATLASCMS_Unfold(argv[1], argv[2], std::stoi(argv[3]), argv[4], argv[5]);
  else if(argc == 7) retVal += toyATLASCMS_Unfold(argv[1], argv[2], std::stoi(argv[3]), argv[4], argv[5], argv[6]);
  else if(argc == 8) retVal += toyATLASCMS_Unfold(argv[1], argv[2], std::stoi(argv[3]), argv[4], argv[5], argv[6], argv[7]);
  else if(argc == 9) retVal += toyATLASCMS_Unfold(argv[1], argv[2], std::stoi(argv[3]), argv[4], argv[5], argv[6], argv[7], argv[8]);

  return retVal;
}
