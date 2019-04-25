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

void setPrior(TH2D* res_p, TH1D* prior_p)
{
  for(Int_t bIY = 0; bIY < res_p->GetYaxis()->GetNbins(); ++bIY){
    Double_t total = 0.0;
    for(Int_t bIX = 0; bIX < res_p->GetXaxis()->GetNbins(); ++bIX){
      total += res_p->GetBinContent(bIX+1, bIY+1);
    }

    if(total <= TMath::Power(10, -20)) continue;

    total = prior_p->GetBinContent(bIY+1)/total;
    for(Int_t bIX = 0; bIX < res_p->GetXaxis()->GetNbins(); ++bIX){
      Double_t val = res_p->GetBinContent(bIX+1, bIY+1)*total;
      Double_t err = res_p->GetBinError(bIX+1, bIY+1)*total;

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

void doUnfold(TH2D* matrix_p, TH1D* prior_p, TH1D* data_p, Int_t nBayes, Int_t bayesVals[], std::string tagStr, TH1D* unfolded_p[], TH2D* newMatrix_p)
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
  setPrior(matrixClone_p, prior_p);
  
  for(Int_t bIX = 0; bIX < matrixClone_p->GetXaxis()->GetNbins(); ++bIX){
    for(Int_t bIY = 0; bIY < matrixClone_p->GetYaxis()->GetNbins(); ++bIY){
      newMatrix_p->SetBinContent(bIX+1, bIY+1, matrixClone_p->GetBinContent(bIX+1, bIY+1));
      newMatrix_p->SetBinError(bIX+1, bIY+1, 0.0);
    }
  }
  
  
  RooUnfoldResponse* rooRes_p = NULL;
  RooUnfoldBayes* rooBayes_p = NULL;
  
  macroHistToSubsetHistX(matrixClone_p, reco_p, true);
  macroHistToSubsetHistY(matrixClone_p, gen_p, true);

  TH1D* tempUnfold_h = NULL;
  /*        
  rooRes_p = new RooUnfoldResponse(reco_p, gen_p, matrixClone_p);
  rooBayes_p = new RooUnfoldBayes(rooRes_p, data_p, 5, false, "temp");
  rooBayes_p->SetVerbose(-1);
  rooBayes_p->SetNToys(1000);
  tempUnfold_h = (TH1D*)rooBayes_p->Hreco(RooUnfold::kCovToy);
  setPrior(matrixClone_p, tempUnfold_h);

  delete rooRes_p;
  delete rooBayes_p;
  delete tempUnfold_h;
  */
  for(Int_t bI = 0; bI < nBayes; ++bI){   
    TH2D* matrixClone2_p = (TH2D*)matrixClone_p->Clone("matrixClone2_p");    
    macroHistToSubsetHistX(matrixClone2_p, reco_p, true);
    macroHistToSubsetHistY(matrixClone2_p, gen_p, true);
        
    rooRes_p = new RooUnfoldResponse(reco_p, gen_p, matrixClone2_p);
    rooBayes_p = new RooUnfoldBayes(rooRes_p, data_p, bayesVals[bI], false, "temp");
    rooBayes_p->SetVerbose(-1);
    rooBayes_p->SetNToys(1000);
    tempUnfold_h = (TH1D*)rooBayes_p->Hreco(RooUnfold::kCovToy);

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

int toyATLASCMS_Unfold(const std::string inFileName, const std::string startTag, bool doClean, const std::string inFileName2 = "", std::string newGenBins = "", std::string newRecoBins = "", std::string newReducedBins = "") 
{
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
  /*
  const Int_t nPowers = 4;
  const Double_t powers[nPowers] = {2.,3.,4.,5.};
  */

  std::cout << "LINE: " << __LINE__ << std::endl;

  const Int_t nPowers = 3;
  Double_t powers[nPowers] = {7., 8., 9.};

  //  if(startTag.find("R10") != std::string::npos) ++powers[0];
  
  TH2D* cmsR10Response_h[nPowers];
  TH2D* cmsR10Response2_h[nPowers];
  TH2D* cmsR10Response_NewPrior_h[nPowers][nPowers];
  TH1D* genSpectraPower_p[nPowers];
  //  TH1D* genSpectraCMSR10TruncPower_p[nPowers];
  //  TH1D* recoCMSR10SpectraPower_p[nPowers];
  TH2D* recoGenSpectraCMSR10TruncPower_p[nPowers];

  
  TH2D* cmsR10Response_Rebin_h[nPowers];
  TH1D* genSpectraPower_Rebin_p[nPowers];
  TH1D* genSpectraCMSR10TruncPower_Rebin_p[nPowers];
  TH1D* recoCMSR10SpectraPower_Rebin_p[nPowers];
  TH2D* recoGenSpectraCMSR10TruncPower_Rebin_p[nPowers];

  std::cout << "LINE: " << __LINE__ << std::endl;
  std::string startTag2 = lowerToUpper(startTag.substr(0, startTag.find("R"))) + startTag.substr(startTag.find("R"), startTag.size());;
  
  for(Int_t pI = 0; pI < nPowers; ++pI){
    std::string powStr = "Power" + prettyString(powers[pI], 1, true);
    cmsR10Response_h[pI] = (TH2D*)inFile_p->Get((startTag + "Response" + powStr + "_h").c_str());      
    std::cout << "LINE: " << __LINE__ << std::endl;
    std::cout << startTag + "Response" + powStr + "_h" << ", " << inFileName << std::endl;
    std::cout << cmsR10Response_h[pI]->GetName() << std::endl;
    cmsR10Response2_h[pI] = NULL;
    if(doFile2) cmsR10Response2_h[pI] = (TH2D*)inFile2_p->Get((startTag + "Response" + powStr + "_h").c_str());      

    genSpectraPower_p[pI] = (TH1D*)inFile_p->Get(("genSpectra" + powStr + "_h").c_str());
    //    genSpectraCMSR10TruncPower_p[pI] = (TH1D*)inFile_p->Get(("genSpectraCMSR10Trunc" + powStr + "_h").c_str());
    //    recoCMSR10SpectraPower_p[pI] = (TH1D*)inFile_p->Get(("recoCMSR10Spectra" + powStr + "_h").c_str());
    recoGenSpectraCMSR10TruncPower_p[pI] = (TH2D*)inFile_p->Get(("recoGenSpectra" + startTag2 + "Trunc" + powStr + "_h").c_str());
  }

  std::cout << "LINE: " << __LINE__ << std::endl;

  const Int_t nMaxBins = 1000;
  Int_t nGenBins = genSpectraPower_p[0]->GetXaxis()->GetNbins();
  Double_t genBins[nMaxBins];
  for(Int_t bI = 0; bI < genSpectraPower_p[0]->GetXaxis()->GetNbins()+1; ++bI){
    genBins[bI] = genSpectraPower_p[0]->GetXaxis()->GetBinLowEdge(bI+1);
  }

  std::cout << "LINE: " << __LINE__ << std::endl;

  Int_t nRecoBins = cmsR10Response_h[0]->GetXaxis()->GetNbins();
  Double_t recoBins[nMaxBins];
  for(Int_t bI = 0; bI < cmsR10Response_h[0]->GetXaxis()->GetNbins()+1; ++bI){
    recoBins[bI] = cmsR10Response_h[0]->GetXaxis()->GetBinLowEdge(bI+1);
  }

  std::cout << "LINE: " << __LINE__ << std::endl;

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

  std::cout << "LINE: " << __LINE__ << std::endl;

  binsVect = binsStrToVect(newRecoBins);
  if(binsVect.size() != 0){
    if(!replaceBins(binsVect, &nRecoBins, recoBins)){
      std::cout << "GIVEN nRecoBins doesnt represent a subset of macrobins. return 1" << std::endl;
      return 1;
    }
  }

  std::cout << "LINE: " << __LINE__ << std::endl;

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

  std::cout << "LINE: " << __LINE__ << std::endl;
  std::cout << "NEW BINS: " << std::endl;
  for(Int_t gI = 0; gI < nGenBins+1; ++gI){
    std::cout << " " << gI << ": " << genBins[gI] << std::endl;
  }

  std::cout << "LINE: " << __LINE__ << std::endl;
  std::cout << "NEW BINS: " << std::endl;
  for(Int_t gI = 0; gI < nRecoBins+1; ++gI){
    std::cout << " " << gI << ": " << recoBins[gI] << std::endl;
  }

  std::cout << "LINE: " << __LINE__ << std::endl;

  Int_t nGeneralBins = 0;
  Double_t generalBins[nMaxBins];
  for(Int_t bI = 0; bI < nGenBins+1; ++bI){    
    if(genBins[bI] < recoBins[0]) continue;
    if(genBins[bI] > recoBins[nRecoBins]) continue;
    
    generalBins[nGeneralBins] = genBins[bI];
    ++nGeneralBins;
  }
  --nGeneralBins;

  std::cout << "LINE: " << __LINE__ << std::endl;

  binsVect = binsStrToVect(newReducedBins);
  if(binsVect.size() != 0){
    if(!replaceBins(binsVect, &nGeneralBins, generalBins)){
      std::cout << "GIVEN nGeneralBins doesnt represent a subset of macrobins. return 1" << std::endl;
      return 1;
    }
  }

  std::cout << "LINE: " << __LINE__ << std::endl;

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

  std::cout << "LINE: " << __LINE__ << std::endl;
  std::cout << "LINE: " << __LINE__ << std::endl;
  
  std::cout << "LINE: " << __LINE__ << std::endl;
  
  for(Int_t pI = 0; pI < nPowers; ++pI){
    std::cout << "LINE: " << __LINE__ << std::endl;
    std::string powStr = "Power" + prettyString(powers[pI], 1, true);
    cmsR10Response_Rebin_h[pI] = new TH2D((startTag + "Response" + powStr + "_Rebin_h").c_str(), "", nRecoBins, recoBins, nGenBins, genBins);      
    genSpectraPower_Rebin_p[pI] = new TH1D(("genSpectra" + powStr + "_Rebin_h").c_str(), "", nGenBins, genBins);
    genSpectraCMSR10TruncPower_Rebin_p[pI] = new TH1D(("genSpectra" + startTag2 + "Trunc" + powStr + "_Rebin_h").c_str(), "", nGenBins, genBins);
    recoGenSpectraCMSR10TruncPower_Rebin_p[pI] = new TH2D(("recoGenSpectra" + startTag2 + "Trunc" + powStr + "_Rebin_h").c_str(), "", nRecoBins, recoBins, nGenBins, genBins);
    recoCMSR10SpectraPower_Rebin_p[pI] = new TH1D(("reco" + startTag2 + "Spectra" + powStr + "_Rebin_h").c_str(), "", nRecoBins, recoBins);

    std::cout << "LINE: " << __LINE__ << ", " << doFile2 << std::endl;

    if(!doFile2) macroHistToSubsetHist(cmsR10Response_h[pI], cmsR10Response_Rebin_h[pI], true);
    else macroHistToSubsetHist(cmsR10Response2_h[pI], cmsR10Response_Rebin_h[pI], true);

    std::cout << "LINE: " << __LINE__ << std::endl;

    std::cout << "LINE: " << __LINE__ << std::endl;
      
    macroHistToSubsetHist(genSpectraPower_p[pI], genSpectraPower_Rebin_p[pI], true);
    
    macroHistToSubsetHist(recoGenSpectraCMSR10TruncPower_p[pI], recoGenSpectraCMSR10TruncPower_Rebin_p[pI], true);

    if(doClean){
      for(Int_t bIX = 0; bIX < nRecoBins; ++bIX){
	for(Int_t bIY = 0; bIY < nGenBins; ++bIY){
	  Double_t val = recoGenSpectraCMSR10TruncPower_Rebin_p[pI]->GetBinContent(bIX+1, bIY+1);
	  Double_t err = recoGenSpectraCMSR10TruncPower_Rebin_p[pI]->GetBinError(bIX+1, bIY+1);
	  if(val <= TMath::Power(10, -20)) continue;

	  if(err/val > 1./TMath::Sqrt(10)){
	    recoGenSpectraCMSR10TruncPower_Rebin_p[pI]->SetBinContent(bIX+1, bIY+1, 0.0);
	    recoGenSpectraCMSR10TruncPower_Rebin_p[pI]->SetBinError(bIX+1, bIY+1, 0.0);
	  }
	}
      }
    }

    macroHistToSubsetHistY(recoGenSpectraCMSR10TruncPower_Rebin_p[pI], genSpectraCMSR10TruncPower_Rebin_p[pI], true);
    macroHistToSubsetHistX(recoGenSpectraCMSR10TruncPower_Rebin_p[pI], recoCMSR10SpectraPower_Rebin_p[pI], true);
  
    //    macroHistToSubsetHist(genSpectraCMSR10TruncPower_p[pI], genSpectraCMSR10TruncPower_Rebin_p[pI], true);
    //    macroHistToSubsetHist(recoCMSR10SpectraPower_p[pI], recoCMSR10SpectraPower_Rebin_p[pI], true);

    ratioPrint(cmsR10Response_Rebin_h[pI], genSpectraCMSR10TruncPower_Rebin_p[pI], true);
    relErrPrint(genSpectraCMSR10TruncPower_Rebin_p[pI]);
    genSpectraCMSR10TruncPower_Rebin_p[pI]->Print("ALL");


    std::cout << "LINE: " << __LINE__ << std::endl;
    /*
    macroHistToSubsetHistY(cmsR10Response_Rebin_h[pI], genSpectraCMSR10TruncPower_Rebin_p[pI], true);       
    genSpectraCMSR10TruncPower_Rebin_p[pI]->Print("ALL");

    std::cout << "LINE: " << __LINE__ << std::endl;

    macroHistToSubsetHist(genSpectraCMSR10TruncPower_p[pI], genSpectraCMSR10TruncPower_Rebin_p[pI], true);
    

    std::cout << "LINE: " << __LINE__ << std::endl;
    */
    /*
    macroHistToSubsetHistX(cmsR10Response_Rebin_h[pI], recoCMSR10SpectraPower_Rebin_p[pI], true);       
    macroHistToSubsetHistY(cmsR10Response_Rebin_h[pI], genSpectraCMSR10TruncPower_Rebin_p[pI], true);       
    */
    std::string powOrigStr = "OrigPower" + prettyString(powers[pI], 1, true);
      
    for(Int_t pI2 = 0; pI2 < nPowers; ++pI2){
      std::string powPriorStr = "PriorPower" + prettyString(powers[pI2], 1, true);
      std::cout << "PRENAME: " << startTag << "Response_NewPrior_" << powOrigStr << "_" << powPriorStr << "_h" << std::endl;
      cmsR10Response_NewPrior_h[pI][pI2] = new TH2D((startTag + "Response_NewPrior_" + powOrigStr + "_" + powPriorStr + "_h").c_str(), "", nRecoBins, recoBins, nGenBins, genBins);
    }
  }

  
  TH2D* subsetRes_p = new TH2D("subsetRes_h", ";Reco Normalized to Unity;Gen", nRecoBins2, recoBins2, nGenBins2, genBins2);
  TH2D* subsetResErr_p = new TH2D("subsetResErr_h", ";Reco;Gen", nRecoBins2, recoBins2, nGenBins2, genBins2);
  std::cout << "LINE: " << __LINE__ << std::endl;
  macroHistToSubsetHist(cmsR10Response_Rebin_h[0], subsetRes_p, true, true, true); 
  macroHistToSubsetHist(cmsR10Response_Rebin_h[0], subsetResErr_p, true, true, true); 

  centerTitles(subsetRes_p);
  centerTitles(subsetResErr_p);
  
  std::cout << "LINE: " << __LINE__ << std::endl;

  std::string saveStr = startTag + "_Gen" + prettyString(genBins[0],1,true) + "to" + prettyString(genBins[nGenBins], 1, true) + "_Reco" + prettyString(recoBins[0],1,true) + "to" + prettyString(recoBins[nRecoBins], 1, true) + "_nFill" + std::to_string(nFill) + "_" + file2Str + "_" + cleanStr + "_" + dateStr + ".pdf";
  std::cout << "LINE: " << __LINE__ << std::endl;
  
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
	  
	  Double_t val = cmsR10Response_Rebin_h[pI]->GetBinContent(bIX+1, bIY+1);
	  Double_t err = cmsR10Response_Rebin_h[pI]->GetBinError(bIX+1, bIY+1);

	  if(val <= TMath::Power(10, -20)) continue;
	  
	  if(err/val > 1./TMath::Sqrt(10)){
	    cmsR10Response_Rebin_h[pI]->SetBinContent(bIX+1, bIY+1, 0.0);
	    cmsR10Response_Rebin_h[pI]->SetBinError(bIX+1, bIY+1, 0.0);
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
  canv_p->SaveAs(("pdfDir/" + dateStr + "/nFill" + std::to_string(nFill) + "/subsetRes_" + saveStr).c_str());
  delete canv_p;
  
  delete subsetRes_p;
  delete subsetResErr_p;

  std::cout << "LINE: " << __LINE__ << std::endl;

  /*  
  for(Int_t pI = 0; pI < nPowers; ++pI){
    std::string powOrigStr = "OrigPower" + prettyString(powers[pI], 1, true);

    for(Int_t pI2 = 0; pI2 < nPowers; ++pI2){
      std::string powPriorStr = "PriorPower" + prettyString(powers[pI2], 1, true);
      cmsR10Response_NewPrior_h[pI][pI2] = new TH2D((startTag + "Response_NewPrior_" + powOrigStr + "_" + powPriorStr + "_h").c_str(), "", nRecoBins, recoBins, nGenBins, genBins);
    }
  }
  */
  
  std::cout << "STARTED UNFOLDING" << std::endl;
  //Unfolding testing
  const Int_t nBayes = 4;
  Int_t bayesVals[nBayes] = {1,3,10,100};

  //ONE FOR UNFOLDING MATRIX YADA YADA
  TH1D* cmsR10PowerUnfolds[nPowers][nPowers][nBayes];
  TH1D* cmsR10PowerUnfoldsGenRats[nPowers][nPowers][nBayes];

  for(Int_t pI = 0; pI < nPowers; ++pI){
    for(Int_t pI2 = 0; pI2 < nPowers; ++pI2){
      doUnfold(cmsR10Response_Rebin_h[pI2], genSpectraCMSR10TruncPower_Rebin_p[pI], recoCMSR10SpectraPower_Rebin_p[pI2], nBayes, bayesVals, startTag + "UnfoldPower" + prettyString(powers[pI], 1, true) + "_SpectPower" + prettyString(powers[pI2], 1, true), cmsR10PowerUnfolds[pI][pI2], cmsR10Response_NewPrior_h[pI][pI2]);    
      
      for(Int_t bI = 0; bI < nBayes; ++bI){
	cmsR10PowerUnfoldsGenRats[pI][pI2][bI] = new TH1D((startTag + "PowerUnfoldsGenRats_UnfoldPow" + prettyString(powers[pI], 1, true) + "_SpectPow" + prettyString(powers[pI2], 1, true) + "_Bayes" + std::to_string(bayesVals[bI]) + "_h").c_str(), ";Jet pT;Unfolded/Gen", nGeneralBins, generalBins);
	centerTitles(cmsR10PowerUnfoldsGenRats[pI][pI2][bI]);
	
	TH1D* tempHist_p = new TH1D("tempHist_p", "", nGeneralBins, generalBins);
	TH1D* tempHist2_p = new TH1D("tempHist2_p", "", nGeneralBins, generalBins);
	macroHistToSubsetHist(cmsR10PowerUnfolds[pI][pI2][bI], tempHist_p);
	macroHistToSubsetHist(genSpectraCMSR10TruncPower_Rebin_p[pI2], tempHist2_p);
	
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
	  
	  cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->SetBinContent(gI+1, val);
	  cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->SetBinError(gI+1, 0.0);
	}

	delete tempHist_p;
	delete tempHist2_p;
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

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  checkMakeDir("pdfDir/" + dateStr + "/nFill" + std::to_string(nFill));

  Double_t max = -1;
  Double_t min = 10;

  for(Int_t pI = 0; pI < nPowers; ++pI){
    for(Int_t pI2 = 0; pI2 < nPowers; ++pI2){

      for(Int_t bI = 0; bI < nBayes; ++bI){
	if(bayesVals[bI] != 100) continue;

	if(max < getMax(cmsR10PowerUnfoldsGenRats[pI][pI2][bI])) max = getMax(cmsR10PowerUnfoldsGenRats[pI][pI2][bI]);
	if(min >= getMin(cmsR10PowerUnfoldsGenRats[pI][pI2][bI])) min = getMin(cmsR10PowerUnfoldsGenRats[pI][pI2][bI]);
      }
    }
  }

  //  Double_t interval = max-min;
  //  max += interval/10.;
  //  min -= interval/10.;  
  max = 1.5;
  min = 0.5;
    
  for(Int_t pI = 0; pI < nPowers; ++pI){
    cmsR10Response_Rebin_h[pI]->Write("", TObject::kOverwrite);
    
    genSpectraPower_Rebin_p[pI]->Write("", TObject::kOverwrite);
    genSpectraCMSR10TruncPower_Rebin_p[pI]->Write("", TObject::kOverwrite);
    recoCMSR10SpectraPower_Rebin_p[pI]->Write("", TObject::kOverwrite);

    if(!doFile2){
      canv_p = new TCanvas("canv_p", "canv_p", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);

      genSpectraPower_Rebin_p[pI]->SetMarkerStyle(24);
      genSpectraPower_Rebin_p[pI]->SetMarkerSize(1);
      genSpectraPower_Rebin_p[pI]->SetMarkerColor(1);
      genSpectraPower_Rebin_p[pI]->SetLineColor(1);
      
      genSpectraCMSR10TruncPower_Rebin_p[pI]->SetMarkerStyle(25);
      genSpectraCMSR10TruncPower_Rebin_p[pI]->SetMarkerSize(1);
      genSpectraCMSR10TruncPower_Rebin_p[pI]->SetMarkerColor(4);
      genSpectraCMSR10TruncPower_Rebin_p[pI]->SetLineColor(4);
      
      recoCMSR10SpectraPower_Rebin_p[pI]->SetMarkerStyle(28);
      recoCMSR10SpectraPower_Rebin_p[pI]->SetMarkerSize(1);
      recoCMSR10SpectraPower_Rebin_p[pI]->SetMarkerColor(2);
      recoCMSR10SpectraPower_Rebin_p[pI]->SetLineColor(2);

      TLegend* leg_p = new TLegend(0.45, 0.75, 0.9, 0.95);
      leg_p->SetBorderSize(0);
      leg_p->SetFillStyle(0);
      leg_p->SetTextFont(43);
      leg_p->SetTextSize(16);

      TLatex* label_p = new TLatex();
      label_p->SetNDC();
      label_p->SetTextFont(43);
      label_p->SetTextSize(16);

      Double_t total = 0.0;
      for(Int_t bIY = 0; bIY < recoGenSpectraCMSR10TruncPower_Rebin_p[pI]->GetNbinsY(); ++bIY){
	total += recoGenSpectraCMSR10TruncPower_Rebin_p[pI]->GetBinContent(1, bIY+1);
      }
      
      Int_t pos = -1;
      Double_t total2 = 0.0;
      for(Int_t bIY = 0; bIY < recoGenSpectraCMSR10TruncPower_Rebin_p[pI]->GetNbinsY(); ++bIY){
	total2 += recoGenSpectraCMSR10TruncPower_Rebin_p[pI]->GetBinContent(1, bIY+1);
	if(total2/total > 0.01){
	  pos = bIY;
	  break;
	}
      }
      
      Double_t center = recoGenSpectraCMSR10TruncPower_Rebin_p[pI]->GetYaxis()->GetBinCenter(pos+1);
      Double_t lineMin = genSpectraPower_Rebin_p[pI]->GetMinimum();
      Double_t lineMax = genSpectraPower_Rebin_p[pI]->GetMaximum();
      
      std::cout << "POS: " << pos << ", " << center << std::endl;
      
      const Int_t nRecoSubsetBins = 1;
      Double_t recoSubsetBins[nRecoSubsetBins+1];
      recoSubsetBins[0] = recoBins[0];
      recoSubsetBins[1] = recoBins[1];
      TH2D* subset2D_p = new TH2D("subset2D_p", "", nRecoSubsetBins, recoSubsetBins, nGenBins, genBins);
      macroHistToSubsetHist(recoGenSpectraCMSR10TruncPower_Rebin_p[pI], subset2D_p, true);
      TH1D* subset1D_p = new TH1D("subset1D_p", "", nGenBins, genBins);
      macroHistToSubsetHistY(subset2D_p, subset1D_p, true);
      
      subset1D_p->SetMarkerColor(4);
      subset1D_p->SetLineColor(4);
      
      subset1D_p->Print("ALL");
      
      genSpectraPower_Rebin_p[pI]->DrawCopy("HIST E1");
      genSpectraCMSR10TruncPower_Rebin_p[pI]->DrawCopy("HIST E1 SAME");
      recoCMSR10SpectraPower_Rebin_p[pI]->DrawCopy("HIST E1 SAME");
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
      leg_p->AddEntry(genSpectraCMSR10TruncPower_Rebin_p[pI], "Gen. w/ Reco.", "P L");
      leg_p->AddEntry(recoCMSR10SpectraPower_Rebin_p[pI], "Reco.", "P L");
      leg_p->AddEntry(subset1D_p, "Gen. From 1st Reco. Bin", "P L");
      
      leg_p->Draw("SAME");
      label_p->DrawLatex(0.5, 0.65, ("Reco pT > " + prettyString(recoBins[0],1,false)).c_str());

      canv_p->SaveAs(("pdfDir/" + dateStr + "/nFill" + std::to_string(nFill) + "/" + genSpectraCMSR10TruncPower_Rebin_p[pI]->GetName() + "_" + saveStr).c_str());
      delete canv_p;

      delete leg_p;
      delete label_p;
      delete subset2D_p;
      delete subset1D_p;
    }
    
    for(Int_t pI2 = 0; pI2 < nPowers; ++pI2){
      cmsR10Response_NewPrior_h[pI][pI2]->Write("", TObject::kOverwrite);
        
      canv_p = new TCanvas("canv_p", "canv_p", 450, 450);
      canv_p->SetTopMargin(0.01);
      canv_p->SetRightMargin(0.01);

      TLegend* leg_p = new TLegend(0.7, 0.15, 0.9, 0.35);
      leg_p->SetBorderSize(0);
      leg_p->SetFillStyle(0);
      
      for(Int_t bI = 0; bI < nBayes; ++bI){
	cmsR10PowerUnfolds[pI][pI2][bI]->Write("", TObject::kOverwrite);
	delete cmsR10PowerUnfolds[pI][pI2][bI];

	
	if(true/*bayesVals[bI] == 100*/){

	  cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->SetMaximum(max);
	  cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->SetMinimum(min);

	  cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->SetMarkerStyle(styles[bI%nStyle]);
	  cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->SetMarkerColor(kPal.getColor(bI%nCol));
	  cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->SetLineColor(kPal.getColor(bI%nCol));
								 
	  if(bI != 0) cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->DrawCopy("HIST E1 P SAME");
	  else cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->DrawCopy("HIST E1 P");

	  if(bayesVals[bI] == 100){
	    Double_t mean = 0.0;
	    
	    for(Int_t bIX = 0; bIX < cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->GetNbinsX(); ++bIX){
	      mean += cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->GetBinContent(bIX+1);
	    }
	    mean /= (Double_t)cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->GetNbinsX();
	    
	    Double_t sigma = 0.0;
	    
	    for(Int_t bIX = 0; bIX < cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->GetNbinsX(); ++bIX){
	      Double_t val = cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->GetBinContent(bIX+1) - mean;
	      sigma += val*val;
	    }
	    
	    sigma /= (Double_t)(cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->GetNbinsX() - 1);
	    sigma = TMath::Sqrt(sigma);
	    
	    TLatex* label_p = new TLatex();
	    label_p->SetNDC();
	    label_p->SetTextFont(43);
	    label_p->SetTextSize(14);
	    label_p->DrawLatex(0.6, 0.95, ("#mu=" + prettyString(mean, 3, false)).c_str());
	    label_p->DrawLatex(0.6, 0.87, ("#sigma=" + prettyString(sigma, 3, false)).c_str());
	    Double_t statsRes = cmsR10Response_Rebin_h[pI2]->Integral()*100;
	    Double_t statsReco = recoCMSR10SpectraPower_Rebin_p[pI2]->Integral()*100;
	    
	    label_p->DrawLatex(0.6, 0.79, ("IntRes=" + prettyStringE(statsRes, 3, false)).c_str());
	    label_p->DrawLatex(0.6, 0.71, ("IntReco=" + prettyStringE(statsReco, 3, false)).c_str());
	    label_p->DrawLatex(0.6, 0.63, ("Reco pT > " + prettyString(recoBins[0], 1, false)).c_str());
	  
	    delete label_p;
	  }
	}

	leg_p->AddEntry(cmsR10PowerUnfoldsGenRats[pI][pI2][bI], ("Bayes=" + std::to_string(bayesVals[bI])).c_str(), "P L");
	cmsR10PowerUnfoldsGenRats[pI][pI2][bI]->Write("", TObject::kOverwrite);
      }

      std::cout << "LINE: " << __LINE__ << std::endl;
      gStyle->SetOptStat(0);
      std::cout << "LINE: " << __LINE__ << std::endl;
      leg_p->Draw("SAME");
      
      canv_p->SaveAs(("pdfDir/" + dateStr + "/nFill" + std::to_string(nFill) + "/" + cmsR10PowerUnfoldsGenRats[pI][pI2][nBayes-1]->GetName() + "_" + saveStr).c_str());	  
      std::cout << "LINE: " << __LINE__ << std::endl;
      delete canv_p;
      delete leg_p;
      
      for(Int_t bI = 0; bI < nBayes; ++bI){
	delete cmsR10PowerUnfoldsGenRats[pI][pI2][bI];
      }
      std::cout << "LINE: " << __LINE__ << std::endl;
    }
  }
 
  for(Int_t pI = 0; pI < nPowers; ++pI){
    delete cmsR10Response_Rebin_h[pI];
    delete genSpectraPower_Rebin_p[pI];
    delete genSpectraCMSR10TruncPower_Rebin_p[pI];
    delete recoCMSR10SpectraPower_Rebin_p[pI];

    for(Int_t pI2 = 0; pI2 < nPowers; ++pI2){
      delete cmsR10Response_NewPrior_h[pI][pI2];
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
  
  return 0;
}


int main(int argc, char* argv[])
{
  if(argc < 4 || argc > 8){
    std::cout << "Usage: ./bin/toyATLASCMS_Unfold.exe <inFileName> <startTag> <doClean> <inFileName2> <newGenBins> <newRecoBins> <newReducedBins>" << std::endl;
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

  return retVal;
}
