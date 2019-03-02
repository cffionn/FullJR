//cpp dependencies
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TCanvas.h"
#include "TDatime.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMatrix.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

//RooUnfold dependencies
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldResponse.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/smallOrLargeR.h"
#include "MainAnalysis/include/texSlideCreator.h"

//Non-local FullJR dependencies (Utility, etc.)
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vanGoghPalette.h"

void correctForFakes(TH1D* rawHist_p, std::vector<double> binEdges, std::vector<double> fakeFactor)
{
  double deltaBinEdge = 0.1;

  std::vector<int> corrPos;
  for(unsigned int bI = 0; bI < binEdges.size()-1; ++bI){
    int tempPos = -1;

    for(Int_t bIX = 0; bIX < rawHist_p->GetNbinsX(); ++bIX){
      if(TMath::Abs(binEdges[bI] - rawHist_p->GetBinLowEdge(bIX+1)) > deltaBinEdge) continue;
      if(TMath::Abs(binEdges[bI+1] - rawHist_p->GetBinLowEdge(bIX+2)) > deltaBinEdge) continue;

      tempPos = bIX+1;
      break;
    }

    if(tempPos == -1){
      std::cout << "Warning! - corrPos = -1, mismatch in binnings." << std::endl;
      std::cout << " Input fake bin: " << binEdges[bI] << "-" << binEdges[bI+1] << std::endl;
      std::cout << " Histogram options: ";
      for(Int_t bIX = 0; bIX < rawHist_p->GetNbinsX(); ++bIX){
	std::cout << rawHist_p->GetBinLowEdge(bIX+1) << ", ";
      }
      std::cout << std::endl;
    }
    corrPos.push_back(tempPos);
  }

  for(unsigned int cI = 0; cI < corrPos.size(); ++cI){
    if(corrPos[cI] == -1) continue;

    double binVal = rawHist_p->GetBinContent(corrPos[cI]);
    binVal -= binVal*fakeFactor[cI];
    double binErr = TMath::Sqrt(binVal);

    rawHist_p->SetBinContent(corrPos[cI], binVal);
    rawHist_p->SetBinError(corrPos[cI], binErr);
  }

  return;
}

void getPearsTMatrix(RooUnfoldBayes* bayes_p, TH2D** covarianceToPlot)
{
  TMatrixD tempCovBayes = (TMatrixD)bayes_p->Ereco(RooUnfold::kCovToy);
  TMatrixD* pearsonCoefsBayes = (TMatrixD*)tempCovBayes.Clone("pearsonCoefsBayes");

  for(Int_t rI = 0; rI < pearsonCoefsBayes->GetNrows(); ++rI){
    for(Int_t cI = 0; cI < pearsonCoefsBayes->GetNcols(); ++cI){
      Double_t valBayes = tempCovBayes(rI, cI);
      bool isGoodDiag = tempCovBayes(rI, rI) > 0.0 && tempCovBayes(cI, cI) > 0.0;
      bool isGoodVal = valBayes > 0.0;
      
      if(isGoodDiag) valBayes /= TMath::Sqrt(tempCovBayes(rI, rI)*tempCovBayes(cI, cI));
      else if(!isGoodVal) valBayes = 0;
      else{
	std::cout << "Warning diag is zero but off diag val is non-zero" << std::endl;
      }
      
      (*(pearsonCoefsBayes))(rI, cI) = valBayes;
    }
  }
  
  (*covarianceToPlot) = new TH2D(*pearsonCoefsBayes);
  
  return;
}

int unfoldRawData(const std::string inDataFileName, const std::string inResponseName, const std::string selectJtAlgo = "")
{
  vanGoghPalette vg;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  const std::string dateStr2 = std::to_string(date->GetYear()) + "." + std::to_string(date->GetMonth()) + "." + std::to_string(date->GetDay());
  delete date;
  
  TFile* responseFile_p = new TFile(inResponseName.c_str(), "READ");
  std::vector<std::string> responseJetDirList = returnRootFileContentsList(responseFile_p, "TDirectoryFile", "JetAnalyzer");
  std::cout << "Printing " << responseJetDirList.size() << " response jets..." << std::endl;

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  checkMakeDir("pdfDir/" + dateStr + "/Unfold");

  //editing
  std::vector<std::vector<std::string> > slideTitlesPerAlgo;
  std::vector<std::vector<std::vector<std::string> > > pdfPerSlidePerAlgo;

  for(unsigned int jI = 0; jI < responseJetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << responseJetDirList.size() << ": " << responseJetDirList[jI] << std::endl;

    std::string tempStr = responseJetDirList[jI];
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    checkMakeDir("pdfDir/" + dateStr + "/Unfold_" + tempStr);
    slideTitlesPerAlgo.push_back({});
    pdfPerSlidePerAlgo.push_back({});
  }

  cutPropagator cutPropResponse;
  cutPropResponse.Clean();
  cutPropResponse.GetAllVarFromFile(responseFile_p);

  responseFile_p->Close();
  delete responseFile_p;

  TFile* dataFile_p = new TFile(inDataFileName.c_str(), "READ");
  std::vector<std::string> dataJetDirList = returnRootFileContentsList(dataFile_p, "TDirectoryFile", "JetAnalyzer");
  std::cout << "Printing " << dataJetDirList.size() << " data jets..." << std::endl;
  for(unsigned int jI = 0; jI < dataJetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << dataJetDirList.size() << ": " << dataJetDirList[jI] << std::endl;
  }

  cutPropagator cutPropData;
  cutPropData.Clean();
  cutPropData.GetAllVarFromFile(dataFile_p);

  dataFile_p->Close();
  delete dataFile_p;

  std::cout << cutPropResponse.GetPthats().size() << std::endl;

  if(!cutPropResponse.CheckPropagatorsMatch(cutPropData, false, true)){
    std::cout << "unfoldRawData - Cuts listed in data file \'" << inDataFileName << "\' and response file \'" << inResponseName << "\' do not match. return 1" << std::endl;
    return 1;
  }

  Int_t valForForLoops = 100000000;
  if(doLocalDebug || doGlobalDebug){
    std::cout << "DOLOCALDEBUG or DOGLOBALDEBUG in unfoldRawData: Setting all histogram array sizes to 1 for faster processing" << std::endl;
    valForForLoops = 1000000;
  }

  const Int_t isDataPP = cutPropData.GetIsPP();
  const Int_t isResponsePP = cutPropResponse.GetIsPP();
  
  const Int_t nMaxSyst = 20;
  const Int_t nSyst = TMath::Min(valForForLoops, cutPropResponse.GetNSyst());
  std::vector<std::string> systStr = cutPropResponse.GetSystStr();
  
  /*
  const Int_t nSyst = 1;
  std::vector<std::string> systStr = {""};
  */

  if(nSyst > nMaxSyst){
    std::cout << "nSyst \'" << nSyst << "\' is greater than nMaxSyst \'" << nMaxSyst << "\'. return 1" << std::endl;
    return 1;
  }

  const Int_t nResponseMod = 1;
  std::vector<double> responseMod = {0.10};
  std::vector<double> jerVarData = {0.10};
  
  /*
  const Int_t nResponseMod = TMath::Min(valForForLoops, cutPropResponse.GetNResponseMod());
  std::vector<double> responseMod = cutPropResponse.GetResponseMod();
  std::vector<double> jerVarData = cutPropResponse.GetJERVarData();
  */

  const Int_t nMaxCentBins = 4;
  const Int_t nCentBins = TMath::Min(valForForLoops, cutPropData.GetNCentBins());
  if(nCentBins > nMaxCentBins){
    std::cout << "nCentBins \'" << nCentBins << "\' is greater than nMaxCentBins \'" << nMaxCentBins << "\'. return 1" << std::endl;
    return 1;
  }

  std::vector<Int_t> centBinsLow = cutPropData.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutPropData.GetCentBinsHi();

  const Int_t nMaxJtPtBins = 50;
  const Int_t nGenJtPtSmallBinsSmallRCent0to10 = cutPropData.GetNGenJtPtSmallBinsSmallRCent0to10();
  const Int_t nGenJtPtLargeBinsSmallRCent0to10 = cutPropData.GetNGenJtPtLargeBinsSmallRCent0to10();
  const Int_t nGenJtPtSmallBinsLargeRCent0to10 = cutPropData.GetNGenJtPtSmallBinsLargeRCent0to10();
  const Int_t nGenJtPtLargeBinsLargeRCent0to10 = cutPropData.GetNGenJtPtLargeBinsLargeRCent0to10();

  const Int_t nGenJtPtSmallBinsSmallRCent10to30 = cutPropData.GetNGenJtPtSmallBinsSmallRCent10to30();
  const Int_t nGenJtPtLargeBinsSmallRCent10to30 = cutPropData.GetNGenJtPtLargeBinsSmallRCent10to30();
  const Int_t nGenJtPtSmallBinsLargeRCent10to30 = cutPropData.GetNGenJtPtSmallBinsLargeRCent10to30();
  const Int_t nGenJtPtLargeBinsLargeRCent10to30 = cutPropData.GetNGenJtPtLargeBinsLargeRCent10to30();

  const Int_t nGenJtPtSmallBinsSmallRCent30to50 = cutPropData.GetNGenJtPtSmallBinsSmallRCent30to50();
  const Int_t nGenJtPtLargeBinsSmallRCent30to50 = cutPropData.GetNGenJtPtLargeBinsSmallRCent30to50();
  const Int_t nGenJtPtSmallBinsLargeRCent30to50 = cutPropData.GetNGenJtPtSmallBinsLargeRCent30to50();
  const Int_t nGenJtPtLargeBinsLargeRCent30to50 = cutPropData.GetNGenJtPtLargeBinsLargeRCent30to50();

  const Int_t nGenJtPtSmallBinsSmallRCent50to90 = cutPropData.GetNGenJtPtSmallBinsSmallRCent50to90();
  const Int_t nGenJtPtLargeBinsSmallRCent50to90 = cutPropData.GetNGenJtPtLargeBinsSmallRCent50to90();
  const Int_t nGenJtPtSmallBinsLargeRCent50to90 = cutPropData.GetNGenJtPtSmallBinsLargeRCent50to90();
  const Int_t nGenJtPtLargeBinsLargeRCent50to90 = cutPropData.GetNGenJtPtLargeBinsLargeRCent50to90();

  std::vector<Double_t> genJtPtSmallBinsSmallRTemp = cutPropData.GetGenJtPtSmallBinsSmallR();
  std::vector<Double_t> genJtPtLargeBinsSmallRTemp = cutPropData.GetGenJtPtLargeBinsSmallR();
  std::vector<Double_t> genJtPtSmallBinsLargeRTemp = cutPropData.GetGenJtPtSmallBinsLargeR();
  std::vector<Double_t> genJtPtLargeBinsLargeRTemp = cutPropData.GetGenJtPtLargeBinsLargeR();
  
  const Int_t nJtAbsEtaBins = 1;
  std::vector<Double_t> jtAbsEtaBinsLowTemp = {0.0};
  std::vector<Double_t> jtAbsEtaBinsHiTemp = {2.0};
  
  /*
  const Int_t nJtAbsEtaBins = TMath::Min(valForForLoops, cutPropData.GetNJtAbsEtaBins());
  std::vector<Double_t> jtAbsEtaBinsLowTemp = cutPropData.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHiTemp = cutPropData.GetJtAbsEtaBinsHi();
  */

  /*
  const Int_t nID = TMath::Min(valForForLoops, cutPropData.GetNID());
  std::vector<std::string> idStr = cutPropData.GetIdStr();
  */

  const Int_t nID = 1;
  std::vector<std::string> idStr = {"LightMUAndCHID"};
  
  Double_t jtAbsEtaBinsLow[nJtAbsEtaBins];
  Double_t jtAbsEtaBinsHi[nJtAbsEtaBins];
  std::cout << "nJtAbsEtaBins: ";
  for(Int_t jI = 0; jI < nJtAbsEtaBins; ++jI){
    jtAbsEtaBinsLow[jI] = jtAbsEtaBinsLowTemp[jI];
    jtAbsEtaBinsHi[jI] = jtAbsEtaBinsHiTemp[jI];
    std::cout << " " << jtAbsEtaBinsLow[jI] << "-" << jtAbsEtaBinsHi[jI] << ",";
  }
  std::cout << std::endl;

  std::cout << "nCentBins: " << nCentBins << std::endl;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::cout << " " << cI << "/" << nCentBins << ": " << centBinsLow[cI] << "-" << centBinsHi[cI] << std::endl;
  }

  std::cout << "Raw data from file: \'" << inDataFileName << "\'" << std::endl;
  std::cout << "Response data from file: \'" << inResponseName << "\'" << std::endl;

  //Reduce to match jet dirs
  unsigned int pos = 0;
  while(pos < dataJetDirList.size()){
    bool isFound = false;
    for(unsigned int i = 0; i < responseJetDirList.size(); ++i){
      if(dataJetDirList[pos].size() != responseJetDirList[i].size()) continue;
      if(dataJetDirList[pos].find(responseJetDirList[i]) == std::string::npos) continue;
      
      isFound = true;
      break;
    }

    if(isFound) ++pos;
    else dataJetDirList.erase(dataJetDirList.begin()+pos);
  }

  pos = 0;
  while(pos < responseJetDirList.size()){
    bool isFound = false;
    for(unsigned int i = 0; i < dataJetDirList.size(); ++i){
      if(responseJetDirList[pos].size() != dataJetDirList[i].size()) continue;
      if(responseJetDirList[pos].find(dataJetDirList[i]) == std::string::npos) continue;
      
      isFound = true;
      break;
    }

    if(isFound) ++pos;
    else responseJetDirList.erase(responseJetDirList.begin()+pos);
  }

  std::cout << "Shared jets to process: " << std::endl;
  for(unsigned int i = 0; i < responseJetDirList.size(); ++i){
    std::cout << " " << i << "/" << responseJetDirList.size() << ": " << responseJetDirList[i] << std::endl;
  }

  if(selectJtAlgo.size() != 0){
    std::cout << "Restricting to \'" << selectJtAlgo << "\'." << std::endl;
    unsigned int pos = 0;
    while(pos < responseJetDirList.size()){
      if(responseJetDirList[pos].find(selectJtAlgo) == std::string::npos) responseJetDirList.erase(responseJetDirList.begin()+pos);
      else ++pos;
    }

    if(responseJetDirList.size() == 0){
      std::cout << "No jet dirs after selection for " << selectJtAlgo << std::endl;
      return 1;
    }
  }


  smallOrLargeR rReader;
  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent0to10(nGenJtPtSmallBinsSmallRCent0to10)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent0to10(nGenJtPtSmallBinsLargeRCent0to10)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent0to10(nGenJtPtLargeBinsSmallRCent0to10)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent0to10(nGenJtPtLargeBinsLargeRCent0to10)) return 1;

  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent10to30(nGenJtPtSmallBinsSmallRCent10to30)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent10to30(nGenJtPtSmallBinsLargeRCent10to30)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent10to30(nGenJtPtLargeBinsSmallRCent10to30)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent10to30(nGenJtPtLargeBinsLargeRCent10to30)) return 1;

  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent30to50(nGenJtPtSmallBinsSmallRCent30to50)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent30to50(nGenJtPtSmallBinsLargeRCent30to50)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent30to50(nGenJtPtLargeBinsSmallRCent30to50)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent30to50(nGenJtPtLargeBinsLargeRCent30to50)) return 1;

  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent50to90(nGenJtPtSmallBinsSmallRCent50to90)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent50to90(nGenJtPtSmallBinsLargeRCent50to90)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent50to90(nGenJtPtLargeBinsSmallRCent50to90)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent50to90(nGenJtPtLargeBinsLargeRCent50to90)) return 1;

  if(!rReader.CheckGenJtPtSmallBinsSmallR(genJtPtSmallBinsSmallRTemp)) return 1;
  if(!rReader.CheckGenJtPtSmallBinsLargeR(genJtPtSmallBinsLargeRTemp)) return 1;
  if(!rReader.CheckGenJtPtLargeBinsSmallR(genJtPtLargeBinsSmallRTemp)) return 1;
  if(!rReader.CheckGenJtPtLargeBinsLargeR(genJtPtLargeBinsLargeRTemp)) return 1;

  const Int_t nSmallLargeBins = 2;
  std::string smallLargeBinsStr[nSmallLargeBins] = {"SmallBins", "LargeBins"};

  const Int_t nMaxDataJet = 10;
  const Int_t nDataJet = responseJetDirList.size();

  if(nDataJet > nMaxDataJet){
    std::cout << "nDataJet \'" << nDataJet << "\' is greater than nMaxDataJet \'" << nMaxDataJet << "\'. return 1" << std::endl;
    return 1;
  }

  std::string outFileName = inDataFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");

  std::string outFileName2 = inResponseName;
  while(outFileName2.find("/") != std::string::npos){outFileName2.replace(0, outFileName2.find("/")+1, "");}
  if(outFileName2.find(".txt") != std::string::npos) outFileName2.replace(outFileName2.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName2.find(".root") != std::string::npos) outFileName2.replace(outFileName2.find(".root"), std::string(".root").size(), "");

  std::string debugStr = "";
  if(doLocalDebug || doGlobalDebug) debugStr = "DEBUG_";

  std::string selectJtAlgoStr = selectJtAlgo;
  if(selectJtAlgoStr.size() != 0) selectJtAlgoStr = selectJtAlgoStr + "_";

  const Int_t sizeToTruncName = 40;
  while(outFileName.size() > sizeToTruncName){outFileName = outFileName.substr(0,outFileName.size()-1);}
  while(outFileName2.size() > sizeToTruncName){outFileName2 = outFileName2.substr(0,outFileName2.size()-1);}

  const Int_t nSuperBayes = 0;
  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);

  outFileName = "output/" + dateStr + "/" + outFileName + "_" + outFileName2 + "_UnfoldRawData_NSuperBayes" + std::to_string(nSuperBayes) + "_" + selectJtAlgoStr + debugStr + dateStr + ".root";

  while(outFileName.find("__") != std::string::npos){outFileName.replace(outFileName.find("__"), 2, "_");}

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);

  const Int_t nBayes = 37;
  const Int_t nBigBayesSymm = 3;
  const Int_t nBayesBig = nBayes - (2*nBigBayesSymm + 1);
  Int_t temp100Pos = -1;
  Int_t bayesVal[nBayes];
  for(Int_t bI = 0; bI < nBayes; ++bI){
    if(bI >= nBayesBig) bayesVal[bI] = 100 + 1 + nBigBayesSymm - (nBayes - bI);
    else bayesVal[bI] = bI+1;

    if(bayesVal[bI] == 100) temp100Pos = bI;
  }


  const Int_t bayes100Pos = temp100Pos;

  std::cout << "Bayes: " << std::endl;
  for(Int_t bI = 0; bI < nBayes; ++bI){
    std::cout << " " << bI << "/" << nBayes << ": " << bayesVal[bI] << std::endl;
  }

  std::cout << "N Super Bayes: " << nSuperBayes << std::endl;

  //  return 0;

  const Int_t nHistDim = nSmallLargeBins*nDataJet*nCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst;
  std::vector<std::string> histTag;
  std::vector<int> histBestBayes;
  
  cutPropData.SetNBayes(nBayes);
  cutPropData.SetNBigBayesSymm(nBigBayesSymm);
  cutPropData.SetBayesVal(nBayes, bayesVal);
  cutPropData.SetNSuperBayes(nSuperBayes);

  cutPropData.SetNResponseMod(nResponseMod);
  cutPropData.SetResponseMod(responseMod);
  cutPropData.SetJERVarData(jerVarData);

  cutPropData.SetNSyst(nSyst);
  cutPropData.SetSystStr(systStr);

  const Int_t nBayesDraw = TMath::Min(valForForLoops, 4);
  TDirectory* dir_p[nMaxDataJet];
  TH1D* jtPtUnfolded_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSmallLargeBins][nBayes];
  RooUnfoldBayes* bayes_p[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSmallLargeBins][nBayes];
  Int_t histTermPos[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSmallLargeBins];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = responseJetDirList[jI];
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dir_p[jI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());

    const Int_t rVal = getRVal(tempStr);
    const bool isSmallR = rReader.GetIsSmallR(rVal);

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      
      const Int_t nGenJtPtSmallBins = rReader.GetSmallOrLargeRNBins(isSmallR, true, true, centStr);
      const Int_t nGenJtPtLargeBins = rReader.GetSmallOrLargeRNBins(isSmallR, true, false, centStr);
      
      Double_t genJtPtSmallBins[nMaxJtPtBins+1];
      Double_t genJtPtLargeBins[nMaxJtPtBins+1];
      rReader.GetSmallOrLargeRBins(isSmallR, true, nGenJtPtSmallBins+1, genJtPtSmallBins, true);
      rReader.GetSmallOrLargeRBins(isSmallR, true, nGenJtPtLargeBins+1, genJtPtLargeBins, false);
      
      const Int_t nGenJtPtBins[nSmallLargeBins] = {nGenJtPtSmallBins, nGenJtPtLargeBins};
      Double_t genJtPtBins[nSmallLargeBins][nMaxJtPtBins+1];
      
      for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
	for(Int_t xI = 0; xI < nGenJtPtBins[binsI]+1; ++xI){
	  if(binsI == 0) genJtPtBins[binsI][xI] = genJtPtSmallBins[xI];
	  else if(binsI == 1) genJtPtBins[binsI][xI] = genJtPtLargeBins[xI];
	}
      }

      if(isDataPP) centStr = "PP_" + centStr;
      else centStr = "PbPb_" + centStr;

      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){
          const std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);

	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

	    for(Int_t sI = 0; sI < nSyst; ++sI){
	      std::string tempSystStr = "_" + systStr[sI] + "_";
	      while(tempSystStr.find("__") != std::string::npos){tempSystStr.replace(tempSystStr.find("__"), 2, "_");}

	      for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
		histTermPos[jI][cI][idI][mI][aI][sI][binsI] = -1;

		for(Int_t bI = 0; bI < nBayes; ++bI){
		  std::string bayesStr = "Bayes" + std::to_string(bayesVal[bI]);
		  
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI][binsI][bI] = new TH1D(("jtPtUnfolded_RecoGenAsymm_" + smallLargeBinsStr[binsI] + "_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + bayesStr + "_h").c_str(), ";Unfolded Jet p_{T};Counts", nGenJtPtBins[binsI], genJtPtBins[binsI]);
		  centerTitles(jtPtUnfolded_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI][binsI][bI]);
		  setSumW2(jtPtUnfolded_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI][binsI][bI]);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  responseFile_p = new TFile(inResponseName.c_str(), "READ");
  //  TH2D* response_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSmallLargeBins];
  RooUnfoldResponse* rooResponse_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSmallLargeBins];
  
  TH1D* recoJtPt_GoodGen_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst];
  TH1D* genJtPt_GoodReco_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSmallLargeBins];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = responseJetDirList[jI];
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      if(isResponsePP) centStr = "PP_" + centStr;
      else centStr = "PbPb_" + centStr;

      for(Int_t iI = 0; iI < nID; ++iI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){
	  const std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);

	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	    
	    for(Int_t sI = 0; sI < nSyst; ++sI){
	      std::string tempSystStr = "_" + systStr[sI] + "_";
	      while(tempSystStr.find("__") != std::string::npos){tempSystStr.replace(tempSystStr.find("__"), 2, "_");}
	     
	      recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI] = (TH1D*)responseFile_p->Get((tempStr + "/recoJtPt_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "GoodGen_h").c_str());

	      for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
		//		response_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI] = (TH2D*)responseFile_p->Get((tempStr + "/response_" + smallLargeBinsStr[binsI] + "_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "RecoGenAsymm_h").c_str());

		if(jI == 0 && cI == 0 && iI == 0 && mI == 0 && aI == 0 && sI == 0 && binsI == 0){
		  std::cout << "NAME: " << tempStr + "/rooResponse_" + smallLargeBinsStr[binsI] + "_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "RecoGenAsymm_h" << std::endl;
		}
	      
		rooResponse_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI] = (RooUnfoldResponse*)responseFile_p->Get((tempStr + "/rooResponse_" + smallLargeBinsStr[binsI] + "_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "RecoGenAsymm_h").c_str());
	      
		if(jI == 0 && cI == 0 && iI == 0 && mI == 0 && aI == 0 && sI == 0 && binsI == 0){
		  std::cout << rooResponse_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI]->GetName() << std::endl;
		}

		genJtPt_GoodReco_h[jI][cI][iI][mI][aI][sI][binsI] = (TH1D*)responseFile_p->Get((tempStr + "/genJtPt_" + smallLargeBinsStr[binsI] + "_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "GoodReco_h").c_str());
	      }
	    }
	  }
	}
      }
    }
  }

  dataFile_p = new TFile(inDataFileName.c_str(), "READ");
  TH1D* jtPtRaw_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nJtAbsEtaBins];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = responseJetDirList[jI];
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      if(isDataPP) centStr = "PP_" + centStr;
      else centStr = "PbPb_" + centStr;

      for(Int_t iI = 0; iI < nID; ++iI){
        for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
          const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

	  const std::string name = tempStr + "/jtPtRaw_RecoGenAsymm_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_h";

	  jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI] = (TH1D*)dataFile_p->Get(name.c_str());
	}
      }
    }
  }

  std::cout << "Start Unfolding..." << std::endl;
  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = responseJetDirList[jI];
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    std::cout << " Unfolding " << jI << "/" << nDataJet << ": " << tempStr << std::endl;

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(!isDataPP) std::cout << "  " << centBinsLow[cI] << "-" << centBinsHi[cI] << "%..." << std::endl;
      
      for(Int_t iI = 0; iI < nID; ++iI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){	 
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    for(Int_t sI = 0; sI < nSyst; ++sI){

 	      for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){

		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		std::string histName = jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][0]->GetName();
		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		bool highLight = false;
		if(histName.find("akCs3PU3PFFlowJetAnalyzer") == std::string::npos) highLight = false;
		else if(histName.find("NoID") == std::string::npos) highLight = false;
		else if(histName.find("ResponseMod0p00") == std::string::npos) highLight = false;
		else if(histName.find("AbsEta0p5to1p0") == std::string::npos) highLight = false;
		else if(histName.find("JECUpMC") == std::string::npos) highLight = false;
		
		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;

		RooUnfoldResponse* rooResSuperClone_p = (RooUnfoldResponse*)rooResponse_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI]->Clone("rooResSuperClone_p");

		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		TH2D* initRes_p = (TH2D*)rooResSuperClone_p->Hresponse()->Clone("initRes_p");

		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		TH1D* initMeas_p = (TH1D*)rooResSuperClone_p->Hmeasured()->Clone("initMeas_p");
		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		TH1D* initTrue_p = (TH1D*)rooResSuperClone_p->Htruth()->Clone("initTrue_p");
		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		if(doLocalDebug || doGlobalDebug) std::cout << jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->GetName() << std::endl;

		TH1D* rawClone_p = (TH1D*)jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->Clone("rawClone_p"); // Using a clone so we can modify for Fake err;
		
		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		
		//NOTE: ERRORS HERE ARE HARD CODED RELATIVE VALUES BASED ON PLOTS IN AN
		if(isStrSame(systStr[sI], "Fake") && !isDataPP){
		  if(tempStr.find("akCs10PU3PFFlow") != std::string::npos){
		    if(centBinsLow[cI] == 0 && centBinsHi[cI] == 10) correctForFakes(rawClone_p, {100, 150, 200, 250}, {0.425, 0.3, 0.21});
		    else if(centBinsLow[cI] == 10 && centBinsHi[cI] == 30) correctForFakes(rawClone_p, {100, 150, 200}, {0.08, 0.075});
		  }
		  else if(tempStr.find("akCs8PU3PFFlow") != std::string::npos){
		    if(centBinsLow[cI] == 0 && centBinsHi[cI] == 10) correctForFakes(rawClone_p, {100, 150, 200, 250}, {0.31, 0.16, 0.04});
		    else if(centBinsLow[cI] == 10 && centBinsHi[cI] == 30) correctForFakes(rawClone_p, {100, 150}, {0.04}); 
		  }
		  else if(tempStr.find("akCs6PU3PFFlow") != std::string::npos){
		    if(centBinsLow[cI] == 0 && centBinsHi[cI] == 10) correctForFakes(rawClone_p, {100, 150, 200}, {0.07, 0.07}); //taking max since the two bins are not descending and statistically compatible
		  }
		  else if(tempStr.find("akCs4PU3PFFlow") != std::string::npos){
		    if(centBinsLow[cI] == 0 && centBinsHi[cI] == 10) correctForFakes(rawClone_p, {100, 150}, {0.07});
		  }
		}
		
		
		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		
		if(highLight){
		  std::cout << "DOING PRE-HIGHLIGHT" << std::endl;
		  std::cout << " Print initTrue: " << std::endl;
		  initTrue_p->Print("ALL");
		  
		  std::cout << " Print initMeas: " << std::endl;
		  initMeas_p->Print("ALL");
		  
		  std::cout << " Print from rooUnfold: " << std::endl;
		  rooResSuperClone_p->Print("ALL");
		  
		  std::cout << " Print from TH2" << std::endl;
		  for(Int_t bIY = 0; bIY < initRes_p->GetYaxis()->GetNbins(); ++bIY){
		    std::string firstString = std::to_string(bIY+1);
		    if(firstString.size() == 0) firstString = " " + firstString;
		    std::cout << firstString;
		    
		    for(Int_t bIX = 0; bIX < initRes_p->GetXaxis()->GetNbins(); ++bIX){
		      std::string value = prettyString(initRes_p->GetBinContent(bIX+1, bIY+1), 6, false);
		      
		      while(value.size() < 10){value = " " + value;}
		      while(value.size() > 10){value.replace(value.size()-2, 1, "");}
		      value = "'" + value + "'";
		      std::cout << value;
		    }
		    std::cout << std::endl;
		  }
		}
	      

		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		
		for(Int_t bsI = 0; bsI < nSuperBayes; ++bsI){
		  RooUnfoldBayes superBayes(rooResSuperClone_p, rawClone_p, 3, false, "name");
		  superBayes.SetVerbose(-1);
		  TH1D* unfold_h = (TH1D*)superBayes.Hreco(RooUnfold::kCovToy);
		  Double_t tot = unfold_h->Integral();
		  unfold_h->Scale(1./tot);
		  
		  for(Int_t bIY = 0; bIY < initRes_p->GetNbinsY(); ++bIY){
		    Double_t sumOverX = 0.0;
		    for(Int_t bIX = 0; bIX < initRes_p->GetNbinsX(); ++bIX){
		      sumOverX += initRes_p->GetBinContent(bIX+1, bIY+1);
		    }
		    
		    Double_t scaleFactor = 1;
		    if(sumOverX > 0) scaleFactor = unfold_h->GetBinContent(bIY+1)/sumOverX;
		    
		    //		  std::cout << "Check bin width match: " << unfold_h->GetBinLowEdge(bIY+1) << "-" << unfold_h->GetBinLowEdge(bIY+2) << ", " << initRes_p->GetYaxis()->GetBinLowEdge(bIY+1) << "-" << initRes_p->GetYaxis()->GetBinLowEdge(bIY+2) << "." << std::endl;

		    for(Int_t bIX = 0; bIX < initRes_p->GetNbinsX(); ++bIX){
		      initRes_p->SetBinContent(bIX+1, bIY+1, initRes_p->GetBinContent(bIX+1, bIY+1)*scaleFactor);
		      initRes_p->SetBinError(bIX+1, bIY+1, initRes_p->GetBinError(bIX+1, bIY+1)*scaleFactor);
		    }
		  }
		  
		  delete rooResSuperClone_p;
		  rooResSuperClone_p = NULL;
		  rooResSuperClone_p = new RooUnfoldResponse(initMeas_p, unfold_h, initRes_p, "rooResSuperClone_p");
		}
		
		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		
		if(highLight){
		  std::cout << "DOING POST-HIGHLIGHT" << std::endl;
		  std::cout << " Print from roo" << std::endl;
		  rooResSuperClone_p->Print("ALL");
		  
		  std::cout << " Print from TH2" << std::endl;
		  for(Int_t bIY = 0; bIY < initRes_p->GetYaxis()->GetNbins(); ++bIY){
		    std::string firstString = std::to_string(bIY+1);
		    if(firstString.size() == 0) firstString = " " + firstString;
		    std::cout << firstString;
		    
		    for(Int_t bIX = 0; bIX < initRes_p->GetXaxis()->GetNbins(); ++bIX){
		      std::string value = prettyString(initRes_p->GetBinContent(bIX+1, bIY+1), 6, false);
		      
		      while(value.size() < 10){value = " " + value;}
		      while(value.size() > 10){value.replace(value.size()-2, 1, "");}
		      value = "'" + value + "'";
		      std::cout << value;
		    }
		    std::cout << std::endl;
		  }
		}
		
		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		
		delete initRes_p;
		delete initMeas_p;
		delete initTrue_p;
		
		for(Int_t bI = 0; bI < nBayes; ++bI){	    		
		  RooUnfoldResponse* rooResClone_p = (RooUnfoldResponse*)rooResSuperClone_p->Clone("rooResClone_p");

		  //		  if(jI == 0 && cI == 0 && iI == 2 && mI == 1 && aI == nJtAbsEtaBins-1 && sI == 0 && binsI == 1){
		  if(jI == 0 && cI == 0 && iI == 0 && mI == 0 && aI == 0 && sI == 0 && binsI == 0 && bI == 4){
		    std::cout << "PRINT UNFOLDING RAW: " << std::endl;
		    rawClone_p->Print("ALL");
		    std::cout << "PRINT UNFOLDING RAW 2: " << std::endl;
		    jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->Print("ALL");
		  }

		  if(cI == 0){
		    std::cout << "PRINTNOW A: " << std::endl;
		    jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->Print("ALL");
		    std::cout << "PRINTNOW B: " << std::endl;
		    rooResClone_p->Print("ALL");
		  }
		  
		  bayes_p[jI][cI][iI][mI][aI][sI][binsI][bI] = new RooUnfoldBayes(rooResClone_p, rawClone_p, bayesVal[bI], false, ("name_" + std::to_string(bI)).c_str());
		  bayes_p[jI][cI][iI][mI][aI][sI][binsI][bI]->SetVerbose(-1);	    
		  TH1D* unfold_h = (TH1D*)bayes_p[jI][cI][iI][mI][aI][sI][binsI][bI]->Hreco(RooUnfold::kCovToy);	  

		  if(cI == 0){
		    std::cout << "PRINTNOW C: " << std::endl;
		    unfold_h->Print("ALL");
		  }
		  
		  for(Int_t bIX = 0; bIX < unfold_h->GetNbinsX(); ++bIX){
		    jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->SetBinContent(bIX+1, unfold_h->GetBinContent(bIX+1));
		    jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->SetBinError(bIX+1, unfold_h->GetBinError(bIX+1));
		  }

		  if(jI == 0 && cI == 0 && iI == 0 && mI == 0 && aI == 0 && sI == 0 && binsI == 0 && bI == 4){
		    std::cout << "PRINT UNFOLDING UNFOLDED 1: " << std::endl;
		    unfold_h->Print("ALL");
		    std::cout << "PRINT UNFOLDING UNFOLDED 2: " << std::endl;
		    jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->Print("ALL");
		  }

		  delete rooResClone_p;
		}
		
		if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " <<  __LINE__ << std::endl;
		
		delete rawClone_p;
		delete rooResSuperClone_p;
	      }
	    }
	  }
	}
      }
    }
  }

  std::cout << "Unfolding complete." << std::endl;

  dataFile_p->Close();
  delete dataFile_p;

  responseFile_p->Close();
  delete responseFile_p;

  outFile_p->cd();

  const Double_t yPadFrac = 0.45;
  const Double_t marg = 0.12;

  const Int_t nStyles = 6;
  const Int_t styles[nStyles] = {24, 25, 27, 28, 46, 44};

  const Int_t nColors = 5;
  const Int_t colors[nColors] = {1, vg.getColor(0), vg.getColor(1), vg.getColor(2), vg.getColor(4)};

  std::vector<std::vector<std::string> > pdfNames;

  Int_t nPDFTotal = 0;

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    outFile_p->cd();
    dir_p[jI]->cd();

    std::string tempStr = responseJetDirList[jI];
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");


    Double_t lowPtTruncVal = 200.;
    Bool_t isBigJet = tempStr.find("ak8") != std::string::npos || tempStr.find("ak10") != std::string::npos || tempStr.find("akCs8") != std::string::npos || tempStr.find("akCs10") != std::string::npos;
    if(isBigJet) lowPtTruncVal = 300.;
    
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      std::string centStr2 = "PP";
      if(!isDataPP){
	centStr = "PbPb_" + centStr;
	centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";
      }
      else centStr = "PP_" + centStr;


      for(Int_t iI = 0; iI < nID; ++iI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){	
	  const std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);

	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	    //	    if(jtAbsEtaStr.find("AbsEta0p0to2p0") == std::string::npos) continue;
	    const std::string jtAbsEtaStr2 = prettyString(jtAbsEtaBinsLow[aI], 1, false) + "<|#eta|<" +  prettyString(jtAbsEtaBinsHi[aI], 1, false);

	    
	    pdfNames.push_back({});

	    for(Int_t sI = 0; sI < nSyst; ++sI){	    	   
              const std::string tempSystStr = systStr[sI] + "_";

	      for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
		TLegend* leg_p = new TLegend(0.15, 0.05, 0.5, 0.6);
		leg_p->SetBorderSize(0.0);
		leg_p->SetFillStyle(0);
		leg_p->SetFillColor(0);
		leg_p->SetTextFont(43);
		leg_p->SetTextSize(14);
		
		TLatex* label_p = new TLatex();
		label_p->SetTextFont(43);
		label_p->SetTextSize(14);
		label_p->SetNDC();                  
		
		const Int_t nPads = 11;
		const Double_t nPadsX = 4;
		const Double_t nPadsY = 2;
		
		Double_t padWidth = 1./nPadsX;
		Double_t padHeight = 1./nPadsY;
		
		
		const Double_t padsXLow[nPads] = {0*padWidth, 
						  0*padWidth, 
						  0*padWidth, 
						  1*padWidth, 
						  1*padWidth, 
						  2*padWidth, 
						  3*padWidth, 
						  0*padWidth, 
						  1*padWidth, 
						  2*padWidth,
						  3*padWidth};
		
		const Double_t padsXHi[nPads] = {1*padWidth, 
						 1*padWidth, 
						 1*padWidth, 
						 2*padWidth, 
						 2*padWidth, 
						 3*padWidth, 
						 4*padWidth, 
						 1*padWidth, 
						 2*padWidth, 
						 3*padWidth, 
						 4*padWidth};
		
		const Double_t padsYLow[nPads] = {0.5 + 0.5*yPadFrac,
						  0.5 + 0.5*yPadFrac - (0.5*yPadFrac - marg*nPadsY/nPadsX)/2, 
						  1*padHeight, 
						  0.5 + 0.5*yPadFrac - (0.5*yPadFrac - marg*nPadsY/nPadsX)/2, 
						  1*padHeight,
						  1*padHeight, 
						  1*padHeight, 
						  0*padHeight, 
						  0*padHeight, 
						  0*padHeight, 
						  0*padHeight};
		
		const Double_t padsYHi[nPads] = {2*padHeight,
						 0.5 + 0.5*yPadFrac, 
						 0.5 + 0.5*yPadFrac - (0.5*yPadFrac - marg*nPadsY/nPadsX)/2, 
						 2*padHeight, 
						 0.5 + 0.5*yPadFrac - (0.5*yPadFrac - marg*nPadsY/nPadsX)/2,
						 2*padHeight,
						 2*padHeight, 
						 1*padHeight, 
						 1*padHeight,
						 1*padHeight, 
						 1*padHeight};
		
		const Double_t bottomMarg[nPads] = {0.001, 0.001, marg*(nPadsY/nPadsX)/(padsYHi[2] - padsYLow[2]), 0.001, marg*(nPadsY/nPadsX)/(padsYHi[4] - padsYLow[4]), marg*(nPadsY/nPadsX)/(padsYHi[5] - padsYLow[5]), marg*(nPadsY/nPadsX)/(padsYHi[6] - padsYLow[6]), marg*(nPadsY/nPadsX)/(padsYHi[7] - padsYLow[7]), marg*(nPadsY/nPadsX)/(padsYHi[8] - padsYLow[8]), marg*(nPadsY/nPadsX)/(padsYHi[9] - padsYLow[9]), marg*(nPadsY/nPadsX)/(padsYHi[10] - padsYLow[10])};
		
		TCanvas* canv_p = new TCanvas("canv_p", "", 450*4, 450*2);
		TPad* pads[nPads];
		for(Int_t padI = 0; padI < nPads; ++padI){
		  canv_p->cd();
		  pads[padI] = new TPad(("pad" + std::to_string(padI)).c_str(), "", padsXLow[padI], padsYLow[padI], padsXHi[padI], padsYHi[padI]);
		  pads[padI]->SetLeftMargin(marg);
		  pads[padI]->SetTopMargin(0.01);
		  pads[padI]->SetBottomMargin(bottomMarg[padI]);
		  pads[padI]->SetRightMargin(0.002);
		  pads[padI]->Draw();
		}
		
		//start here for terminalpos purposes
		canv_p->cd();
		pads[2]->cd();
		
		int terminalPos = -1;
		int terminalPos5 = -1;
		int terminalPos5AndPears = -1;
		double terminalSumMin = 999999;
	      				
		TH1D* bandValLow_p = (TH1D*)jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayes100Pos]->Clone("bandValLow");
		TH1D* bandValHi_p = (TH1D*)jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayes100Pos]->Clone("bandValHi");
	   
		for(Int_t bI = bayes100Pos-nBigBayesSymm; bI <= bayes100Pos+nBigBayesSymm; ++bI){
		  for(Int_t bIX = 0; bIX < bandValLow_p->GetNbinsX(); ++bIX){
		    if(jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetBinContent(bIX+1) < bandValLow_p->GetBinContent(bIX+1)){
		      bandValLow_p->SetBinContent(bIX+1, jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetBinContent(bIX+1));
		    }
		    if(jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetBinContent(bIX+1) > bandValHi_p->GetBinContent(bIX+1)){
		      bandValHi_p->SetBinContent(bIX+1, jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetBinContent(bIX+1));
		    }
		  }
		}
		
		TH1D* clones_p[nBayes];
		
		Double_t max = -1;
		Double_t min = 100;
		for(Int_t bI = 0; bI < nBayes; ++bI){
		  clones_p[bI] = (TH1D*)jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->Clone(("clone_" + std::to_string(bI)).c_str());

		  bool sub1Perc = true;
		  bool sub5Perc = true;
		  for(Int_t bIX = 0; bIX < clones_p[bI]->GetNbinsX(); ++bIX){
		    double binCenter = (clones_p[bI]->GetBinLowEdge(bIX+1) + clones_p[bI]->GetBinLowEdge(bIX+2))/2.;
		    double binContent = clones_p[bI]->GetBinContent(bIX+1);
		    double bandLowContent = bandValLow_p->GetBinContent(bIX+1);
		    double bandHiContent = bandValHi_p->GetBinContent(bIX+1);
		    
		    if(binCenter > lowPtTruncVal && binCenter < 1000.){
		      if(sub1Perc){
			if(binContent/bandHiContent > 1.01) sub1Perc = false;
			else if(binContent/bandLowContent < .99) sub1Perc = false;
		      }
		      if(sub5Perc){
			if(binContent/bandHiContent > 1.05) sub5Perc = false;
			else if(binContent/bandLowContent < .95) sub5Perc = false;
		      }
		    }
		  }
		  
		  if(sub1Perc && terminalPos < 0 && bI < nBayesBig){
		    terminalPos = bI;
		    //		    histTermPos[jI][cI][iI][mI][aI][sI][binsI] = terminalPos;
		  }
		  if(sub5Perc && terminalPos5 < 0 && bI < nBayesBig) terminalPos5 = bI;
		  if(sub5Perc && bI < nBayesBig){
		    TH2D* covarianceToPlot = NULL;
		    getPearsTMatrix(bayes_p[jI][cI][iI][mI][aI][sI][binsI][bI], &covarianceToPlot);

		    bool keepMinOppDiag = true;
		    Double_t sum = 0;
		    for(Int_t bIX = 0; bIX < covarianceToPlot->GetNbinsX(); ++bIX){
		      for(Int_t bIY = 0; bIY < covarianceToPlot->GetNbinsY(); ++bIY){
			if(bIX == bIY) continue;//continue on diagonal

			Double_t val = covarianceToPlot->GetBinContent(bIX+1, bIY+1);
			if(bIX == 0 && bIY == covarianceToPlot->GetNbinsY()-1 && TMath::Abs(val) > 0.15) keepMinOppDiag = false;

			sum += TMath::Abs(val);
		      }
		    }

		    if(sum < terminalSumMin && keepMinOppDiag){
		      terminalSumMin = sum;
		      terminalPos5AndPears = bI;
		      histTermPos[jI][cI][iI][mI][aI][sI][binsI] = terminalPos5AndPears;
		    }

		    delete covarianceToPlot;
		  }
	      		
		  clones_p[bI]->Divide(jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayes100Pos]);
		  Double_t tempMax = getMax(clones_p[bI]);
		  Double_t tempMin = getMinGTZero(clones_p[bI]);
		  if(tempMax > max) max = tempMax;
		  if(tempMin < min) min = tempMin;
		}
		
		if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
		
		delete bandValLow_p;
		delete bandValHi_p;
		
		if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
		
		Double_t interval = (max - min)/10.;
		max += interval;
		min -= interval;
		
		const Int_t nTempBins = 10;
		Double_t tempBins[nTempBins+1];
		getLinBins(min, max, nTempBins, tempBins);
		
		if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

		if(max > 1.15){
		  max = 1.15;
		  interval = (max - min)/10.;
		}
		if(min < .85){
		  min = 0.85;
		  interval = (max - min)/10.;
		}
			       		
		for(Int_t bI = 1; bI < nBayes; ++bI){
		  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << bI << std::endl;
		  
		  clones_p[bI]->SetMaximum(max);
		  clones_p[bI]->SetMinimum(min);
		  
		  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << bI << std::endl;
		  
		  clones_p[bI]->GetYaxis()->SetTitle("Ratio w/ Bayes 100");
		  clones_p[bI]->GetYaxis()->SetTitleSize(9);
		  clones_p[bI]->GetYaxis()->SetLabelSize(9);
		  clones_p[bI]->GetYaxis()->SetTitleOffset(clones_p[bI]->GetYaxis()->GetTitleOffset()*1.5);
		  clones_p[bI]->GetYaxis()->SetNdivisions(505);
		  
		  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << bI << std::endl;
		  
		  if(bI < nBayesDraw){
		    Int_t bayesPos = bI;
		    if(bI == nBayesDraw-1 && histTermPos[jI][cI][iI][mI][aI][sI][binsI] >= nBayesDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI][binsI];
		    
		    clones_p[bayesPos]->SetMarkerStyle(styles[bayesPos%nStyles]);
		    clones_p[bayesPos]->SetMarkerColor(colors[bayesPos%nColors]);
		    clones_p[bayesPos]->SetLineColor(colors[bayesPos%nColors]);
		    clones_p[bayesPos]->SetMarkerSize(1.0);
		    
		    clones_p[bayesPos]->GetXaxis()->SetTitleFont(43);
		    clones_p[bayesPos]->GetYaxis()->SetTitleFont(43);
		    clones_p[bayesPos]->GetXaxis()->SetLabelFont(43);
		    clones_p[bayesPos]->GetYaxis()->SetLabelFont(43);
		    
		    clones_p[bayesPos]->GetXaxis()->SetTitleSize(14);
		    clones_p[bayesPos]->GetYaxis()->SetTitleSize(14);
		    clones_p[bayesPos]->GetXaxis()->SetLabelSize(14);
		    clones_p[bayesPos]->GetYaxis()->SetLabelSize(14);
		    
		    clones_p[bayesPos]->GetXaxis()->SetTitleOffset(5.);
		    clones_p[bayesPos]->GetYaxis()->SetTitleOffset(1.5);
		  		    
		    if(bayesPos == 1) clones_p[bayesPos]->DrawCopy("HIST E1 P");
		    else clones_p[bayesPos]->DrawCopy("HIST E1 P SAME");
		  }
		}
		
		bool doLogX = false;
		if(jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][0]->GetBinWidth(1)*3 < jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][0]->GetBinWidth(jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][0]->GetNbinsX()-1)) doLogX = true;
		
		gStyle->SetOptStat(0);
		if(doLogX) gPad->SetLogx();
		
		for(Int_t bI = 0; bI < nBayes; ++bI){
		  delete clones_p[bI];
		  clones_p[bI] = NULL;
		}
		
		drawWhiteBox(900, 1100, .00, min*.999);

		canv_p->cd();	       
		pads[0]->cd();
	      
		min = 1000000000;
		max = -1;
		
		for(Int_t bI = 0; bI < nBayesDraw; ++bI){
		  Int_t bayesPos = bI;
		  if(bI == nBayesDraw-1 && histTermPos[jI][cI][iI][mI][aI][sI][binsI] >= nBayesDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI][binsI];
		  
		  Double_t tempMax = getMax(jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayesPos]);
		  Double_t tempMin = getMinGTZero(jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayesPos]);
		  if(tempMax > max) max = tempMax;
		  if(tempMin < min) min = tempMin;
		}
		
		Double_t globalMin = getNearestFactor10Down(min, 3);
		Double_t globalMax = getNearestFactor10Up(max, 2);
		
		for(Int_t bI = 0; bI < nBayes; ++bI){
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->SetMaximum(globalMax);
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->SetMinimum(globalMin);
		  
		  
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetXaxis()->SetTitleFont(43);
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetYaxis()->SetTitleFont(43);
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetXaxis()->SetLabelFont(43);
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetYaxis()->SetLabelFont(43);
		  
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetXaxis()->SetTitleSize(14);
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetYaxis()->SetTitleSize(14);
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetXaxis()->SetLabelSize(14);
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetYaxis()->SetLabelSize(14);

		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetXaxis()->SetTitleOffset(5.);
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->GetYaxis()->SetTitleOffset(1.5);
		  
		  if(bI < nBayesDraw){
		    Int_t bayesPos = bI;
		    if(bI == nBayesDraw-1 && histTermPos[jI][cI][iI][mI][aI][sI][binsI] >= nBayesDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI][binsI];
		
		    jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayesPos]->SetMarkerStyle(styles[bayesPos%nStyles]);
		    jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayesPos]->SetMarkerColor(colors[bayesPos%nColors]);
		    jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayesPos]->SetLineColor(colors[bayesPos%nColors]);
		    jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayesPos]->SetMarkerSize(1.0);
		    
		    if(bayesPos == 0){
		      jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayesPos]->DrawCopy("HIST E1 P");
		      
		      jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetMarkerStyle(1);
		      jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetMarkerSize(0.001);
		      jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetMarkerColor(0);
		      jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetLineColor(1);
		      jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetLineWidth(2);
		      jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->DrawCopy("HIST E1 SAME");
		      
		      Double_t rescale = jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->Integral()/recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->Integral();
		      
		      recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->Scale(rescale);
		      recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetMarkerStyle(1);
		      recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetMarkerSize(0.001);
		      recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetMarkerColor(0);
		      recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetLineColor(4);
		      recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetLineStyle(2);
		      recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetLineWidth(2);
		      recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->DrawCopy("HIST E1 SAME");
		      
		      leg_p->AddEntry(jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI], "Folded", "L");
		      leg_p->AddEntry(recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI], "Folded MC", "L");
		    }
		    else jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayesPos]->DrawCopy("HIST E1 P SAME"); 		
		    leg_p->AddEntry(jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayesPos], ("Bayes=" + std::to_string(bayesVal[bayesPos])).c_str(), "P L");
		    
		  }
		}
		
		canv_p->cd();
		pads[0]->cd();
		gStyle->SetOptStat(0);
		gPad->SetLogy();
		gPad->SetTicks(1,2);
		
		if(doLogX) gPad->SetLogx();
		
		leg_p->Draw("SAME");
		
		label_p->DrawLatex(0.55, 0.94, tempStr.c_str());
		label_p->DrawLatex(0.55, 0.88, centStr2.c_str());
		label_p->DrawLatex(0.55, 0.82, jtAbsEtaStr2.c_str());
		label_p->DrawLatex(0.55, 0.76, resStr.c_str());
		label_p->DrawLatex(0.55, 0.70, idStr[iI].c_str());
		label_p->DrawLatex(0.55, 0.64, systStr[sI].c_str());
		
		canv_p->cd();
		pads[1]->cd();
				
		max = -1;
		min = 100;
		
		for(Int_t bI = 0; bI < nBayes; ++bI){
		  clones_p[bI] = (TH1D*)jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->Clone(("clone_" + std::to_string(bI)).c_str());
		  clones_p[bI]->Divide(jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][0]);
		  Double_t tempMax = getMax(clones_p[bI]);
		  Double_t tempMin = getMinGTZero(clones_p[bI]);
		  if(max < tempMax) max = tempMax;
		  if(min > tempMin) min = tempMin;
		}
		
		interval = (max - min)/10.;
		max += interval;
		min -= interval;
	      
		if(max > 1.15) max = 1.15;
		if(min < .85) min = 0.85;
		
		for(Int_t bI = 0; bI < nBayes; ++bI){
		  clones_p[bI]->SetMaximum(max);
		  clones_p[bI]->SetMinimum(min);
		  
		  clones_p[bI]->GetYaxis()->SetTitle("Ratio w/ Bayes0");
		  clones_p[bI]->GetYaxis()->SetTitleSize(9);
		  clones_p[bI]->GetYaxis()->SetLabelSize(9);
		  clones_p[bI]->GetYaxis()->SetTitleOffset(clones_p[bI]->GetYaxis()->GetTitleOffset()*1.5);
		  clones_p[bI]->GetYaxis()->SetNdivisions(505);	       		
		  
		  if(bI < nBayesDraw){
		    Int_t bayesPos = bI;
		    if(bI == nBayesDraw-1 && histTermPos[jI][cI][iI][mI][aI][sI][binsI] >= nBayesDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI][binsI];
		    
		    if(bayesPos == 0) clones_p[bayesPos]->DrawCopy("HIST E1 P");
		    else clones_p[bayesPos]->DrawCopy("HIST E1 P SAME");
		  }
		}
		
		gStyle->SetOptStat(0);
		if(doLogX) gPad->SetLogx();
		
		for(Int_t bI = 0; bI < nBayes; ++bI){
		  delete clones_p[bI];
		  clones_p[bI] = NULL;
		}
		
		
		canv_p->cd();
		pads[3]->cd();
		
		//Refold EDITING HERE		
		jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetMaximum(globalMax*20.);
		jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetMinimum(globalMin/20.);
		jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->DrawCopy("HIST E1");
		
		label_p->DrawLatex(0.35, .94, "Refolding + Comp. to data");
		
		for(Int_t bI = 0; bI < nBayesDraw; ++bI){
		  Int_t bayesPos = bI;
		  if(bI == nBayesDraw-1 && histTermPos[jI][cI][iI][mI][aI][sI][binsI] >= nBayesDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI][binsI];
		  
		  canv_p->cd();
		  pads[3]->cd();
		  
		  TH1D* tempClone = (TH1D*)rooResponse_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI]->ApplyToTruth(jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bayesPos], "tempClone");

		  centerTitles(tempClone);
		  tempClone->SetMarkerStyle(styles[bayesPos%nStyles]);
		  tempClone->SetMarkerColor(colors[bayesPos%nColors]);
		  tempClone->SetLineColor(colors[bayesPos%nColors]);
		  tempClone->SetMarkerSize(1.0);
		  
		  tempClone->GetXaxis()->SetTitleFont(43);
		  tempClone->GetYaxis()->SetTitleFont(43);
		  tempClone->GetXaxis()->SetLabelFont(43);
		  tempClone->GetYaxis()->SetLabelFont(43);
		  
		  tempClone->GetXaxis()->SetTitleSize(14);
		  tempClone->GetYaxis()->SetTitleSize(14);
		  tempClone->GetXaxis()->SetLabelSize(14);
		  tempClone->GetYaxis()->SetLabelSize(14);
		  
		  tempClone->GetXaxis()->SetTitleOffset(5.);
		  tempClone->GetYaxis()->SetTitleOffset(1.5);
		  
		  tempClone->SetTitle("");
		  
		  tempClone->DrawCopy("HIST E1 P SAME");
		  
		  canv_p->cd();
		  pads[4]->cd();

		  tempClone->Divide(jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]);
		  
		  if(bayesPos == 0){
		    tempClone->SetMaximum(1.25);
		    tempClone->SetMinimum(0.75);
		    tempClone->GetYaxis()->SetNdivisions(505);
		    tempClone->GetYaxis()->SetTitle("Ratio");
		    tempClone->DrawCopy("HIST E1 P");
		  }
		  else tempClone->DrawCopy("HIST E1 P SAME");
		  
		  delete tempClone;
		}
		
		drawWhiteBox(900, 1100, .00, min*.999);
		
		label_p->SetNDC(0);
		//	      label_p->DrawLatex(100, min - interval*3, "100");
		if(!isBigJet) label_p->DrawLatex(200, .77, "200");
		label_p->DrawLatex(400, .77, "400");
		label_p->DrawLatex(600, .77, "600");
		label_p->DrawLatex(1000, .77, "1000");
		label_p->SetNDC(1);
	      		
		const Int_t nPearsToDraw = 6;
		for(Int_t bI = 0; bI < nPearsToDraw; ++bI){
		  Int_t bayesPos = bI;
		  if(bI == nPearsToDraw-1){
		    if(histTermPos[jI][cI][iI][mI][aI][sI][binsI] >= nPearsToDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI][binsI];
		    else continue;
		  }

		  canv_p->cd();
		  pads[5+bI]->cd();

		  TH2D* covarianceToPlot = NULL;
		  getPearsTMatrix(bayes_p[jI][cI][iI][mI][aI][sI][binsI][bayesPos], &covarianceToPlot);
		    
		  covarianceToPlot->SetMarkerSize(1.75);		  
		  gStyle->SetPaintTextFormat("1.3f");
		  
		  std::vector<Float_t> tempVals;
		  for(Int_t bIX = 0; bIX < covarianceToPlot->GetNbinsX(); ++bIX){
		    for(Int_t bIY = 0; bIY < covarianceToPlot->GetNbinsY(); ++bIY){
		      Double_t val = covarianceToPlot->GetBinContent(bIX+1, bIY+1);
		      tempVals.push_back(val);
		      covarianceToPlot->SetBinContent(bIX+1, bIY+1, TMath::Abs(val));
		    }
		  }

		  covarianceToPlot->SetMaximum(1.05);
		  covarianceToPlot->SetMinimum(-0.05);
		  covarianceToPlot->DrawCopy("COL");
		  label_p->SetNDC(0);
		  unsigned int tempPos = 0;
		  for(Int_t bIX = 0; bIX < covarianceToPlot->GetNbinsX(); ++bIX){
		    Double_t xLow = covarianceToPlot->GetXaxis()->GetBinLowEdge(bIX+1);
		    Double_t xCent = covarianceToPlot->GetXaxis()->GetBinCenter(bIX+1);
                    for(Int_t bIY = 0; bIY < covarianceToPlot->GetNbinsY(); ++bIY){
		      Double_t yCent = covarianceToPlot->GetYaxis()->GetBinCenter(bIY+1);
		      if(TMath::Abs(tempVals[tempPos]) > 0.15) label_p->DrawLatex((xCent + xLow)/2., yCent, ("#color[2]{" + prettyString(tempVals[tempPos], 3, false) + "}").c_str());
		      else label_p->DrawLatex((xCent + xLow)/2., yCent, prettyString(tempVals[tempPos], 3, false).c_str());
		      ++tempPos;
		    }
		  }

		  label_p->SetNDC();
		  label_p->DrawLatex(0.35, .94, ("#bf{Pearson, Iteration " + std::to_string(bayesVal[bayesPos])+ "}").c_str());
		  
		  //	      gPad->SetLogx();
		  //	      gPad->SetLogy();
		  
		  drawWhiteBox(-0.6, covarianceToPlot->GetXaxis()->GetBinLowEdge(covarianceToPlot->GetNbinsX()+1)+.25, -0.6, -.01);
		  drawWhiteBox(-0.6, -.01, -0.6, covarianceToPlot->GetXaxis()->GetBinLowEdge(covarianceToPlot->GetNbinsX()+1)+.25);
		  
		  label_p->SetNDC(0);
		  
		  Int_t binPos200 = -1;
		  Int_t binPos300 = -1;
		  Int_t binPos1000 = -1;
		  
		  for(Int_t bIX = 0; bIX < genJtPt_GoodReco_h[jI][cI][iI][mI][aI][sI][binsI]->GetNbinsX(); ++bIX){
		    Int_t binLowEdge = genJtPt_GoodReco_h[jI][cI][iI][mI][aI][sI][binsI]->GetBinLowEdge(bIX+1);
		    if(binLowEdge > 199 && binLowEdge < 201) binPos200 = bIX;
		    if(binLowEdge > 299 && binLowEdge < 301) binPos300 = bIX;
		    if(binLowEdge > 999 && binLowEdge < 1001) binPos1000 = bIX;
		    
		    label_p->DrawLatex(-0.7, 0.05 + bIX, std::to_string(int(genJtPt_GoodReco_h[jI][cI][iI][mI][aI][sI][binsI]->GetBinLowEdge(bIX+1))).c_str());
		    label_p->DrawLatex(0.05 + bIX, -0.3, std::to_string(int(genJtPt_GoodReco_h[jI][cI][iI][mI][aI][sI][binsI]->GetBinLowEdge(bIX+1))).c_str());
		  }
		  label_p->SetNDC(1);
		  
		  TLine* line_p = new TLine();
		  line_p->SetLineStyle(2);
		  line_p->SetLineColor(2);
		  line_p->SetLineWidth(line_p->GetLineWidth()*2);
		  
		  Int_t binLow = binPos200;
		  if(isBigJet) binLow = binPos300;
		  
		  line_p->DrawLine(binLow, binLow, binLow, binPos1000);
		  line_p->DrawLine(binLow, binPos1000, binPos1000, binPos1000);
		  line_p->DrawLine(binLow, binLow, binPos1000, binLow);
		  line_p->DrawLine(binPos1000, binLow, binPos1000, binPos1000);
		
		  delete line_p;
		  delete covarianceToPlot;
		}
		
		canv_p->cd();
		pads[3]->cd();
		gPad->SetLogy();
		if(doLogX) gPad->SetLogx();
		
		canv_p->cd();
		pads[4]->cd();
		if(doLogX) gPad->SetLogx();
		
		
		//	      label_p->SetNDC(0);
		canv_p->cd();
		pads[0]->cd();
		
		if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	      
		const std::string tempHistTag = tempStr + "_" + smallLargeBinsStr[binsI] + "_" + centStr + "_" + idStr[iI] + "_ResponseMod" + prettyString(responseMod[mI], 2, true) + "_" + jtAbsEtaStr + "_" + systStr[sI];
		histTag.push_back(tempHistTag);

		if(terminalPos5AndPears >= 0) histBestBayes.push_back(bayesVal[terminalPos5AndPears]);
		else histBestBayes.push_back(-1);
		/*
		if(terminalPos >= 0) histBestBayes.push_back(bayesVal[terminalPos]);
		else histBestBayes.push_back(-1);
		*/
	      
		if(terminalPos5 >= 0) label_p->DrawLatex(0.4, 0.33, ("Term. 5% at Bayes=" + std::to_string(bayesVal[terminalPos5])).c_str());
		else label_p->DrawLatex(0.4, 0.33, "Doesn't term. at 5% level");

		if(terminalPos5AndPears >= 0) label_p->DrawLatex(0.4, 0.26, ("Term. 5%+PearsMin at Bayes=" + std::to_string(bayesVal[terminalPos5AndPears])).c_str());
		else label_p->DrawLatex(0.4, 0.26, "Doesn't term. at 5% level+Pears");
		
		if(terminalPos >= 0){
		  //		label_p->DrawLatex(110, tempBins[8], ("Term. at Bayes=" + std::to_string(terminalPos+1)).c_str());
		  label_p->DrawLatex(0.4, 0.19, ("Term. 1% at Bayes=" + std::to_string(bayesVal[terminalPos])).c_str());
		  
		  Double_t maxDeltaCenter = -1;
		  Double_t maxDelta = 0;
		  Double_t lastDeltaCenter = -1;
		  Double_t lastDelta = 0;
		  for(Int_t bIX = 0; bIX < jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][terminalPos]->GetNbinsX(); ++bIX){
		    double center = (jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][terminalPos]->GetBinLowEdge(bIX+1) + jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][terminalPos]->GetBinLowEdge(bIX+2))/2.;
		    
		    if(center < lowPtTruncVal) continue;
		    if(center > 1000.) continue;
		    
		    double content1 = jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][terminalPos]->GetBinContent(bIX+1);
		    double content1Max = jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][nBayes-1]->GetBinContent(bIX+1);
		    double content1MaxMin1 = jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][nBayes-2]->GetBinContent(bIX+1);
		    
		    if(content1 == 0 && content1Max != 0){
		      maxDelta = 100;
		      maxDeltaCenter = center;
		    }
		    else if(TMath::Abs(content1 - content1Max)/content1 > maxDelta){
		      maxDelta = TMath::Abs(content1 - content1Max)/content1;
		      maxDeltaCenter = center;
		    }
		    
		    if(content1MaxMin1 == 0 && content1Max != 0){
		      lastDelta = 100;
		      lastDeltaCenter = center;
		    }
		    else if(TMath::Abs(content1MaxMin1 - content1Max)/content1MaxMin1 > lastDelta){
		      lastDelta = TMath::Abs(content1MaxMin1 - content1Max)/content1MaxMin1;
		      maxDeltaCenter = center;
		    }
		  }
		  
		  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
		  
		  std::cout << "MaxDelta: " << maxDelta << ", " << maxDeltaCenter << std::endl;
		  std::cout << "LastDelta: " << lastDelta << ", " << lastDeltaCenter << std::endl;
		  
		  //		label_p->DrawLatex(110, tempBins[6], ("MaxDelta: " + prettyString(maxDelta*100, 2, false) + "%").c_str());
		  //		label_p->DrawLatex(110, tempBins[4], ("LastDelta: " + prettyString(lastDelta*100, 2, false) + "%").c_str());
		  label_p->DrawLatex(0.4, 0.12, ("MaxDelta: " + prettyString(maxDelta*100, 2, false) + "%").c_str());
		  label_p->DrawLatex(0.4, 0.05, ("LastDelta: " + prettyString(lastDelta*100, 2, false) + "%").c_str());
		}
		else label_p->DrawLatex(0.4, 0.19, ("Doesn't term. for nBayes=" + std::to_string(nBayes+1)).c_str());
		//	      else label_p->DrawLatex(110, (max + min)/2., ("Doesn't term. for nBayes=" + std::to_string(nBayes+1)).c_str());
		
		canv_p->cd();
		pads[2]->cd();
		
		label_p->SetNDC(0);
		
		label_p->DrawLatex(100, min - interval*2, "100");
		label_p->DrawLatex(200, min - interval*2, "200");
		label_p->DrawLatex(400, min - interval*2, "400");
		label_p->DrawLatex(600, min - interval*2, "600");
		label_p->DrawLatex(1000, min - interval*2, "1000");

		for(Int_t pI = 0; pI < nPads; ++pI){
		  canv_p->cd();
		  pads[pI]->cd();
		  gPad->RedrawAxis();
		  gPad->SetTicks(1, 2);
		}
	    
		const std::string saveName = "jtPtUnfolded_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr +  "AllBayes_RecoGenAsymm_" + debugStr + dateStr + ".pdf";
	   
			
		std::string titleStr = tempStr + ", " + centStr2 + ", " + idStr[iI] + ", $" + jtAbsEtaStr2 + "$, " + tempSystStr;
		while(titleStr.find("#") != std::string::npos){titleStr.replace(titleStr.find("#"), 1, "\\");}
		slideTitlesPerAlgo[jI].push_back(titleStr);
		pdfPerSlidePerAlgo[jI].push_back({});
		pdfPerSlidePerAlgo[jI][pdfPerSlidePerAlgo[jI].size()-1].push_back(saveName);
		
		const std::string finalSaveName = "pdfDir/" + dateStr + "/Unfold_" + tempStr + "/" + saveName;
		quietSaveAs(canv_p, finalSaveName);
		++nPDFTotal;

		for(Int_t pI = 0; pI < nPads; ++pI){
		  delete pads[pI];
		}
		
		delete canv_p;
		delete leg_p;
		delete label_p;
	      }
	    }
	  }

	  
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    for(Int_t sI = 0; sI < nSyst; ++sI){
	      for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
		for(Int_t bI = 0; bI < nBayes; ++bI){		
		  jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI]->Write("", TObject::kOverwrite);
		  delete jtPtUnfolded_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][binsI][bI];
		  delete bayes_p[jI][cI][iI][mI][aI][sI][binsI][bI];
		}
	      }
	    }
	  }
	}
      }
    }
  }

  /*
  //TEX FILE
  const std::string textWidth = "0.24";
  std::string texFileName = outFileName;
  texFileName.replace(texFileName.find(".root"), std::string(".root").size(), ".tex");
  texFileName.replace(0, std::string("output").size(), "pdfDir/" + dateStr + "/Unfold");

  std::ofstream texFile(texFileName.c_str());

  texFile << "\\RequirePackage{xspace}" << std::endl;
  texFile << "\\RequirePackage{amsmath}" << std::endl;
  texFile << std::endl;

  texFile << "\\documentclass[xcolor=dvipsnames]{beamer}" << std::endl;
  texFile << "\\usetheme{Warsaw}" << std::endl;
  texFile << "\\setbeamercolor{structure}{fg=NavyBlue!90!NavyBlue}" << std::endl;
  texFile << "\\setbeamercolor{footlinecolor}{fg=white,bg=lightgray}" << std::endl;
  texFile << std::endl;

  texFile << "\\newcommand{\\pt}{\\ensuremath{p_{\\mathrm{T}}}\\xspace}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamersize{text margin left=3pt,text margin right=3pt}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamerfont{frametitle}{size=\\tiny}" << std::endl;
  texFile << "\\setbeamertemplate{frametitle}" << std::endl;
  texFile << "{" << std::endl;
  texFile << "    \\nointerlineskip" << std::endl;
  texFile << "    \\begin{beamercolorbox}[sep=0.05cm, ht=1.0em, wd=\\paperwidth]{frametitle}" << std::endl;
  texFile << "        \\vbox{}\\vskip-2ex%" << std::endl;
  texFile << "        \\strut\\insertframetitle\\strut" << std::endl;
  texFile << "        \\vskip-0.8ex%" << std::endl;
  texFile << "    \\end{beamercolorbox}" << std::endl;
  texFile << "}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamertemplate{footline}{%" << std::endl;
  texFile << "  \\begin{beamercolorbox}[sep=.4em,wd=\\paperwidth,leftskip=0.5cm,rightskip=0.5cm]{footlinecolor}" << std::endl;
  texFile << "    \\hspace{0.075cm}%" << std::endl;
  texFile << "    \\hfill\\insertauthor \\hfill\\insertpagenumber" << std::endl;
  texFile << "  \\end{beamercolorbox}%" << std::endl;
  texFile << "}" << std::endl;
  texFile << "\\setbeamertemplate{navigation symbols}{}" << std::endl;
  texFile << std::endl;
  
  texFile << "\\setbeamertemplate{itemize item}[circle]" << std::endl;
  texFile << "\\setbeamertemplate{itemize subitem}[circle]" << std::endl;
  texFile << "\\setbeamertemplate{itemize subsubitem}[circle]" << std::endl;
  texFile << "\\setbeamercolor{itemize item}{fg=black}" << std::endl;
  texFile << "\\setbeamercolor{itemize subitem}{fg=black}" << std::endl;
  texFile << "\\setbeamercolor{itemize subsubitem}{fg=black}" << std::endl;
  texFile << std::endl;

  texFile << "\\definecolor{links}{HTML}{00BFFF}" << std::endl;
  texFile << "\\hypersetup{colorlinks,linkcolor=,urlcolor=links}" << std::endl;
  texFile << std::endl;

  texFile << "\\author[CM]{Chris McGinn}" << std::endl;
  texFile << std::endl;

  texFile << "\\begin{document}" << std::endl;
  texFile << "\\begin{frame}" << std::endl;
  texFile << "\\frametitle{\\centerline{Unfolding termination (" << dateStr2 << ")}}" << std::endl;
  texFile << " \\begin{itemize}" << std::endl;
  texFile << "  \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "  \\item{Placeholder}" << std::endl;
  texFile << "  \\begin{itemize}" << std::endl;
  texFile << "   \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "   \\item{Placeholder}" << std::endl;
  texFile << "  \\end{itemize}" << std::endl;
  texFile << " \\end{itemize}" << std::endl;
  texFile << "\\end{frame}" << std::endl;


  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int spI = 0; spI < pdfNames.size(); ++spI){
    if(pdfNames.at(spI).at(0).find("AbsEta0p0to2p0") == std::string::npos) continue;

    std::string centStr = "PP";

    if(!isDataPP){
      centStr = pdfNames.at(spI).at(0).substr(pdfNames.at(spI).at(0).find("Cent"), pdfNames.at(spI).at(0).size());
      centStr.replace(centStr.find("_"), centStr.size(), "");
    }

    std::string jtStr = pdfNames.at(spI).at(0).substr(pdfNames.at(spI).at(0).find("_ak")+1, pdfNames.at(spI).at(0).size());
    jtStr.replace(jtStr.find("_"), jtStr.size(), "");

    std::string resStr = pdfNames.at(spI).at(0).substr(pdfNames.at(spI).at(0).find("_Response")+1, pdfNames.at(spI).at(0).size());
    resStr.replace(resStr.find("_"), resStr.size(), "");

    std::string idStr = "";
    if(pdfNames.at(spI).at(0).find("_NoID") != std::string::npos) idStr = "NoID";
    else if(pdfNames.at(spI).at(0).find("_LightMUID") != std::string::npos) idStr = "LightMUID";
    else if(pdfNames.at(spI).at(0).find("_LightMUAndCHID") != std::string::npos) idStr = "LightMUAndCHID";
    else if(pdfNames.at(spI).at(0).find("_FullLight") != std::string::npos) idStr = "FullLight";
    else if(pdfNames.at(spI).at(0).find("_FullTight") != std::string::npos) idStr = "FullTight";

    if(idStr.find("LightMUAndCHID") == std::string::npos) continue;
    if(resStr.find("0p10") == std::string::npos) continue;
    
    texFile << "\\begin{frame}" << std::endl;
    texFile << "\\frametitle{\\centerline{" << jtStr << ", " << centStr << ", " << resStr << ", " << idStr << "}}" << std::endl;

    for(unsigned int mI = 0; mI < pdfNames.at(spI).size(); ++mI){
      texFile << "\\includegraphics[width=" << textWidth << "\\textwidth]{" << pdfNames.at(spI).at(mI) << "}";
      if(mI == 3 || mI == 7) texFile << "\\\\";
      texFile << std::endl;
    }

    texFile << "\\begin{itemize}" << std::endl;
    texFile << "\\fontsize{8}{8}\\selectfont" << std::endl;
    texFile << "\\item{test}" << std::endl;
    texFile << "\\end{itemize}" << std::endl;
    texFile << "\\end{frame}" << std::endl;
  }


  texFile << "\\end{document}" << std::endl;
  texFile << std::endl;

  texFile.close();
  */

  std::cout << "Total nPDF: " << nPDFTotal << std::endl;
  
  for(Int_t jI = 0; jI < nDataJet; ++jI){
    texSlideCreator tex;

    std::string tempStr = responseJetDirList[jI];
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    tex.Clean();
    tex.Init(outFileName);
    tex.InitDir("pdfDir/" + dateStr + "/Unfold_" + tempStr + "/");
    tex.SetAuthor("Chris McGinn");

    tex.SetSlideTitles(slideTitlesPerAlgo[jI]);
    tex.SetSlidePdfs(pdfPerSlidePerAlgo[jI]);
    if(!(tex.CreateTexSlides())){
      std::cout << "Warning: .tex slide creation failed" << std::endl;
    }
  }

  outFile_p->cd();
  TDirectory* cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");
  TDirectory* subDir_p = (TDirectory*)cutDir_p->mkdir("subDir");
  TDirectory* unfoldDir_p = (TDirectory*)cutDir_p->mkdir("unfoldDir");

  if(nHistDim != (int)histTag.size()) std::cout << "WARNING: nHistDim (" << nHistDim << ") != histTag.size() (" << histTag.size() << ")" << std::endl;
  if(nHistDim != (int)histBestBayes.size()) std::cout << "WARNING: nHistDim (" << nHistDim << ") != histBestBayes.size() (" << histBestBayes.size() << ")" << std::endl;

  cutPropData.SetNHistDim(nHistDim);
  cutPropData.SetHistTag(histTag);
  cutPropData.SetHistBestBayes(histBestBayes);

  if(!cutPropData.WriteAllVarToFile(outFile_p, cutDir_p, subDir_p, unfoldDir_p)) std::cout << "Warning: Cut writing has failed" << std::endl;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3 && argc != 4){
    std::cout << "Usage: ./bin/unfoldRawData.exe <inDataFileName> <inResponseName> <selectJtAlgo-Opt>" << std::endl;
    return 1;
  }

  int retVal = 0;

  if(argc == 3) retVal += unfoldRawData(argv[1], argv[2]);
  else if(argc == 4) retVal += unfoldRawData(argv[1], argv[2], argv[3]);

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
