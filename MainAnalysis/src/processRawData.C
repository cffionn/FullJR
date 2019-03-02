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
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/doLocalDebug.h"
#include "MainAnalysis/include/smallOrLargeR.h"

//Non-local FullJR dependencies (Utility, etc.)
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/etaPhiFunc.h"
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vanGoghPalette.h"

int processRawData(const std::string inDataFileName, const std::string inResponseName, bool isDataPP = false, const std::string tagStr = "")
{
  cppWatch timer;
  timer.start();

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  const std::string dateStr2 = std::to_string(date->GetYear()) + "." + std::to_string(date->GetMonth()) + "." + std::to_string(date->GetDay());
  delete date;

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  checkMakeDir("pdfDir/" + dateStr + "/Process");

  vanGoghPalette vg;
  const Int_t nStyles = 5;
  const Int_t styles[nStyles] = {21, 24, 34, 33, 25};
  const Int_t colors[nStyles] = {1, vg.getColor(0), vg.getColor(1), vg.getColor(2), vg.getColor(3)};
  const Double_t yPadFrac = 0.4;
  const Double_t marg = 0.12;

  TFile* responseFile_p = new TFile(inResponseName.c_str(), "READ");
  std::vector<std::string> jetDirList = returnRootFileContentsList(responseFile_p, "TDirectoryFile", "JetAnalyzer");

  std::cout << "Validating " << jetDirList.size() << " jets..." << std::endl;
  for(unsigned int jI = 0; jI < jetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << jetDirList.size() << ": " << jetDirList[jI] << std::endl;
  }

  cutPropagator cutProp;
  cutProp.Clean();
  cutProp.GetAllVarFromFile(responseFile_p);
  Int_t nCentBinsTemp = cutProp.GetNCentBins();
  std::vector<Int_t> centBinsLow = cutProp.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutProp.GetCentBinsHi();

  Float_t jtAbsEtaMaxTemp = cutProp.GetJtAbsEtaMax();

  const Int_t nJtAlgos = cutProp.GetNJtAlgos();
  std::vector<std::string> jtAlgos = cutProp.GetJtAlgos();
  std::vector<double> minJtPtCutTemp = cutProp.GetMinJtPtCut();
  std::vector<double> multiJtPtCutTemp = cutProp.GetMultiJtPtCut();
  std::vector<int> recoTruncPosTemp = cutProp.GetRecoTruncPos();

  Int_t nRecoJtPtBinsSmallRCent0to10Temp = cutProp.GetNRecoJtPtBinsSmallRCent0to10();
  Int_t nRecoJtPtBinsLargeRCent0to10Temp = cutProp.GetNRecoJtPtBinsLargeRCent0to10();
  Int_t nGenJtPtSmallBinsSmallRCent0to10Temp = cutProp.GetNGenJtPtSmallBinsSmallRCent0to10();
  Int_t nGenJtPtLargeBinsSmallRCent0to10Temp = cutProp.GetNGenJtPtLargeBinsSmallRCent0to10();
  Int_t nGenJtPtSmallBinsLargeRCent0to10Temp = cutProp.GetNGenJtPtSmallBinsLargeRCent0to10();
  Int_t nGenJtPtLargeBinsLargeRCent0to10Temp = cutProp.GetNGenJtPtLargeBinsLargeRCent0to10();

  Int_t nRecoJtPtBinsSmallRCent10to30Temp = cutProp.GetNRecoJtPtBinsSmallRCent10to30();
  Int_t nRecoJtPtBinsLargeRCent10to30Temp = cutProp.GetNRecoJtPtBinsLargeRCent10to30();
  Int_t nGenJtPtSmallBinsSmallRCent10to30Temp = cutProp.GetNGenJtPtSmallBinsSmallRCent10to30();
  Int_t nGenJtPtLargeBinsSmallRCent10to30Temp = cutProp.GetNGenJtPtLargeBinsSmallRCent10to30();
  Int_t nGenJtPtSmallBinsLargeRCent10to30Temp = cutProp.GetNGenJtPtSmallBinsLargeRCent10to30();
  Int_t nGenJtPtLargeBinsLargeRCent10to30Temp = cutProp.GetNGenJtPtLargeBinsLargeRCent10to30();

  Int_t nRecoJtPtBinsSmallRCent30to50Temp = cutProp.GetNRecoJtPtBinsSmallRCent30to50();
  Int_t nRecoJtPtBinsLargeRCent30to50Temp = cutProp.GetNRecoJtPtBinsLargeRCent30to50();
  Int_t nGenJtPtSmallBinsSmallRCent30to50Temp = cutProp.GetNGenJtPtSmallBinsSmallRCent30to50();
  Int_t nGenJtPtLargeBinsSmallRCent30to50Temp = cutProp.GetNGenJtPtLargeBinsSmallRCent30to50();
  Int_t nGenJtPtSmallBinsLargeRCent30to50Temp = cutProp.GetNGenJtPtSmallBinsLargeRCent30to50();
  Int_t nGenJtPtLargeBinsLargeRCent30to50Temp = cutProp.GetNGenJtPtLargeBinsLargeRCent30to50();

  Int_t nRecoJtPtBinsSmallRCent50to90Temp = cutProp.GetNRecoJtPtBinsSmallRCent50to90();
  Int_t nRecoJtPtBinsLargeRCent50to90Temp = cutProp.GetNRecoJtPtBinsLargeRCent50to90();
  Int_t nGenJtPtSmallBinsSmallRCent50to90Temp = cutProp.GetNGenJtPtSmallBinsSmallRCent50to90();
  Int_t nGenJtPtLargeBinsSmallRCent50to90Temp = cutProp.GetNGenJtPtLargeBinsSmallRCent50to90();
  Int_t nGenJtPtSmallBinsLargeRCent50to90Temp = cutProp.GetNGenJtPtSmallBinsLargeRCent50to90();
  Int_t nGenJtPtLargeBinsLargeRCent50to90Temp = cutProp.GetNGenJtPtLargeBinsLargeRCent50to90();

  std::vector<Double_t> recoJtPtBinsSmallRTemp = cutProp.GetRecoJtPtBinsSmallR();
  std::vector<Double_t> recoJtPtBinsLargeRTemp = cutProp.GetRecoJtPtBinsLargeR();
  std::vector<Double_t> genJtPtSmallBinsSmallRTemp = cutProp.GetGenJtPtSmallBinsSmallR();
  std::vector<Double_t> genJtPtLargeBinsSmallRTemp = cutProp.GetGenJtPtLargeBinsSmallR();
  std::vector<Double_t> genJtPtSmallBinsLargeRTemp = cutProp.GetGenJtPtSmallBinsLargeR();
  std::vector<Double_t> genJtPtLargeBinsLargeRTemp = cutProp.GetGenJtPtLargeBinsLargeR();

  Int_t nJtAbsEtaBinsTemp = cutProp.GetNJtAbsEtaBins();
  std::vector<Double_t> jtAbsEtaBinsLowTemp = cutProp.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHiTemp = cutProp.GetJtAbsEtaBinsHi();

  bool isResponsePP = cutProp.GetIsPP();

  Int_t nIDTemp = cutProp.GetNID();
  std::vector<std::string> idStr = cutProp.GetIdStr();
  std::vector<double> jtPfCHMFCutLow = cutProp.GetJtPfCHMFCutLow();
  std::vector<double> jtPfCHMFCutHi = cutProp.GetJtPfCHMFCutHi();
  std::vector<double> jtPfMUMFCutLow = cutProp.GetJtPfMUMFCutLow();
  std::vector<double> jtPfMUMFCutHi = cutProp.GetJtPfMUMFCutHi();
  std::vector<double> jtPfNHFCutLow = cutProp.GetJtPfNHFCutLow();
  std::vector<double> jtPfNHFCutHi = cutProp.GetJtPfNHFCutHi();
  std::vector<double> jtPfNEFCutLow = cutProp.GetJtPfNEFCutLow();
  std::vector<double> jtPfNEFCutHi = cutProp.GetJtPfNEFCutHi();
  std::vector<double> jtPfMUFCutLow = cutProp.GetJtPfMUFCutLow();
  std::vector<double> jtPfMUFCutHi = cutProp.GetJtPfMUFCutHi();
  std::vector<double> jtPfCHFCutLow = cutProp.GetJtPfCHFCutLow();
  std::vector<double> jtPfCHFCutHi = cutProp.GetJtPfCHFCutHi();
  std::vector<double> jtPfCEFCutLow = cutProp.GetJtPfCEFCutLow();
  std::vector<double> jtPfCEFCutHi = cutProp.GetJtPfCEFCutHi();
  std::vector<int> jtPfMinMult = cutProp.GetJtPfMinMult();
  std::vector<int> jtPfMinChgMult = cutProp.GetJtPfMinChgMult();

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(nCentBinsTemp < 0) std::cout << "nCentBins less than 0. please check input file. return 1" << std::endl;
  if(nIDTemp < 0) std::cout << "nID less than 0. please check input file. return 1" << std::endl;
  if(nRecoJtPtBinsSmallRCent0to10Temp < 0) std::cout << "nRecoJtPtBinsSmallRCent0to10Temp less than 0. please check input file. return 1" << std::endl;
  if(nRecoJtPtBinsLargeRCent0to10Temp < 0) std::cout << "nRecoJtPtBinsLargeRCent0to10Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtSmallBinsSmallRCent0to10Temp < 0) std::cout << "nGenJtPtSmallBinsSmallRCent0to10Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtLargeBinsSmallRCent0to10Temp < 0) std::cout << "nGenJtPtLargeBinsSmallRCent0to10Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtSmallBinsLargeRCent0to10Temp < 0) std::cout << "nGenJtPtSmallBinsLargeRCent0to10Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtLargeBinsLargeRCent0to10Temp < 0) std::cout << "nGenJtPtLargeBinsLargeRCent0to10Temp less than 0. please check input file. return 1" << std::endl;
  if(nRecoJtPtBinsSmallRCent10to30Temp < 0) std::cout << "nRecoJtPtBinsSmallRCent10to30Temp less than 0. please check input file. return 1" << std::endl;
  if(nRecoJtPtBinsLargeRCent10to30Temp < 0) std::cout << "nRecoJtPtBinsLargeRCent10to30Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtSmallBinsSmallRCent10to30Temp < 0) std::cout << "nGenJtPtSmallBinsSmallRCent10to30Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtLargeBinsSmallRCent10to30Temp < 0) std::cout << "nGenJtPtLargeBinsSmallRCent10to30Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtSmallBinsLargeRCent10to30Temp < 0) std::cout << "nGenJtPtSmallBinsLargeRCent10to30Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtLargeBinsLargeRCent10to30Temp < 0) std::cout << "nGenJtPtLargeBinsLargeRCent10to30Temp less than 0. please check input file. return 1" << std::endl;
  if(nRecoJtPtBinsSmallRCent30to50Temp < 0) std::cout << "nRecoJtPtBinsSmallRCent30to50Temp less than 0. please check input file. return 1" << std::endl;
  if(nRecoJtPtBinsLargeRCent30to50Temp < 0) std::cout << "nRecoJtPtBinsLargeRCent30to50Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtSmallBinsSmallRCent30to50Temp < 0) std::cout << "nGenJtPtSmallBinsSmallRCent30to50Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtLargeBinsSmallRCent30to50Temp < 0) std::cout << "nGenJtPtLargeBinsSmallRCent30to50Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtSmallBinsLargeRCent30to50Temp < 0) std::cout << "nGenJtPtSmallBinsLargeRCent30to50Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtLargeBinsLargeRCent30to50Temp < 0) std::cout << "nGenJtPtLargeBinsLargeRCent30to50Temp less than 0. please check input file. return 1" << std::endl;
  if(nRecoJtPtBinsSmallRCent50to90Temp < 0) std::cout << "nRecoJtPtBinsSmallRCent50to90Temp less than 0. please check input file. return 1" << std::endl;
  if(nRecoJtPtBinsLargeRCent50to90Temp < 0) std::cout << "nRecoJtPtBinsLargeRCent50to90Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtSmallBinsSmallRCent50to90Temp < 0) std::cout << "nGenJtPtSmallBinsSmallRCent50to90Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtLargeBinsSmallRCent50to90Temp < 0) std::cout << "nGenJtPtLargeBinsSmallRCent50to90Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtSmallBinsLargeRCent50to90Temp < 0) std::cout << "nGenJtPtSmallBinsLargeRCent50to90Temp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtLargeBinsLargeRCent50to90Temp < 0) std::cout << "nGenJtPtLargeBinsLargeRCent50to90Temp less than 0. please check input file. return 1" << std::endl;
  if(nJtAbsEtaBinsTemp < 0) std::cout << "nJtAbsEtaBinsTemp less than 0. please check input file. return 1" << std::endl;
  if(jtAbsEtaMaxTemp < -98) std::cout << "jtAbsEtaMaxTemp less than -98. please check input file. return 1" << std::endl;

  if(nCentBinsTemp < 0 || nRecoJtPtBinsSmallRCent0to10Temp < 0 || nRecoJtPtBinsLargeRCent0to10Temp < 0 || nRecoJtPtBinsSmallRCent10to30Temp < 0 || nRecoJtPtBinsLargeRCent10to30Temp < 0 || nRecoJtPtBinsSmallRCent30to50Temp < 0 || nRecoJtPtBinsLargeRCent30to50Temp < 0 || nRecoJtPtBinsSmallRCent50to90Temp < 0 || nRecoJtPtBinsLargeRCent50to90Temp < 0 || nJtAbsEtaBinsTemp < 0 || jtAbsEtaMaxTemp < -98 || nIDTemp < 0){
    responseFile_p->Close();
    delete responseFile_p;
    return 1;
  }

  smallOrLargeR rReader;
  if(!rReader.CheckNRecoJtPtBinsSmallRCent0to10(nRecoJtPtBinsSmallRCent0to10Temp)) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent0to10(nRecoJtPtBinsLargeRCent0to10Temp)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent0to10(nGenJtPtSmallBinsSmallRCent0to10Temp)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent0to10(nGenJtPtLargeBinsSmallRCent0to10Temp)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent0to10(nGenJtPtSmallBinsLargeRCent0to10Temp)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent0to10(nGenJtPtLargeBinsLargeRCent0to10Temp)) return 1;

  if(!rReader.CheckNRecoJtPtBinsSmallRCent10to30(nRecoJtPtBinsSmallRCent10to30Temp)) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent10to30(nRecoJtPtBinsLargeRCent10to30Temp)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent10to30(nGenJtPtSmallBinsSmallRCent10to30Temp)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent10to30(nGenJtPtLargeBinsSmallRCent10to30Temp)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent10to30(nGenJtPtSmallBinsLargeRCent10to30Temp)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent10to30(nGenJtPtLargeBinsLargeRCent10to30Temp)) return 1;

  if(!rReader.CheckNRecoJtPtBinsSmallRCent30to50(nRecoJtPtBinsSmallRCent30to50Temp)) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent30to50(nRecoJtPtBinsLargeRCent30to50Temp)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent30to50(nGenJtPtSmallBinsSmallRCent30to50Temp)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent30to50(nGenJtPtLargeBinsSmallRCent30to50Temp)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent30to50(nGenJtPtSmallBinsLargeRCent30to50Temp)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent30to50(nGenJtPtLargeBinsLargeRCent30to50Temp)) return 1;

  if(!rReader.CheckNRecoJtPtBinsSmallRCent50to90(nRecoJtPtBinsSmallRCent50to90Temp)) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent50to90(nRecoJtPtBinsLargeRCent50to90Temp)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent50to90(nGenJtPtSmallBinsSmallRCent50to90Temp)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent50to90(nGenJtPtLargeBinsSmallRCent50to90Temp)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent50to90(nGenJtPtSmallBinsLargeRCent50to90Temp)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent50to90(nGenJtPtLargeBinsLargeRCent50to90Temp)) return 1;

  if(!rReader.CheckRecoJtPtBinsSmallR(recoJtPtBinsSmallRTemp)) return 1;
  if(!rReader.CheckRecoJtPtBinsLargeR(recoJtPtBinsLargeRTemp)) return 1;
  if(!rReader.CheckGenJtPtSmallBinsSmallR(genJtPtSmallBinsSmallRTemp)) return 1;
  if(!rReader.CheckGenJtPtLargeBinsSmallR(genJtPtLargeBinsSmallRTemp)) return 1;
  if(!rReader.CheckGenJtPtSmallBinsLargeR(genJtPtSmallBinsLargeRTemp)) return 1;
  if(!rReader.CheckGenJtPtLargeBinsLargeR(genJtPtLargeBinsLargeRTemp)) return 1;

  const Float_t jtAbsEtaMax = jtAbsEtaMaxTemp;
  const Int_t nMaxCentBins = 4;
  const Int_t nCentBins = nCentBinsTemp;

  if(nCentBins > nMaxCentBins){
    std::cout << "nCentBins \'" << nCentBins << "\' is greater than nMaxCentBins \'" << nMaxCentBins << "\'. return 1" << std::endl;
    return 1;
  }

  const Int_t nSmallLargeBins = 2;
  const std::string smallLargeBinsStr[nSmallLargeBins] = {"SmallBins", "LargeBins"};
  const Int_t nMaxJtPtBins = 50;
  const Int_t nRecoJtPtBinsArr = TMath::Max(TMath::Max(TMath::Max(nRecoJtPtBinsSmallRCent0to10Temp, nRecoJtPtBinsLargeRCent0to10Temp), TMath::Max(nRecoJtPtBinsSmallRCent10to30Temp, nRecoJtPtBinsLargeRCent10to30Temp)), TMath::Max(TMath::Max(nRecoJtPtBinsSmallRCent30to50Temp, nRecoJtPtBinsLargeRCent30to50Temp), TMath::Max(nRecoJtPtBinsSmallRCent50to90Temp, nRecoJtPtBinsLargeRCent50to90Temp)));
  
  if(nMaxJtPtBins < nRecoJtPtBinsArr){
    std::cout << "nRecoJtPtBinsArr \'" << nRecoJtPtBinsArr << "\' is less than nMaxJtPtBins \'" << nMaxJtPtBins << "\'. return 1" << std::endl;
    return 1;
  }

  const Int_t nMaxJtAbsEtaBins = 6;
  const Int_t nJtAbsEtaBins = nJtAbsEtaBinsTemp;

  if(nMaxJtAbsEtaBins < nJtAbsEtaBins){
    std::cout << "nJtAbsEtaBins \'" << nJtAbsEtaBins << "\' is less than nMaxJtAbsEtaBins \'" << nMaxJtAbsEtaBins << "\'. return 1" << std::endl;
    return 1;
  }

  Double_t jtAbsEtaBinsLow[nMaxJtAbsEtaBins];
  Double_t jtAbsEtaBinsHi[nMaxJtAbsEtaBins];
  std::cout << "nJtAbsEtaBins: ";
  for(Int_t jI = 0; jI < nJtAbsEtaBins; ++jI){
    jtAbsEtaBinsLow[jI] = jtAbsEtaBinsLowTemp[jI];
    jtAbsEtaBinsHi[jI] = jtAbsEtaBinsHiTemp[jI];
    std::cout << " " << jtAbsEtaBinsLow[jI] << "-" << jtAbsEtaBinsHi[jI] << ",";
  }
  std::cout << std::endl;

  const Int_t nMaxID = 6;
  const Int_t nID = nIDTemp; 

  if(nMaxID < nID){
    std::cout << "nID \'" << nID << "\' is less than nMaxID\'" << nMaxID << "\'. return 1" << std::endl;
    return 1;
  }


  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  std::cout << "nCentBins: " << nCentBins << std::endl;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::cout << " " << cI << "/" << nCentBins << ": " << centBinsLow[cI] << "-" << centBinsHi[cI] << std::endl;
  }

  std::cout << "Raw data from file: \'" << inDataFileName << "\'" << std::endl;

  std::vector<std::string> fileList;

  if(inDataFileName.find(".root") != std::string::npos){
    fileList.push_back(inDataFileName);
  }
  else if(inDataFileName.find(".txt") != std::string::npos){
    std::ifstream file(inDataFileName.c_str());
    std::string tempStr;

    while(std::getline(file, tempStr)){
      while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
      else{
        if(tempStr.substr(0, std::string("ISPP=").size()).find("ISPP=") != std::string::npos){
          tempStr.replace(0,tempStr.find("=")+1, "");
          isDataPP = std::stoi(tempStr);
        }
        else std::cout << "WARNING: Line in \'" << inDataFileName << "\', \'" << tempStr << "\' is invalid. check input" << std::endl;
      }
    }

    file.close();
  }
  else{
    std::cout << "Given inDataFileName \'" << inDataFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  if(fileList.size() == 0){
    std::cout << "Given inDataFileName \'" << inDataFileName << "\' is gives no valid root files. return 1" << std::endl;
    return 1;
  }
  else if(isDataPP != isResponsePP){
    std::cout << "Response isDataPP == " << isResponsePP << ", data isDataPP == " << isDataPP << " aren't equivalent. return 1" << std::endl;
    return 1;
  }

  std::cout << "Checking for matched inputs..." << std::endl;
  TFile* inDataFile_p = TFile::Open(mntToXRootdFileString(fileList[0]).c_str(), "READ");
  std::vector<std::string> dataTreeList = returnRootFileContentsList(inDataFile_p, "TTree", "JetAna");

  unsigned int pos = 0;
  while(dataTreeList.size() > pos){   
    bool isFound = false;
    std::string tempJet = dataTreeList[pos];
    tempJet.replace(tempJet.find("/"), tempJet.size() - tempJet.find("/"), "");

    for(unsigned int jI = 0; jI < jetDirList.size(); ++jI){
      if(jetDirList[jI].find(tempJet) != std::string::npos && tempJet.size() == jetDirList[jI].size()){
	isFound = true;
	break;
      }
    }

    if(isFound) ++pos;
    else dataTreeList.erase(dataTreeList.begin() + pos);
  }

  inDataFile_p->Close();
  delete inDataFile_p;
  inDataFile_p = NULL;

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nMaxDataJet = 10;
  const Int_t nDataJet = dataTreeList.size();

  if(nMaxDataJet < nDataJet){
    std::cout << "nDataJet \'" << nDataJet << "\' is less than nMaxDataJet\'" << nMaxDataJet << "\'. return 1" << std::endl;
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

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  
  const Int_t sizeToTruncName = 40;
  while(outFileName.size() > sizeToTruncName){outFileName = outFileName.substr(0,outFileName.size()-1);}
  while(outFileName2.size() > sizeToTruncName){outFileName2 = outFileName2.substr(0,outFileName2.size()-1);}
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + outFileName2 + "_ProcessRawData_";
  if(tagStr.size() != 0) outFileName = outFileName + tagStr + "_";
  outFileName = outFileName + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);
  TDirectory* dir_p[nMaxDataJet];
  TH1D* jtPtRaw_RecoGenSymm_h[nMaxDataJet][nMaxCentBins][nMaxID][nMaxJtAbsEtaBins][nSmallLargeBins];
  TH1D* jtPtRaw_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nMaxID][nMaxJtAbsEtaBins];

  TH1D* multijetAJ_All_h[nMaxDataJet][nMaxCentBins][nMaxID][nMaxJtAbsEtaBins][nMaxJtPtBins];
  TH1D* multijetAJ_Pass_h[nMaxDataJet][nMaxCentBins][nMaxID][nMaxJtAbsEtaBins][nMaxJtPtBins];
  TH1D* multijetAJ_Fail_h[nMaxDataJet][nMaxCentBins][nMaxID][nMaxJtAbsEtaBins][nMaxJtPtBins];

  Double_t minJtPtCut[nMaxDataJet];
  Double_t multiJtPtCut[nMaxDataJet];
  Int_t recoTruncPos[nMaxDataJet];
  for(Int_t jI = 0; jI < nDataJet; ++jI){
    int pos = -1;
    std::string treeStr1 = dataTreeList[jI];
    while(treeStr1.find("/") != std::string::npos){treeStr1.replace(treeStr1.find("/"), treeStr1.size(), "");}

    for(Int_t jI2 = 0; jI2 < nJtAlgos; ++jI2){
      std::string treeStr2 = jtAlgos[jI2];
      while(treeStr2.find("/") != std::string::npos){treeStr2.replace(treeStr2.find("/"), treeStr2.size(), "");}

      if(isStrSame(treeStr1, treeStr2)){
	pos = jI2;
	break;
      }
    }

    multiJtPtCut[jI] = multiJtPtCutTemp[pos];
    minJtPtCut[jI] = minJtPtCutTemp[pos];
    recoTruncPos[jI] = recoTruncPosTemp[pos];
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  Double_t minRecoJtPt = 99999999.;
  Double_t minGenJtPt = 99999999.;

  Int_t rValI[nDataJet];
  bool isSmallR[nDataJet];
  Int_t nRecoJtPtBins[nMaxDataJet][nMaxCentBins];
  Int_t nGenJtPtSmallBins[nMaxDataJet][nMaxCentBins];
  Int_t nGenJtPtLargeBins[nMaxDataJet][nMaxCentBins];
  Int_t nGenJtPtBins[nMaxDataJet][nMaxCentBins][nSmallLargeBins];

  Double_t recoJtPtBins[nMaxDataJet][nMaxCentBins][nMaxJtPtBins+1];
  Double_t genJtPtSmallBins[nMaxDataJet][nMaxCentBins][nMaxJtPtBins+1];
  Double_t genJtPtLargeBins[nMaxDataJet][nMaxCentBins][nMaxJtPtBins+1];
  Double_t genJtPtBins[nMaxDataJet][nMaxCentBins][nSmallLargeBins][nMaxJtPtBins+1];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = dataTreeList[jI];
    tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    rValI[jI] = getRVal(tempStr);
    isSmallR[jI] = rReader.GetIsSmallR(rValI[jI]);

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

      nRecoJtPtBins[jI][cI] = rReader.GetSmallOrLargeRNBins(isSmallR[jI], false, false, centStr);
      nGenJtPtSmallBins[jI][cI] = rReader.GetSmallOrLargeRNBins(isSmallR[jI], true, true, centStr);
      nGenJtPtLargeBins[jI][cI] = rReader.GetSmallOrLargeRNBins(isSmallR[jI], true, false, centStr);

      rReader.GetSmallOrLargeRBins(isSmallR[jI], false, nRecoJtPtBins[jI][cI]+1, recoJtPtBins[jI][cI], false);
      rReader.GetSmallOrLargeRBins(isSmallR[jI], true, nGenJtPtSmallBins[jI][cI]+1, genJtPtSmallBins[jI][cI], true);
      rReader.GetSmallOrLargeRBins(isSmallR[jI], true, nGenJtPtLargeBins[jI][cI]+1, genJtPtLargeBins[jI][cI], false);

      for(Int_t bIX = 0; bIX < nRecoJtPtBins[jI][cI]; ++bIX){
	if(minRecoJtPt > recoJtPtBins[jI][cI][bIX]) minRecoJtPt = recoJtPtBins[jI][cI][bIX];
      }

      nGenJtPtBins[jI][cI][0] = nGenJtPtSmallBins[jI][cI];
      nGenJtPtBins[jI][cI][1] = nGenJtPtLargeBins[jI][cI];

      for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
        for(Int_t bIX = 0; bIX < nGenJtPtBins[jI][cI][binsI]+1; ++bIX){
          if(binsI == 0) genJtPtBins[jI][cI][binsI][bIX] = genJtPtSmallBins[jI][cI][bIX];
          else if(binsI == 1) genJtPtBins[jI][cI][binsI][bIX] = genJtPtLargeBins[jI][cI][bIX];

	  if(minGenJtPt > genJtPtBins[jI][cI][binsI][bIX]) minGenJtPt = genJtPtBins[jI][cI][binsI][bIX];
        }
      }

    }
  }

  const Int_t nMultijetAJ = 9;
  const Double_t multijetAJ[nMultijetAJ+1] = {-0.5, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.1};
  //  Double_t multijetAJFactors[nMultijetAJ];
  //  for(Int_t mI = 0; mI < nMultijetAJ; ++mI){
  //    multijetAJFactors[mI] = (multijetAJ[mI+1] - multijetAJ[mI])/0.10;
  //  }
  
  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = dataTreeList[jI];
    tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dir_p[jI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());
    std::cout << " " << jI << "/" << nDataJet << ": " << tempStr << std::endl;
    
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);      
      if(isDataPP) centStr = "PP_" + centStr;
      else centStr = "PbPb_" + centStr;
      
      for(Int_t idI = 0; idI < nID; ++idI){
   	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
   	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	  
   	  for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
   	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI] = new TH1D(("jtPtRaw_RecoGenSymm_" + smallLargeBinsStr[binsI] + "_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + jtAbsEtaStr + "_h").c_str(), ";Raw Jet p_{T};Counts", nGenJtPtBins[jI][cI][binsI], genJtPtBins[jI][cI][binsI]);
   	    std::vector<TH1*> tempVect = {jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]};
   	    setSumW2(tempVect);
   	    centerTitles(tempVect);
   	  }
	  
   	  jtPtRaw_RecoGenAsymm_h[jI][cI][idI][aI] = new TH1D(("jtPtRaw_RecoGenAsymm_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + jtAbsEtaStr + "_h").c_str(), ";Raw Jet p_{T};Counts", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI]);
	  
	  for(Int_t jIter = 0; jIter < nRecoJtPtBins[jI][cI]; ++jIter){
	    if(jIter < nRecoJtPtBins[jI][cI]){
	      const std::string jtPtStr = "JtPt" + prettyString(recoJtPtBins[jI][cI][jIter], 1, true) + "to" + prettyString(recoJtPtBins[jI][cI][jIter+1], 1, true);
	      	      
	      multijetAJ_All_h[jI][cI][idI][aI][jIter] = new TH1D(("multijetAJ_All_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + jtAbsEtaStr + "_" + jtPtStr + "_h").c_str(), ";Multijet A_{J};Counts", nMultijetAJ, multijetAJ);
	      multijetAJ_Pass_h[jI][cI][idI][aI][jIter] = new TH1D(("multijetAJ_Pass_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + jtAbsEtaStr + "_" + jtPtStr + "_h").c_str(), ";Multijet A_{J};Counts", nMultijetAJ, multijetAJ);
	      multijetAJ_Fail_h[jI][cI][idI][aI][jIter] = new TH1D(("multijetAJ_Fail_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + jtAbsEtaStr + "_" + jtPtStr + "_h").c_str(), ";Multijet A_{J};Counts", nMultijetAJ, multijetAJ);
	      
	      std::vector<TH1*> tempVect = {multijetAJ_All_h[jI][cI][idI][aI][jIter], multijetAJ_Pass_h[jI][cI][idI][aI][jIter], multijetAJ_Fail_h[jI][cI][idI][aI][jIter]};
	      setSumW2(tempVect);
	      centerTitles(tempVect);
	    }
	    else{
	      multijetAJ_All_h[jI][cI][idI][aI][jIter] = NULL;
	      multijetAJ_Pass_h[jI][cI][idI][aI][jIter] = NULL;
	      multijetAJ_Fail_h[jI][cI][idI][aI][jIter] = NULL;
	    }
	  }
	  
	  std::vector<TH1*> tempVect = {jtPtRaw_RecoGenAsymm_h[jI][cI][idI][aI]};
	  setSumW2(tempVect);
	  centerTitles(tempVect);
	}
      }
    }
  }
  
  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::cout << "Processing " << fileList.size() << " files..., isDataPP == " << isDataPP << std::endl;
  std::cout << "Jets: ";
  for(unsigned int jI = 0; jI < dataTreeList.size(); ++jI){
    std::cout << dataTreeList[jI] << ", ";
  }
  std::cout << std::endl;
  
  
  //Jet Var  
  const Int_t nMaxJet = 500;
  Int_t nref_[nMaxDataJet];
  Float_t jtpt_[nMaxDataJet][nMaxJet];
  Float_t jteta_[nMaxDataJet][nMaxJet];
  Float_t jtphi_[nMaxDataJet][nMaxJet];
  Float_t jtPfCHMF_[nMaxDataJet][nMaxJet];
  Float_t jtPfMUMF_[nMaxDataJet][nMaxJet];
  Float_t jtPfNHF_[nMaxDataJet][nMaxJet];
  Float_t jtPfNEF_[nMaxDataJet][nMaxJet];
  Float_t jtPfMUF_[nMaxDataJet][nMaxJet];
  Float_t jtPfCHF_[nMaxDataJet][nMaxJet];
  Float_t jtPfCEF_[nMaxDataJet][nMaxJet];
  Int_t jtPfNHM_[nMaxDataJet][nMaxJet];
  Int_t jtPfNEM_[nMaxDataJet][nMaxJet];
  Int_t jtPfMUM_[nMaxDataJet][nMaxJet];
  Int_t jtPfCHM_[nMaxDataJet][nMaxJet];
  Int_t jtPfCEM_[nMaxDataJet][nMaxJet];
  
  Float_t vz_;
  Float_t hiHF_;
  Int_t hiBin_;
  unsigned int run_, lumi_;
  unsigned long long evt_;
  
  Int_t HBHENoiseFilterResultRun2Loose_ = -1;
  Int_t pprimaryVertexFilter_ = -1;
  Int_t pBeamScrapingFilter_ = -1;
  Int_t phfCoincFilter3_ = -1;
  Int_t pclusterCompatibilityFilter_ = -1;
  
  goodGlobalSelection globalSel;
  globalSel.setIsPbPb(!isDataPP);
  
  std::vector<Int_t> centPos;
  if(isDataPP){for(Int_t cI = 0; cI < nCentBins; ++cI){centPos.push_back(cI);}}
  else centPos.push_back(-1);

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "File " << fI << "/" << fileList.size() << ": " << fileList[fI] << std::endl;
    
    TFile* inDataFile_p = TFile::Open(mntToXRootdFileString(fileList[fI]).c_str(), "READ");
    TTree* jetTrees_p[nMaxDataJet];
    TTree* hiTree_p = (TTree*)inDataFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* skimTree_p = (TTree*)inDataFile_p->Get("skimanalysis/HltTree");
    
    for(Int_t jI = 0; jI < nDataJet; ++jI){
      jetTrees_p[jI] = NULL;
      jetTrees_p[jI] = (TTree*)inDataFile_p->Get(dataTreeList[jI].c_str());
      
      jetTrees_p[jI]->SetBranchStatus("*", 0);
      jetTrees_p[jI]->SetBranchStatus("nref", 1);
      jetTrees_p[jI]->SetBranchStatus("jtpt", 1);
      jetTrees_p[jI]->SetBranchStatus("jtphi", 1);
      jetTrees_p[jI]->SetBranchStatus("jteta", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfCHMF", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfMUMF", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfNHF", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfNEF", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfMUF", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfCHF", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfCEF", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfNHM", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfNEM", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfMUM", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfCHM", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfCEM", 1);

      jetTrees_p[jI]->SetBranchAddress("nref", &(nref_[jI]));
      jetTrees_p[jI]->SetBranchAddress("jtpt", jtpt_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtphi", jtphi_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jteta", jteta_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfCHMF", jtPfCHMF_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfMUMF", jtPfMUMF_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfNHF", jtPfNHF_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfNEF", jtPfNEF_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfMUF", jtPfMUF_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfCHF", jtPfCHF_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfCEF", jtPfCEF_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfNHM", jtPfNHM_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfNEM", jtPfNEM_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfMUM", jtPfMUM_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfCHM", jtPfCHM_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfCEM", jtPfCEM_[jI]);
    }
    
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("hiHF", 1);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("hiHF", &hiHF_);
    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("evt", &evt_);    
    
    skimTree_p->SetBranchStatus("*", 0);
    skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
    
    if(!isDataPP){
      skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
      skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
      skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);
      
      skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
      skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
      skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);
    }
    else{
      skimTree_p->SetBranchStatus("pBeamScrapingFilter", 1);
      skimTree_p->SetBranchStatus("pPAprimaryVertexFilter", 1);
      
      skimTree_p->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter_);
      skimTree_p->SetBranchAddress("pPAprimaryVertexFilter", &pprimaryVertexFilter_);
    }
    
    const Int_t nEntries = TMath::Min(100000000, (Int_t)hiTree_p->GetEntries());
    const Int_t printInterval = TMath::Max(1, nEntries/20);
  
    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(nEntries >= 50000 && entry%printInterval == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;
      
      //      std::cout << __LINE__ << std::endl;
      hiTree_p->GetEntry(entry);
      skimTree_p->GetEntry(entry);

      globalSel.setVz(vz_);
      globalSel.setHiHF(hiHF_);
      globalSel.setPprimaryVertexFilter(pprimaryVertexFilter_);
      globalSel.setPBeamScrapingFilter(pBeamScrapingFilter_);
      globalSel.setPhfCoincFilter3(phfCoincFilter3_);
      globalSel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
      globalSel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);

      //      std::cout << __LINE__ << std::endl;

      if(!globalSel.isGood()) continue;
    
      if(!isDataPP){
	centPos[0] = -1;
	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  if(hiBin_/2 >= centBinsLow[cI] && hiBin_/2 < centBinsHi[cI]){
	    centPos[0] = cI;
	    break;
	  }
	}
	if(centPos[0] < 0) continue;
      }

      //      std::cout << __LINE__ << std::endl;

      for(Int_t tI = 0; tI < nDataJet; ++tI){jetTrees_p[tI]->GetEntry(entry);}

      for(Int_t tI = 0; tI < nDataJet; ++tI){
	Double_t tempLeadingPt_ = -999;
	Double_t tempLeadingPhi_ = -999;
	std::vector<Double_t> tempSubleadingPt_;
	std::vector<Double_t> tempSubleadingPhi_;
	Int_t tempLeadingPos_ = -1;

	//	std::cout << __LINE__ << std::endl;

	for(Int_t jI = 0; jI < nref_[tI]; ++jI){
	  if(TMath::Abs(jteta_[tI][jI]) > jtAbsEtaMax) continue;
	  
	  if(jtpt_[tI][jI] > tempLeadingPt_){
	    tempLeadingPt_ = jtpt_[tI][jI];
	    tempLeadingPhi_ = jtphi_[tI][jI];
	    tempLeadingPos_ = jI;
	  }
	}

	for(auto const & cent : centPos){
	  Int_t tempLeadingFillPos_ = -1;
	
	  //	  std::cout << __LINE__ << std::endl;
	  	       
 	  if(/*tempLeadingPt_ > jtPtBins[0] && */tempLeadingPt_ > minJtPtCut[tI]){
	    for(Int_t jI = 0; jI < nRecoJtPtBins[tI][cent]; ++jI){
	      if(recoJtPtBins[tI][cent][jI] <= tempLeadingPt_ && tempLeadingPt_ < recoJtPtBins[tI][cent][jI+1]){
		tempLeadingFillPos_ = jI;
		break;
	      }
	    }
	    if(tempLeadingFillPos_ < 0) tempLeadingFillPos_ = nRecoJtPtBins[tI][cent]-1;
	    
	    for(Int_t jI = 0; jI < nref_[tI]; ++jI){
	      if(jI == tempLeadingPos_) continue;
	      if(TMath::Abs(jteta_[tI][jI]) > jtAbsEtaMax) continue;
	      if(jtpt_[tI][jI] < multiJtPtCut[jI]) continue;
	      if(TMath::Abs(getDPHI(tempLeadingPhi_, jtphi_[tI][jI])) < 3.*TMath::Pi()/4.) continue;
	      //	    if(TMath::Abs(getDPHI(tempLeadingPhi_, jtphi_[tI][jI])) < 5.*TMath::Pi()/8. && TMath::Abs(getDPHI(tempLeadingPhi_, jtphi_[tI][jI])) > 3.*TMath::Pi()/8.) continue;
	      
	      tempSubleadingPt_.push_back(jtpt_[tI][jI]);
	      tempSubleadingPhi_.push_back(jtphi_[tI][jI]);
	    }
	    
	    if(tempSubleadingPt_.size() == 0){
	      tempSubleadingPt_.push_back(0.);
	      if(tempLeadingPhi_ > 0) tempSubleadingPhi_.push_back(tempLeadingPhi_ - TMath::Pi());
	      else tempSubleadingPhi_.push_back(tempLeadingPhi_ + TMath::Pi());
	      
	      /*
		if(tempSubleadingPhi_.at(tempSubleadingPhi_.size()-1) > 0) tempSubleadingPhi_.at(tempSubleadingPhi_.size()-1) += TMath::Pi()*4.9/6.;
		else tempSubleadingPhi_.at(tempSubleadingPhi_.size()-1) += TMath::Pi()*4.9/6.;
	      */
	    }
	  }
	  else tempLeadingPos_ = -1;

	  //	  std::cout << __LINE__ << std::endl;
	  
	  Int_t leadID[nMaxID];
	  for(Int_t idI = 0; idI < nID; ++idI){leadID[idI] = false;}

	  for(Int_t jI = 0; jI < nref_[tI]; ++jI){
	    if(TMath::Abs(jteta_[tI][jI]) > jtAbsEtaMax) continue;
	    if(jtpt_[tI][jI] < minRecoJtPt && jtpt_[tI][jI] < minGenJtPt) continue;

	    std::vector<int> jtAbsEtaPoses;
	    for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	      if(TMath::Abs(jteta_[tI][jI]) >= jtAbsEtaBinsLow[aI] && TMath::Abs(jteta_[tI][jI]) < jtAbsEtaBinsHi[aI]){
		jtAbsEtaPoses.push_back(aI);
	      }
	    }
	    
	    std::vector<int> idPoses;
	    //	  bool passesLast = false;
	    
	    for(Int_t idI = 0; idI < nID; ++idI){
	      if(jtPfCHMFCutLow[idI] > jtPfCHMF_[tI][jI]) continue;
	      if(jtPfCHMFCutHi[idI] < jtPfCHMF_[tI][jI]) continue;
	      if(jtPfMUMFCutLow[idI] > jtPfMUMF_[tI][jI]) continue;
	      if(jtPfMUMFCutHi[idI] < jtPfMUMF_[tI][jI]) continue;
	      if(jtPfNHF_[tI][jI] < jtPfNHFCutLow[idI]) continue;
	      if(jtPfNHF_[tI][jI] > jtPfNHFCutHi[idI]) continue;
	      if(jtPfNEF_[tI][jI] < jtPfNEFCutLow[idI]) continue;
	      if(jtPfNEF_[tI][jI] > jtPfNEFCutHi[idI]) continue;
	      if(jtPfMUF_[tI][jI] < jtPfMUFCutLow[idI]) continue;
	      if(jtPfMUF_[tI][jI] > jtPfMUFCutHi[idI]) continue;
	      if(jtPfCHF_[tI][jI] < jtPfCHFCutLow[idI]) continue;
	      if(jtPfCHF_[tI][jI] > jtPfCHFCutHi[idI]) continue;
	      if(jtPfCEF_[tI][jI] < jtPfCEFCutLow[idI]) continue;
	      if(jtPfCEF_[tI][jI] > jtPfCEFCutHi[idI]) continue;
	      
	      //std::cout << jtPfMinMult.size() << ", " << jtPfMinChgMult.size() << std::endl;

	      if(jtPfCEM_[tI][jI] + jtPfNEM_[tI][jI] + jtPfCHM_[tI][jI] + jtPfNHM_[tI][jI] + jtPfMUM_[tI][jI] < jtPfMinMult[idI]) continue;
	      if(jtPfCHM_[tI][jI] < jtPfMinChgMult[idI]) continue;
	      
	      //	    if(idI == nID-1) passesLast = true;

	      if(tempLeadingPos_ == jI) leadID[idI] = true;
	      
	      idPoses.push_back(idI);
	    }

	    //	    std::cout << __LINE__ << std::endl;

	    //	  if(dataTreeList[tI].find("ak4PF") != std::string::npos && !passesLast && jtpt_[tI][jI] > 600.){
	    //std::cout << "Jet fail: " << entry << ", " << jtpt_[tI][jI] << std::endl;
	    //	  }
	  
	    for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI){
	      for(unsigned int idI = 0; idI < idPoses.size(); ++idI){
		bool goodReco[nSmallLargeBins] = {(jtpt_[tI][jI] >= genJtPtBins[tI][cent][0][0] && jtpt_[tI][jI] < genJtPtBins[tI][cent][0][nGenJtPtBins[tI][cent][0]]), (jtpt_[tI][jI] >= genJtPtBins[tI][cent][1][0] && jtpt_[tI][jI] < genJtPtBins[tI][cent][1][nGenJtPtBins[tI][cent][1]])};
		bool goodRecoTrunc = (jtpt_[tI][jI] >= recoJtPtBins[tI][cent][0] && jtpt_[tI][jI] < recoJtPtBins[tI][cent][nRecoJtPtBins[tI][cent]]);//since we do rDep binning keep same as goodreco for now
	      
		if(goodReco[0]) jtPtRaw_RecoGenSymm_h[tI][cent][idPoses[idI]][jtAbsEtaPoses[aI]][0]->Fill(jtpt_[tI][jI]);
		if(goodReco[1]) jtPtRaw_RecoGenSymm_h[tI][cent][idPoses[idI]][jtAbsEtaPoses[aI]][1]->Fill(jtpt_[tI][jI]);
		if(goodRecoTrunc) jtPtRaw_RecoGenAsymm_h[tI][cent][idPoses[idI]][jtAbsEtaPoses[aI]]->Fill(jtpt_[tI][jI]);
	      }

	      if(tempLeadingPos_ == jI){
		double axis = tempLeadingPhi_;
		if(axis > 0) axis -= TMath::Pi();
		else axis += TMath::Pi();
		
		for(unsigned int pI = 0; pI < tempSubleadingPt_.size()-1; ++pI){
		  for(unsigned int pI2 = pI+1; pI2 < tempSubleadingPt_.size(); ++pI2){
		    
		    if(tempSubleadingPt_[pI]*TMath::Cos(TMath::Abs(getDPHI(axis, tempSubleadingPhi_[pI]))) < tempSubleadingPt_[pI2]*TMath::Cos(TMath::Abs(getDPHI(axis, tempSubleadingPhi_[pI2]))) ){
		      Float_t tempPt = tempSubleadingPt_[pI];
		      Float_t tempPhi = tempSubleadingPhi_[pI];
		      
		      tempSubleadingPt_[pI] = tempSubleadingPt_[pI2];
		      tempSubleadingPhi_[pI] = tempSubleadingPhi_[pI2];
		      
		      tempSubleadingPt_[pI2] = tempPt;
		      tempSubleadingPhi_[pI2] = tempPhi;
		    }
		  }
		}

		//		std::cout << __LINE__ << std::endl;
		
		double projSub = 0;
		//		std::cout << __LINE__ << std::endl;
		for(unsigned int pI = 0; pI < TMath::Min((unsigned int)2, (unsigned int)tempSubleadingPt_.size()); ++pI){
		  //		std::cout << __LINE__ << std::endl;
		  projSub += tempSubleadingPt_[pI]*TMath::Cos(TMath::Abs(getDPHI(axis, tempSubleadingPhi_[pI])));
		  //		std::cout << __LINE__ << std::endl;
		}

		//		std::cout << __LINE__ << std::endl;

		double aj = (tempLeadingPt_ - projSub)/(tempLeadingPt_ + projSub);

		//		std::cout << __LINE__ << std::endl;
		
		for(Int_t idI = 0; idI < nID; ++idI){
		  
		  //		  std::cout << __LINE__ << std::endl;

		  if(aj < multijetAJ[0]) aj = (multijetAJ[0] + multijetAJ[1])/2.;
		  else if(aj > multijetAJ[nMultijetAJ]) aj = (multijetAJ[nMultijetAJ - 1] + multijetAJ[nMultijetAJ])/2.;
	
		  //		std::cout << __LINE__ << std::endl;

		  //		std::cout << nRecoJtPtBins[tI][cent] << std::endl;

		  //		std::cout << tI << ", " << cent << ", " << idI << ", " << jtAbsEtaPoses[aI] << ", " << tempLeadingFillPos_ << std::endl;

		  multijetAJ_All_h[tI][cent][idI][jtAbsEtaPoses[aI]][tempLeadingFillPos_]->Fill(aj);
		  //		std::cout << __LINE__ << std::endl;

		  if(leadID[idI]) multijetAJ_Pass_h[tI][cent][idI][jtAbsEtaPoses[aI]][tempLeadingFillPos_]->Fill(aj);
		  else multijetAJ_Fail_h[tI][cent][idI][jtAbsEtaPoses[aI]][tempLeadingFillPos_]->Fill(aj);
		  //		std::cout << __LINE__ << std::endl;
		}
		//std::cout << __LINE__ << std::endl;
	      }
		//std::cout << __LINE__ << std::endl;
	    }
		//std::cout << __LINE__ << std::endl;
	  }      
		//std::cout << __LINE__ << std::endl;
	}
		//std::cout << __LINE__ << std::endl;
      }
		//std::cout << __LINE__ << std::endl;
    }

    //std::cout << __LINE__ << std::endl;

    inDataFile_p->Close();
    delete inDataFile_p;
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  outFile_p->cd();


  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<std::string> spectraNames;
  std::vector<std::vector<std::string> > multijetNames;

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = dataTreeList[jI];
    tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    for(Int_t cI = 0; cI < nCentBins; ++cI){
	std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
	std::string centStr2 =  std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";

      if(isDataPP){
 	centStr = "PP_" + centStr;
	centStr2 = "PP " + centStr2;
      }
      else{
 	centStr = "PbPb_" + centStr;
	centStr2 = "PbPb " + centStr2;
      }

      for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	  const std::string jtAbsEtaStr2 = prettyString(jtAbsEtaBinsLow[aI], 1, true) + "<|#eta|<" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	  
	  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
	  canv_p->SetTopMargin(0.01);
	  canv_p->SetRightMargin(0.002);
	  canv_p->SetLeftMargin(0.01);
	  canv_p->SetBottomMargin(0.01);
	  
	  TPad* pads[3];
	  canv_p->cd();
	  pads[0] = new TPad("pad0", "", 0.0, yPadFrac, 1.0, 1.0);
	  pads[0]->SetLeftMargin(marg);
	  pads[0]->SetTopMargin(0.01);
	  pads[0]->SetBottomMargin(0.001);
	  pads[0]->SetRightMargin(0.002);
	  pads[0]->Draw();
	  
	  canv_p->cd();
	  pads[1] = new TPad("pad1", "", 0.0, yPadFrac - (yPadFrac - marg)/2., 1.0, yPadFrac);
	  pads[1]->Draw();
	  pads[1]->SetLeftMargin(marg);
	  pads[1]->SetTopMargin(0.001);
	  pads[1]->SetBottomMargin(0.001);
	  pads[1]->SetRightMargin(0.002);
	  
	  canv_p->cd();
	  pads[2] = new TPad("pad2", "", 0.0, 0.0, 1.0, yPadFrac - (yPadFrac - marg)/2.);
	  pads[2]->Draw();
	  pads[2]->SetLeftMargin(marg);
	  pads[2]->SetTopMargin(0.001);
	  pads[2]->SetBottomMargin(marg/(yPadFrac - (yPadFrac - marg)/2.));
	  pads[2]->SetRightMargin(0.002);
	  
	  canv_p->cd();
	  pads[0]->cd();
	  
	  double max = -1;
	  double min = 10000000;
	  
	  double maxRat = -1;
	  double minRat = 10000000;
	  
	  for(Int_t idI = 0; idI < nID; ++idI){
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->SetMarkerColor(colors[idI]);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->SetLineColor(colors[idI]);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->SetMarkerStyle(styles[idI]);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->SetMarkerSize(1.);	  
	    
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetXaxis()->SetTitleFont(43);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetXaxis()->SetTitleSize(14);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetYaxis()->SetTitleFont(43);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetYaxis()->SetTitleSize(14);
	    
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetXaxis()->SetLabelFont(43);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetXaxis()->SetLabelSize(14);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetYaxis()->SetLabelFont(43);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetYaxis()->SetLabelSize(14);
	    
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetXaxis()->SetNdivisions(505);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetYaxis()->SetNdivisions(505);

	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetXaxis()->SetTitleOffset(5.0);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetYaxis()->SetTitleOffset(2.0);
	    
	    for(Int_t bIX = 0; bIX < jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetNbinsX(); ++bIX){
	      if(jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetBinContent(bIX+1) > max) max = jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetBinContent(bIX+1);
	      if(jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetBinContent(bIX+1) < min && jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetBinContent(bIX+1) > 0) min = jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetBinContent(bIX+1);
	      
	      if(jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetBinContent(bIX+1)/jtPtRaw_RecoGenSymm_h[jI][cI][0][aI][binsI]->GetBinContent(bIX+1) > maxRat) maxRat = jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetBinContent(bIX+1)/jtPtRaw_RecoGenSymm_h[jI][cI][0][aI][binsI]->GetBinContent(bIX+1);
	      if(jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetBinContent(bIX+1)/jtPtRaw_RecoGenSymm_h[jI][cI][0][aI][binsI]->GetBinContent(bIX+1) < minRat) minRat = jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetBinContent(bIX+1)/jtPtRaw_RecoGenSymm_h[jI][cI][0][aI][binsI]->GetBinContent(bIX+1);
	    }
	  }

	  TLegend* leg_p = new TLegend(0.10, 0.1, 0.5, 0.5);
	  leg_p->SetBorderSize(0.0);
	  leg_p->SetFillStyle(0);
	  leg_p->SetFillColor(0);
	  leg_p->SetTextFont(43);
	  leg_p->SetTextSize(14);

	  for(Int_t idI = 0; idI < nID; ++idI){
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->SetMaximum(max*20.);
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->SetMinimum(min/20.);
	    
	    std::string id = idStr[idI] + ", N=" + std::to_string((int)(jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->GetEntries()));
	    leg_p->AddEntry(jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI], id.c_str(), "P L");

	    if(idI == 0) jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->DrawCopy("HIST E1 P");
	    else jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->DrawCopy("HIST E1 P SAME");
	  }
	  
	  TLatex* label_p = new TLatex();
	  label_p->SetTextFont(43);
	  label_p->SetTextSize(14);
	  label_p->SetNDC();
	  
	  label_p->DrawLatex(0.5, 0.94, tempStr.c_str());
	  label_p->DrawLatex(0.5, 0.88, centStr2.c_str());
	  label_p->DrawLatex(0.5, 0.82, jtAbsEtaStr2.c_str());
	  
	  leg_p->Draw("SAME");
	  
	  bool doLogX = false;
	  if(jtPtRaw_RecoGenSymm_h[jI][cI][0][aI][binsI]->GetBinWidth(1)*3 < jtPtRaw_RecoGenSymm_h[jI][cI][0][aI][binsI]->GetBinWidth(jtPtRaw_RecoGenSymm_h[jI][cI][0][aI][binsI]->GetNbinsX()-1)) doLogX = true;
	  
	  gPad->SetLogy();
	  if(doLogX) gPad->SetLogx();
	  gStyle->SetOptStat(0);
	  
	  gPad->SetTicks(1, 2);
	  gPad->RedrawAxis();
	  
	  canv_p->cd();
	  pads[1]->cd();
	  
	  for(Int_t idI = 1; idI < nID; ++idI){
	    TH1D* clone_p = (TH1D*)jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->Clone("temp");
	    clone_p->Divide(jtPtRaw_RecoGenSymm_h[jI][cI][0][aI][binsI]);
	    
	    if(maxRat < 1.05) maxRat = 1.05;
	    if(minRat > 0.95) minRat = 0.95;
	    
	    double interval = maxRat - minRat;
	    maxRat += interval/5.;
	    minRat -= interval/5.;
	    
	    clone_p->SetMaximum(maxRat);
	    clone_p->SetMinimum(minRat);
	    
	    clone_p->GetYaxis()->SetTitle("Ratio");
	    
	    if(idI == 1) clone_p->DrawCopy("HIST E1 P");
	    else clone_p->DrawCopy("HIST E1 P SAME");
	    
	    delete clone_p;
	  }
	  
	  if(doLogX) gPad->SetLogx();
	  gStyle->SetOptStat(0);
	  
	  gPad->SetTicks(1, 2);
	  gPad->RedrawAxis();
	  
	  
	  canv_p->cd();
	  pads[2]->cd();
	  
	  for(Int_t idI = 1; idI < nID; ++idI){
	    TH1D* clone_p = (TH1D*)jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->Clone("temp");
	    clone_p->Divide(jtPtRaw_RecoGenSymm_h[jI][cI][0][aI][binsI]);
	    clone_p->SetMaximum(1.015);
	    clone_p->SetMinimum(0.965);
	    
	    clone_p->GetYaxis()->SetTitle("Ratio (Zoom)");
	    
	    if(idI == 1) clone_p->DrawCopy("HIST E1 P");
	    else clone_p->DrawCopy("HIST E1 P SAME");
	    delete clone_p;
	  }
	  
	  drawWhiteBox(900, 1100, .93, .9649);

	  label_p->SetNDC(0);
	  label_p->DrawLatex(100, .955, "100");
	  label_p->DrawLatex(200, .955, "200");
	  label_p->DrawLatex(400, .955, "400");
	  label_p->DrawLatex(600, .955, "600");
	  label_p->DrawLatex(1000, .955, "1000");
	  label_p->DrawLatex(2000, .955, "2000");
	  
	  if(doLogX) gPad->SetLogx();
	  gStyle->SetOptStat(0);
	  
	  gPad->SetTicks(1, 2);
	  gPad->RedrawAxis();
       
	  std::string spectraName = "jtPtRaw_" + smallLargeBinsStr[binsI] + "_" + tempStr + "_" + centStr + "_" + jtAbsEtaStr + "_" + dateStr + ".pdf";
	  spectraNames.push_back(spectraName);
	  
	  const std::string saveName = "pdfDir/" + dateStr + "/Process/" + spectraName;
	  quietSaveAs(canv_p, saveName);
	  
	  delete pads[0];
	  delete pads[1];
	  delete pads[2];
	  delete canv_p;
	  delete leg_p;
	  delete label_p;
	}
      }
    
      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	const std::string jtAbsEtaStr2 = prettyString(jtAbsEtaBinsLow[aI], 1, false) + "<|#eta|<" + prettyString(jtAbsEtaBinsHi[aI], 1, false);
      
	multijetNames.push_back({});
	multijetNames.push_back({});
            
	for(Int_t jetI = 0; jetI < nRecoJtPtBins[jI][cI]; ++jetI){
	  if(multijetAJ_All_h[jI][cI][0][aI][jetI]->GetEntries() == 0) continue;

	  const std::string jtPtStr = "JtPt" + prettyString(recoJtPtBins[jI][cI][jetI], 1, true) + "to" + prettyString(recoJtPtBins[jI][cI][jetI+1], 1, true);
	  const std::string jtPtStr2 =  prettyString(recoJtPtBins[jI][cI][jetI], 1, false) + "<p_{T}<" + prettyString(recoJtPtBins[jI][cI][jetI+1], 1, false);
	
	  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
	  canv_p->SetTopMargin(0.01);
	  canv_p->SetRightMargin(0.002);
	  canv_p->SetLeftMargin(0.01);
	  canv_p->SetBottomMargin(0.01);
	  
	  TPad* pads[2];
	  canv_p->cd();
	  pads[0] = new TPad("pad0", "", 0.0, yPadFrac, 1.0, 1.0);
	  pads[0]->SetLeftMargin(marg);
	  pads[0]->SetTopMargin(0.01);
	  pads[0]->SetBottomMargin(0.001);
	  pads[0]->SetRightMargin(0.002);
	  pads[0]->Draw();
	  
	  canv_p->cd();
	  pads[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, yPadFrac);
	  pads[1]->Draw();
	  pads[1]->SetLeftMargin(marg);
	  pads[1]->SetTopMargin(0.001);
	  pads[1]->SetBottomMargin(marg/yPadFrac);
	  pads[1]->SetRightMargin(0.002);
	  
	  canv_p->cd();
	  pads[0]->cd();

	  double max = -1;
	  double min = 10000000;
	  
	  for(Int_t idI = 0; idI < nID; ++idI){
	    Int_t passColor = vg.getColor(1);
	    if(idI == 0) passColor = 1;
	    Int_t failColor = vg.getColor(2);

	    for(Int_t bIX = 0; bIX < multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetNbinsX(); ++bIX){
	      /*
	      multijetAJ_All_h[jI][cI][idI][aI][jetI]->SetBinContent(bIX+1, multijetAJ_All_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1)/multijetAJFactors[bIX]);
	      multijetAJ_All_h[jI][cI][idI][aI][jetI]->SetBinError(bIX+1, multijetAJ_All_h[jI][cI][idI][aI][jetI]->GetBinError(bIX+1)/multijetAJFactors[bIX]);

	      multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->SetBinContent(bIX+1, multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1)/multijetAJFactors[bIX]);
	      multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->SetBinError(bIX+1, multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetBinError(bIX+1)/multijetAJFactors[bIX]);

	      multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->SetBinContent(bIX+1, multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1)/multijetAJFactors[bIX]);
	      multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->SetBinError(bIX+1, multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetBinError(bIX+1)/multijetAJFactors[bIX]);
	      */
	    }

	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->SetMarkerColor(passColor);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->SetLineColor(passColor);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->SetMarkerStyle(styles[idI]);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->SetMarkerSize(1.);	  
	    
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetXaxis()->SetTitleFont(43);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetXaxis()->SetTitleSize(14);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetYaxis()->SetTitleFont(43);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetYaxis()->SetTitleSize(14);
	    
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetXaxis()->SetLabelFont(43);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetXaxis()->SetLabelSize(14);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetYaxis()->SetLabelFont(43);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetYaxis()->SetLabelSize(14);
	    
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetXaxis()->SetNdivisions(505);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetYaxis()->SetNdivisions(505);
	    
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetXaxis()->SetTitleOffset(3.0);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetYaxis()->SetTitleOffset(2.0);

	    multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->SetMarkerColor(failColor);
            multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->SetLineColor(failColor);
            multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->SetMarkerStyle(styles[idI]);
            multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->SetMarkerSize(1.);

	    
	    for(Int_t bIX = 0; bIX < multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetNbinsX(); ++bIX){
	      if(multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1) > max) max = multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1);
	      if(multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1) < min && multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1) > 0) min = multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1);
	      
	      if(multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1) > max) max = multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1);
	      if(multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1) < min && multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1) > 0) min = multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1);
	    }
	  }
	  
	  TLegend* leg_p = new TLegend(0.6, 0.75, 0.9, 0.95);
	  leg_p->SetBorderSize(0.0);
	  leg_p->SetFillStyle(0);
	  leg_p->SetFillColor(0);
	  leg_p->SetTextFont(43);
	  leg_p->SetTextSize(14);

	  for(Int_t idI = 0; idI < nID; ++idI){
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->SetMaximum(max*20.);
	    multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->SetMinimum(min/20.);

	    std::string id = idStr[idI];
	    leg_p->AddEntry(multijetAJ_Pass_h[jI][cI][idI][aI][jetI], id.c_str(), "P L");

	    if(idI == 0) multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->DrawCopy("HIST E1 P");
	    else{
	      multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->DrawCopy("HIST E1 P SAME");
	      multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->DrawCopy("HIST E1 P SAME");
	    }
	  }

	  TLatex* label_p = new TLatex();
	  label_p->SetTextFont(43);
	  label_p->SetTextSize(14);
	  label_p->SetNDC();

	  label_p->DrawLatex(0.2, 0.94, tempStr.c_str());
	  label_p->DrawLatex(0.2, 0.88, centStr2.c_str());
	  label_p->DrawLatex(0.2, 0.82, jtAbsEtaStr2.c_str());
	  label_p->DrawLatex(0.2, 0.76, jtPtStr2.c_str());

	  leg_p->Draw("SAME");

	  gPad->SetLogy();
	  gStyle->SetOptStat(0);
	  gPad->SetTicks(1, 2);
	  gPad->RedrawAxis();
	  
	  canv_p->cd();
	  pads[1]->cd();
   	  
	  label_p->SetTextSize(12);

	  for(Int_t idI = 1; idI < nID; ++idI){
	    TH1D* clone_p = (TH1D*)multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->Clone("temp");
	    clone_p->Divide(multijetAJ_All_h[jI][cI][0][aI][jetI]);
	    double totalFail = 0;
	    for(Int_t bIX = 0; bIX < multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetNbinsX(); ++bIX){
	      totalFail += multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1);
	    }
	    double highFail = multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetBinContent(multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetNbinsX());
	    //	    std::cout << "High fail 1: " << highFail << std::endl;
	    //	    highFail *= multijetAJFactors[multijetAJ_Fail_h[jI][cI][idI][aI][jetI]->GetNbinsX()-1];
	    //	    std::cout << " " << highFail << std::endl;

	    for(Int_t bIX = 0; bIX < multijetAJ_All_h[jI][cI][0][aI][jetI]->GetNbinsX(); ++bIX){
	      if(multijetAJ_All_h[jI][cI][0][aI][jetI]->GetBinContent(bIX+1) > 0 && multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1) == 0) clone_p->SetBinContent(bIX+1, 0.0000001);
	      else if(multijetAJ_Pass_h[jI][cI][idI][aI][jetI]->GetBinContent(bIX+1) == 0) clone_p->SetBinContent(bIX+1, -100.);
	      clone_p->SetBinError(bIX+1, 0.0000001);
	    }

	    clone_p->SetMaximum(1.1);
	    clone_p->SetMinimum(-0.1);

	    clone_p->GetYaxis()->SetTitle("Pass Ratio");
	    
	    if(idI == 1) clone_p->DrawCopy("HIST E1 P");
	    else clone_p->DrawCopy("HIST E1 P SAME");

	    std::string labelStr = idStr[idI] + ": " + std::to_string((int)highFail) + "/" + std::to_string((int)totalFail);
	    if(((int)totalFail) > 0) labelStr = labelStr + ", " + prettyString(highFail/totalFail, 3, false);

	    label_p->DrawLatex(0.15, 0.8 - 0.09*(idI-1), labelStr.c_str());
	    
	    delete clone_p;
	  }

	  gStyle->SetOptStat(0);
	  gPad->SetTicks(1, 2);
	  gPad->RedrawAxis();
	
	  std::string multijetName = "multijetAJ_" + tempStr + "_" + centStr + "_" + jtAbsEtaStr + "_" + jtPtStr + "_" + dateStr + ".pdf";
	  multijetNames[multijetNames.size()-2].push_back(multijetName);
	  multijetNames[multijetNames.size()-1].push_back(multijetName);
	  const std::string saveName = "pdfDir/" + dateStr + "/Process/" + multijetName;
       	  quietSaveAs(canv_p, saveName);

	  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	  
	  delete pads[0];
	  delete pads[1];
	  delete canv_p;
	  delete leg_p;	  		  
	  delete label_p;

	  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	}
      }
    }
  }  

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  //  spectraNames
  //  multijetNames

  const std::string textWidth = "0.24";
  std::string texFileName = outFileName;
  texFileName.replace(texFileName.find(".root"), std::string(".root").size(), ".tex");
  texFileName.replace(0, std::string("output").size(), ("pdfDir/" + dateStr + "/Process").c_str());
  
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

  texFile << "\\setbeamertemplate{frametitle}" << std::endl;
  texFile << "{" << std::endl;
  texFile << "    \\nointerlineskip" << std::endl;
  texFile << "    \\begin{beamercolorbox}[sep=0.1cm, ht=1.0em, wd=\\paperwidth]{frametitle}" << std::endl;
  texFile << "        \\vbox{}\\vskip-2ex%" << std::endl;
  texFile << "        \\strut\\insertframetitle\\strut" << std::endl;
  texFile << "        \\vskip-0.8ex%" << std::endl;
  texFile << "    \\end{beamercolorbox}" << std::endl;
  texFile << "}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamertemplate{footline}{%" << std::endl;
  texFile << "  \\begin{beamercolorbox}[sep=.4em,wd=\\paperwidth,leftskip=0.5cm,rightskip=0.5cm]{footlinecolor}" << std::endl;
  texFile << "    \\hspace{0.15cm}%" << std::endl;
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

  texFile << "\\author[CM]{Placeholder}" << std::endl;
  texFile << std::endl;

  texFile << "\\begin{document}" << std::endl;
  texFile << "\\begin{frame}" << std::endl;
  texFile << "\\frametitle{\\centerline{JEC Validation (" << dateStr2 << ")}}" << std::endl;
  texFile << " \\begin{itemize}" << std::endl;
  texFile << "  \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "  \\item{Placeholder}" << std::endl;
  texFile << "  \\begin{itemize}" << std::endl;
  texFile << "   \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "   \\item{Placeholder}" << std::endl;
  texFile << "  \\end{itemize}" << std::endl;
  texFile << " \\end{itemize}" << std::endl;
  texFile << "\\end{frame}" << std::endl;

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int spI = 0; spI < spectraNames.size(); ++spI){
    if(spectraNames[spI].find("AbsEta0p0to2p0") == std::string::npos) continue;
    
    std::string centStr = "PP";
    
    if(!isDataPP){
      centStr = spectraNames[spI].substr(spectraNames[spI].find("Cent"), spectraNames[spI].size());
      centStr.replace(centStr.find("_"), centStr.size(), "");
    }

    std::string jtStr = spectraNames[spI].substr(spectraNames[spI].find("_ak")+1, spectraNames[spI].size());
    jtStr.replace(jtStr.find("_"), jtStr.size(), "");

    texFile << "\\begin{frame}" << std::endl;
    texFile << "\\frametitle{\\centerline{" << jtStr << ", " << centStr << "}}" << std::endl;
    texFile << "\\includegraphics[width=" << textWidth << "\\textwidth]{" << spectraNames[spI] << "}" << std::endl;
    
    for(unsigned int mI = 0; mI < multijetNames[spI].size(); ++mI){    
      texFile << "\\includegraphics[width=" << textWidth << "\\textwidth]{" << multijetNames[spI][mI] << "}";
      if(mI == 2 || mI == 6) texFile << "\\\\";
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

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    outFile_p->cd();
    dir_p[jI]->cd();

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){

	  for(Int_t binsI = 0; binsI < nSmallLargeBins; ++binsI){
	    jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI]->Write("", TObject::kOverwrite);
	    delete jtPtRaw_RecoGenSymm_h[jI][cI][idI][aI][binsI];
	  }
	  
	  jtPtRaw_RecoGenAsymm_h[jI][cI][idI][aI]->Write("", TObject::kOverwrite);
	  delete jtPtRaw_RecoGenAsymm_h[jI][cI][idI][aI];

	  for(Int_t jIter = 0; jIter < nRecoJtPtBins[jI][cI]; ++jIter){
	    multijetAJ_All_h[jI][cI][idI][aI][jIter]->Write("", TObject::kOverwrite);
	    delete multijetAJ_All_h[jI][cI][idI][aI][jIter];
	    
	    multijetAJ_Pass_h[jI][cI][idI][aI][jIter]->Write("", TObject::kOverwrite);
	    delete multijetAJ_Pass_h[jI][cI][idI][aI][jIter];
	    
	    multijetAJ_Fail_h[jI][cI][idI][aI][jIter]->Write("", TObject::kOverwrite);
	    delete multijetAJ_Fail_h[jI][cI][idI][aI][jIter];
	  }
	}
      }
    }
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  outFile_p->cd();
  TDirectory* cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");
  TDirectory* subFileDir_p = (TDirectory*)outFile_p->mkdir("cutDir/fullFiles");

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  cutPropagator cutPropOut;
  cutPropOut.Clean();
  cutPropOut.SetIsPP(isDataPP);
  cutPropOut.SetJtAbsEtaMax(jtAbsEtaMax);
  cutPropOut.SetNJtAlgos(nDataJet);
  cutPropOut.SetJtAlgos(dataTreeList);
  cutPropOut.SetMinJtPtCut(nDataJet, minJtPtCut);
  cutPropOut.SetMultiJtPtCut(nDataJet, multiJtPtCut);
  cutPropOut.SetRecoTruncPos(nDataJet, recoTruncPos);
  cutPropOut.SetNSmallR(cutProp.GetNSmallR());
  cutPropOut.SetSmallRVals(cutProp.GetSmallRVals());
  cutPropOut.SetNLargeR(cutProp.GetNLargeR());
  cutPropOut.SetLargeRVals(cutProp.GetLargeRVals());

  cutPropOut.SetNRecoJtPtBinsSmallRCent0to10(nRecoJtPtBinsSmallRCent0to10Temp);
  cutPropOut.SetNRecoJtPtBinsLargeRCent0to10(nRecoJtPtBinsLargeRCent0to10Temp);
  cutPropOut.SetNGenJtPtSmallBinsSmallRCent0to10(cutProp.GetNGenJtPtSmallBinsSmallRCent0to10());
  cutPropOut.SetNGenJtPtLargeBinsSmallRCent0to10(cutProp.GetNGenJtPtLargeBinsSmallRCent0to10());
  cutPropOut.SetNGenJtPtSmallBinsLargeRCent0to10(cutProp.GetNGenJtPtSmallBinsLargeRCent0to10());
  cutPropOut.SetNGenJtPtLargeBinsLargeRCent0to10(cutProp.GetNGenJtPtLargeBinsLargeRCent0to10());

  cutPropOut.SetNRecoJtPtBinsSmallRCent10to30(nRecoJtPtBinsSmallRCent10to30Temp);
  cutPropOut.SetNRecoJtPtBinsLargeRCent10to30(nRecoJtPtBinsLargeRCent10to30Temp);
  cutPropOut.SetNGenJtPtSmallBinsSmallRCent10to30(cutProp.GetNGenJtPtSmallBinsSmallRCent10to30());
  cutPropOut.SetNGenJtPtLargeBinsSmallRCent10to30(cutProp.GetNGenJtPtLargeBinsSmallRCent10to30());
  cutPropOut.SetNGenJtPtSmallBinsLargeRCent10to30(cutProp.GetNGenJtPtSmallBinsLargeRCent10to30());
  cutPropOut.SetNGenJtPtLargeBinsLargeRCent10to30(cutProp.GetNGenJtPtLargeBinsLargeRCent10to30());

  cutPropOut.SetNRecoJtPtBinsSmallRCent30to50(nRecoJtPtBinsSmallRCent30to50Temp);
  cutPropOut.SetNRecoJtPtBinsLargeRCent30to50(nRecoJtPtBinsLargeRCent30to50Temp);
  cutPropOut.SetNGenJtPtSmallBinsSmallRCent30to50(cutProp.GetNGenJtPtSmallBinsSmallRCent30to50());
  cutPropOut.SetNGenJtPtLargeBinsSmallRCent30to50(cutProp.GetNGenJtPtLargeBinsSmallRCent30to50());
  cutPropOut.SetNGenJtPtSmallBinsLargeRCent30to50(cutProp.GetNGenJtPtSmallBinsLargeRCent30to50());
  cutPropOut.SetNGenJtPtLargeBinsLargeRCent30to50(cutProp.GetNGenJtPtLargeBinsLargeRCent30to50());

  cutPropOut.SetNRecoJtPtBinsSmallRCent50to90(nRecoJtPtBinsSmallRCent50to90Temp);
  cutPropOut.SetNRecoJtPtBinsLargeRCent50to90(nRecoJtPtBinsLargeRCent50to90Temp);
  cutPropOut.SetNGenJtPtSmallBinsSmallRCent50to90(cutProp.GetNGenJtPtSmallBinsSmallRCent50to90());
  cutPropOut.SetNGenJtPtLargeBinsSmallRCent50to90(cutProp.GetNGenJtPtLargeBinsSmallRCent50to90());
  cutPropOut.SetNGenJtPtSmallBinsLargeRCent50to90(cutProp.GetNGenJtPtSmallBinsLargeRCent50to90());
  cutPropOut.SetNGenJtPtLargeBinsLargeRCent50to90(cutProp.GetNGenJtPtLargeBinsLargeRCent50to90());

  cutPropOut.SetRecoJtPtBinsSmallR(recoJtPtBinsSmallRTemp);
  cutPropOut.SetRecoJtPtBinsLargeR(recoJtPtBinsLargeRTemp);
  cutPropOut.SetGenJtPtSmallBinsSmallR(cutProp.GetGenJtPtSmallBinsSmallR());
  cutPropOut.SetGenJtPtLargeBinsSmallR(cutProp.GetGenJtPtLargeBinsSmallR());
  cutPropOut.SetGenJtPtSmallBinsLargeR(cutProp.GetGenJtPtSmallBinsLargeR());
  cutPropOut.SetGenJtPtLargeBinsLargeR(cutProp.GetGenJtPtLargeBinsLargeR());
  cutPropOut.SetNJtAbsEtaBins(nJtAbsEtaBins);
  cutPropOut.SetJtAbsEtaBinsLow(nJtAbsEtaBins, jtAbsEtaBinsLow);
  cutPropOut.SetJtAbsEtaBinsHi(nJtAbsEtaBins, jtAbsEtaBinsHi);
  cutPropOut.SetNPthats(0);
  cutPropOut.SetPthats({});
  cutPropOut.SetPthatWeights({});
  cutPropOut.SetNCentBins(nCentBins);
  cutPropOut.SetCentBinsLow(centBinsLow);
  cutPropOut.SetCentBinsHi(centBinsHi);
  cutPropOut.SetNID(nID);
  cutPropOut.SetIdStr(idStr);
  cutPropOut.SetJtPfCHMFCutLow(jtPfCHMFCutLow);
  cutPropOut.SetJtPfCHMFCutHi(jtPfCHMFCutHi);
  cutPropOut.SetJtPfMUMFCutLow(jtPfMUMFCutLow);
  cutPropOut.SetJtPfMUMFCutHi(jtPfMUMFCutHi);
  cutPropOut.SetJtPfNHFCutLow(jtPfNHFCutLow);
  cutPropOut.SetJtPfNHFCutHi(jtPfNHFCutHi);
  cutPropOut.SetJtPfNEFCutLow(jtPfNEFCutLow);
  cutPropOut.SetJtPfNEFCutHi(jtPfNEFCutHi);
  cutPropOut.SetJtPfMUFCutLow(jtPfMUFCutLow);
  cutPropOut.SetJtPfMUFCutHi(jtPfMUFCutHi);
  cutPropOut.SetJtPfCHFCutLow(jtPfCHFCutLow);
  cutPropOut.SetJtPfCHFCutHi(jtPfCHFCutHi);
  cutPropOut.SetJtPfCEFCutLow(jtPfCEFCutLow);
  cutPropOut.SetJtPfCEFCutHi(jtPfCEFCutHi);
  cutPropOut.SetJtPfMinMult(jtPfMinMult);
  cutPropOut.SetJtPfMinChgMult(jtPfMinChgMult);

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(!cutPropOut.WriteAllVarToFile(outFile_p, cutDir_p, subFileDir_p)) std::cout << "Warning: Cut writing has failed" << std::endl;

  outFile_p->Close();
  delete outFile_p;

  responseFile_p->Close();
  delete responseFile_p;

  timer.stop();
  std::cout << "Total run time wall: " << timer.totalWall() << std::endl;
  std::cout << "Total run time cpu: " << timer.totalCPU() << std::endl; 
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3 && argc != 4 && argc != 5){
    std::cout << "Usage: ./bin/processRawData.exe <inDataFileName> <inResponseName> <isPP-opt> <tagStr-opt>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 3) retVal += processRawData(argv[1], argv[2]);
  else if(argc == 4) retVal += processRawData(argv[1], argv[2], std::stoi(argv[3]));
  else if(argc == 5) retVal += processRawData(argv[1], argv[2], std::stoi(argv[3]), argv[4]);

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
