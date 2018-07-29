#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDatime.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/etaPhiFunc.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/ncollFunctions_5TeV.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/smearingFuncs.h"
#include "Utility/include/specialHYDJETEventExclude.h"
#include "Utility/include/vanGoghPalette.h"

int checkGoodJetBadJetPF(const std::string inName, bool isPP = false)
{
  std::vector<std::string> fileList;
  std::vector<double> pthats;
  std::vector<double> pthatWeights;

  if(inName.find(".root") != std::string::npos){
    fileList.push_back(inName);
    pthats.push_back(1.);
    pthatWeights.push_back(1.);
  }
  else if(inName.find(".txt") != std::string::npos){
    std::ifstream file(inName.c_str());
    std::string tempStr;

    while(std::getline(file, tempStr)){
      while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
      else{
        if(tempStr.substr(0, std::string("PTHAT=").size()).find("PTHAT=") != std::string::npos){
          tempStr.replace(0, std::string("PTHAT=").size(), "");
          while(tempStr.find(",") != std::string::npos){
            pthats.push_back(std::stod(tempStr.substr(0, tempStr.find(","))));
            tempStr.replace(0, tempStr.find(",")+1, "");
          }
          if(tempStr.size() != 0) pthats.push_back(std::stod(tempStr));
        }
        else if(tempStr.substr(0, std::string("PTHATWEIGHTS=").size()).find("PTHATWEIGHTS=") != std::string::npos){
          tempStr.replace(0, std::string("PTHATWEIGHTS=").size(), "");
          while(tempStr.find(",") != std::string::npos){
            pthatWeights.push_back(std::stod(tempStr.substr(0, tempStr.find(","))));
            tempStr.replace(0, tempStr.find(",")+1, "");
          }
          if(tempStr.size() != 0) pthatWeights.push_back(std::stod(tempStr));
        }
        else if(tempStr.substr(0, std::string("ISPP=").size()).find("ISPP=") != std::string::npos){
          tempStr.replace(0,tempStr.find("=")+1, "");
          isPP = std::stoi(tempStr);
        }
        else std::cout << "WARNING: Line in \'" << inName << "\', \'" << tempStr << "\' is invalid. check input" << std::endl;
      }
    }

    file.close();
  }
  else{
    std::cout << "Given inName \'" << inName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  if(fileList.size() == 0){
    std::cout << "Given inName \'" << inName << "\' is gives no valid root files. return 1" << std::endl;
    return 1;
  }
  else if(pthats.size() == 0){
    std::cout << "Given inName \'" << inName << "\' contains no pthat list. return 1" << std::endl;
    return 1;
  }
  else if(pthatWeights.size() == 0){
    std::cout << "Given inName \'" << inName << "\' contains no pthatWeights list. return 1" << std::endl;
    return 1;
  }

  std::cout << "Pthats and weights: " << std::endl;
  for(unsigned int pI = 0; pI < pthats.size(); ++pI){
    std::cout << " " << pI << "/" << pthats.size() << ": " << pthats.at(pI) << ", " << pthatWeights.at(pI) << std::endl;
  }

  std::cout << "IsPP: " << isPP << std::endl;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  std::string outFileName = inName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
  else if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  outFileName = "output/" + outFileName + "_CheckGoodJetBadJetPF_" + dateStr + ".root";
  checkMakeDir("output");

  const Int_t nCentBinsPerma = 4;
  const Int_t centBinsLowPerma[nCentBinsPerma] = {0, 10, 30, 50};
  const Int_t centBinsHiPerma[nCentBinsPerma] = {10, 30, 50, 90};

  Int_t nCentBinsTemp = nCentBinsPerma;
  if(isPP) nCentBinsTemp = 1;
  
  const Int_t nCentBins = nCentBinsTemp;
  std::vector<Int_t> centBinsLow;
  std::vector<Int_t> centBinsHi;

  if(isPP){
    centBinsLow.push_back(0);
    centBinsHi.push_back(100);
  }
  else{
    for(Int_t cI = 0; cI < nCentBinsPerma; ++cI){
      centBinsLow.push_back(centBinsLowPerma[cI]);
      centBinsHi.push_back(centBinsHiPerma[cI]);
    }
  }

  const Double_t jtRecoPtMin = 100.;
  const Double_t jtRecoPtMax = 1100.;
  const Int_t nJtRecoPtBins = 40;
  Double_t jtRecoPtBins[nJtRecoPtBins+1];
  getLinBins(jtRecoPtMin, jtRecoPtMax, nJtRecoPtBins, jtRecoPtBins);

  const Int_t nJtRecoPtFineBins = 200;
  Double_t jtRecoPtFineBins[nJtRecoPtFineBins+1];
  getLinBins(jtRecoPtMin, jtRecoPtMax, nJtRecoPtFineBins, jtRecoPtFineBins);

  //  const Double_t jtRefGoodPtMin = 350.;
  //  const Double_t jtRefBadPtMax = 200.;
  const Double_t jtAbsEtaMax = 2.;

  TFile* inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(0)).c_str(), "READ");
  std::vector<std::string> jetList = returnRootFileContentsList(inFile_p, "TTree", "JetAnalyzer");
  
  unsigned int pos = 0;
  while(jetList.size() > pos){
    if(jetList.at(pos).find("akPu") != std::string::npos) jetList.erase(jetList.begin()+pos);
    //    else if(jetList.at(pos).find("akCs4PU3PFFlow") != std::string::npos) ++pos;
    else ++pos;
  }
  

  inFile_p->Close();
  delete inFile_p;
  inFile_p = NULL;

  const Int_t nJetAlgos = jetList.size();
  std::vector<std::string> jetAlgos;
  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    jetAlgos.push_back(jetList.at(jI).substr(0, jetList.at(jI).find("JetAna")));
    std::cout << jetAlgos.at(jI) << std::endl;
  }

  Int_t posR4Temp = -1;
  for(int jI = 0; jI < nJetAlgos; ++jI){
    if(isPP){
      if(jetList.at(jI).find("ak4") != std::string::npos && posR4Temp < 0) posR4Temp = jI;
    }
    else{
      if(jetList.at(jI).find("akCs4PU3PFFlow") != std::string::npos && posR4Temp < 0) posR4Temp = jI;
    }
  }
  const Int_t posR4 = posR4Temp;

  std::cout << "PosR4: " << posR4 << std::endl;
  if(posR4 < 0){
    return 1;
  }

  const Int_t nJtRecoPtBins2 = 2;
  const Int_t jtRecoPtBins2Low[nJtRecoPtBins2] ={100, 300};
  const Int_t jtRecoPtBins2Hi[nJtRecoPtBins2] ={300, 1100};

  const Int_t nWeighted = 2;
  const std::string weightStr[nWeighted] = {"Unweighted", "Weighted"};

  //FULL ID taken from here https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016, FullLight and FullTight correspond to previous levels of cuts + the Loose and TightLepVeto versions, respectively                                                                                                                
    const Int_t nID = 5;
  const std::string idStr[nID] = {"NoID", "LightMUID", "LightMUAndCHID", "FullLight", "FullTight"};
  const Double_t jtPfCHMFCutLow[nID] = {0.0, 0.0, 0.00, 0.00, 0.00};
  const Double_t jtPfCHMFCutHi[nID] = {1.0, 1.0, 0.90, 0.90, 0.90};
  const Double_t jtPfMUMFCutLow[nID] = {0.0, 0.0, 0.00, 0.00, 0.00};
  const Double_t jtPfMUMFCutHi[nID] = {1.0, 0.60, 0.60, 0.60, 0.60};
  const Double_t jtPfNHFCutLow[nID] = {0.0, 0.0, 0.0, 0.0, 0.0};
  const Double_t jtPfNHFCutHi[nID] = {1.0, 1.0, 1.0, 0.99, 0.90};
  const Double_t jtPfNEFCutLow[nID] = {0.0, 0.0, 0.0, 0.0, 0.0};
  const Double_t jtPfNEFCutHi[nID] = {1.0, 1.0, 1.0, 0.99, 0.90};
  const Int_t jtPfMinMult[nID] = {1, 1, 1, 2, 2};
  const Double_t jtPfMUFCutLow[nID] = {0.0, 0.0, 0.0, 0.0, 0.0};
  const Double_t jtPfMUFCutHi[nID] = {1.0, 1.0, 1.0, 1.0, 0.8};
  const Double_t jtPfCHFCutLow[nID] = {0.0, 0.0, 0.0, 0.0000001, 0.0000001};
  const Double_t jtPfCHFCutHi[nID] = {1.0, 1.0, 1.0, 1.0, 1.0};
  const Int_t jtPfMinChgMult[nID] = {0, 0, 0, 1, 1};
  const Double_t jtPfCEFCutLow[nID] = {0.0, 0.0, 0.0, 0.0, 0.0};
  const Double_t jtPfCEFCutHi[nID] = {1.0, 1.0, 1.0, 0.99, 0.90};
  

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);

  TH1D* jtPt_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH1D* jtPt_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH1D* jtPt_FineTot_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH1D* jtPfCHF_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfCEF_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNHF_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNEF_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfMUF_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];

  TH2D* jtPfCHF_jtPfCEF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHF_jtPfNHF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHF_jtPfNEF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHF_jtPfMUF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfCEF_jtPfNHF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCEF_jtPfNEF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCEF_jtPfMUF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfNHF_jtPfNEF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfNHF_jtPfMUF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfNEF_jtPfMUF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH1D* jtPfCHMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfCEMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNHMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNEMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfMUMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];

  TH2D* jtPfCHMF_jtPfCEMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHMF_jtPfNHMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHMF_jtPfNEMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHMF_jtPfMUMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfCEMF_jtPfNHMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCEMF_jtPfNEMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCEMF_jtPfMUMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfNHMF_jtPfNEMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfNHMF_jtPfMUMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfNEMF_jtPfMUMF_Good_h[nJetAlgos][nCentBins][nWeighted][nID];


  TH1D* jtPfCHM_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfCEM_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNHM_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNEM_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfMUM_Good_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];

  TH1D* jtPfCHF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfCEF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNHF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNEF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfMUF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];

  TH2D* jtPfCHF_jtPfCEF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHF_jtPfNHF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHF_jtPfNEF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHF_jtPfMUF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfCEF_jtPfNHF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCEF_jtPfNEF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCEF_jtPfMUF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfNHF_jtPfNEF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfNHF_jtPfMUF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfNEF_jtPfMUF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH1D* jtPfCHMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfCEMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNHMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNEMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfMUMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];

  TH2D* jtPfCHMF_jtPfCEMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHMF_jtPfNHMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHMF_jtPfNEMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCHMF_jtPfMUMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfCEMF_jtPfNHMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCEMF_jtPfNEMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfCEMF_jtPfMUMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfNHMF_jtPfNEMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];
  TH2D* jtPfNHMF_jtPfMUMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH2D* jtPfNEMF_jtPfMUMF_Bad_h[nJetAlgos][nCentBins][nWeighted][nID];

  TH1D* jtPfCHM_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfCEM_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNHM_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfNEM_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];
  TH1D* jtPfMUM_Bad_h[nJetAlgos][nCentBins][nWeighted][nID][nJtRecoPtBins2];

  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
      if(isPP) centStr = "PP";
      
      for(Int_t wI = 0; wI < nWeighted; ++wI){

	for(Int_t idI = 0; idI < nID; ++idI){
	  jtPt_Good_h[jI][cI][wI][idI] = new TH1D(("jtPt_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";Jet p_{T} (Good);Counts", nJtRecoPtBins, jtRecoPtBins);
	  jtPt_Bad_h[jI][cI][wI][idI] = new TH1D(("jtPt_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";Jet p_{T} (Bad);Counts", nJtRecoPtBins, jtRecoPtBins);
	  jtPt_FineTot_h[jI][cI][wI][idI] = new TH1D(("jtPt_FineTot_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";Jet p_{T} (Fine Tot);Counts", nJtRecoPtFineBins, jtRecoPtFineBins);

	  for(Int_t rI = 0; rI < nJtRecoPtBins2; ++rI){
	    std::string recoPtStr = "RecoJtPt" + std::to_string(jtRecoPtBins2Low[rI]) + "to" + std::to_string(jtRecoPtBins2Hi[rI]);

	    jtPfCHF_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCHF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCHF;Counts", 10, 0.0, 1.0);
	    jtPfCEF_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCEF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCEF;Counts", 10, 0.0, 1.0);
	    jtPfNHF_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNHF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNHF;Counts", 10, 0.0, 1.0);
	    jtPfNEF_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNEF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNEF;Counts", 10, 0.0, 1.0);
	    jtPfMUF_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfMUF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfMUF;Counts", 10, 0.0, 1.0);

	    jtPfCHMF_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCHMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCHMF;Counts", 10, 0.0, 1.0);
	    jtPfCEMF_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCEMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCEMF;Counts", 10, 0.0, 1.0);
	    jtPfNHMF_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNHMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNHMF;Counts", 10, 0.0, 1.0);
	    jtPfNEMF_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNEMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNEMF;Counts", 10, 0.0, 1.0);
	    jtPfMUMF_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfMUMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfMUMF;Counts", 10, 0.0, 1.0);

	    jtPfCHM_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCHM_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCHM;Counts", 20, 0.0, 40.0);
	    jtPfCEM_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCEM_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCEM;Counts", 20, 0.0, 40.0);
	    jtPfNHM_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNHM_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNHM;Counts", 20, 0.0, 40.0);
	    jtPfNEM_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNEM_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNEM;Counts", 20, 0.0, 40.0);
	    jtPfMUM_Good_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfMUM_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfMUM;Counts", 20, 0.0, 40.0);
	    
	    jtPfCHF_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCHF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCHF;Counts", 10, 0.0, 1.0);
	    jtPfCEF_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCEF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCEF;Counts", 10, 0.0, 1.0);
	    jtPfNHF_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNHF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNHF;Counts", 10, 0.0, 1.0);
	    jtPfNEF_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNEF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNEF;Counts", 10, 0.0, 1.0);
	    jtPfMUF_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfMUF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfMUF;Counts", 10, 0.0, 1.0);

	    jtPfCHMF_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCHMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCHMF;Counts", 10, 0.0, 1.0);
	    jtPfCEMF_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCEMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCEMF;Counts", 10, 0.0, 1.0);
	    jtPfNHMF_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNHMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNHMF;Counts", 10, 0.0, 1.0);
	    jtPfNEMF_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNEMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNEMF;Counts", 10, 0.0, 1.0);
	    jtPfMUMF_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfMUMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfMUMF;Counts", 10, 0.0, 1.0);

	    jtPfCHM_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCHM_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCHM;Counts", 20, 0.0, 40.0);
	    jtPfCEM_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfCEM_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfCEM;Counts", 20, 0.0, 40.0);
	    jtPfNHM_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNHM_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNHM;Counts", 20, 0.0, 40.0);
	    jtPfNEM_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfNEM_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfNEM;Counts", 20, 0.0, 40.0);
	    jtPfMUM_Bad_h[jI][cI][wI][idI][rI] = new TH1D(("jtPfMUM_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_" + recoPtStr + "_h").c_str(), ";jtPfMUM;Counts", 20, 0.0, 40.0);

	    std::vector<TH1*> tempVect = {jtPfCHF_Good_h[jI][cI][wI][idI][rI], jtPfCEF_Good_h[jI][cI][wI][idI][rI], jtPfNHF_Good_h[jI][cI][wI][idI][rI], jtPfNEF_Good_h[jI][cI][wI][idI][rI], jtPfMUF_Good_h[jI][cI][wI][idI][rI], jtPfCHM_Good_h[jI][cI][wI][idI][rI], jtPfCEM_Good_h[jI][cI][wI][idI][rI], jtPfNHM_Good_h[jI][cI][wI][idI][rI], jtPfNEM_Good_h[jI][cI][wI][idI][rI], jtPfMUM_Good_h[jI][cI][wI][idI][rI], jtPfCHMF_Good_h[jI][cI][wI][idI][rI], jtPfCEMF_Good_h[jI][cI][wI][idI][rI], jtPfNHMF_Good_h[jI][cI][wI][idI][rI], jtPfNEMF_Good_h[jI][cI][wI][idI][rI], jtPfMUMF_Good_h[jI][cI][wI][idI][rI], jtPfCHF_Bad_h[jI][cI][wI][idI][rI], jtPfCEF_Bad_h[jI][cI][wI][idI][rI], jtPfNHF_Bad_h[jI][cI][wI][idI][rI], jtPfNEF_Bad_h[jI][cI][wI][idI][rI], jtPfMUF_Bad_h[jI][cI][wI][idI][rI], jtPfCHM_Bad_h[jI][cI][wI][idI][rI], jtPfCEM_Bad_h[jI][cI][wI][idI][rI], jtPfNHM_Bad_h[jI][cI][wI][idI][rI], jtPfNEM_Bad_h[jI][cI][wI][idI][rI], jtPfMUM_Bad_h[jI][cI][wI][idI][rI], jtPfCHMF_Bad_h[jI][cI][wI][idI][rI], jtPfCEMF_Bad_h[jI][cI][wI][idI][rI], jtPfNHMF_Bad_h[jI][cI][wI][idI][rI], jtPfNEMF_Bad_h[jI][cI][wI][idI][rI], jtPfMUMF_Bad_h[jI][cI][wI][idI][rI]};

	    centerTitles(tempVect);
	    setSumW2(tempVect);
	  }
	  
	  jtPfCHF_jtPfCEF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCHF_jtPfCEF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHF;jtPfCEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHF_jtPfNHF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCHF_jtPfNHF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHF;jtPfNHF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHF_jtPfNEF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCHF_jtPfNEF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHF_jtPfMUF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCHF_jtPfMUF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  jtPfCEF_jtPfNHF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCEF_jtPfNHF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEF;jtPfNHF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCEF_jtPfNEF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCEF_jtPfNEF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCEF_jtPfMUF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCEF_jtPfMUF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  jtPfNHF_jtPfNEF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfNHF_jtPfNEF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNHF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfNHF_jtPfMUF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfNHF_jtPfMUF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNHF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  jtPfNEF_jtPfMUF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfNEF_jtPfMUF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNEF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  
	  jtPfCHMF_jtPfCEMF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCHMF_jtPfCEMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHMF;jtPfCEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHMF_jtPfNHMF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCHMF_jtPfNHMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHMF;jtPfNHMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHMF_jtPfNEMF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCHMF_jtPfNEMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHMF_jtPfMUMF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCHMF_jtPfMUMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  jtPfCEMF_jtPfNHMF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCEMF_jtPfNHMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEMF;jtPfNHMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCEMF_jtPfNEMF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCEMF_jtPfNEMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCEMF_jtPfMUMF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfCEMF_jtPfMUMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  jtPfNHMF_jtPfNEMF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfNHMF_jtPfNEMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNHMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfNHMF_jtPfMUMF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfNHMF_jtPfMUMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNHMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  jtPfNEMF_jtPfMUMF_Good_h[jI][cI][wI][idI] = new TH2D(("jtPfNEMF_jtPfMUMF_Good_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNEMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  
	  
	  jtPfCHF_jtPfCEF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCHF_jtPfCEF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHF;jtPfCEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHF_jtPfNHF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCHF_jtPfNHF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHF;jtPfNHF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHF_jtPfNEF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCHF_jtPfNEF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHF_jtPfMUF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCHF_jtPfMUF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  jtPfCEF_jtPfNHF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCEF_jtPfNHF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEF;jtPfNHF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCEF_jtPfNEF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCEF_jtPfNEF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCEF_jtPfMUF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCEF_jtPfMUF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	
	  jtPfNHF_jtPfNEF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfNHF_jtPfNEF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNHF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfNHF_jtPfMUF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfNHF_jtPfMUF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNHF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  jtPfNEF_jtPfMUF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfNEF_jtPfMUF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNEF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  	  
	  jtPfCHMF_jtPfCEMF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCHMF_jtPfCEMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHMF;jtPfCEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHMF_jtPfNHMF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCHMF_jtPfNHMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHMF;jtPfNHMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHMF_jtPfNEMF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCHMF_jtPfNEMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCHMF_jtPfMUMF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCHMF_jtPfMUMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCHMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  jtPfCEMF_jtPfNHMF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCEMF_jtPfNHMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEMF;jtPfNHMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCEMF_jtPfNEMF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCEMF_jtPfNEMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfCEMF_jtPfMUMF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfCEMF_jtPfMUMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfCEMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  jtPfNHMF_jtPfNEMF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfNHMF_jtPfNEMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNHMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  jtPfNHMF_jtPfMUMF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfNHMF_jtPfMUMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNHMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
	  
	  jtPfNEMF_jtPfMUMF_Bad_h[jI][cI][wI][idI] = new TH2D(("jtPfNEMF_jtPfMUMF_Bad_" + weightStr[wI] + "_" + jetAlgos.at(jI) + "_" + centStr + "_" + idStr[idI] + "_h").c_str(), ";jtPfNEMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);   
		
	  std::vector<TH1*> tempVect = {jtPt_Good_h[jI][cI][wI][idI], jtPt_Bad_h[jI][cI][wI][idI], jtPt_FineTot_h[jI][cI][wI][idI], jtPfCHF_jtPfCEF_Good_h[jI][cI][wI][idI], jtPfCHF_jtPfNHF_Good_h[jI][cI][wI][idI], jtPfCHF_jtPfNEF_Good_h[jI][cI][wI][idI], jtPfCHF_jtPfMUF_Good_h[jI][cI][wI][idI], jtPfCEF_jtPfNHF_Good_h[jI][cI][wI][idI], jtPfCEF_jtPfNEF_Good_h[jI][cI][wI][idI], jtPfCEF_jtPfMUF_Good_h[jI][cI][wI][idI], jtPfNHF_jtPfNEF_Good_h[jI][cI][wI][idI], jtPfNHF_jtPfMUF_Good_h[jI][cI][wI][idI], jtPfNEF_jtPfMUF_Good_h[jI][cI][wI][idI], jtPfCHF_jtPfCEF_Bad_h[jI][cI][wI][idI], jtPfCHF_jtPfNHF_Bad_h[jI][cI][wI][idI], jtPfCHF_jtPfNEF_Bad_h[jI][cI][wI][idI], jtPfCHF_jtPfMUF_Bad_h[jI][cI][wI][idI], jtPfCEF_jtPfNHF_Bad_h[jI][cI][wI][idI], jtPfCEF_jtPfNEF_Bad_h[jI][cI][wI][idI], jtPfCEF_jtPfMUF_Bad_h[jI][cI][wI][idI], jtPfNHF_jtPfNEF_Bad_h[jI][cI][wI][idI], jtPfNHF_jtPfMUF_Bad_h[jI][cI][wI][idI], jtPfNEF_jtPfMUF_Bad_h[jI][cI][wI][idI], jtPfCHMF_jtPfCEMF_Good_h[jI][cI][wI][idI], jtPfCHMF_jtPfNHMF_Good_h[jI][cI][wI][idI], jtPfCHMF_jtPfNEMF_Good_h[jI][cI][wI][idI], jtPfCHMF_jtPfMUMF_Good_h[jI][cI][wI][idI], jtPfCEMF_jtPfNHMF_Good_h[jI][cI][wI][idI], jtPfCEMF_jtPfNEMF_Good_h[jI][cI][wI][idI], jtPfCEMF_jtPfMUMF_Good_h[jI][cI][wI][idI], jtPfNHMF_jtPfNEMF_Good_h[jI][cI][wI][idI], jtPfNHMF_jtPfMUMF_Good_h[jI][cI][wI][idI], jtPfNEMF_jtPfMUMF_Good_h[jI][cI][wI][idI], jtPfCHMF_jtPfCEMF_Bad_h[jI][cI][wI][idI], jtPfCHMF_jtPfNHMF_Bad_h[jI][cI][wI][idI], jtPfCHMF_jtPfNEMF_Bad_h[jI][cI][wI][idI], jtPfCHMF_jtPfMUMF_Bad_h[jI][cI][wI][idI], jtPfCEMF_jtPfNHMF_Bad_h[jI][cI][wI][idI], jtPfCEMF_jtPfNEMF_Bad_h[jI][cI][wI][idI], jtPfCEMF_jtPfMUMF_Bad_h[jI][cI][wI][idI], jtPfNHMF_jtPfNEMF_Bad_h[jI][cI][wI][idI], jtPfNHMF_jtPfMUMF_Bad_h[jI][cI][wI][idI], jtPfNEMF_jtPfMUMF_Bad_h[jI][cI][wI][idI]};

	  centerTitles(tempVect);     
	  setSumW2(tempVect);
	}
      }
    }
  }

  const Int_t nMaxJet_ = 500;
  Float_t pthat_[nJetAlgos];
  Int_t nref_[nJetAlgos];
  Float_t jtpt_[nJetAlgos][nMaxJet_];
  Float_t jteta_[nJetAlgos][nMaxJet_];
  Float_t jtphi_[nJetAlgos][nMaxJet_];
  Float_t refpt_[nJetAlgos][nMaxJet_];
  Float_t refphi_[nJetAlgos][nMaxJet_];
  Float_t refeta_[nJetAlgos][nMaxJet_];
  Float_t jtPfCHF_[nJetAlgos][nMaxJet_];
  Float_t jtPfCEF_[nJetAlgos][nMaxJet_];
  Float_t jtPfNHF_[nJetAlgos][nMaxJet_];
  Float_t jtPfNEF_[nJetAlgos][nMaxJet_];
  Float_t jtPfMUF_[nJetAlgos][nMaxJet_];

  Float_t jtPfCHMF_[nJetAlgos][nMaxJet_];
  Float_t jtPfCEMF_[nJetAlgos][nMaxJet_];
  Float_t jtPfNHMF_[nJetAlgos][nMaxJet_];
  Float_t jtPfNEMF_[nJetAlgos][nMaxJet_];
  Float_t jtPfMUMF_[nJetAlgos][nMaxJet_];

  Int_t jtPfCHM_[nJetAlgos][nMaxJet_];
  Int_t jtPfCEM_[nJetAlgos][nMaxJet_];
  Int_t jtPfNHM_[nJetAlgos][nMaxJet_];
  Int_t jtPfNEM_[nJetAlgos][nMaxJet_];
  Int_t jtPfMUM_[nJetAlgos][nMaxJet_];

  Int_t ngen_[nJetAlgos];
  Float_t genpt_[nJetAlgos][nMaxJet_];
  Float_t geneta_[nJetAlgos][nMaxJet_];
  Float_t genphi_[nJetAlgos][nMaxJet_];
  Int_t gensubid_[nJetAlgos][nMaxJet_];

  unsigned int run_, lumi_;
  unsigned long long evt_;
  Int_t hiBin_;
  Float_t vz_;
  Float_t hiHF_;

  Int_t HBHENoiseFilterResultRun2Loose_ = -1;
  Int_t pprimaryVertexFilter_ = -1;
  Int_t pBeamScrapingFilter_ = -1;
  Int_t phfCoincFilter3_ = -1;
  Int_t pclusterCompatibilityFilter_ = -1;

  goodGlobalSelection globalSel;
  globalSel.setIsPbPb(!isPP);

  specialHYDJETEventExclude specialSel;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "Processing file " << fI << "/" << fileList.size() << ": \'" << fileList.at(fI) << "\'" << std::endl;

    inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(fI)).c_str(), "READ");
    TTree* jetTree_p[nJetAlgos];

    for(Int_t jI = 0; jI < nJetAlgos; ++jI){
      jetTree_p[jI] = (TTree*)inFile_p->Get(jetList.at(jI).c_str());

      jetTree_p[jI]->SetBranchStatus("*", 0);
      jetTree_p[jI]->SetBranchStatus("pthat", 1);
      jetTree_p[jI]->SetBranchStatus("nref", 1);
      jetTree_p[jI]->SetBranchStatus("jtpt", 1);
      jetTree_p[jI]->SetBranchStatus("jteta", 1);
      jetTree_p[jI]->SetBranchStatus("jtphi", 1);
      jetTree_p[jI]->SetBranchStatus("refpt", 1);
      jetTree_p[jI]->SetBranchStatus("refeta", 1);
      jetTree_p[jI]->SetBranchStatus("refphi", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfCHF", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfCEF", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfNHF", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfNEF", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfMUF", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfCHMF", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfCEMF", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfNHMF", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfNEMF", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfMUMF", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfCHM", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfCEM", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfNHM", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfNEM", 1);
      jetTree_p[jI]->SetBranchStatus("jtPfMUM", 1);
      jetTree_p[jI]->SetBranchStatus("ngen", 1);
      jetTree_p[jI]->SetBranchStatus("genpt", 1);
      jetTree_p[jI]->SetBranchStatus("geneta", 1);
      jetTree_p[jI]->SetBranchStatus("genphi", 1);
      jetTree_p[jI]->SetBranchStatus("gensubid", 1);
      
      jetTree_p[jI]->SetBranchAddress("pthat", &(pthat_[jI]));
      jetTree_p[jI]->SetBranchAddress("nref", &(nref_[jI]));
      jetTree_p[jI]->SetBranchAddress("jtpt", jtpt_[jI]);
      jetTree_p[jI]->SetBranchAddress("jteta", jteta_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtphi", jtphi_[jI]);
      jetTree_p[jI]->SetBranchAddress("refpt", refpt_[jI]);
      jetTree_p[jI]->SetBranchAddress("refeta", refeta_[jI]);
      jetTree_p[jI]->SetBranchAddress("refphi", refphi_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfCHF", jtPfCHF_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfCEF", jtPfCEF_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfNHF", jtPfNHF_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfNEF", jtPfNEF_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfMUF", jtPfMUF_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfCHMF", jtPfCHMF_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfCEMF", jtPfCEMF_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfNHMF", jtPfNHMF_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfNEMF", jtPfNEMF_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfMUMF", jtPfMUMF_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfCHM", jtPfCHM_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfCEM", jtPfCEM_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfNHM", jtPfNHM_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfNEM", jtPfNEM_[jI]);
      jetTree_p[jI]->SetBranchAddress("jtPfMUM", jtPfMUM_[jI]);
      jetTree_p[jI]->SetBranchAddress("ngen", &(ngen_[jI]));
      jetTree_p[jI]->SetBranchAddress("genpt", genpt_[jI]);
      jetTree_p[jI]->SetBranchAddress("geneta", geneta_[jI]);
      jetTree_p[jI]->SetBranchAddress("genphi", genphi_[jI]);
      jetTree_p[jI]->SetBranchAddress("gensubid", gensubid_[jI]);
    }


    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");

    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("hiHF", 1);
    hiTree_p->SetBranchStatus("vz", 1);

    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("evt", &evt_);
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("hiHF", &hiHF_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    
    TTree* skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");

    skimTree_p->SetBranchStatus("*", 0);
    skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);

    if(!isPP){
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

    const Int_t nEntries = TMath::Min(1000000000, (Int_t)jetTree_p[0]->GetEntries());
    const Int_t printInterval = TMath::Max(1, nEntries/20);

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(nEntries >= 50000 && entry%printInterval == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;

      for(Int_t jI = 0; jI < nJetAlgos; ++jI){
	jetTree_p[jI]->GetEntry(entry);

	const Int_t tempNMatch = ngen_[jI];
	Bool_t isMatched[tempNMatch];

	for(Int_t gI = 0; gI < tempNMatch; ++gI){
	  isMatched[gI] = false;

	  for(Int_t gI2 = gI+1; gI2 < tempNMatch; ++gI2){
	    if(genpt_[jI][gI] < genpt_[jI][gI2]){
	      Float_t tempGenPt = genpt_[jI][gI];
	      Float_t tempGenEta = geneta_[jI][gI];
	      Float_t tempGenPhi = genphi_[jI][gI];
	      Int_t tempGenSubId = gensubid_[jI][gI];

	      genpt_[jI][gI] = genpt_[jI][gI2];
	      genphi_[jI][gI] = genphi_[jI][gI2];
	      geneta_[jI][gI] = geneta_[jI][gI2];
	      gensubid_[jI][gI] = gensubid_[jI][gI2];

	      genpt_[jI][gI2] = tempGenPt;
	      genphi_[jI][gI2] = tempGenPhi;
	      geneta_[jI][gI2] = tempGenEta;
	      gensubid_[jI][gI2] = tempGenSubId;
	    }
	  }
	}	

	for(Int_t rI = 0; rI < nref_[jI]; ++rI){
	  if(refpt_[jI][rI] < 0 || (jtpt_[jI][rI] > 200 && refpt_[jI][rI]/jtpt_[jI][rI] < 0.10)){
	    
	    for(Int_t gI = 0; gI < tempNMatch; ++gI){
	      if(isMatched[gI]) continue;
	      if(gensubid_[gI] == 0) continue;
	      if(getDR(geneta_[jI][gI], genphi_[jI][gI], jteta_[jI][rI], jtphi_[jI][rI]) < 0.2){
		refpt_[jI][rI] = genpt_[jI][gI];
		refphi_[jI][rI] = genphi_[jI][gI];
		refeta_[jI][rI] = geneta_[jI][gI];
		isMatched[gI] = true;
		break;
	      }
	    }	    
	  }
	}
      }

      hiTree_p->GetEntry(entry);
      skimTree_p->GetEntry(entry);

      globalSel.setVz(vz_);
      globalSel.setHiHF(hiHF_);
      globalSel.setPprimaryVertexFilter(pprimaryVertexFilter_);
      globalSel.setPBeamScrapingFilter(pBeamScrapingFilter_);
      globalSel.setPhfCoincFilter3(phfCoincFilter3_);
      globalSel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
      globalSel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);

      if(!globalSel.isGood()) continue;

      Int_t centPos = -1;
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(hiBin_/2 >= centBinsLow.at(cI) && hiBin_/2 < centBinsHi.at(cI)){
	  centPos = cI;
	  break;
	}
      }

      if(centPos < 0) continue;
    
      bool badJetSpecialSel = specialSel.CheckEventBadJet(ngen_[posR4], genpt_[posR4], genphi_[posR4], geneta_[posR4], gensubid_[posR4], entry);
    
      Double_t pthatWeight_ = -1;
      for(unsigned int pI = 0; pI < pthats.size()-1; ++pI){
        if(pthats.at(pI) <= pthat_[0] && pthat_[0] < pthats.at(pI+1)){
          pthatWeight_ = pthatWeights.at(pI);
          break;
        }
      }
      if(pthat_[0] > pthats.at(pthats.size()-1)) pthatWeight_ = pthatWeights.at(pthatWeights.size()-1);

      if(pthatWeight_ < 0){
	std::cout << "WARNING - NO WEIGHT FOR pthat \'" << pthat_[0] << "\'. Set to 1" << std::endl;
        pthatWeight_ = 1.;
      }

      Double_t fullWeight_ = pthatWeight_;
      if(!isPP){
	Double_t ncollWeight_ = findNcoll_Renorm(hiBin_);
	fullWeight_ *= ncollWeight_;
      }

      for(Int_t tI = 0; tI < nJetAlgos; ++tI){
	for(Int_t jI = 0; jI < nref_[tI]; ++jI){
	  Bool_t passesID[nID];

	  std::vector<Int_t> recoPtPos;

	  for(Int_t rI = 0; rI < nJtRecoPtBins2; ++rI){
	    if(jtpt_[tI][jI] >= jtRecoPtBins2Low[rI] && jtpt_[tI][jI] < jtRecoPtBins2Hi[rI]){
	      recoPtPos.push_back(rI);
	    }
	  }

	  for(Int_t idI = 0; idI < nID; ++idI){
	    passesID[idI] = false;
	  }

	  for(Int_t idI = 0; idI < nID; ++idI){
 	    if(jtPfCHMF_[tI][jI] < jtPfCHMFCutLow[idI]) continue;
 	    if(jtPfCHMF_[tI][jI] > jtPfCHMFCutHi[idI]) continue;
 	    if(jtPfMUMF_[tI][jI] < jtPfMUMFCutLow[idI]) continue;
 	    if(jtPfMUMF_[tI][jI] > jtPfMUMFCutHi[idI]) continue;
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

	    if(jtPfCEM_[tI][jI] + jtPfNEM_[tI][jI] + jtPfCHM_[tI][jI] + jtPfNHM_[tI][jI] + jtPfMUM_[tI][jI] < jtPfMinMult[idI]) continue;
	    if(jtPfCHM_[tI][jI] < jtPfMinChgMult[idI]) continue;

	    passesID[idI] = true;
	  }

	  if(TMath::Abs(jteta_[tI][jI]) > jtAbsEtaMax) continue;
	  if(jtpt_[tI][jI] < jtRecoPtMin) continue;
	  
	  //	  Double_t res = getResForPtEtaCentAlgo(jtpt_[tI][jI], jteta_[tI][jI], hiBin_/2, "akCs4PU3PFFlow");	
	  Double_t res = getResForPtEtaCentAlgo(jtpt_[tI][jI], jteta_[tI][jI], hiBin_/2, jetAlgos.at(tI));	

	  bool isGood = refpt_[tI][jI] > TMath::Max(0., jtpt_[tI][jI] - jtpt_[tI][jI]*3.*res);
	  bool isBad = refpt_[tI][jI] < TMath::Max(0., jtpt_[tI][jI] - jtpt_[tI][jI]*5.*res);
	  	  	  
	  if(!isPP && badJetSpecialSel) continue;
	  
	  if(isGood){
	    for(Int_t wI = 0; wI < nWeighted; ++wI){
	      Double_t weight_ = 1.;
	      if(wI == 1) weight_ = fullWeight_;
	      
	      for(Int_t idI = 0; idI < nID; ++idI){
		if(passesID[idI]){
		  jtPt_Good_h[tI][centPos][wI][idI]->Fill(jtpt_[tI][jI], weight_);
		  jtPt_FineTot_h[tI][centPos][wI][idI]->Fill(jtpt_[tI][jI], weight_);

				
		  if(jtpt_[tI][jI] >= 247.5 && jtpt_[tI][jI] <= 262.5 && pthatWeight_ > 0.01){
		    if(tI == posR4 && idI == 0 && wI == 0) std::cout << "Bad jet weight between  247.5-262.5, entry: " << entry << std::endl;
		  }
		  

		  for(unsigned int rI = 0; rI < recoPtPos.size(); ++rI){
		    jtPfCHF_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCHF_[tI][jI], weight_);
		    jtPfCEF_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCEF_[tI][jI], weight_);
		    jtPfNHF_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNHF_[tI][jI], weight_);
		    jtPfNEF_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNEF_[tI][jI], weight_);
		    jtPfMUF_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfMUF_[tI][jI], weight_);
		    
		    jtPfCHM_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCHM_[tI][jI], weight_);
		    jtPfCEM_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCEM_[tI][jI], weight_);
		    jtPfNHM_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNHM_[tI][jI], weight_);
		    jtPfNEM_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNEM_[tI][jI], weight_);
		    jtPfMUM_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfMUM_[tI][jI], weight_);

		    jtPfCHMF_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCHMF_[tI][jI], weight_);
		    jtPfCEMF_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCEMF_[tI][jI], weight_);
		    jtPfNHMF_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNHMF_[tI][jI], weight_);
		    jtPfNEMF_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNEMF_[tI][jI], weight_);
		    jtPfMUMF_Good_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfMUMF_[tI][jI], weight_);
		  }

		  jtPfCHF_jtPfCEF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCHF_[tI][jI], jtPfCEF_[tI][jI], weight_);
		  jtPfCHF_jtPfNHF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCHF_[tI][jI], jtPfNHF_[tI][jI], weight_);
		  jtPfCHF_jtPfNEF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCHF_[tI][jI], jtPfNEF_[tI][jI], weight_);
		  jtPfCHF_jtPfMUF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCHF_[tI][jI], jtPfMUF_[tI][jI], weight_);
		  
		  jtPfCEF_jtPfNHF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCEF_[tI][jI], jtPfNHF_[tI][jI], weight_);
		  jtPfCEF_jtPfNEF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCEF_[tI][jI], jtPfNEF_[tI][jI], weight_);
		  jtPfCEF_jtPfMUF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCEF_[tI][jI], jtPfMUF_[tI][jI], weight_);
		  
		  jtPfNHF_jtPfNEF_Good_h[tI][centPos][wI][idI]->Fill(jtPfNHF_[tI][jI], jtPfNEF_[tI][jI], weight_);
		  jtPfNHF_jtPfMUF_Good_h[tI][centPos][wI][idI]->Fill(jtPfNHF_[tI][jI], jtPfMUF_[tI][jI], weight_);
		  
		  jtPfNEF_jtPfMUF_Good_h[tI][centPos][wI][idI]->Fill(jtPfNEF_[tI][jI], jtPfMUF_[tI][jI], weight_);
		  
		  
		  jtPfCHMF_jtPfCEMF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCHMF_[tI][jI], jtPfCEMF_[tI][jI], weight_);
		  jtPfCHMF_jtPfNHMF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCHMF_[tI][jI], jtPfNHMF_[tI][jI], weight_);
		  jtPfCHMF_jtPfNEMF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCHMF_[tI][jI], jtPfNEMF_[tI][jI], weight_);
		  jtPfCHMF_jtPfMUMF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCHMF_[tI][jI], jtPfMUMF_[tI][jI], weight_);
		  
		  jtPfCEMF_jtPfNHMF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCEMF_[tI][jI], jtPfNHMF_[tI][jI], weight_);
		  jtPfCEMF_jtPfNEMF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCEMF_[tI][jI], jtPfNEMF_[tI][jI], weight_);
		  jtPfCEMF_jtPfMUMF_Good_h[tI][centPos][wI][idI]->Fill(jtPfCEMF_[tI][jI], jtPfMUMF_[tI][jI], weight_);
		  
		  jtPfNHMF_jtPfNEMF_Good_h[tI][centPos][wI][idI]->Fill(jtPfNHMF_[tI][jI], jtPfNEMF_[tI][jI], weight_);
		  jtPfNHMF_jtPfMUMF_Good_h[tI][centPos][wI][idI]->Fill(jtPfNHMF_[tI][jI], jtPfMUMF_[tI][jI], weight_);
		  
		  jtPfNEMF_jtPfMUMF_Good_h[tI][centPos][wI][idI]->Fill(jtPfNEMF_[tI][jI], jtPfMUMF_[tI][jI], weight_);
		}
	      }
	    }
	  }	 
	  else if(isBad){
	    for(Int_t wI = 0; wI < nWeighted; ++wI){
	      Double_t weight_ = 1.;
	      if(wI == 1) weight_ = fullWeight_;
	      
	      for(Int_t idI = 0; idI < nID; ++idI){
		if(passesID[idI]){
		  jtPt_Bad_h[tI][centPos][wI][idI]->Fill(jtpt_[tI][jI], weight_);
		  jtPt_FineTot_h[tI][centPos][wI][idI]->Fill(jtpt_[tI][jI], weight_);

		  
		  if(jtpt_[tI][jI] >= 247.5 && jtpt_[tI][jI] <= 262.5 && pthatWeight_ > 0.01){
		    if(tI == posR4 && idI == 0 && wI == 0) std::cout << "Bad jet weight between  247.5-262.5, entry: " << entry << std::endl;
		  }
		  

		  if(centPos == 2 && jtpt_[tI][jI] > 300 && jtpt_[tI][jI] < 400 && idI == nID - 1 && wI == 1){
		    //		    std::cout << "Warning bad jet: " << jtpt_[tI][jI] << ", weight==" << weight_ << ", entry==" << entry << std::endl;
		    //		    std::cout << " " << refpt_[tI][jI] << "= refpt" << std::endl;
		  }

		  for(unsigned int rI = 0; rI < recoPtPos.size(); ++rI){
		    jtPfCHF_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCHF_[tI][jI], weight_);
		    jtPfCEF_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCEF_[tI][jI], weight_);
		    jtPfNHF_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNHF_[tI][jI], weight_);
		    jtPfNEF_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNEF_[tI][jI], weight_);
		    jtPfMUF_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfMUF_[tI][jI], weight_);
		    
		    jtPfCHM_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCHM_[tI][jI], weight_);
		    jtPfCEM_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCEM_[tI][jI], weight_);
		    jtPfNHM_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNHM_[tI][jI], weight_);
		    jtPfNEM_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNEM_[tI][jI], weight_);
		    jtPfMUM_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfMUM_[tI][jI], weight_);
		    
		    jtPfCHMF_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCHMF_[tI][jI], weight_);
		    jtPfCEMF_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfCEMF_[tI][jI], weight_);
		    jtPfNHMF_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNHMF_[tI][jI], weight_);
		    jtPfNEMF_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfNEMF_[tI][jI], weight_);
		    jtPfMUMF_Bad_h[tI][centPos][wI][idI][recoPtPos.at(rI)]->Fill(jtPfMUMF_[tI][jI], weight_);
		  }
		  
		  jtPfCHF_jtPfCEF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCHF_[tI][jI], jtPfCEF_[tI][jI], weight_);
		  jtPfCHF_jtPfNHF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCHF_[tI][jI], jtPfNHF_[tI][jI], weight_);
		  jtPfCHF_jtPfNEF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCHF_[tI][jI], jtPfNEF_[tI][jI], weight_);
		  jtPfCHF_jtPfMUF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCHF_[tI][jI], jtPfMUF_[tI][jI], weight_);
		  
		  jtPfCEF_jtPfNHF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCEF_[tI][jI], jtPfNHF_[tI][jI], weight_);
		  jtPfCEF_jtPfNEF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCEF_[tI][jI], jtPfNEF_[tI][jI], weight_);
		  jtPfCEF_jtPfMUF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCEF_[tI][jI], jtPfMUF_[tI][jI], weight_);
		  
		  jtPfNHF_jtPfNEF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfNHF_[tI][jI], jtPfNEF_[tI][jI], weight_);
		  jtPfNHF_jtPfMUF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfNHF_[tI][jI], jtPfMUF_[tI][jI], weight_);
		  
		  jtPfNEF_jtPfMUF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfNEF_[tI][jI], jtPfMUF_[tI][jI], weight_);	  
		  
		  
		  jtPfCHMF_jtPfCEMF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCHMF_[tI][jI], jtPfCEMF_[tI][jI], weight_);
		  jtPfCHMF_jtPfNHMF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCHMF_[tI][jI], jtPfNHMF_[tI][jI], weight_);
		  jtPfCHMF_jtPfNEMF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCHMF_[tI][jI], jtPfNEMF_[tI][jI], weight_);
		  jtPfCHMF_jtPfMUMF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCHMF_[tI][jI], jtPfMUMF_[tI][jI], weight_);
		  
		  jtPfCEMF_jtPfNHMF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCEMF_[tI][jI], jtPfNHMF_[tI][jI], weight_);
		  jtPfCEMF_jtPfNEMF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCEMF_[tI][jI], jtPfNEMF_[tI][jI], weight_);
		  jtPfCEMF_jtPfMUMF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfCEMF_[tI][jI], jtPfMUMF_[tI][jI], weight_);
		  
		  jtPfNHMF_jtPfNEMF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfNHMF_[tI][jI], jtPfNEMF_[tI][jI], weight_);
		  jtPfNHMF_jtPfMUMF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfNHMF_[tI][jI], jtPfMUMF_[tI][jI], weight_);
		  
		  jtPfNEMF_jtPfMUMF_Bad_h[tI][centPos][wI][idI]->Fill(jtPfNEMF_[tI][jI], jtPfMUMF_[tI][jI], weight_);	  
		}
	      }
	    }
	  }
	}
      }
    }

    inFile_p->Close();
    delete inFile_p;
  }

 
  std::vector<TH1*> goodHists;
  std::vector<TH1*> badHists;
  std::vector<TH1*> fineHists;

  std::vector<TH2*> goodHistsTH2;
  std::vector<TH2*> badHistsTH2;

  goodHists.reserve(41*nCentBins);
  badHists.reserve(41*nCentBins);

  goodHistsTH2.reserve(41*nCentBins);
  badHistsTH2.reserve(41*nCentBins);

  for(Int_t tI = 0; tI < nJetAlgos; ++tI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      for(Int_t wI = 0; wI < nWeighted; ++wI){
	for(Int_t idI = 0; idI < nID; ++idI){

	  goodHists.push_back(jtPt_Good_h[tI][cI][wI][idI]);
	  fineHists.push_back(jtPt_FineTot_h[tI][cI][wI][idI]);
	  
	  for(Int_t rI = 0; rI < nJtRecoPtBins2; ++rI){
	    goodHists.push_back(jtPfCHF_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfCEF_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfNHF_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfNEF_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfMUF_Good_h[tI][cI][wI][idI][rI]);
	    
	    goodHists.push_back(jtPfCHM_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfCEM_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfNHM_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfNEM_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfMUM_Good_h[tI][cI][wI][idI][rI]);

	    goodHists.push_back(jtPfCHMF_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfCEMF_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfNHMF_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfNEMF_Good_h[tI][cI][wI][idI][rI]);
	    goodHists.push_back(jtPfMUMF_Good_h[tI][cI][wI][idI][rI]);
	  }

	  goodHistsTH2.push_back(jtPfCHF_jtPfCEF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfCHF_jtPfNHF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfCHF_jtPfNEF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfCHF_jtPfMUF_Good_h[tI][cI][wI][idI]);
	  
	  goodHistsTH2.push_back(jtPfCEF_jtPfNHF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfCEF_jtPfNEF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfCEF_jtPfMUF_Good_h[tI][cI][wI][idI]);
      
	  goodHistsTH2.push_back(jtPfNHF_jtPfNEF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfNHF_jtPfMUF_Good_h[tI][cI][wI][idI]);
	  
	  goodHistsTH2.push_back(jtPfNEF_jtPfMUF_Good_h[tI][cI][wI][idI]);
	  
	  
	  goodHistsTH2.push_back(jtPfCHMF_jtPfCEMF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfCHMF_jtPfNHMF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfCHMF_jtPfNEMF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfCHMF_jtPfMUMF_Good_h[tI][cI][wI][idI]);
	  
	  goodHistsTH2.push_back(jtPfCEMF_jtPfNHMF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfCEMF_jtPfNEMF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfCEMF_jtPfMUMF_Good_h[tI][cI][wI][idI]);
	  
	  goodHistsTH2.push_back(jtPfNHMF_jtPfNEMF_Good_h[tI][cI][wI][idI]);
	  goodHistsTH2.push_back(jtPfNHMF_jtPfMUMF_Good_h[tI][cI][wI][idI]);
	  
	  goodHistsTH2.push_back(jtPfNEMF_jtPfMUMF_Good_h[tI][cI][wI][idI]);
	  
	  badHists.push_back(jtPt_Bad_h[tI][cI][wI][idI]);
	
          for(Int_t rI = 0; rI < nJtRecoPtBins2; ++rI){
	    badHists.push_back(jtPfCHF_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfCEF_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfNHF_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfNEF_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfMUF_Bad_h[tI][cI][wI][idI][rI]);
	    
	    badHists.push_back(jtPfCHM_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfCEM_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfNHM_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfNEM_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfMUM_Bad_h[tI][cI][wI][idI][rI]);

	    badHists.push_back(jtPfCHMF_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfCEMF_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfNHMF_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfNEMF_Bad_h[tI][cI][wI][idI][rI]);
	    badHists.push_back(jtPfMUMF_Bad_h[tI][cI][wI][idI][rI]);
	  }

	  badHistsTH2.push_back(jtPfCHF_jtPfCEF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfCHF_jtPfNHF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfCHF_jtPfNEF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfCHF_jtPfMUF_Bad_h[tI][cI][wI][idI]);
	
	  badHistsTH2.push_back(jtPfCEF_jtPfNHF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfCEF_jtPfNEF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfCEF_jtPfMUF_Bad_h[tI][cI][wI][idI]);
	  
	  badHistsTH2.push_back(jtPfNHF_jtPfNEF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfNHF_jtPfMUF_Bad_h[tI][cI][wI][idI]);
	  
	  badHistsTH2.push_back(jtPfNEF_jtPfMUF_Bad_h[tI][cI][wI][idI]);
	
	  badHistsTH2.push_back(jtPfCHMF_jtPfCEMF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfCHMF_jtPfNHMF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfCHMF_jtPfNEMF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfCHMF_jtPfMUMF_Bad_h[tI][cI][wI][idI]);
	  
	  badHistsTH2.push_back(jtPfCEMF_jtPfNHMF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfCEMF_jtPfNEMF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfCEMF_jtPfMUMF_Bad_h[tI][cI][wI][idI]);
	  
	  badHistsTH2.push_back(jtPfNHMF_jtPfNEMF_Bad_h[tI][cI][wI][idI]);
	  badHistsTH2.push_back(jtPfNHMF_jtPfMUMF_Bad_h[tI][cI][wI][idI]);
	  
	  badHistsTH2.push_back(jtPfNEMF_jtPfMUMF_Bad_h[tI][cI][wI][idI]);
	}
      }
    }
  }

  //std::cout << __LINE__ << std::endl;
  for(unsigned int gI = 0; gI < goodHists.size(); ++gI){
    std::cout << " Good " << gI << "/" << goodHists.size() << ": " << goodHists.at(gI) << ", " << goodHists.at(gI)->GetName() << std::endl;
  }

  specialSel.PrintExcludedNumbers();

  outFile_p->cd();

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    checkMakeDir("pdfDir/" + dateStr + "/" + jetAlgos.at(jI));
  }


  vanGoghPalette vg;

  const int nPlots = 2;
  const int styles[nPlots] = {20, 21};
  const int colors[nPlots] = {vg.getColor(0), vg.getColor(1)};

  double maxValPt = -1;
  double minValPt = 100000;

  for(unsigned int i = 0; i < goodHists.size(); ++i){
    std::string canvStr = goodHists.at(i)->GetName();
    if(canvStr.find("jtPt_") == std::string::npos) continue;

    for(Int_t bIX = 0; bIX < goodHists.at(i)->GetNbinsX(); ++bIX){
      if(maxValPt < goodHists.at(i)->GetBinContent(bIX+1)) maxValPt = goodHists.at(i)->GetBinContent(bIX+1);
      if(minValPt > goodHists.at(i)->GetBinContent(bIX+1) && goodHists.at(i)->GetBinContent(bIX+1) > 0) minValPt = goodHists.at(i)->GetBinContent(bIX+1);

      if(maxValPt < badHists.at(i)->GetBinContent(bIX+1)) maxValPt = badHists.at(i)->GetBinContent(bIX+1);
      if(minValPt > badHists.at(i)->GetBinContent(bIX+1) && badHists.at(i)->GetBinContent(bIX+1) > 0) minValPt = badHists.at(i)->GetBinContent(bIX+1);
    }
  }

  maxValPt *= 5.;
  minValPt /= 5.;

  const Float_t yPadFrac = 0.35 + .125;

  //std::cout << __LINE__ << std::endl;

  TLegend* leg_p = new TLegend(0.5, 0.5, 0.9, 0.9);
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(12);

  leg_p->AddEntry(goodHists.at(0), "3#sigma Good", "P L");
  leg_p->AddEntry(badHists.at(0), "5#sigma Bad", "P L");

  for(unsigned int i = 0; i < goodHists.size(); ++i){
    std::cout << "Processing " << goodHists.at(i) << std::endl;

    std::string canvStr = goodHists.at(i)->GetName();
    canvStr.replace(canvStr.find("_Good"), std::string("_Good").size(), "");

    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(gPad->GetBottomMargin()*1.3);
    gPad->SetLeftMargin(gPad->GetBottomMargin());

    //std::cout << __LINE__ << std::endl;

    const Int_t nSpect = 4;
    TPad* spectRat[nSpect];
    for(Int_t sI = 0; sI < nSpect; ++sI){
      spectRat[sI] = NULL;
    }
    
    //std::cout << __LINE__ << std::endl;

    std::string type = canvStr.substr(0, canvStr.find("_")+1);

    if(false /*canvStr.find("jtPt_") == std::string::npos*/){
      if(goodHists.at(i)->Integral() > 0) goodHists.at(i)->Scale(1./goodHists.at(i)->Integral());
      if(badHists.at(i)->Integral() > 0) badHists.at(i)->Scale(1./badHists.at(i)->Integral());

      //std::cout << __LINE__ << std::endl;

      goodHists.at(i)->GetYaxis()->SetTitle("Counts (Norm. To 1)");
      badHists.at(i)->GetYaxis()->SetTitle("Counts (Norm. To 1)");
    }
    else{
      //std::cout << __LINE__ << std::endl;


      gPad->SetBottomMargin(0.01);
      gPad->SetLeftMargin(0.01);

      goodHists.at(i)->SetMaximum(maxValPt);
      goodHists.at(i)->SetMinimum(minValPt);

      badHists.at(i)->SetMaximum(maxValPt);
      badHists.at(i)->SetMinimum(minValPt);

      canv_p->cd();
      spectRat[0] = new TPad("pad0", "", 0.0, yPadFrac, 1.0, 1.0);
      spectRat[0]->SetBottomMargin(0.001);
      spectRat[0]->SetTopMargin(0.01);
      spectRat[0]->SetRightMargin(0.01);
      spectRat[0]->SetLeftMargin(0.10);
      spectRat[0]->Draw("SAME");

      canv_p->cd();
      spectRat[1] = new TPad("pad1", "", 0.0, yPadFrac - .125, 1.0, yPadFrac);
      //      spectRat[1]->SetBottomMargin(spectRat[1]->GetLeftMargin()*3.);
      spectRat[1]->SetBottomMargin(0.001);
      spectRat[1]->SetTopMargin(0.001);
      spectRat[1]->SetRightMargin(0.01);
      spectRat[1]->SetLeftMargin(0.10);
      spectRat[1]->Draw("SAME");

      canv_p->cd();
      spectRat[2] = new TPad("pad1", "", 0.0, yPadFrac - .25, 1.0, yPadFrac - 0.125);
      //      spectRat[1]->SetBottomMargin(spectRat[1]->GetLeftMargin()*3.);
      spectRat[2]->SetBottomMargin(0.001);
      spectRat[2]->SetTopMargin(0.001);
      spectRat[2]->SetRightMargin(0.01);
      spectRat[2]->SetLeftMargin(0.10);
      spectRat[2]->Draw("SAME");

      canv_p->cd();
      spectRat[3] = new TPad("pad1", "", 0.0, 0.0, 1.0, yPadFrac - 0.25);
      //      spectRat[1]->SetBottomMargin(spectRat[1]->GetLeftMargin()*3.);
      spectRat[3]->SetBottomMargin(0.10*1./(yPadFrac - 0.125));
      spectRat[3]->SetTopMargin(0.001);
      spectRat[3]->SetRightMargin(0.01);
      spectRat[3]->SetLeftMargin(0.10);
      spectRat[3]->Draw("SAME");

      canv_p->cd();
      spectRat[0]->cd();
    }

    //std::cout << __LINE__ << std::endl;

    if(canvStr.find("jtPt_") == std::string::npos){
      double maxVal = -1;
      double minVal = 9999;
      
      for(Int_t bIX = 0; bIX < goodHists.at(i)->GetNbinsX(); ++bIX){
	if(maxVal < goodHists.at(i)->GetBinContent(bIX+1)) maxVal = goodHists.at(i)->GetBinContent(bIX+1);
	if(maxVal < badHists.at(i)->GetBinContent(bIX+1)) maxVal = badHists.at(i)->GetBinContent(bIX+1);
	
	if(minVal > goodHists.at(i)->GetBinContent(bIX+1) && goodHists.at(i)->GetBinContent(bIX+1) > 0) minVal = goodHists.at(i)->GetBinContent(bIX+1);
	if(minVal > badHists.at(i)->GetBinContent(bIX+1) && badHists.at(i)->GetBinContent(bIX+1) > 0) minVal = badHists.at(i)->GetBinContent(bIX+1);
      }
      
      goodHists.at(i)->SetMaximum(maxVal*5.);
      goodHists.at(i)->SetMinimum(minVal/5.);
    }

    goodHists.at(i)->SetMarkerColor(colors[0]);
    goodHists.at(i)->SetMarkerStyle(styles[0]);
    goodHists.at(i)->SetMarkerSize(0.8);
    goodHists.at(i)->SetLineColor(colors[0]);

    badHists.at(i)->SetMarkerColor(colors[1]);
    badHists.at(i)->SetMarkerStyle(styles[1]);
    badHists.at(i)->SetMarkerSize(0.6);
    badHists.at(i)->SetLineColor(colors[1]);
    
    //std::cout << __LINE__ << std::endl;
    //    canv_p->cd();

    goodHists.at(i)->GetXaxis()->SetTitleFont(43);
    goodHists.at(i)->GetYaxis()->SetTitleFont(43);
    goodHists.at(i)->GetXaxis()->SetLabelFont(43);
    goodHists.at(i)->GetYaxis()->SetLabelFont(43);

    goodHists.at(i)->GetXaxis()->SetTitleSize(12);
    goodHists.at(i)->GetYaxis()->SetTitleSize(12);
    goodHists.at(i)->GetXaxis()->SetLabelSize(12);
    goodHists.at(i)->GetYaxis()->SetLabelSize(12);

    //std::cout << __LINE__ << std::endl;

    goodHists.at(i)->GetYaxis()->SetTitleOffset(1.6);
    goodHists.at(i)->GetXaxis()->SetTitleOffset(3.);
      
    goodHists.at(i)->DrawCopy("HIST E1 P");
    badHists.at(i)->DrawCopy("HIST E1 P SAME");
    gStyle->SetOptStat(0);
    leg_p->Draw("SAME");

    //std::cout << __LINE__ << std::endl;

    if(true/*canvStr.find("jtPt_") != std::string::npos*/){
      TLatex* label_p = new TLatex();
      //      label_p->SetNDC();
      label_p->SetTextFont(43);
      label_p->SetTextSize(16);
      
      std::string centStr = "PP";
      std::string centStr2 = "PP";
      int cent = 0;
      std::string algoStr = "";

      for(unsigned int jI = 0; jI < jetAlgos.size(); ++jI){
	if(canvStr.find("_" + jetAlgos.at(jI) + "_") != std::string::npos){
	  algoStr = jetAlgos.at(jI);
	  break;
	}
      }

      if(algoStr.size() == 0){
	std::cout << "WARNING NO JET ALGO CHOSEN" << std::endl;
      }

      //std::cout << __LINE__ << std::endl;

      if(!isPP){
	//std::cout << __LINE__ << std::endl;
      
	//	algoStr = "akCs4PU3PFFlow";
	centStr = canvStr.substr(canvStr.find("Cent")+4, canvStr.size());
	centStr.replace(centStr.find("_"), centStr.size(), "");
	cent = std::stoi(centStr.substr(0, centStr.find("to")));

	//std::cout << __LINE__ << std::endl;

	centStr2 = centStr;
	centStr2.replace(centStr2.find("to"), 2, "-");
	centStr2 = centStr2 + "%";

	//std::cout << __LINE__ << std::endl;
      }

      //std::cout << __LINE__ << std::endl;

      Double_t xVal1 = (jtRecoPtBins[0] + jtRecoPtBins[1])/2.;
      Double_t xVal2 = (jtRecoPtBins[nJtRecoPtBins/2 - 1] + jtRecoPtBins[nJtRecoPtBins/2])/2.;
      Double_t xVal3 = jtRecoPtBins[nJtRecoPtBins - 2];

      Double_t res1Good = TMath::Max(0., xVal1 - 3.*xVal1*getResForPtEtaCentAlgo(xVal1, 0, cent, algoStr));
      Double_t res1Bad = TMath::Max(0., xVal1 - 5.*xVal1*getResForPtEtaCentAlgo(xVal1, 0, cent, algoStr));

      Double_t res2Good = TMath::Max(0., xVal2 - 3.*xVal2*getResForPtEtaCentAlgo(xVal2, 0, cent, algoStr));
      Double_t res2Bad = TMath::Max(0., xVal2 - 5.*xVal2*getResForPtEtaCentAlgo(xVal2, 0, cent, algoStr));


      if(type.find("jtPt_") != std::string::npos){
	label_p->DrawLatex(xVal1, maxValPt/4, ("At: " + prettyString(xVal1, 1, false)).c_str());
	label_p->DrawLatex(xVal1, maxValPt/16, ("3#sigma Good: " + prettyString(res1Good, 1, false)).c_str());
	label_p->DrawLatex(xVal1, maxValPt/(16*4), ("5#sigma Bad: " + prettyString(res1Bad, 1, false)).c_str());
	
	label_p->DrawLatex(xVal2, maxValPt/4, ("At: " + prettyString(xVal2, 1, false)).c_str());
	label_p->DrawLatex(xVal2, maxValPt/16, ("3#sigma Good: " + prettyString(res2Good, 1, false)).c_str());
	label_p->DrawLatex(xVal2, maxValPt/(16*4), ("5#sigma Bad: " + prettyString(res2Bad, 1, false)).c_str());
	label_p->DrawLatex(xVal3, maxValPt/4, centStr2.c_str());
      }
      else label_p->DrawLatex(0.8, maxValPt/4, centStr2.c_str());


      //std::cout << __LINE__ << std::endl;

      delete label_p;
    }

    gPad->SetLogy();
  
    std::string algoStr = "";
    std::string algoStr2 = "";
    for(unsigned int aI = 0; aI < jetAlgos.size(); ++aI){
      std::string test = "_" + jetAlgos.at(aI) + "_";
      if(canvStr.find(test) != std::string::npos){
	algoStr = "_" + jetAlgos.at(aI) + "_";
	algoStr2 = jetAlgos.at(aI);
	break;
      }
    }
    
    //std::cout << __LINE__ << std::endl;
    //okokok
    if(true/*canvStr.find("jtPt_") != std::string::npos*/){
      //std::cout << __LINE__ << std::endl;
      spectRat[1]->cd();

      //std::cout << __LINE__ << std::endl;

      std::string centStr = "PP";
      if(!isPP){
	centStr = canvStr.substr(canvStr.find("Cent"), canvStr.size());
	centStr.replace(centStr.find("_"), centStr.size(), "");
      }
      //std::cout << __LINE__ << std::endl;

      TH1D* tempHist_p = (TH1D*)goodHists.at(i)->Clone("temp");
      Int_t denomPos = -1;
   
      std::string recoPtStr = "";
      if(canvStr.find("RecoJtPt") != std::string::npos){
	recoPtStr = canvStr.substr(canvStr.find("RecoJtPt"), canvStr.size());
	recoPtStr.replace(recoPtStr.find("_"), recoPtStr.size(), "");
      }

      //std::cout << __LINE__ << std::endl;
    

      

      for(unsigned int iter = 0; iter < goodHists.size(); ++iter){
	std::string canvStr2 = goodHists.at(iter)->GetName();
	if(canvStr2.find(type) == std::string::npos) continue;
	if(canvStr2.find(centStr) == std::string::npos) continue;

	if(recoPtStr.size() != 0){
	  if(canvStr2.find(recoPtStr) == std::string::npos) continue;
	}

	//	std::cout << "hist \'" << canvStr << "\' gets here \'" << canvStr2 << "\'" << std::endl;            
	if(canvStr.find("Weighted") != std::string::npos){
	  if(canvStr2.find("Weighted") != std::string::npos && canvStr2.find(algoStr) != std::string::npos && canvStr2.find("NoID") != std::string::npos){
	    denomPos = iter;
	    break;
	  }
	}
	else{
	  if(canvStr2.find("Unweighted") != std::string::npos && canvStr2.find(algoStr) != std::string::npos && canvStr2.find("NoID") != std::string::npos){
	    denomPos = iter;
	    break;
	  }
	}
      }
    
      //std::cout << __LINE__ << std::endl;

      std::cout << "Found denom for hist \'" << canvStr << "\': " << denomPos << ", " << goodHists.at(denomPos)->GetName() << std::endl; 
     
      tempHist_p->GetYaxis()->SetTitle("Cut/Raw");
      tempHist_p->Divide(goodHists.at(denomPos));
      std::cout << "Divide: " << goodHists.at(i)->GetName() << "/" << goodHists.at(denomPos)->GetName() << std::endl;
      tempHist_p->SetMaximum(1.03);
      tempHist_p->SetMinimum(0.97);
      tempHist_p->GetXaxis()->SetTitleFont(43);
      tempHist_p->GetYaxis()->SetTitleFont(43);
      tempHist_p->GetXaxis()->SetLabelFont(43);
      tempHist_p->GetYaxis()->SetLabelFont(43);

      tempHist_p->GetXaxis()->SetTitleSize(12);
      tempHist_p->GetYaxis()->SetTitleSize(12);
      tempHist_p->GetXaxis()->SetLabelSize(12);
      tempHist_p->GetYaxis()->SetLabelSize(12);

      tempHist_p->GetYaxis()->SetNdivisions(505);
      tempHist_p->GetXaxis()->SetTitleOffset(3.0);
      tempHist_p->GetYaxis()->SetTitleOffset(1.6);    
      tempHist_p->DrawCopy("HIST E1 P");

      TLine* line_p = new TLine();
      line_p->SetLineStyle(2);
      line_p->DrawLine(tempHist_p->GetBinLowEdge(1), 1., tempHist_p->GetBinLowEdge(tempHist_p->GetNbinsX()), 1.);
      delete line_p;
      
      delete tempHist_p;

      spectRat[2]->cd();

      tempHist_p = (TH1D*)badHists.at(i)->Clone("temp");
      TH1D* tempHist2_p = (TH1D*)badHists.at(i)->Clone("temp2");
      denomPos = -1;
      Int_t denomPos2 = -1;
      
      for(unsigned int iter = 0; iter < badHists.size(); ++iter){
	std::string canvStr2 = badHists.at(iter)->GetName();
	if(canvStr2.find(type) == std::string::npos) continue;
	if(canvStr2.find(centStr) == std::string::npos) continue;

	if(recoPtStr.size() != 0){
	  if(canvStr2.find(recoPtStr) == std::string::npos) continue;
	}

	if(canvStr.find("Weighted") != std::string::npos){
	  if(canvStr2.find("Weighted") != std::string::npos && canvStr2.find(algoStr) != std::string::npos && canvStr2.find("NoID") != std::string::npos){
	    denomPos = iter;
	    //	    break;
	  }
	  else if(canvStr2.find("Weighted") != std::string::npos && canvStr2.find(algoStr) != std::string::npos && canvStr2.find("LightMUID") != std::string::npos){
	    denomPos2 = iter;
	  }
	}
	else{
	  if(canvStr2.find("Unweighted") != std::string::npos && canvStr2.find(algoStr) != std::string::npos && canvStr2.find("NoID") != std::string::npos){
	    denomPos = iter;
	    //	    break;
	  }
	  else if(canvStr2.find("Weighted") != std::string::npos && canvStr2.find(algoStr) != std::string::npos && canvStr2.find("LightMUID") != std::string::npos){
	    denomPos2 = iter;
	  }
	}
      }
    
      //std::cout << __LINE__ << std::endl;

      std::cout << "Found denom for hist \'" << canvStr << "\': " << denomPos << std::endl; 
     
      tempHist_p->GetYaxis()->SetTitle("Cut/Raw");
      tempHist_p->Divide(badHists.at(denomPos));
      tempHist2_p->Divide(badHists.at(denomPos2));
      std::cout << "Divide: " << badHists.at(i)->GetName() << "/" << badHists.at(denomPos)->GetName() << " (" << badHists.at(denomPos2)->GetName() << "), " << recoPtStr << std::endl;
      tempHist_p->SetMaximum(1.05);
      tempHist_p->SetMinimum(0.00);
      tempHist_p->GetXaxis()->SetTitleFont(43);
      tempHist_p->GetYaxis()->SetTitleFont(43);
      tempHist_p->GetXaxis()->SetLabelFont(43);
      tempHist_p->GetYaxis()->SetLabelFont(43);

      tempHist_p->GetXaxis()->SetTitleSize(12);
      tempHist_p->GetYaxis()->SetTitleSize(12);
      tempHist_p->GetXaxis()->SetLabelSize(12);
      tempHist_p->GetYaxis()->SetLabelSize(12);

      tempHist_p->GetYaxis()->SetNdivisions(505);
      tempHist_p->GetXaxis()->SetTitleOffset(4.0);
      tempHist_p->GetYaxis()->SetTitleOffset(1.6);

      Double_t min = 5.;

      for(Int_t bI = 0; bI < tempHist_p->GetNbinsX(); ++bI){
	if(tempHist_p->GetBinContent(bI+1) > 0 && tempHist_p->GetBinContent(bI+1) < min) min = tempHist_p->GetBinContent(bI+1);
	if(tempHist2_p->GetBinContent(bI+1) > 0 && tempHist2_p->GetBinContent(bI+1) < min) min = tempHist2_p->GetBinContent(bI+1);
      }

      min /= 10.;

      tempHist_p->SetMaximum(4.);
      tempHist_p->SetMinimum(min);

      tempHist_p->GetXaxis()->SetTitle("Reco. Jet p_{T}");
      tempHist_p->DrawCopy("HIST E1 P");

      tempHist2_p->SetMarkerStyle(25);
      if(canvStr.find("NoID") == std::string::npos) tempHist2_p->DrawCopy("HIST E1 P SAME");

      gPad->SetLogy();

      delete tempHist_p;
      delete tempHist2_p;

      spectRat[3]->cd();

      tempHist_p = (TH1D*)badHists.at(i)->Clone("temp");
      denomPos = -1;
      
      for(unsigned int iter = 0; iter < goodHists.size(); ++iter){
	std::string canvStr2 = goodHists.at(iter)->GetName();
	if(canvStr2.find(type) == std::string::npos) continue;
	if(canvStr2.find(centStr) == std::string::npos) continue;

	if(recoPtStr.size() != 0){
	  if(canvStr2.find(recoPtStr) == std::string::npos) continue;
	}

	if(canvStr.find("Weighted") != std::string::npos){
	  if(canvStr2.find("Weighted") != std::string::npos && canvStr2.find(algoStr) != std::string::npos && canvStr2.find("NoID") != std::string::npos){
	    denomPos = iter;
	    //	    break;
	  }
	  else if(canvStr2.find("Weighted") != std::string::npos && canvStr2.find(algoStr) != std::string::npos && canvStr2.find("LightMUID") != std::string::npos){
	    denomPos2 = iter;
	  }
	}
	else{
	  if(canvStr2.find("Unweighted") != std::string::npos && canvStr2.find(algoStr) != std::string::npos && canvStr2.find("NoID") != std::string::npos){
	    denomPos = iter;
	    //	    break;
	  }
	  else if(canvStr2.find("Weighted") != std::string::npos && canvStr2.find(algoStr) != std::string::npos && canvStr2.find("LightMUID") != std::string::npos){
	    denomPos2 = iter;
	  }
	}
      }
    
      //std::cout << __LINE__ << std::endl;

      std::cout << "Found denom for hist \'" << canvStr << "\': " << denomPos << std::endl; 
     
      tempHist_p->GetYaxis()->SetTitle("Bad/Good");
      tempHist_p->Divide(goodHists.at(denomPos));
      tempHist_p->GetXaxis()->SetTitleFont(43);
      tempHist_p->GetYaxis()->SetTitleFont(43);
      tempHist_p->GetXaxis()->SetLabelFont(43);
      tempHist_p->GetYaxis()->SetLabelFont(43);

      tempHist_p->GetXaxis()->SetTitleSize(12);
      tempHist_p->GetYaxis()->SetTitleSize(12);
      tempHist_p->GetXaxis()->SetLabelSize(12);
      tempHist_p->GetYaxis()->SetLabelSize(12);

      tempHist_p->GetYaxis()->SetNdivisions(505);
      tempHist_p->GetXaxis()->SetTitleOffset(4.0);
      tempHist_p->GetYaxis()->SetTitleOffset(1.6);

      min = 5.;

      for(Int_t bI = 0; bI < tempHist_p->GetNbinsX(); ++bI){
	if(tempHist_p->GetBinContent(bI+1) > 0 && tempHist_p->GetBinContent(bI+1) < min) min = tempHist_p->GetBinContent(bI+1);
      }

      min /= 10.;

      tempHist_p->SetMaximum(4.);
      tempHist_p->SetMinimum(min);

      tempHist_p->GetXaxis()->SetTitle("Reco. Jet p_{T}");
      tempHist_p->DrawCopy("HIST E1 P");

      gPad->SetLogy();

      delete tempHist_p;
    }
    
    //std::cout << __LINE__ << std::endl;

    canv_p->SaveAs(("pdfDir/" + dateStr + "/" + algoStr2 + "/" + canvStr + "_" + dateStr + ".pdf").c_str());

    if(true/*canvStr.find("jtPt_") != std::string::npos*/){
      for(Int_t sI = 0; sI < nSpect; ++sI){
	delete spectRat[sI];
      }
    }
    delete canv_p;

    goodHists.at(i)->Write("", TObject::kOverwrite);
    badHists.at(i)->Write("", TObject::kOverwrite);
  }

  for(unsigned int i = 0; i < fineHists.size(); ++i){
    fineHists.at(i)->Write("", TObject::kOverwrite);
    delete fineHists.at(i);
  }

  //std::cout << __LINE__ << std::endl;

  for(unsigned int i = 0; i < goodHists.size(); ++i){
    delete goodHists.at(i);
    delete badHists.at(i);
  }


  for(unsigned int i = 0; i < goodHistsTH2.size(); ++i){
    std::string canvStr = goodHistsTH2.at(i)->GetName();
    canvStr.replace(canvStr.find("_Good"), std::string("_Good").size(), "");

    TCanvas* canv_p = new TCanvas("canv_p", "", 450*2, 450);
    
    if(goodHistsTH2.at(i)->Integral() > 0) goodHistsTH2.at(i)->Scale(1./goodHistsTH2.at(i)->Integral());
    if(badHistsTH2.at(i)->Integral() > 0) badHistsTH2.at(i)->Scale(1./badHistsTH2.at(i)->Integral());

    double maxVal = -1;
    double minVal = 9999;

    for(Int_t bIX = 0; bIX < goodHistsTH2.at(i)->GetNbinsX(); ++bIX){
      for(Int_t bIY = 0; bIY < goodHistsTH2.at(i)->GetNbinsY(); ++bIY){
	if(maxVal < goodHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1)) maxVal = goodHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1);
	if(maxVal < badHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1)) maxVal = badHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1);
	
	if(minVal > goodHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1) && goodHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1) > 0) minVal = goodHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1);
	if(minVal > badHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1) && badHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1) > 0) minVal = badHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1);
      }
    }

    goodHistsTH2.at(i)->SetMaximum(maxVal*5.);
    goodHistsTH2.at(i)->SetMinimum(minVal/5.);

    canv_p->Divide(2, 1);
    canv_p->cd();
    canv_p->cd(1);

    goodHistsTH2.at(i)->SetTitle("GOOD");
    goodHistsTH2.at(i)->DrawCopy("COLZ");
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    canv_p->cd();
    canv_p->cd(2);

    badHistsTH2.at(i)->SetTitle("BAD");
    badHistsTH2.at(i)->DrawCopy("COLZ");
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    //    canv_p->SaveAs(("pdfDir/" + dateStr + "/" + algoStr2 + "/" + canvStr + "_" + dateStr + ".pdf").c_str());
    delete canv_p;

    goodHistsTH2.at(i)->Write("", TObject::kOverwrite);
    badHistsTH2.at(i)->Write("", TObject::kOverwrite);

    delete goodHistsTH2.at(i);
    delete badHistsTH2.at(i);
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage: ./bin/checkGoodJetBadJetPF.exe <inFileName> <isPP-Opt>" << std::endl;
    return 1;
  }

  int retVal = 0; 
  if(argc == 2) retVal += checkGoodJetBadJetPF(argv[1]);
  else if(argc == 3) retVal += checkGoodJetBadJetPF(argv[1], std::stoi(argv[2]));

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
