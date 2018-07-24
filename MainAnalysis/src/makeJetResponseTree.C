//cpp dependencies
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>

//ROOT dependencies
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNamed.h"
#include "TMath.h"
#include "TRandom3.h"

//RooUnfold dependencies
#include "src/RooUnfoldResponse.h"

//Local FullJR (MainAnalysis) dependencies
#include "MainAnalysis/include/cutPropagator.h"

//Non-local FullJR (Utility, etc.) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/ncollFunctions_5TeV.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/specialHYDJETEventExclude.h"


template <typename T>
std::string to_string_with_precision(const T a_value, const int n)
{
  std::ostringstream out;
  out << std::setprecision(n) << a_value;
  return out.str();
}

int makeJetResponseTree(const std::string inName, bool isPP = false)
{
  TRandom3* randGen_p = new TRandom3(0);

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

  TFile* inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(0)).c_str(), "READ");
  std::vector<std::string> responseTrees = returnRootFileContentsList(inFile_p, "TTree", "JetAna");

  inFile_p->Close();
  delete inFile_p;
  inFile_p = NULL;

  Int_t posR4Temp = -1;
  const Int_t nTrees = responseTrees.size();
  std::cout << "Making response matrices for the following " << nTrees << " jet trees: " << std::endl;
  for(int jI = 0; jI < nTrees; ++jI){
    std::cout << " " << jI << "/" << nTrees << ": " << responseTrees.at(jI) << std::endl;

    if(responseTrees.at(jI).find("akCs4") != std::string::npos && posR4Temp < 0) posR4Temp = jI;
  }
  const Int_t posR4 = posR4Temp;

  std::cout << "PosR4: " << posR4 << std::endl;

  const Int_t nCentBinsPerma = 4;
  const Int_t centBinsLowPerma[nCentBinsPerma] = {0, 10, 30, 50};
  const Int_t centBinsHiPerma[nCentBinsPerma] = {10, 30, 50, 90};

  Int_t nCentBinsTemp = 1;
  if(!isPP) nCentBinsTemp = nCentBinsPerma;

  const Int_t nCentBins =nCentBinsTemp;
  std::vector<Int_t> centBinsLow, centBinsHi;
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

  const Int_t nPthatBins = 100;
  const Float_t pthatLow = pthats.at(0);
  const Float_t pthatHi = pthats.at(pthats.size()-1)*2.;
  Double_t pthatBins[nPthatBins+1];
  getLinBins(pthatLow, pthatHi, nPthatBins, pthatBins);

  const Int_t nCentBins2 = 100;
  const Float_t centBinsLow2 = 0;
  const Float_t centBinHi2 = 100;
  Double_t centBins2[nCentBins2+1];
  getLinBins(centBinsLow2, centBinHi2, nCentBins2, centBins2);

  const double fracParaFills = 0.1;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  std::string outFileName = inName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
  else if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  outFileName = "output/" + outFileName + "_JetResponse_" + dateStr + ".root";

  checkMakeDir("output");

  const Double_t jtAbsEtaMax = 2.;
  Int_t anomolousJetCount = 0;
  const Int_t maxAnomolousJet = 1000;

  const Int_t nJtAbsEtaBins = 5;
  const Double_t jtAbsEtaBinsLow[nJtAbsEtaBins] = {0.0, 0.5, 1.0, 1.5, 0.0};
  const Double_t jtAbsEtaBinsHi[nJtAbsEtaBins] = {0.5, 1.0, 1.5, 2.0, 2.0};

  const Int_t nJtPtBins = 10;
  const Double_t jtPtBins[nJtPtBins+1] = {100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1100.};

  const Int_t nID = 3;
  const std::string idStr[nID] = {"NoID", "LightMUID", "LightMUAndCHID"};
  const Double_t jtPfCHMFCutLow[nID] = {0.0, 0.0, 0.00};
  const Double_t jtPfCHMFCutHi[nID] = {1.0, 1.0, 0.90};
  const Double_t jtPfMUMFCutLow[nID] = {0.0, 0.0, 0.00};
  const Double_t jtPfMUMFCutHi[nID] = {1.0, 0.60, 0.60};
  //LightMUAndCHID == jtPfCHMF < 0.9 && jtPfMUMF < 0.6

  const Int_t nResponseBins = 300;
  Double_t responseBins[nResponseBins+1];
  getLinBins(0., 9., nResponseBins, responseBins);

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TDirectory* generalDir_p = (TDirectory*)outFile_p->mkdir("generalHistDir");
  generalDir_p->cd();
  TH1D* pthat_h = new TH1D("pthat_h", ";p_{T} Hat;Counts (Unweighted)", nPthatBins, pthatBins);
  TH1D* pthatWeighted_h = new TH1D("pthatWeighted_h", ";p_{T} Hat;Counts (Weighted)", nPthatBins, pthatBins);
  TH1D* pthatFullWeighted_h = new TH1D("pthatFullWeighted_h", ";p_{T} Hat;Counts (Full Weighted)", nPthatBins, pthatBins);
  TH1D* pthatFullRatio_h = new TH1D("pthatFullRatio_h", ";p_{T} Hat;Ratio", nPthatBins, pthatBins);

  centerTitles({pthat_h, pthatWeighted_h, pthatFullWeighted_h, pthatFullRatio_h});
  setSumW2({pthat_h, pthatWeighted_h, pthatFullWeighted_h, pthatFullRatio_h});

  TH1D* centrality_h = NULL;
  TH1D* centralityWeighted_h = NULL;
  TH1D* centralityFullWeighted_h = NULL;
  TH1D* centralityFullRatio_h = NULL;

  if(!isPP){
    centrality_h = new TH1D("centrality_h", ";Centrality (%);Counts (Unweighted)", nCentBins2, centBins2);
    centralityWeighted_h = new TH1D("centralityWeighted_h", ";Centrality (%);Counts (Weighted)", nCentBins2, centBins2);
    centralityFullWeighted_h = new TH1D("centralityFullWeighted_h", ";Centrality (%);Counts (Full Weighted)", nCentBins2, centBins2);
    centralityFullRatio_h = new TH1D("centralityFullRatio_h", ";Centrality (%);Ratio", nCentBins2, centBins2);
    
    centerTitles({centrality_h, centralityWeighted_h, centralityFullWeighted_h, centralityFullRatio_h});
    setSumW2({centrality_h, centralityWeighted_h, centralityFullWeighted_h, centralityFullRatio_h});
  }

  TDirectory* dir_p[nTrees] = {NULL};
  TH1D* recoJtPt_h[nTrees][nCentBins][nID][nJtAbsEtaBins];
  TH1D* recoJtPt_RecoTrunc_h[nTrees][nCentBins][nID][nJtAbsEtaBins];
  TH1D* genJtPt_h[nTrees][nCentBins][nID][nJtAbsEtaBins];
  TH1D* recoJtPtPerGenPtBin_h[nTrees][nCentBins][nID][nJtAbsEtaBins][nJtPtBins];
  TH1D* genJtPtPerRecoPtBin_h[nTrees][nCentBins][nID][nJtAbsEtaBins][nJtPtBins];
  TH1D* recoJtPtPerGenPtBinWeighted_h[nTrees][nCentBins][nID][nJtAbsEtaBins][nJtPtBins];
  TH1D* genJtPtPerRecoPtBinWeighted_h[nTrees][nCentBins][nID][nJtAbsEtaBins][nJtPtBins];

  TH1D* recoJtPt_ParaFills_h[nTrees][nCentBins][nID][nJtAbsEtaBins];
  TH1D* recoJtPt_RecoTrunc_ParaFills_h[nTrees][nCentBins][nID][nJtAbsEtaBins];
  TH1D* genJtPt_ParaFills_h[nTrees][nCentBins][nID][nJtAbsEtaBins];

  TH2D* response_h[nTrees][nCentBins][nID][nJtAbsEtaBins];
  TH2D* response_RecoTrunc_h[nTrees][nCentBins][nID][nJtAbsEtaBins];
  RooUnfoldResponse* rooResponse_RecoTrunc_h[nTrees][nCentBins][nID][nJtAbsEtaBins];

  for(Int_t dI = 0; dI < nTrees; ++dI){
    outFile_p->cd();
    std::string dirName = responseTrees.at(dI);
    dirName = dirName.substr(0, dirName.find("/"));

    dir_p[dI] = (TDirectory*)outFile_p->mkdir(dirName.c_str());

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
      if(isPP) centStr = "PP";

      for(Int_t iI = 0; iI < nID; ++iI){
	
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	  
	  recoJtPt_h[dI][cI][iI][aI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
	  recoJtPt_RecoTrunc_h[dI][cI][iI][aI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_RecoTrunc_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
	  genJtPt_h[dI][cI][iI][aI] = new TH1D(("genJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
	  
	  recoJtPt_ParaFills_h[dI][cI][iI][aI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_ParaFills_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
	  recoJtPt_RecoTrunc_ParaFills_h[dI][cI][iI][aI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_RecoTrunc_ParaFills_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
	  genJtPt_ParaFills_h[dI][cI][iI][aI] = new TH1D(("genJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_ParaFills_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
	  
	  response_h[dI][cI][iI][aI] = new TH2D(("response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_h").c_str(), ";Reco. Jet p_{T};Gen. Jet p_{T}", nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);
	  response_RecoTrunc_h[dI][cI][iI][aI] = new TH2D(("response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_RecoTrunc_h").c_str(), ";Reco. Jet p_{T};Gen. Jet p_{T}", nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);
	  
	  rooResponse_RecoTrunc_h[dI][cI][iI][aI] = new RooUnfoldResponse(("rooResponse_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_RecoTrunc_h").c_str(), "");
	  rooResponse_RecoTrunc_h[dI][cI][iI][aI]->Setup(recoJtPt_RecoTrunc_h[dI][cI][iI][aI], genJtPt_h[dI][cI][iI][aI]);
	  
	  centerTitles({recoJtPt_RecoTrunc_h[dI][cI][iI][aI], recoJtPt_h[dI][cI][iI][aI], genJtPt_h[dI][cI][iI][aI], recoJtPt_RecoTrunc_ParaFills_h[dI][cI][iI][aI], recoJtPt_ParaFills_h[dI][cI][iI][aI], genJtPt_ParaFills_h[dI][cI][iI][aI], response_h[dI][cI][iI][aI], response_RecoTrunc_h[dI][cI][iI][aI]});
	  setSumW2({recoJtPt_RecoTrunc_h[dI][cI][iI][aI], recoJtPt_h[dI][cI][iI][aI], genJtPt_h[dI][cI][iI][aI], recoJtPt_RecoTrunc_ParaFills_h[dI][cI][iI][aI], recoJtPt_ParaFills_h[dI][cI][iI][aI], genJtPt_ParaFills_h[dI][cI][iI][aI], response_h[dI][cI][iI][aI], response_RecoTrunc_h[dI][cI][iI][aI]});
	  
	  
	  for(Int_t jI = 0; jI < nJtPtBins; ++jI){
	    const std::string jtPtStr = "Pt" + prettyString(jtPtBins[jI], 1, true) + "to" + prettyString(jtPtBins[jI+1], 1, true);
	    const std::string jtPtStr2 = prettyString(jtPtBins[jI], 1, false) + "< p_{T,Gen} <" + prettyString(jtPtBins[jI+1], 1, false);
	    const std::string jtPtStr3 = prettyString(jtPtBins[jI], 1, false) + "< p_{T,Reco} <" + prettyString(jtPtBins[jI+1], 1, false);
	  
	    recoJtPtPerGenPtBin_h[dI][cI][iI][aI][jI] = new TH1D(("recoJtPtPerGenPtBin" + dirName + "_" + centStr + "_" + idStr[iI] + "_Gen" + jtPtStr + "_" + jtAbsEtaStr + "_h").c_str(), (";Reco p_{T}/Gen p_{T};Counts (" + jtPtStr2 + ")").c_str(), nResponseBins, responseBins);
	    
	    genJtPtPerRecoPtBin_h[dI][cI][iI][aI][jI] = new TH1D(("genJtPtPerRecoPtBin" + dirName + "_" + centStr + "_" + idStr[iI] + "_Reco" + jtPtStr + "_" + jtAbsEtaStr + "_h").c_str(), (";Reco p_{T}/Gen p_{T};Counts (" + jtPtStr3 + ")").c_str(), nResponseBins, responseBins);
	    
	    recoJtPtPerGenPtBinWeighted_h[dI][cI][iI][aI][jI] = new TH1D(("recoJtPtPerGenPtBin" + dirName + "_" + centStr + "_" + idStr[iI] + "_Gen" + jtPtStr + "_" + jtAbsEtaStr + "_Weighted_h").c_str(), (";Reco p_{T}/Gen p_{T};Counts (" + jtPtStr2 + ")").c_str(), nResponseBins, responseBins);
	    
	    genJtPtPerRecoPtBinWeighted_h[dI][cI][iI][aI][jI] = new TH1D(("genJtPtPerRecoPtBin" + dirName + "_" + centStr + "_" + idStr[iI] + "_Reco" + jtPtStr + "_" + jtAbsEtaStr + "_Weighted_h").c_str(), (";Reco p_{T}/Gen p_{T};Counts (" + jtPtStr3 + ")").c_str(), nResponseBins, responseBins);
	    
	    centerTitles({recoJtPtPerGenPtBin_h[dI][cI][iI][aI][jI], genJtPtPerRecoPtBin_h[dI][cI][iI][aI][jI], recoJtPtPerGenPtBinWeighted_h[dI][cI][iI][aI][jI], genJtPtPerRecoPtBinWeighted_h[dI][cI][iI][aI][jI]});
	  }
	}
      }
    }
  }


  const Int_t nMaxJet = 500;
  goodGlobalSelection globalSel;
  globalSel.setIsPbPb(!isPP);

  specialHYDJETEventExclude specialSel;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "Processing file " << fI << "/" << fileList.size() << ": \'" << fileList.at(fI) << "\'" << std::endl;

    inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(fI)).c_str(), "READ");
    TTree* jetTrees_p[nTrees] = {NULL};
    Int_t nref_[nTrees];
    Float_t jtpt_[nTrees][nMaxJet];
    Float_t rawpt_[nTrees][nMaxJet];
    Float_t jteta_[nTrees][nMaxJet];
    Float_t refpt_[nTrees][nMaxJet];
    Float_t jtPfCHF_[nTrees][nMaxJet];
    Float_t jtPfCEF_[nTrees][nMaxJet];
    Float_t jtPfNHF_[nTrees][nMaxJet];
    Float_t jtPfNEF_[nTrees][nMaxJet];
    Float_t jtPfMUF_[nTrees][nMaxJet];
    Float_t jtPfCHMF_[nTrees][nMaxJet];
    Float_t jtPfCEMF_[nTrees][nMaxJet];
    Float_t jtPfNHMF_[nTrees][nMaxJet];
    Float_t jtPfNEMF_[nTrees][nMaxJet];
    Float_t jtPfMUMF_[nTrees][nMaxJet];
    Int_t jtPfCHM_[nTrees][nMaxJet];
    Int_t jtPfCEM_[nTrees][nMaxJet];
    Int_t jtPfNHM_[nTrees][nMaxJet];
    Int_t jtPfNEM_[nTrees][nMaxJet];
    Int_t jtPfMUM_[nTrees][nMaxJet];

    Int_t ngen_[nTrees];
    Float_t genpt_[nTrees][nMaxJet];
    Float_t genphi_[nTrees][nMaxJet];
    Float_t geneta_[nTrees][nMaxJet];
    Int_t gensubid_[nTrees][nMaxJet];

    for(Int_t tI = 0; tI < nTrees; ++tI){
      jetTrees_p[tI] = (TTree*)inFile_p->Get(responseTrees.at(tI).c_str());
      jetTrees_p[tI]->SetBranchStatus("*", 0);
      jetTrees_p[tI]->SetBranchStatus("nref", 1);
      jetTrees_p[tI]->SetBranchStatus("jtpt", 1);
      jetTrees_p[tI]->SetBranchStatus("rawpt", 1);
      jetTrees_p[tI]->SetBranchStatus("jteta", 1);
      jetTrees_p[tI]->SetBranchStatus("refpt", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfCHF", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfCEF", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfNHF", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfNEF", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfMUF", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfCHMF", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfCEMF", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfNHMF", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfNEMF", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfMUMF", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfCHM", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfCEM", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfNHM", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfNEM", 1);
      jetTrees_p[tI]->SetBranchStatus("jtPfMUM", 1);
      jetTrees_p[tI]->SetBranchStatus("ngen", 1);
      jetTrees_p[tI]->SetBranchStatus("genpt", 1);
      jetTrees_p[tI]->SetBranchStatus("geneta", 1);
      jetTrees_p[tI]->SetBranchStatus("genphi", 1);
      jetTrees_p[tI]->SetBranchStatus("gensubid", 1);

      jetTrees_p[tI]->SetBranchAddress("nref", &(nref_[tI]));
      jetTrees_p[tI]->SetBranchAddress("jtpt", jtpt_[tI]);
      jetTrees_p[tI]->SetBranchAddress("rawpt", rawpt_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jteta", jteta_[tI]);
      jetTrees_p[tI]->SetBranchAddress("refpt", refpt_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfCHF", jtPfCHF_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfCEF", jtPfCEF_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfNHF", jtPfNHF_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfNEF", jtPfNEF_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfMUF", jtPfMUF_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfCHF", jtPfCHMF_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfCEF", jtPfCEMF_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfNHF", jtPfNHMF_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfNEF", jtPfNEMF_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfMUF", jtPfMUMF_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfCHM", jtPfCHM_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfCEM", jtPfCEM_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfNHM", jtPfNHM_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfNEM", jtPfNEM_[tI]);
      jetTrees_p[tI]->SetBranchAddress("jtPfMUM", jtPfMUM_[tI]);
      jetTrees_p[tI]->SetBranchAddress("ngen", &(ngen_[tI]));
      jetTrees_p[tI]->SetBranchAddress("genpt", genpt_[tI]);
      jetTrees_p[tI]->SetBranchAddress("geneta", geneta_[tI]);
      jetTrees_p[tI]->SetBranchAddress("genphi", genphi_[tI]);
      jetTrees_p[tI]->SetBranchAddress("gensubid", gensubid_[tI]);
    }
    
    Float_t pthat_;

    jetTrees_p[0]->SetBranchStatus("pthat", 1);
    jetTrees_p[0]->SetBranchAddress("pthat", &pthat_);

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

    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");

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

    const Int_t nEntries = jetTrees_p[0]->GetEntries();
    const Int_t printInterval = TMath::Max(1, nEntries/20);

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(nEntries >= 50000 && entry%printInterval == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;

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

      for(Int_t tI = 0; tI < nTrees; ++tI){jetTrees_p[tI]->GetEntry(entry);}
    
      Int_t centPos = -1;

      if(isPP) centPos = 0;
      else{
	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  if(centBinsLow.at(cI)*2 <= hiBin_ && hiBin_ < centBinsHi.at(cI)*2){
	    centPos = cI;
	    break;
	  }
	}
	if(centPos < 0) continue;
      }
    
      bool badJetSpecialSel = specialSel.CheckEventBadJet(ngen_[posR4], genpt_[posR4], genphi_[posR4], geneta_[posR4], gensubid_[posR4], entry);
      if(badJetSpecialSel) continue;

      Double_t pthatWeight_ = -1;
      for(unsigned int pI = 0; pI < pthats.size()-1; ++pI){
	if(pthats.at(pI) <= pthat_ && pthat_ < pthats.at(pI+1)){
	  pthatWeight_ = pthatWeights.at(pI);
	  break;
	}
      }
      if(pthat_ > pthats.at(pthats.size()-1)) pthatWeight_ = pthatWeights.at(pthatWeights.size()-1);

      if(pthatWeight_ < 0){
	std::cout << "WARNING - NO WEIGHT FOR pthat \'" << pthat_ << "\'. Set to 1" << std::endl;
	pthatWeight_ = 1.;
      }

      Double_t ncollWeight_ = 1.;
      if(!isPP) findNcoll_Renorm(hiBin_);
      Double_t fullWeight_ = ncollWeight_*pthatWeight_;      

      pthat_h->Fill(pthat_);
      pthatWeighted_h->Fill(pthat_, pthatWeight_);
      pthatFullWeighted_h->Fill(pthat_, fullWeight_);

      if(!isPP){
	centrality_h->Fill(hiBin_/2.);
	centralityWeighted_h->Fill(hiBin_/2., ncollWeight_);
	centralityFullWeighted_h->Fill(hiBin_/2., fullWeight_);
      }

      Bool_t isPara = randGen_p->Uniform(0., 1.) < fracParaFills;

      for(Int_t tI = 0; tI < nTrees; ++tI){
	for(Int_t jI = 0; jI < nref_[tI]; ++jI){
	  if(TMath::Abs(jteta_[tI][jI]) > jtAbsEtaMax) continue;

	  bool goodReco = (jtpt_[tI][jI] >= jtPtBins[0] && jtpt_[tI][jI] < jtPtBins[nJtPtBins]);
	  bool goodTruth = (refpt_[tI][jI] >= jtPtBins[0] && refpt_[tI][jI] < jtPtBins[nJtPtBins]);

	  if(!goodReco && !goodTruth) continue;
      	
	  //via https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
	  bool passesTightLepVeto = jtPfCHF_[tI][jI] < 0.9 && jtPfCEF_[tI][jI] < 0.9;
	  passesTightLepVeto = passesTightLepVeto && jtPfCHM_[tI][jI] + jtPfCEM_[tI][jI] + jtPfNHM_[tI][jI] + jtPfNEM_[tI][jI] + jtPfMUM_[tI][jI] > 1;
	  passesTightLepVeto = passesTightLepVeto && jtPfMUF_[tI][jI] < 0.8;
	  passesTightLepVeto = passesTightLepVeto && jtPfCHF_[tI][jI] > 0.0;
	  passesTightLepVeto = passesTightLepVeto && jtPfCHM_[tI][jI] > 0;
	  passesTightLepVeto = passesTightLepVeto && jtPfCEF_[tI][jI] < 0.9;

	  //	  if(!passesTightLepVeto) continue;
	  std::vector<bool> passesID;
	  for(Int_t iI = 0; iI < nID; ++iI){
	    bool pass = jtPfCHMFCutLow[iI] <= jtPfCHMF_[tI][jI] && jtPfCHMF_[tI][jI] <= jtPfCHMFCutHi[iI];
	    pass = pass && jtPfMUMFCutLow[iI] <= jtPfMUMF_[tI][jI] && jtPfMUMF_[tI][jI] <= jtPfMUMFCutHi[iI];
	    passesID.push_back(pass);
	  }
	  
	  bool goodRecoTrunc = (jtpt_[tI][jI] >= jtPtBins[1] && jtpt_[tI][jI] < jtPtBins[nJtPtBins-1]);
	  Int_t recoJtPos = -1;
	  Int_t genJtPos = -1;

	  for(Int_t posI = 0; posI < nJtPtBins; ++posI){
	    if(jtPtBins[posI] <= jtpt_[tI][jI] && jtPtBins[posI+1] > jtpt_[tI][jI]){
	      recoJtPos = posI;
	      break;
	    }
	  }

	  for(Int_t posI = 0; posI < nJtPtBins; ++posI){
	    if(jtPtBins[posI] <= refpt_[tI][jI] && jtPtBins[posI+1] > refpt_[tI][jI]){
	      genJtPos = posI;
	      break;
	    }
	  }

	  //	  if(recoJtPos < 0) std::cout << "WARNING: recoJtPos -1 for jtpt==" << jtpt_[tI][jI] << std::endl;
	  //	  if(genJtPos < 0) std::cout << "WARNING: genJtPos -1 for refpt==" << refpt_[tI][jI] << std::endl;
	
	  std::vector<int> jtAbsEtaPoses;
          for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
            if(TMath::Abs(jteta_[tI][jI]) >= jtAbsEtaBinsLow[aI] && TMath::Abs(jteta_[tI][jI]) < jtAbsEtaBinsHi[aI]){
              jtAbsEtaPoses.push_back(aI);
            }
          }

	  for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI){
	    for(unsigned int iI = 0; iI < passesID.size(); ++iI){
	      if(!passesID.at(iI)) continue;

	      if(!isPara){
		if(goodTruth && goodReco){
		  recoJtPt_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Fill(jtpt_[tI][jI], fullWeight_);
		  genJtPt_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Fill(refpt_[tI][jI], fullWeight_);
		  response_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Fill(jtpt_[tI][jI], refpt_[tI][jI], fullWeight_);
		  
		  if(genJtPos >= 0){
		    //		std::cout << __LINE__ << ", " << genJtPos << ", " << jtpt_[tI][jI] << ", " << refpt_[tI][jI] << ", " << fullWeight_ << std::endl;
		    //		std::cout << recoJtPtPerGenPtBin_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)][genJtPos]->GetName() << std::endl;
		    //		std::cout << recoJtPtPerGenPtBinWeighted_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)][genJtPos]->GetName() << std::endl;
		    
		    if(anomolousJetCount < maxAnomolousJet && jtpt_[tI][jI] > 800. && refpt_[tI][jI] < 200. && passesID.at(1) && responseTrees.at(tI).find("akPu") == std::string::npos){
		      std::cout << "Anomolous jet \'" << anomolousJetCount << "\':" << std::endl;
		      std::cout << " File: " << fileList.at(fI) << std::endl;
		      std::cout << " Algo, entry: " << responseTrees.at(tI) << ", " << entry << std::endl;
		      std::cout << " reco,gen: " << jtpt_[tI][jI] << ", " << refpt_[tI][jI] << ", " << jteta_[tI][jI] << std::endl;
		      std::cout << " run,lumi,evt: " << run_ << ", " << lumi_ << ", " << evt_ << std::endl;
		      
		      anomolousJetCount++;
		    }
		    
		    recoJtPtPerGenPtBin_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)][genJtPos]->Fill(jtpt_[tI][jI]/refpt_[tI][jI]);
		    recoJtPtPerGenPtBinWeighted_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)][genJtPos]->Fill(jtpt_[tI][jI]/refpt_[tI][jI], fullWeight_);
		    //		std::cout << __LINE__ << std::endl;
		  }
		  
		  if(recoJtPos >= 0 && refpt_[tI][jI] > 0){
		    //		std::cout << __LINE__ << std::endl;
		    genJtPtPerRecoPtBin_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)][recoJtPos]->Fill(jtpt_[tI][jI]/refpt_[tI][jI]);
		    genJtPtPerRecoPtBinWeighted_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)][recoJtPos]->Fill(jtpt_[tI][jI]/refpt_[tI][jI], fullWeight_);
		    //		std::cout << __LINE__ << std::endl;
		  }
		  
		  if(goodRecoTrunc){
		    recoJtPt_RecoTrunc_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Fill(jtpt_[tI][jI], fullWeight_);
		    response_RecoTrunc_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Fill(jtpt_[tI][jI], refpt_[tI][jI], fullWeight_);	    
		    rooResponse_RecoTrunc_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Fill(jtpt_[tI][jI], refpt_[tI][jI], fullWeight_);
		  }
		  else rooResponse_RecoTrunc_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Miss(refpt_[tI][jI], fullWeight_);
		}
		else if(goodTruth){
		  genJtPt_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Fill(refpt_[tI][jI], fullWeight_);
		  rooResponse_RecoTrunc_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Miss(refpt_[tI][jI], fullWeight_);
		}
	      }
	      else{
		if(goodTruth && goodReco){
		  recoJtPt_ParaFills_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Fill(jtpt_[tI][jI], fullWeight_);
		  genJtPt_ParaFills_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Fill(refpt_[tI][jI], fullWeight_);
		  
		  if(goodRecoTrunc) recoJtPt_RecoTrunc_ParaFills_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Fill(jtpt_[tI][jI], fullWeight_);
		}
		else if(goodTruth) genJtPt_ParaFills_h[tI][centPos][iI][jtAbsEtaPoses.at(aI)]->Fill(refpt_[tI][jI], fullWeight_);
	      }
	    }
	  }
	}
      }
    }

    inFile_p->Close();
    delete inFile_p;
    inFile_p = NULL;
  }

  outFile_p->cd();
  generalDir_p->cd();
  pthat_h->Write("", TObject::kOverwrite);
  pthatWeighted_h->Write("", TObject::kOverwrite);
  pthatFullWeighted_h->Write("", TObject::kOverwrite);
  pthatFullRatio_h->Divide(pthatFullWeighted_h, pthatWeighted_h);
  pthatFullRatio_h->Write("", TObject::kOverwrite);

  if(!isPP){
    centrality_h->Write("", TObject::kOverwrite);
    centralityWeighted_h->Write("", TObject::kOverwrite);
    centralityFullWeighted_h->Write("", TObject::kOverwrite);
    centralityFullRatio_h->Divide(centralityFullWeighted_h, centralityWeighted_h);
    centralityFullRatio_h->Write("", TObject::kOverwrite);
  }

  for(Int_t dI = 0; dI < nTrees; ++dI){
    outFile_p->cd();
    dir_p[dI]->cd();

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      for(Int_t iI = 0; iI < nID; ++iI){

	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  recoJtPt_h[dI][cI][iI][aI]->Write("", TObject::kOverwrite);
	  recoJtPt_RecoTrunc_h[dI][cI][iI][aI]->Write("", TObject::kOverwrite);
	  genJtPt_h[dI][cI][iI][aI]->Write("", TObject::kOverwrite);
	
	  recoJtPt_ParaFills_h[dI][cI][iI][aI]->Write("", TObject::kOverwrite);
	  recoJtPt_RecoTrunc_ParaFills_h[dI][cI][iI][aI]->Write("", TObject::kOverwrite);
	  genJtPt_ParaFills_h[dI][cI][iI][aI]->Write("", TObject::kOverwrite);
	  
	  response_h[dI][cI][iI][aI]->Write("", TObject::kOverwrite);
	  response_RecoTrunc_h[dI][cI][iI][aI]->Write("", TObject::kOverwrite);
	  rooResponse_RecoTrunc_h[dI][cI][iI][aI]->Write("", TObject::kOverwrite);
	  
	  for(Int_t jI = 0; jI < nJtPtBins; ++jI){
	    recoJtPtPerGenPtBin_h[dI][cI][iI][aI][jI]->Write("", TObject::kOverwrite);
	    genJtPtPerRecoPtBin_h[dI][cI][iI][aI][jI]->Write("", TObject::kOverwrite);
	    
	    recoJtPtPerGenPtBinWeighted_h[dI][cI][iI][aI][jI]->Write("", TObject::kOverwrite);
	    genJtPtPerRecoPtBinWeighted_h[dI][cI][iI][aI][jI]->Write("", TObject::kOverwrite);
	  }
	}
      }
    }
  }

  delete pthat_h;
  delete pthatWeighted_h;
  delete pthatFullWeighted_h;
  delete pthatFullRatio_h;

  if(!isPP){
    delete centrality_h;
    delete centralityWeighted_h;
    delete centralityFullWeighted_h;
    delete centralityFullRatio_h;
  }

  for(Int_t dI = 0; dI < nTrees; ++dI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      for(Int_t iI = 0; iI < nID; ++iI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  
	  delete recoJtPt_RecoTrunc_h[dI][cI][iI][aI];
	  delete recoJtPt_h[dI][cI][iI][aI];
	  delete genJtPt_h[dI][cI][iI][aI];
	  
	  delete recoJtPt_RecoTrunc_ParaFills_h[dI][cI][iI][aI];
	  delete recoJtPt_ParaFills_h[dI][cI][iI][aI];
	  delete genJtPt_ParaFills_h[dI][cI][iI][aI];
	  
	  delete response_h[dI][cI][iI][aI];
	  delete response_RecoTrunc_h[dI][cI][iI][aI];
	  delete rooResponse_RecoTrunc_h[dI][cI][iI][aI];
	  
	  for(Int_t jI = 0; jI < nJtPtBins; ++jI){
	    delete recoJtPtPerGenPtBin_h[dI][cI][iI][aI][jI];
	    delete genJtPtPerRecoPtBin_h[dI][cI][iI][aI][jI];
	    
	    delete recoJtPtPerGenPtBinWeighted_h[dI][cI][iI][aI][jI];
	    delete genJtPtPerRecoPtBinWeighted_h[dI][cI][iI][aI][jI];
	  }
	}
      }
    }
  }    

  outFile_p->cd();

  TDirectory* cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");

  cutPropagator cutProp;
  cutProp.SetIsPP(isPP);
  cutProp.SetJtAbsEtaMax(jtAbsEtaMax);
  cutProp.SetNJtPtBins(nJtPtBins);
  cutProp.SetJtPtBins(nJtPtBins+1, jtPtBins);
  cutProp.SetNJtAbsEtaBins(nJtAbsEtaBins);
  cutProp.SetJtAbsEtaBinsLow(nJtAbsEtaBins, jtAbsEtaBinsLow);
  cutProp.SetJtAbsEtaBinsHi(nJtAbsEtaBins, jtAbsEtaBinsHi);
  cutProp.SetNPthats(pthats.size());
  cutProp.SetPthats(pthats);
  cutProp.SetPthatWeights(pthatWeights);
  cutProp.SetNCentBins(nCentBins);
  cutProp.SetCentBinsLow(centBinsLow);
  cutProp.SetCentBinsHi(centBinsHi);
  cutProp.SetNID(nID);
  cutProp.SetIdStr(nID, idStr);
  cutProp.SetJtPfCHMFCutLow(nID, jtPfCHMFCutLow);
  cutProp.SetJtPfCHMFCutHi(nID, jtPfCHMFCutHi);
  cutProp.SetJtPfMUMFCutLow(nID, jtPfMUMFCutLow);
  cutProp.SetJtPfMUMFCutHi(nID, jtPfMUMFCutHi);

  if(!cutProp.WriteAllVarToFile(outFile_p, cutDir_p)) std::cout << "Warning: Cut writing has failed" << std::endl;

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage ./bin/makeJetResponseTree.exe <inName> <isPP-Opt>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  if(argc == 2) retVal += makeJetResponseTree(argv[1]);
  else if(argc == 3) retVal += makeJetResponseTree(argv[1], std::stoi(argv[2]));
  return retVal;
}
