//cpp dependencies
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

//ROOT dependencies
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNamed.h"
#include "TMath.h"
#include "TRandom3.h"

//Local FullJR (MainAnalysis) dependencies
#include "MainAnalysis/include/cutPropagator.h"

//Non-local FullJR (Utility, etc.) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/etaPhiFunc.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/ncollFunctions_5TeV.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/scaleErrorTool.h"
#include "Utility/include/specialHYDJETEventExclude.h"
#include "Utility/include/stringUtil.h"

template <typename T>
std::string to_string_with_precision(const T a_value, const int n)
{
  std::ostringstream out;
  out << std::setprecision(n) << a_value;
  return out.str();
}

int makeDeriveFlatResponse(const std::string inName, bool isPP = false)
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

  TFile* inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(0)).c_str(), "READ");
  std::vector<std::string> responseTrees = returnRootFileContentsList(inFile_p, "TTree", "JetAna");  
  
  unsigned int pos = 0;
  while(responseTrees.size() > pos){
    //For testing, uncomment to exclude all but R=0.4 trees
    //    if(responseTrees.at(pos).find("akCs4") == std::string::npos) responseTrees.erase(responseTrees.begin()+pos);

    if(isPP){
      if(responseTrees.at(pos).find("akCs") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
      else if(responseTrees.at(pos).find("akPu") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
      else ++pos;
    }
    else{
      if(responseTrees.at(pos).find("akCs3PF") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
      else if(responseTrees.at(pos).find("akCs4PF") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos); 
      else if(responseTrees.at(pos).find("akPu3PF") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
      else if(responseTrees.at(pos).find("akPu4PF") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
      else if(responseTrees.at(pos).find("akCs3PU3PFJet") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
      else if(responseTrees.at(pos).find("akCs4PU3PFJet") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
      else if(responseTrees.at(pos).find("akCs6PU3PFJet") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
      else if(responseTrees.at(pos).find("akCs8PU3PFJet") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
      else if(responseTrees.at(pos).find("akCs10PU3PFJet") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
      else ++pos;
    }
  }
  
  inFile_p->Close();
  delete inFile_p;
  inFile_p = NULL;

  Int_t posR4Temp = -1;
  const Int_t nTrees = responseTrees.size();

  std::cout << "Making response matrices for the following " << nTrees << " jet trees: " << std::endl;
  for(int jI = 0; jI < nTrees; ++jI){
    std::cout << " " << jI << "/" << nTrees << ": " << responseTrees.at(jI) << std::endl;

    if(responseTrees.at(jI).find("akCs4") != std::string::npos && posR4Temp < 0) posR4Temp = jI;
    else if(responseTrees.at(jI).find("ak4") != std::string::npos && posR4Temp < 0) posR4Temp = jI;
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

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate()) + "_" + std::to_string(date->GetHour());
  delete date;

  const std::string fullPath = std::getenv("FULLJRDIR");

  std::string outFileName = inName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
  else if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  outFileName = "output/" + outFileName + "_FlatGenJetResponse_" + dateStr + ".root";

  checkMakeDir("output");

  const Double_t jtAbsEtaMax = 2.;

  const Int_t nJtAbsEtaBins = 5;
  const Double_t jtAbsEtaBinsLow[nJtAbsEtaBins] = {0.0, 0.5, 1.0, 1.5, 0.0};
  const Double_t jtAbsEtaBinsHi[nJtAbsEtaBins] = {0.5, 1.0, 1.5, 2.0, 2.0};


  const Int_t nJtPtBins = 50;
  Double_t jtPtBins[nJtPtBins+1];
  getLinBins(100, 1100, nJtPtBins, jtPtBins);

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);

  TDirectory* dir_p[nTrees] = {NULL};

  TH1D* genJtPt_h[nTrees][nCentBins][nJtAbsEtaBins];
  TH1D* genJtPt_Flat_h[nTrees][nCentBins][nJtAbsEtaBins];

  for(Int_t dI = 0; dI < nTrees; ++dI){
    outFile_p->cd();
    std::string dirName = responseTrees.at(dI);
    dirName = dirName.substr(0, dirName.find("/"));

    dir_p[dI] = (TDirectory*)outFile_p->mkdir(dirName.c_str());

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
      if(isPP) centStr = "PP";

      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
    
	genJtPt_h[dI][cI][aI] = new TH1D(("genJtPt_" + dirName + "_" + centStr + "_" + jtAbsEtaStr + "_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
	genJtPt_Flat_h[dI][cI][aI] = new TH1D(("genJtPt_Flat_" + dirName + "_" + centStr + "_" + jtAbsEtaStr + "_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
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

    Int_t ngen_[nTrees];
    Float_t genpt_[nTrees][nMaxJet];
    Float_t genphi_[nTrees][nMaxJet];
    Float_t geneta_[nTrees][nMaxJet];
    Int_t gensubid_[nTrees][nMaxJet];

    for(Int_t tI = 0; tI < nTrees; ++tI){
      jetTrees_p[tI] = (TTree*)inFile_p->Get(responseTrees.at(tI).c_str());
      jetTrees_p[tI]->SetBranchStatus("*", 0);
      jetTrees_p[tI]->SetBranchStatus("ngen", 1);
      jetTrees_p[tI]->SetBranchStatus("genpt", 1);
      jetTrees_p[tI]->SetBranchStatus("geneta", 1);
      jetTrees_p[tI]->SetBranchStatus("genphi", 1);
      jetTrees_p[tI]->SetBranchStatus("gensubid", 1);

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

    const Int_t nEntries = TMath::Min((Int_t)1000000000, (Int_t)jetTrees_p[0]->GetEntries());
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
    
      for(Int_t tI = 0; tI < nTrees; ++tI){

	for(Int_t jI = 0; jI < ngen_[tI]; ++jI){
	  if(TMath::Abs(geneta_[tI][jI]) >= jtAbsEtaMax) continue;
	  if(genpt_[tI][jI] < jtPtBins[0]) continue;
	  std::vector<int> jtAbsEtaPoses;
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){

	    if(TMath::Abs(geneta_[tI][jI]) >= jtAbsEtaBinsLow[aI] && TMath::Abs(geneta_[tI][jI]) < jtAbsEtaBinsHi[aI]){
	      jtAbsEtaPoses.push_back(aI);
	    }
	  }

	  for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI){
	    genJtPt_h[tI][centPos][jtAbsEtaPoses.at(aI)]->Fill(genpt_[tI][jI], fullWeight_);
	  }
	}
      }
    }

    inFile_p->Close();
    delete inFile_p;
    inFile_p = NULL;
  }


  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "Processing file " << fI << "/" << fileList.size() << ": \'" << fileList.at(fI) << "\'" << std::endl;

    inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(fI)).c_str(), "READ");
    TTree* jetTrees_p[nTrees] = {NULL};

    Int_t ngen_[nTrees];
    Float_t genpt_[nTrees][nMaxJet];
    Float_t genphi_[nTrees][nMaxJet];
    Float_t geneta_[nTrees][nMaxJet];
    Int_t gensubid_[nTrees][nMaxJet];

    for(Int_t tI = 0; tI < nTrees; ++tI){
      jetTrees_p[tI] = (TTree*)inFile_p->Get(responseTrees.at(tI).c_str());
      jetTrees_p[tI]->SetBranchStatus("*", 0);
      jetTrees_p[tI]->SetBranchStatus("ngen", 1);
      jetTrees_p[tI]->SetBranchStatus("genpt", 1);
      jetTrees_p[tI]->SetBranchStatus("geneta", 1);
      jetTrees_p[tI]->SetBranchStatus("genphi", 1);
      jetTrees_p[tI]->SetBranchStatus("gensubid", 1);

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

    const Int_t nEntries = TMath::Min((Int_t)1000000000, (Int_t)jetTrees_p[0]->GetEntries());
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
    
      for(Int_t tI = 0; tI < nTrees; ++tI){

	for(Int_t jI = 0; jI < ngen_[tI]; ++jI){
	  if(TMath::Abs(geneta_[tI][jI]) >= jtAbsEtaMax) continue;
	  if(genpt_[tI][jI] < jtPtBins[0]) continue;
	  std::vector<int> jtAbsEtaPoses;
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    if(TMath::Abs(geneta_[tI][jI]) >= jtAbsEtaBinsLow[aI] && TMath::Abs(geneta_[tI][jI]) < jtAbsEtaBinsHi[aI]){
	      jtAbsEtaPoses.push_back(aI);
	    }
	  }

	  Int_t genBinPos = -1;
	  for(Int_t bIX = 0; bIX < nJtPtBins; ++bIX){
	    if(genpt_[tI][jI] >= jtPtBins[bIX] && genpt_[tI][jI] < jtPtBins[bIX+1]){
	      genBinPos = bIX+1;
	      break;
	    }
	  }
	  if(genBinPos < 0){
	    if(genpt_[tI][jI] >= jtPtBins[nJtPtBins]) genBinPos = nJtPtBins;
	    else{
	      std::cout << "WARNING: pt " << genpt_[tI][jI] << " undefined" << std::endl;
	    }
	  }
	  

	  for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI){
	    double ptWeight = ((Double_t)(genJtPt_h[tI][centPos][jtAbsEtaPoses.at(aI)]->GetBinContent(1)))/(Double_t)(genJtPt_h[tI][centPos][jtAbsEtaPoses.at(aI)]->GetBinContent(genBinPos));
	    genJtPt_Flat_h[tI][centPos][jtAbsEtaPoses.at(aI)]->Fill(genpt_[tI][jI], fullWeight_*ptWeight);
	  }
	}
      }
    }

    inFile_p->Close();
    delete inFile_p;
    inFile_p = NULL;
  }

  outFile_p->cd();

  for(Int_t dI = 0; dI < nTrees; ++dI){
    outFile_p->cd();
    dir_p[dI]->cd();

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	genJtPt_h[dI][cI][aI]->Write("", TObject::kOverwrite);
	genJtPt_Flat_h[dI][cI][aI]->Write("", TObject::kOverwrite);
      }
    }
  }


  outFile_p->cd();

  for(Int_t dI = 0; dI < nTrees; ++dI){
    outFile_p->cd();
    dir_p[dI]->cd();

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	delete genJtPt_h[dI][cI][aI];
	delete genJtPt_Flat_h[dI][cI][aI];
      }
    }
  }    

  TDirectory* cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");
  TDirectory* subDir_p = (TDirectory*)cutDir_p->mkdir("subDir");

  cutPropagator cutProp;
  cutProp.Clean();
  cutProp.SetInFileNames({inName});
  cutProp.SetInFullFileNames(fileList);
  cutProp.SetIsPP(isPP);
  cutProp.SetJtAbsEtaMax(jtAbsEtaMax);
  cutProp.SetNJtAlgos(nTrees);
  cutProp.SetJtAlgos(responseTrees);
  cutProp.SetNJtAbsEtaBins(nJtAbsEtaBins);
  cutProp.SetJtAbsEtaBinsLow(nJtAbsEtaBins, jtAbsEtaBinsLow);
  cutProp.SetJtAbsEtaBinsHi(nJtAbsEtaBins, jtAbsEtaBinsHi);
  cutProp.SetNPthats(pthats.size());
  cutProp.SetPthats(pthats);
  cutProp.SetPthatWeights(pthatWeights);
  cutProp.SetNCentBins(nCentBins);
  cutProp.SetCentBinsLow(centBinsLow);
  cutProp.SetCentBinsHi(centBinsHi);

  if(!cutProp.WriteAllVarToFile(outFile_p, cutDir_p, subDir_p)) std::cout << "Warning: Cut writing has failed" << std::endl;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage ./bin/makeDeriveFlatResponse.exe <inName> <isPP-Opt>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  if(argc == 2) retVal += makeDeriveFlatResponse(argv[1]);
  else if(argc == 3) retVal += makeDeriveFlatResponse(argv[1], std::stoi(argv[2]));

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
