#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TNamed.h"

#include "Utility/include/returnFileList.h"
#include "Utility/include/doGlobalDebug.h"
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/inToOutFileString.h"

int makeRawRAATree(const std::string inFileName, std::string outFileName = "")
{
  std::vector<std::string> fileList;

  if(checkFile(inFileName) && inFileName.find(".root") != std::string::npos) fileList.push_back(inFileName);
  else if(checkFile(inFileName) && inFileName.find(".txt") != std::string::npos){
    std::ifstream file(inFileName.c_str());
    std::string tempStr;
    while(std::getline(file,tempStr)){
      if(tempStr.find(".root") == std::string::npos) continue;
      fileList.push_back(tempStr);
    }
    file.close();
  }
  else if(checkDir(inFileName)){
    fileList = returnFileList(inFileName, "HiForest");
    unsigned int filePos = 0;
    while(filePos < fileList.size()){
      if(fileList.at(filePos).find("/merged/") != std::string::npos) fileList.erase(fileList.begin()+filePos);
      else ++filePos;
    }
  }

  unsigned int filePos = 0;
  while(filePos < fileList.size()){
    if(fileList.at(filePos).find("/failed/") != std::string::npos) fileList.erase(fileList.begin()+filePos);
    else ++filePos;
  }

  if(fileList.size() == 0){
    std::cout << "inFileName \'" << inFileName << "\' is invalid return 1" << std::endl;
    return 1;
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TFile* inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(0)).c_str(), "READ");
  std::vector<std::string> inJetTreeList_ = returnRootFileContentsList(inFile_p, "TTree", "JetAnalyzer");
  std::vector<std::string> inHITreeList_ = returnRootFileContentsList(inFile_p, "TTree", "hiEvtAnalyzer");
  TTree* hiTree_p = (TTree*)inFile_p->Get(inHITreeList_.at(0).c_str());
  bool isPP = hiTree_p->GetMaximum("hiBin") == hiTree_p->GetMinimum("hiBin");
  unsigned int pos = 0;
  while(inJetTreeList_.size() > pos){
    if(inJetTreeList_.at(pos).find("SD") != std::string::npos) inJetTreeList_.erase(inJetTreeList_.begin()+pos);
    else if(inJetTreeList_.at(pos).find("ak1PF") != std::string::npos) inJetTreeList_.erase(inJetTreeList_.begin()+pos);
    else if(inJetTreeList_.at(pos).find("ak2PF") != std::string::npos) inJetTreeList_.erase(inJetTreeList_.begin()+pos);
    else if(inJetTreeList_.at(pos).find("ak5PF") != std::string::npos) inJetTreeList_.erase(inJetTreeList_.begin()+pos);
    else ++pos;
    //    else if(inJetTreeList_.at(pos).find("ak4PF") == std::string::npos && isPP) inJetTreeList_.erase(inJetTreeList_.begin()+pos);
    //    else if(inJetTreeList_.at(pos).find("akCs4PU3PFFlow") == std::string::npos && !isPP) inJetTreeList_.erase(inJetTreeList_.begin()+pos);
    //    else ++pos;
  }
  inFile_p->Close();
  delete inFile_p;

  if(outFileName.size() == 0) outFileName = "output/" + inToOutFileString(inFileName, "RawRAAHIST_NJet" + std::to_string(inJetTreeList_.size()));
  else outFileName = "output/" + inToOutFileString(outFileName, "RawRAAHIST_NJet" + std::to_string(inJetTreeList_.size()));
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);

  std::map<ULong64_t, Int_t> runLumiMap;
  
  inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(0)).c_str(), "READ");
  bool hasHLT = returnRootFileContentsList(inFile_p, "TTree", "hltanalysis").size() > 0;


  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(inHITreeList_.size() != 1){
    std::cout << "inHITreeList not equal to 1 in size. return 1" << std::endl;
    inFile_p->Close();
    delete inFile_p;
    return 1;
  }
  inFile_p->Close();
  delete inFile_p;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  const int nJetAlgos = inJetTreeList_.size();
  std::vector<int> algoPos;
  for(Int_t i = 0; i < nJetAlgos; ++i){
    algoPos.push_back(i);
  }

  int tempNJetAlgosPlus = nJetAlgos;
  const int nJetAlgosPlus = tempNJetAlgosPlus;

  Int_t hiBin_;
  Float_t vz_;

  const Int_t nMaxJets = 500;
  Int_t nref_[nJetAlgos];
  Float_t jtpt_[nJetAlgos][nMaxJets];
  Float_t rawpt_[nJetAlgos][nMaxJets];
  Float_t jteta_[nJetAlgos][nMaxJets];
  Float_t jtphi_[nJetAlgos][nMaxJets]; 
  Float_t jtPfCHF_[nJetAlgos][nMaxJets];
  Float_t jtPfCEF_[nJetAlgos][nMaxJets];
  Float_t jtPfNHF_[nJetAlgos][nMaxJets];
  Float_t jtPfNEF_[nJetAlgos][nMaxJets];
  Float_t jtPfMUF_[nJetAlgos][nMaxJets];

  Int_t HBHENoiseFilterResultRun2Loose_;
  Int_t pprimaryVertexFilter_;
  Int_t pBeamScrapingFilter_;
  Int_t phfCoincFilter3_;
  Int_t pclusterCompatibilityFilter_;

  TTree* jetTrees_p[nJetAlgos];
  TTree* skimTree_p = NULL;
  TTree* hltTree_p = NULL;
    
  const Float_t jtPtLowPP = 110.;
  const Float_t jtPtLowPbPb = 140.;
  Float_t jtPtLowTemp = jtPtLowPP;
  if(!isPP) jtPtLowTemp = jtPtLowPbPb;
  const Float_t jtPtLow = jtPtLowTemp;

  const Double_t jtEtaMax = 2.0;
  
  outFile_p->cd();

  std::map<ULong64_t, Int_t> runLumi_MBCounts;
  std::map<ULong64_t, Int_t> runLumi_JetCounts;
  std::map<ULong64_t, Int_t> runLumi_MBAndJetCounts;
  Int_t Run_;
  Int_t LumiBlock_;
  Int_t HLT_HIL1MinimumBiasHF2AND_v;
  Int_t HLT_HIL1MinimumBiasHF2AND_part1_v;
  Int_t HLT_HIPuAK4CaloJet100_Eta5p1_v;

  UInt_t run_, lumi_;


  Int_t hiBinOut_[nJetAlgos];
  Int_t nrefOut_[nJetAlgos];
  Float_t jtptOut_[nJetAlgos][nMaxJets];
  Float_t rawptOut_[nJetAlgos][nMaxJets];
  Float_t jtetaOut_[nJetAlgos][nMaxJets];
  Float_t jtphiOut_[nJetAlgos][nMaxJets];
  Float_t jtPfCHFOut_[nJetAlgos][nMaxJets];
  Float_t jtPfCEFOut_[nJetAlgos][nMaxJets];
  Float_t jtPfNHFOut_[nJetAlgos][nMaxJets];
  Float_t jtPfNEFOut_[nJetAlgos][nMaxJets];
  Float_t jtPfMUFOut_[nJetAlgos][nMaxJets];

  TTree* jetTreesOut_p[nJetAlgosPlus];

  TDirectory* algoDir_p[nJetAlgosPlus];

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(Int_t jI = 0; jI < nJetAlgosPlus; ++jI){
    std::string algoString = inJetTreeList_.at(jI);
    algoString = algoString.substr(0, algoString.find("JetAn"));

    outFile_p->cd();
    algoDir_p[jI] = outFile_p->GetDirectory(algoString.c_str());
    if(algoDir_p[jI]) algoDir_p[jI]->cd();
    else{
      algoDir_p[jI] = outFile_p->mkdir(algoString.c_str());
      algoDir_p[jI]->cd();
    }

    std::string jetTreeName = inJetTreeList_.at(jI);
    jetTreeName.replace(jetTreeName.find("/"), jetTreeName.size() - jetTreeName.find("/"), "");
    jetTreeName = jetTreeName + "OUT";
    std::cout << "JetTreeName: " << jetTreeName << std::endl;
    jetTreesOut_p[jI] = NULL;
    jetTreesOut_p[jI] = new TTree(jetTreeName.c_str(), jetTreeName.c_str());
    
    jetTreesOut_p[jI]->Branch("hiBin", &(hiBinOut_[jI]), "hiBin/I");
    jetTreesOut_p[jI]->Branch("nref", &(nrefOut_[jI]), "nref/I");
    jetTreesOut_p[jI]->Branch("jtpt", (jtptOut_[jI]), "jtpt[nref]/F");
    jetTreesOut_p[jI]->Branch("rawpt", (rawptOut_[jI]), "rawpt[nref]/F");
    jetTreesOut_p[jI]->Branch("jtphi", (jtphiOut_[jI]), "jtphi[nref]/F");
    jetTreesOut_p[jI]->Branch("jteta", (jtetaOut_[jI]), "jteta[nref]/F");
    jetTreesOut_p[jI]->Branch("jtPfCHF", (jtPfCHFOut_[jI]), "jtPfCHF[nref]/F");
    jetTreesOut_p[jI]->Branch("jtPfCEF", (jtPfCEFOut_[jI]), "jtPfCEF[nref]/F");
    jetTreesOut_p[jI]->Branch("jtPfNHF", (jtPfNHFOut_[jI]), "jtPfNHF[nref]/F");
    jetTreesOut_p[jI]->Branch("jtPfNEF", (jtPfNEFOut_[jI]), "jtPfNEF[nref]/F");
    jetTreesOut_p[jI]->Branch("jtPfMUF", (jtPfMUFOut_[jI]), "jtPfMUF[nref]/F");
  }

  std::cout << "Procesing files..." << std::endl;
  for(unsigned int fileIter = 0; fileIter < fileList.size(); ++fileIter){
    if(fileIter%100 == 0){
      std::cout << " File: " << fileIter << "/" << fileList.size() << std::endl;
      std::cout << "  \'" << fileList.at(fileIter) << "\'" << std::endl;
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    inFile_p = NULL;
    inFile_p = TFile::Open(fileList.at(fileIter).c_str(), "READ");
    if(inFile_p == NULL){
      std::cout << "Continuing on null \'" << fileList.at(fileIter) << "\'" << std::endl;
      continue;
    }

    if(inFile_p->IsZombie()){
      std::cout << "Continuing on zombie \'" << fileList.at(fileIter) << "\'" << std::endl;

      inFile_p->Close();
      delete inFile_p;
      continue;
    }

    skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");
    if(!isPP && hasHLT) hltTree_p = (TTree*)inFile_p->Get("hltanalysis/HltTree");
    hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
  
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    skimTree_p->SetBranchStatus("*", 0);
    skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    if(!isPP){
      skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
      skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
      skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);
    }
    else{
      skimTree_p->SetBranchStatus("pBeamScrapingFilter", 1);
      skimTree_p->SetBranchStatus("pPAprimaryVertexFilter", 1);
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
    if(!isPP){
      skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
      skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
      skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);
    }
    else{
      skimTree_p->SetBranchAddress("pPAprimaryVertexFilter", &pprimaryVertexFilter_);
      skimTree_p->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter_);
    }


    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(!isPP && hasHLT){
      hltTree_p->SetBranchStatus("*", 0);
      hltTree_p->SetBranchStatus("Run", 1);
      hltTree_p->SetBranchStatus("LumiBlock", 1);
      
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
      if(hltTree_p->GetListOfBranches()->FindObject("HLT_HIL1MinimumBiasHF2AND_v1") == 0){
	std::cout << "Missing branch \'" << fileList.at(fileIter) << "\'" << std::endl;
      }
      
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      
      if(!isPP && hasHLT){
	hltTree_p->SetBranchStatus("HLT_HIL1MinimumBiasHF2AND_v1", 1);
	hltTree_p->SetBranchStatus("HLT_HIL1MinimumBiasHF2AND_part1_v1", 1);
	hltTree_p->SetBranchStatus("HLT_HIPuAK4CaloJet100_Eta5p1_v1", 1);
      }
      
      hltTree_p->SetBranchAddress("Run", &Run_);
      hltTree_p->SetBranchAddress("LumiBlock", &LumiBlock_);
      
      if(!isPP && hasHLT){
	hltTree_p->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_v1", &HLT_HIL1MinimumBiasHF2AND_v);
	hltTree_p->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part1_v1", &HLT_HIL1MinimumBiasHF2AND_part1_v);
	hltTree_p->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1", &HLT_HIPuAK4CaloJet100_Eta5p1_v);
      }
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("lumi", 1);

    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);

    for(int treeIter = 0; treeIter < nJetAlgos; ++treeIter){
      //      std::cout << "Tree " << treeIter << "/" << inJetTreeList_.size() << ": " << inJetTreeList_.at(treeIter) << std::endl;
      jetTrees_p[treeIter] = (TTree*)inFile_p->Get(inJetTreeList_.at(treeIter).c_str());
      //      std::cout << "Tree " << treeIter << "/" << nJetAlgos << ": " << inJetTreeList_.at(treeIter) << std::endl;

      jetTrees_p[treeIter]->SetBranchStatus("*", 0);
      jetTrees_p[treeIter]->SetBranchStatus("nref", 1);
      jetTrees_p[treeIter]->SetBranchStatus("jtpt", 1);
      jetTrees_p[treeIter]->SetBranchStatus("rawpt", 1);
      jetTrees_p[treeIter]->SetBranchStatus("jteta", 1);
      jetTrees_p[treeIter]->SetBranchStatus("jtphi", 1);
      jetTrees_p[treeIter]->SetBranchStatus("jtPfCHF", 1);
      jetTrees_p[treeIter]->SetBranchStatus("jtPfCEF", 1);
      jetTrees_p[treeIter]->SetBranchStatus("jtPfNHF", 1);
      jetTrees_p[treeIter]->SetBranchStatus("jtPfNEF", 1);
      jetTrees_p[treeIter]->SetBranchStatus("jtPfMUF", 1);
      
      jetTrees_p[treeIter]->SetBranchAddress("nref", &nref_[treeIter]);
      jetTrees_p[treeIter]->SetBranchAddress("jtpt", jtpt_[treeIter]);
      jetTrees_p[treeIter]->SetBranchAddress("rawpt", rawpt_[treeIter]);
      jetTrees_p[treeIter]->SetBranchAddress("jteta", jteta_[treeIter]);
      jetTrees_p[treeIter]->SetBranchAddress("jtphi", jtphi_[treeIter]);
      jetTrees_p[treeIter]->SetBranchAddress("jtPfCHF", jtPfCHF_[treeIter]);
      jetTrees_p[treeIter]->SetBranchAddress("jtPfCEF", jtPfCEF_[treeIter]);
      jetTrees_p[treeIter]->SetBranchAddress("jtPfNHF", jtPfNHF_[treeIter]);
      jetTrees_p[treeIter]->SetBranchAddress("jtPfNEF", jtPfNEF_[treeIter]);
      jetTrees_p[treeIter]->SetBranchAddress("jtPfMUF", jtPfMUF_[treeIter]);
    }

    inFile_p->cd();
    const Int_t nEntries = hiTree_p->GetEntries();
  
    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(entry%10000 == 0 && fileList.size() == 1) std::cout << "Entry " << entry << "/" << nEntries << std::endl;
      
      hiTree_p->GetEntry(entry);
      skimTree_p->GetEntry(entry);
      if(!isPP && hasHLT) hltTree_p->GetEntry(entry);

      if(!pprimaryVertexFilter_) continue;
      if(!HBHENoiseFilterResultRun2Loose_) continue;
      if(TMath::Abs(vz_) > 15.) continue;
      
      if(isPP && !pBeamScrapingFilter_) continue;
      if(!isPP && !phfCoincFilter3_) continue;
      if(!isPP && !pclusterCompatibilityFilter_) continue;
      
    
      if(!isPP && hasHLT){
	ULong64_t keyRun = 1000000*(ULong64_t(Run_));
	ULong64_t keyLumi = (ULong64_t(LumiBlock_));
	ULong64_t key = keyRun + keyLumi;
	
	if(Run_ <= 263153){
	  if(runLumi_MBCounts.count(key) == 0){
	    runLumi_MBCounts[key] = 0;
	    runLumi_JetCounts[key] = 0;
	    runLumi_MBAndJetCounts[key] = 0;
	  }
	  if(!isPP){
	    if(HLT_HIL1MinimumBiasHF2AND_v) runLumi_MBCounts[key] += 1;
	    if(HLT_HIPuAK4CaloJet100_Eta5p1_v) runLumi_JetCounts[key] += 1;
	    if(HLT_HIL1MinimumBiasHF2AND_v && HLT_HIPuAK4CaloJet100_Eta5p1_v) runLumi_MBAndJetCounts[key] += 1;
	  }
	}
	else if(Run_ <= 263797){
	  if(runLumi_MBCounts.count(key) == 0){
	    runLumi_MBCounts[key] = 0;
	    runLumi_JetCounts[key] = 0;
	    runLumi_MBAndJetCounts[key] = 0;
	  }
	  if(!isPP){
	    if(HLT_HIL1MinimumBiasHF2AND_part1_v) runLumi_MBCounts[key] += 1;
	    if(HLT_HIPuAK4CaloJet100_Eta5p1_v) runLumi_JetCounts[key] += 1;
	    if(HLT_HIL1MinimumBiasHF2AND_part1_v && HLT_HIPuAK4CaloJet100_Eta5p1_v) runLumi_MBAndJetCounts[key] += 1;
	  }
	}

	runLumiMap[key] = 0;	
      }
      else if(!isPP && !hasHLT){
	ULong64_t keyRun = 1000000*(ULong64_t(run_));
        ULong64_t keyLumi = (ULong64_t(lumi_));
        ULong64_t key = keyRun + keyLumi;
	
	runLumiMap[key] = 0;
      }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      
      for(Int_t jI = 0; jI < nJetAlgosPlus; ++jI){
	std::string algoString = inJetTreeList_.at(jI);
	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	algoString = algoString.substr(0, algoString.find("JetAn"));

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	//	if(algoString.find("4") == std::string::npos) continue;
	
      	if(jI < nJetAlgos) jetTrees_p[jI]->GetEntry(entry);

	hiBinOut_[jI] = hiBin_;
	nrefOut_[jI] = 0;

	if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	
	for(Int_t jtIter = 0; jtIter < nref_[algoPos.at(jI)]; ++jtIter){
	  if(TMath::Abs(jteta_[algoPos.at(jI)][jtIter]) > jtEtaMax) continue;
	  
	  std::vector<int> absEtaPoses;

	  if(jtpt_[algoPos.at(jI)][jtIter] >= jtPtLow){
	
	    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
	    jtptOut_[jI][nrefOut_[jI]] = jtpt_[algoPos.at(jI)][jtIter];
	    rawptOut_[jI][nrefOut_[jI]] = rawpt_[algoPos.at(jI)][jtIter];
	    jtetaOut_[jI][nrefOut_[jI]] = jteta_[algoPos.at(jI)][jtIter];
	    jtphiOut_[jI][nrefOut_[jI]] = jtphi_[algoPos.at(jI)][jtIter];
	    jtPfCHFOut_[jI][nrefOut_[jI]] = jtPfCHF_[algoPos.at(jI)][jtIter];
	    jtPfCEFOut_[jI][nrefOut_[jI]] = jtPfCEF_[algoPos.at(jI)][jtIter];
	    jtPfNHFOut_[jI][nrefOut_[jI]] = jtPfNHF_[algoPos.at(jI)][jtIter];
	    jtPfNEFOut_[jI][nrefOut_[jI]] = jtPfNEF_[algoPos.at(jI)][jtIter];
	    jtPfMUFOut_[jI][nrefOut_[jI]] = jtPfMUF_[algoPos.at(jI)][jtIter];
	    nrefOut_[jI] += 1;
	  }
	}	       

	jetTreesOut_p[jI]->Fill();
      }
    }

    inFile_p->Close();
    delete inFile_p;
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


  outFile_p->cd();

  std::string jtEtaMaxStr = std::to_string(jtEtaMax);

  TNamed jtEtaMaxName("jtEtaMax", jtEtaMaxStr.c_str());
  jtEtaMaxName.Write("", TObject::kOverwrite);

  for(Int_t jI = 0; jI < nJetAlgosPlus; ++jI){
    outFile_p->cd();
    algoDir_p[jI]->cd();
      
    jetTreesOut_p[jI]->Write("", TObject::kOverwrite);
    delete jetTreesOut_p[jI];
  }
  
  outFile_p->cd();

  if(!isPP){
    Int_t mbCounts = 0;
    Int_t jetCounts = 0;
    Int_t mbAndJetCounts = 0;
    
    TTree* countTree_p = new TTree("countTree", "countTree");
    countTree_p->Branch("Run", &Run_, "Run/I");
    countTree_p->Branch("LumiBlock", &LumiBlock_, "LumiBlock/I");

    if(hasHLT){
      countTree_p->Branch("mbCounts", &mbCounts, "mbCounts/I");
      countTree_p->Branch("jetCounts", &jetCounts, "jetCounts/I");
      countTree_p->Branch("mbAndJetCounts", &mbAndJetCounts, "mbAndJetCounts/I");
    }

    std::map<ULong64_t, Int_t>::iterator it;
    for(it = runLumiMap.begin(); it != runLumiMap.end(); it++){
      //    std::cout << "Key: " << it->first << std::endl;
      
      Run_ = it->first/1000000;
      LumiBlock_ = it->first%1000000;

      if(hasHLT){
	mbCounts = runLumi_MBCounts[it->first];
	jetCounts = runLumi_JetCounts[it->first];
	mbAndJetCounts = runLumi_MBAndJetCounts[it->first];
      }

      countTree_p->Fill();
    }
    
    countTree_p->Write("", TObject::kOverwrite);
    delete countTree_p;
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage: ./makeRawRAATree.exe <inFileName> <outFileName-optional>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += makeRawRAATree(argv[1]);
  else if(argc == 3) retVal += makeRawRAATree(argv[1], argv[2]);

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
