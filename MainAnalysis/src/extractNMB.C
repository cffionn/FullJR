//cpp dependencies
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/returnRootFileContentsList.h"

//https://twiki.cern.ch/twiki/pub/CMS/HINUpsilonRaa2016/Jason_MinBiasCounting_2017-02-02.pdf

ULong64_t getRunLumiKey(ULong64_t run, ULong64_t lumi){return run*100000 + lumi;}
ULong64_t runLumiKeyToRun(ULong64_t key){return key/100000;}
ULong64_t runLumiKeyToLumi(ULong64_t key){return key%100000;}

bool fileValid(std::string fileStr, std::string fileNameStr, std::string extStr)
{
  if(!checkFile(fileNameStr) || fileNameStr.find(extStr) == std::string::npos){
    std::cout << fileStr << " \'" << fileNameStr << "\' is not valid. return 1" << std::endl;
    return false;
  }
  return true;
}

int extractNMB(const std::string inJSONFile, const std::string inROOTMB1File, const std::string inROOTMB2File, const std::string inROOTHPFile)
{
  if(!fileValid("inJSONFile", inJSONFile, ".txt")) return 1;
  if(!fileValid("inROOTMB1File", inROOTMB1File, ".root")) return 1;
  if(!fileValid("inROOTMB2File", inROOTMB2File, ".root")) return 1;
  if(!fileValid("inROOTHPFile", inROOTHPFile, ".root")) return 1;
    
  std::vector<std::string> rootFileNames = {inROOTMB1File, inROOTMB2File, inROOTHPFile};
  const std::string hltTreeStr = "hltanalysis/HltTree";
  const std::string hiTreeStr = "hiEvtAnalyzer/HiTree";
  const std::string skimTreeStr = "skimanalysis/HltTree";
  ULong64_t runSplit = 263153; //all runs before this use HLT_HIL1MinimumBiasHF2AND_v1; all runs after use HLT_HIL1MinimumBiasHF2AND_part1_v1
  std::vector<ULong64_t> earlyRunMBCount = {0, 0, 0};
  std::vector<ULong64_t> lateRunMBCount = {0, 0, 0};
  std::vector<ULong64_t> fullRunJet100Count = {0, 0, 0};
  std::vector<ULong64_t> fullRunJet100AndMBCount = {0, 0, 0};

  const UInt_t nCentBins = 6;
  const Int_t centBinsLow[nCentBins] = {0, 0, 0, 10, 30, 50};
  const Int_t centBinsHigh[nCentBins] = {100, 90, 10, 30, 50, 90};

  for(UInt_t fI = 0; fI < rootFileNames.size(); ++fI){
    for(UInt_t cI = 0; cI < nCentBins; ++cI){
      earlyRunMBCount.push_back(0);
      lateRunMBCount.push_back(0);
      fullRunJet100Count.push_back(0);
      fullRunJet100AndMBCount.push_back(0);
    }
  }

  for(auto const & name : rootFileNames){
    TFile* inFile_p = new TFile(name.c_str(), "READ");
    std::vector<std::string> treeList = returnRootFileContentsList(inFile_p, "TTree", ""); 
    inFile_p->Close();
    delete inFile_p;

    bool hasHltTree = false;
    bool hasHiTree = false;
    bool hasSkimTree = false;
    for(auto const & val : treeList){
      if(isStrSame(val, hltTreeStr)) hasHltTree = true;
      else if(isStrSame(val, hiTreeStr)) hasHiTree = true;
      else if(isStrSame(val, skimTreeStr)) hasSkimTree = true;
    }

    if(!hasHltTree){
      std::cout << "File \'" << name << "\' is missing ttree \'" << hltTreeStr << "\'. return 1" << std::endl;
      return 1;
    }
    else if(!hasHiTree){
      std::cout << "File \'" << name << "\' is missing ttree \'" << hiTreeStr << "\'. return 1" << std::endl;
      return 1;
    }
    else if(!hasSkimTree){
      std::cout << "File \'" << name << "\' is missing ttree \'" << skimTreeStr << "\'. return 1" << std::endl;
      return 1;
    }
  }

  std::map<ULong64_t, ULong64_t> jsonMapRunToLumiVect;

  std::ifstream file(inJSONFile.c_str());
  std::string tempStr;
  std::vector<std::string> strToRep = {"\"", "{", "}", " "};
  
  while(std::getline(file, tempStr)){
    for(auto const & val : strToRep){
      while(tempStr.find(val) != std::string::npos){tempStr.replace(tempStr.find(val), 1, "");}
    }

    while(tempStr.find("]]") != std::string::npos){
      std::string tempStr2 = tempStr.substr(0, tempStr.find("]]")+2);
      tempStr.replace(0, tempStr.find("]]")+3, "");

      ULong64_t runNum = std::stol(tempStr2.substr(0, tempStr2.find(":")));
      tempStr2.replace(0, tempStr2.find(":")+1, "");
      
      while(tempStr2.find("[") != std::string::npos){
	while(tempStr2.substr(0,1).find("[") != std::string::npos){tempStr2.replace(0,1,"");}
	ULong64_t num1 = std::stol(tempStr2.substr(0, tempStr2.find(",")));
	tempStr2.replace(0, tempStr2.find(",")+1, "");
	ULong64_t num2  = std::stol(tempStr2.substr(0, tempStr2.find("]")));

	for(ULong64_t i = num1; i <= num2; ++i){
	  ULong64_t key = getRunLumiKey(runNum, i);
	  jsonMapRunToLumiVect[key] = 0;
	}

	if(tempStr2.find("[") != std::string::npos){
	  while(tempStr2.substr(0,1).find("[") == std::string::npos){
	    tempStr2.replace(0, 1, "");
	  }
	}
      }
    }
  }
  file.close();

  Int_t RunHLT, LumiBlockHLT, HLT_HIL1MinimumBiasHF2AND_v1, HLT_HIL1MinimumBiasHF2AND_part1_v1, HLT_HIPuAK4CaloJet100_Eta5p1_v1, hiBin_, HBHENoiseFilterResultRun2Loose_, pprimaryVertexFilter_, phfCoincFilter3_, pclusterCompatibilityFilter_;
  Int_t pBeamScrapingFilter_ = -1;

  Float_t vz_, hiHF_; 
  
  goodGlobalSelection globalSel;
  globalSel.setIsPbPb(true);
  
  for(ULong64_t fI = 0; fI < rootFileNames.size(); ++fI){
    std::string name = rootFileNames[fI];
    TFile* inFile_p = new TFile(name.c_str(), "READ");
    TTree* hltTree_p = (TTree*)inFile_p->Get(hltTreeStr.c_str());
    TTree* hiTree_p = (TTree*)inFile_p->Get(hiTreeStr.c_str());
    TTree* skimTree_p = (TTree*)inFile_p->Get(skimTreeStr.c_str());
    
    hltTree_p->SetBranchStatus("*", 0);
    hltTree_p->SetBranchStatus("LumiBlock", 1);
    hltTree_p->SetBranchStatus("Run", 1);
    hltTree_p->SetBranchStatus("HLT_HIL1MinimumBiasHF2AND_v1", 1);
    hltTree_p->SetBranchStatus("HLT_HIL1MinimumBiasHF2AND_part1_v1", 1);
    hltTree_p->SetBranchStatus("HLT_HIPuAK4CaloJet100_Eta5p1_v1", 1);
    
    hltTree_p->SetBranchAddress("LumiBlock", &LumiBlockHLT);
    hltTree_p->SetBranchAddress("Run", &RunHLT);
    hltTree_p->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_v1", &HLT_HIL1MinimumBiasHF2AND_v1);
    hltTree_p->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_part1_v1", &HLT_HIL1MinimumBiasHF2AND_part1_v1);
    hltTree_p->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1", &HLT_HIPuAK4CaloJet100_Eta5p1_v1);

    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("hiHF", 1);
    hiTree_p->SetBranchStatus("hiBin", 1);

    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("hiHF", &hiHF_);
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    
    skimTree_p->SetBranchStatus("*", 0);
    skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
    skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
    skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);

    skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
    skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
    skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
    skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);

    const ULong64_t nEntries = hltTree_p->GetEntries();
    const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);
    
    std::cout << "Processing " << nEntries << " in file \'" << name << "\'..." << std::endl;
    for(ULong64_t entry = 0; entry < nEntries; ++entry){
      if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
      hltTree_p->GetEntry(entry);
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

      ULong64_t runMain = (ULong64_t)RunHLT;
      ULong64_t lumiMain = (ULong64_t)LumiBlockHLT;
      ULong64_t key = getRunLumiKey(runMain, lumiMain);
      if(jsonMapRunToLumiVect.count(key) == 0){
	std::cout << "WARNING: file \'" << name << "\' is using run/lumi (" << runMain << "/" << lumiMain << ") not found in json..." << std::endl;
      }

      std::vector<UInt_t> centPos;
      for(UInt_t cI = 0; cI < nCentBins; ++cI){
	if(hiBin_/2 >= centBinsLow[cI] && hiBin_/2 < centBinsHigh[cI]){
	  centPos.push_back(cI);
	}
      }
      
      for(auto const & cent : centPos){
	if(runMain <= runSplit){
	  if(HLT_HIL1MinimumBiasHF2AND_v1) ++(earlyRunMBCount[fI*nCentBins + cent]);
	}
	else{
	  if(HLT_HIL1MinimumBiasHF2AND_part1_v1) ++(lateRunMBCount[fI*nCentBins + cent]);
	}
	
	if(HLT_HIPuAK4CaloJet100_Eta5p1_v1){
	  ++(fullRunJet100Count[fI*nCentBins + cent]);
	  
	  if(runMain <= runSplit){
	    if(HLT_HIL1MinimumBiasHF2AND_v1) ++(fullRunJet100AndMBCount[fI*nCentBins + cent]);
	  }
	  else{
	    if(HLT_HIL1MinimumBiasHF2AND_part1_v1) ++(fullRunJet100AndMBCount[fI*nCentBins + cent]);
	  }
	}
      }
    }
    
    inFile_p->Close();
    delete inFile_p;    
  }

  for(ULong64_t fI = 0; fI < rootFileNames.size(); ++fI){
    std::cout << "Counts from \'" << rootFileNames[fI] << "\': " << std::endl;

    for(ULong64_t cI = 0; cI < nCentBins; ++cI){
      std::cout << " Cent bins: " << centBinsLow[cI] << "-" << centBinsHigh[cI] << std::endl;
      std::cout << "  Early run MB Counts: " << earlyRunMBCount[fI*nCentBins + cI] << std::endl;
      std::cout << "  Late run MB Counts: " << lateRunMBCount[fI*nCentBins + cI] << std::endl;
      std::cout << "  Full run Jet100 Counts: " << fullRunJet100Count[fI*nCentBins + cI] << std::endl;
      std::cout << "  Full run Jet100 And MB Counts: " << fullRunJet100AndMBCount[fI*nCentBins + cI] << std::endl;  
    }
  }  

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 5){
    std::cout << "Usage: ./bin/extractNMB.exe <inJSONFile> <inROOTMB1File> <inROOTMB2File> <inROOTHPFile>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += extractNMB(argv[1], argv[2], argv[3], argv[4]);
  return retVal;
}
