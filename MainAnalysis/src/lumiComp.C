//cpp dependencies
#include <algorithm>
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
#include "Utility/include/fileUtilities.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"

//https://twiki.cern.ch/twiki/pub/CMS/HINUpsilonRaa2016/Jason_MinBiasCounting_2017-02-02.pdf

ULong64_t getRunLumiKey(ULong64_t run, ULong64_t lumi){return run*100000 + lumi;}
ULong64_t runLumiKeyToRun(ULong64_t key){return key/100000;}
ULong64_t runLumiKeyToLumi(ULong64_t key){return key%100000;}

int lumiComp(const std::string inJSONFile, const std::string inROOTFile)
{
  if(!fileIsGood(inJSONFile, ".txt")) return 1;
  if(!fileIsGood(inROOTFile, ".root")) return 1;

  TFile* inFile_p = new TFile(inROOTFile.c_str(), "READ");
  std::vector<std::string> treeList = returnRootFileContentsList(inFile_p, "TTree", ""); 
  inFile_p->Close();
  delete inFile_p;

  bool hasHiTree = false;
  const std::string hiTreeStr = "hiEvtAnalyzer/HiTree";
  bool hasHiTreeAlt = false;
  const std::string hiTreeStrAlt = "HiEvt";
  bool hasHltTree = false;
  const std::string hltTreeStr = "hltanalysis/HltTree";

  for(auto const & val : treeList){
    //    std::cout << val << std::endl;
    if(isStrSame(val, hiTreeStr)) hasHiTree = true;
    else if(isStrSame(val, hiTreeStrAlt)) hasHiTreeAlt = true;
    else if(isStrSame(val, hltTreeStr)) hasHltTree = true;
  }

  if(!hasHiTree && !hasHiTreeAlt && !hasHltTree){
    std::cout << "Missing both \'" << hiTreeStr << "\' and \'" << hiTreeStrAlt << "\' and \'" << hltTreeStr << "\'. return 1" << std::endl;
    return 1;
  }

  std::map<ULong64_t, ULong64_t> jsonMapRunToLumiVect;
  std::map<ULong64_t, ULong64_t> jsonMapRunToTotalLumi;
  std::map<ULong64_t, std::vector<ULong64_t> > jsonMapRunToAllLumi;
  std::map<ULong64_t, std::vector<ULong64_t> > jsonMapRunToFoundLumi;
  std::map<ULong64_t, std::vector<ULong64_t> > jsonMapRunToMissingLumi;

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
      jsonMapRunToAllLumi[runNum] = {};

      while(tempStr2.find("[") != std::string::npos){
	while(tempStr2.substr(0,1).find("[") != std::string::npos){tempStr2.replace(0,1,"");}
	ULong64_t num1 = std::stol(tempStr2.substr(0, tempStr2.find(",")));
	tempStr2.replace(0, tempStr2.find(",")+1, "");
	ULong64_t num2  = std::stol(tempStr2.substr(0, tempStr2.find("]")));
	
	for(ULong64_t i = num1; i <= num2; ++i){
	  ULong64_t key = getRunLumiKey(runNum, i);
	  jsonMapRunToLumiVect[key] = 0;
	  ++(jsonMapRunToTotalLumi[runNum]);
	  jsonMapRunToAllLumi[runNum].push_back(i);
	}

	jsonMapRunToFoundLumi[runNum] = {};
	jsonMapRunToMissingLumi[runNum] = {};

	if(tempStr2.find("[") != std::string::npos){
	  while(tempStr2.substr(0,1).find("[") == std::string::npos){
	    tempStr2.replace(0, 1, "");
	  }
	}
      }
    }
  }
  file.close();

  UInt_t runHI;
  UInt_t lumiHI;

  Int_t LumiBlockHLT;
  Int_t RunHLT;
  
  inFile_p = new TFile(inROOTFile.c_str(), "READ");
  TTree* inTree_p = NULL;
  if(hasHiTree || hasHiTreeAlt){
    if(hasHiTree) inTree_p = (TTree*)inFile_p->Get(hiTreeStr.c_str());
    else if(hasHiTreeAlt) inTree_p = (TTree*)inFile_p->Get(hiTreeStrAlt.c_str());

    inTree_p->SetBranchStatus("*", 0);
    inTree_p->SetBranchStatus("run", 1);
    inTree_p->SetBranchStatus("lumi", 1);

    inTree_p->SetBranchAddress("run", &runHI);
    inTree_p->SetBranchAddress("lumi", &lumiHI);
  }
  else if(hasHltTree){
    inTree_p = (TTree*)inFile_p->Get(hltTreeStr.c_str());

    inTree_p->SetBranchStatus("*", 0);
    inTree_p->SetBranchStatus("LumiBlock", 1);
    inTree_p->SetBranchStatus("Run", 1);

    inTree_p->SetBranchAddress("LumiBlock", &LumiBlockHLT);
    inTree_p->SetBranchAddress("Run", &RunHLT);
  }

  const ULong64_t nEntries = inTree_p->GetEntries();
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

  std::cout << "Processing " << nEntries << "..." << std::endl;
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    inTree_p->GetEntry(entry);

    ULong64_t runMain, lumiMain;

    if(hasHiTree || hasHiTreeAlt){
      runMain = (ULong64_t)runHI;
      lumiMain = (ULong64_t)lumiHI;      
    }
    else if(hasHltTree){
      runMain = (ULong64_t)RunHLT;
      lumiMain = (ULong64_t)LumiBlockHLT;
    }

    ULong64_t key = getRunLumiKey(runMain, lumiMain);
    ++(jsonMapRunToLumiVect[key]);

    std::vector<ULong64_t>* tempVect_p = &(jsonMapRunToFoundLumi[runMain]);

    bool lumiFound = false;
    for(unsigned int tI = 0; tI < tempVect_p->size(); ++tI){
      if((*tempVect_p)[tI] == lumiMain){
	lumiFound = true;
	break;
      }
    }
    if(!lumiFound) jsonMapRunToFoundLumi[runMain].push_back(lumiMain);
  }

  inFile_p->Close();
  delete inFile_p;

  Double_t none = 0;
  Double_t total = 0;

  for(auto const & val : jsonMapRunToLumiVect){
    if(val.second == 0) ++none;
    ++total;
  }

  for(auto const & val : jsonMapRunToAllLumi){
    std::vector<ULong64_t>* temp = &(jsonMapRunToAllLumi[val.first]);
    std::sort(temp->begin(), temp->end());
  }

  for(auto const & val : jsonMapRunToFoundLumi){
    std::vector<ULong64_t>* temp = &(jsonMapRunToFoundLumi[val.first]);
    std::sort(temp->begin(), temp->end());
  }

  for(auto const & val : jsonMapRunToFoundLumi){
    for(auto const & a : jsonMapRunToAllLumi[val.first]){
      bool isFound = false;
      for(auto const & f : val.second){
	if(f == a){
	  isFound = true;
	  break;
	}
      }
      if(!isFound) jsonMapRunToMissingLumi[val.first].push_back(a);
    }
  }

  std::cout << "Fraction: " << ((ULong64_t)none) << "/" << ((ULong64_t)total) << "=" << none/total << std::endl;
  for(auto const & val : jsonMapRunToTotalLumi){
    std::cout << " " << val.first << ": " << jsonMapRunToFoundLumi[val.first].size() << "/" << val.second << "=" << prettyString(((double)jsonMapRunToFoundLumi[val.first].size())/(double)val.second, 3 , false) << std::endl;
  }

  std::string jsonMissStr = "";
  for(auto const & val : jsonMapRunToMissingLumi){
    std::string lumiStr = "";
    ULong64_t frontNum = 0;
    ULong64_t prevNum = 0;
    ULong64_t currNum = 0;
    
    for(auto const & lumi : val.second){
      currNum = lumi;

      if(frontNum == 0){
	frontNum = currNum;	
	prevNum = currNum;
      }
      else{
	if(currNum - prevNum != 1 ){
	  lumiStr = lumiStr + "[" + std::to_string(frontNum) + ", " + std::to_string(prevNum) + "], ";
	  frontNum = currNum;
	}
	prevNum = currNum;
      }
    }

    if(frontNum != 0){
      lumiStr = "[" + lumiStr + "[" + std::to_string(frontNum) + ", " + std::to_string(currNum) + "]]";
      jsonMissStr = jsonMissStr + "\"" + std::to_string(val.first) + "\": " + lumiStr + ", ";
    }

    std::cout << val.first << ": " << lumiStr << std::endl;    
  }
  if(jsonMissStr.size() != 0) jsonMissStr = "{" + jsonMissStr.substr(0, jsonMissStr.rfind(",")) + "}";

  std::cout << jsonMissStr << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/lumiComp.exe <inJSONFile> <inROOTFile>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += lumiComp(argv[1], argv[2]);
  return retVal;
}
