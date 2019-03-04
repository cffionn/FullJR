//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TDatime.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TTree.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"

int splitFiles(const std::string inFileName, const int splitLevel)
{ 
  if(!checkFile(inFileName) || inFileName.find(".root") == std::string::npos){
    std::cout << "inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> listOfTrees = returnRootFileContentsList(inFile_p, "TTree", "");
  std::vector<std::string> listOfDirs = returnRootFileContentsList(inFile_p, "TDirectory", "");
  std::vector<std::string> listOfDirs2 = returnRootFileContentsList(inFile_p, "TDirectoryFile", "");
  listOfDirs.insert(std::end(listOfDirs), std::begin(listOfDirs2), std::end(listOfDirs2));
  
  for(unsigned int dI = 0; dI < listOfDirs.size(); ++dI){
    for(unsigned int dI2 = dI+1; dI2 < listOfDirs.size(); ++dI2){
      if(listOfDirs[dI].size() > listOfDirs[dI2].size()){
	std::string tempStr = listOfDirs[dI];
	listOfDirs[dI] = listOfDirs[dI2];
	listOfDirs[dI2] = tempStr;
      }
    }
  }
  
  if(listOfTrees.size() == 0){
    std::cout << "Input file \'" << inFileName << "\' has no TTrees. return 1" << std::endl;
    return 1;
  }

  for(unsigned int tI = 0; tI < listOfTrees.size(); ++tI){
    std::cout << listOfTrees[tI] << std::endl;
  }

  for(unsigned int tI = 0; tI < listOfDirs.size(); ++tI){
    std::cout << listOfDirs[tI] << std::endl;
  }

  std::vector<TTree*> trees_p;
  std::vector<TTree*> outTrees_p;
  std::vector<TDirectoryFile*> outDirs_p;

  TFile* outFile_p = NULL;
  trees_p.reserve(listOfTrees.size());
  outTrees_p.reserve(listOfTrees.size());
  outDirs_p.reserve(listOfDirs.size());

  for(unsigned int dI = 0; dI < listOfDirs.size(); ++dI){outDirs_p.push_back(NULL);}
  
  ULong64_t nEntries = 0;
  
  for(unsigned int lI = 0; lI < listOfTrees.size(); ++lI){
    trees_p.push_back((TTree*)inFile_p->Get(listOfTrees[lI].c_str()));    
    outTrees_p.push_back(NULL);
    
    if(lI == 0) nEntries = trees_p[lI]->GetEntries();
    else if(((ULong64_t)trees_p[lI]->GetEntries()) != nEntries){
      std::cout << "nEntries " << nEntries << " not matching to current tree " << trees_p[lI]->GetEntries() << ". return 1" << std::endl;
      return 1;
    }
  }

  if(splitLevel <= 0 || ((ULong64_t)splitLevel) > nEntries){
    std::cout << "Split level \'" << splitLevel << "\' is invalid (nEntries=" << splitLevel << "). return 1." << std::endl;
    return 1;
  }

  std::string outFileName = inFileName;
  outFileName.replace(outFileName.find(".root"), 5, "");
  outFileName = outFileName + "_SPLIT";

  ULong64_t nSplitEntries = nEntries/splitLevel;
  ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);
  
  std::cout << "Processing trees...." << std::endl;
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;

    if(entry%nSplitEntries == 0 && entry/nSplitEntries != ((ULong64_t)splitLevel)){
      if(entry != 0){
	outFile_p->cd();

	for(unsigned int oI = 0; oI < outTrees_p.size(); ++oI){
	  outFile_p->cd();
	  for(unsigned int dI = 0; dI < listOfDirs.size(); ++dI){
	    if(listOfTrees[oI].find(listOfDirs[dI]) != std::string::npos){
	      outDirs_p[dI]->cd();	      
	      outFile_p->cd(listOfDirs[dI].c_str());
	      break;
	    }
	  }
	  
	  outTrees_p[oI]->Write("", TObject::kOverwrite);
	  delete outTrees_p[oI];
	  outTrees_p[oI] = NULL;
	}

	for(unsigned int dI = 0; dI < listOfDirs.size(); ++dI){outDirs_p[listOfDirs.size()-1-dI] = NULL;}

	outFile_p->Close();
	delete outFile_p;
	outFile_p = NULL;
      }

      std::string tempFileName = outFileName + std::to_string(entry/nSplitEntries) + "OutOf" + std::to_string(splitLevel) + ".root";
      std::cout << "New File: " << tempFileName << std::endl;
      outFile_p = new TFile(tempFileName.c_str(), "RECREATE");

      for(unsigned int dI = 0; dI < outDirs_p.size(); ++dI){
	outDirs_p[dI] = (TDirectoryFile*)outFile_p->mkdir(listOfDirs[dI].c_str());
      }
      
      
      for(unsigned int oI = 0; oI < outTrees_p.size(); ++oI){
	outFile_p->cd();
	for(unsigned int dI = 0; dI < listOfDirs.size(); ++dI){
	  if(listOfTrees[oI].find(listOfDirs[dI]) != std::string::npos){
	    outDirs_p[dI]->cd();
	    outFile_p->cd(listOfDirs[dI].c_str());
	    break;
	  }
	}	
	
	outTrees_p[oI] = trees_p[oI]->CloneTree(0);	
      }
    }

    for(unsigned int tI = 0; tI < trees_p.size(); ++tI){
      trees_p[tI]->GetEntry(entry);
      outTrees_p[tI]->Fill();
    }    
  }

  outFile_p->cd();

  for(unsigned int oI = 0; oI < outTrees_p.size(); ++oI){
    outFile_p->cd();
    for(unsigned int dI = 0; dI < listOfDirs.size(); ++dI){
      if(listOfTrees[oI].find(listOfDirs[dI]) != std::string::npos){
	outDirs_p[dI]->cd();	      
	outFile_p->cd(listOfDirs[dI].c_str());
	break;
      }
    }
    
    outTrees_p[oI]->Write("", TObject::kOverwrite);
    delete outTrees_p[oI];
  }
  
  outFile_p->Close();
  delete outFile_p;
    
  inFile_p->Close();
  delete inFile_p;

  std::cout << "splitFiles Complete!" << std::endl;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/splitFiles.exe <inFileName> <nFinalFiles>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += splitFiles(argv[1], std::stoi(argv[2]));
  return retVal;
}
