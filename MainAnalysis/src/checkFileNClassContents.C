//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/stringUtil.h"

int checkFileNClassContents(const std::string inFileName, const std::string inFileName2 = "")
{
  if(!checkFile(inFileName) || inFileName.find(".root") == std::string::npos){
    std::cout << "Warning: inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> classes;
  std::vector<std::string> contents = returnRootFileContentsList(inFile_p, "", "", -1, &(classes));
  inFile_p->Close();
  delete inFile_p;

  Int_t th1ICounter = 0;
  Int_t th1FCounter = 0;
  Int_t th1DCounter = 0;
  Int_t th2ICounter = 0;
  Int_t th2FCounter = 0;
  Int_t th2DCounter = 0;
  Int_t rooUnfResCounter = 0;

  std::vector<std::string> th1IVect;
  std::vector<std::string> th1FVect;
  std::vector<std::string> th1DVect;
  std::vector<std::string> th2IVect;
  std::vector<std::string> th2FVect;
  std::vector<std::string> th2DVect;
  std::vector<std::string> rooUnfResVect;
  std::vector<std::string> generalHistContents;
  std::vector<std::string> generalHistClasses;

  for(unsigned int lI = 0; lI < classes.size(); ++lI){
    if(isStrSame("TH1I", classes.at(lI))){
      ++th1ICounter;
      th1IVect.push_back(contents.at(lI));
    }
    else if(isStrSame("TH1F", classes.at(lI))){
      ++th1FCounter;
      th1FVect.push_back(contents.at(lI));
    }
    else if(isStrSame("TH1D", classes.at(lI))){
      ++th1DCounter;
      th1DVect.push_back(contents.at(lI));
    }
    else if(isStrSame("TH2I", classes.at(lI))){
      ++th2ICounter;
      th2IVect.push_back(contents.at(lI));
    }
    else if(isStrSame("TH2F", classes.at(lI))){
      ++th2FCounter;
      th2FVect.push_back(contents.at(lI));
    }
    else if(isStrSame("TH2D", classes.at(lI))){
      ++th2DCounter;
      th2DVect.push_back(contents.at(lI));
    }
    else if(isStrSame("RooUnfoldResponse", classes.at(lI))){
      ++rooUnfResCounter;
      rooUnfResVect.push_back(contents.at(lI));
    }

    if(contents.at(lI).find("generalHistDir") != std::string::npos){
      generalHistContents.push_back(contents.at(lI));
      generalHistClasses.push_back(classes.at(lI));
    }
  }

  std::cout << "Contents of \'" << inFileName << "\': " << std::endl;
  std::cout << " TH1I: " << th1ICounter << std::endl;
  std::cout << " TH1F: " << th1FCounter << std::endl;
  std::cout << " TH1D: " << th1DCounter << std::endl;
  std::cout << " TH2I: " << th2ICounter << std::endl;
  std::cout << " TH2F: " << th2FCounter << std::endl;
  std::cout << " TH2D: " << th2DCounter << std::endl;
  std::cout << " rooUnfRes: " << rooUnfResCounter << std::endl;
  
  if(inFileName2.size() != 0 && checkFile(inFileName2) && inFileName2.find(".root") != std::string::npos){
    inFile_p = new TFile(inFileName2.c_str(), "READ");
    std::vector<std::string> classes2;
    std::vector<std::string> contents2 = returnRootFileContentsList(inFile_p, "", "", -1, &(classes2));
    inFile_p->Close();
    delete inFile_p;
    
    Int_t th1ICounter2 = 0;
    Int_t th1FCounter2 = 0;
    Int_t th1DCounter2 = 0;
    Int_t th2ICounter2 = 0;
    Int_t th2FCounter2 = 0;
    Int_t th2DCounter2 = 0;
    Int_t rooUnfResCounter2 = 0;
    
    std::vector<std::string> th1IVect2;
    std::vector<std::string> th1FVect2;
    std::vector<std::string> th1DVect2;
    std::vector<std::string> th2IVect2;
    std::vector<std::string> th2FVect2;
    std::vector<std::string> th2DVect2;
    std::vector<std::string> rooUnfResVect2;
    std::vector<std::string> generalHistContents2;
    std::vector<std::string> generalHistClasses2;

    for(unsigned int lI = 0; lI < classes2.size(); ++lI){
      if(isStrSame("TH1I", classes2.at(lI))){
	++th1ICounter2;
	th1IVect2.push_back(contents2.at(lI));
      }
      else if(isStrSame("TH1F", classes2.at(lI))){
	++th1FCounter2;
	th1FVect2.push_back(contents2.at(lI));
      }
      else if(isStrSame("TH1D", classes2.at(lI))){
	++th1DCounter2;
	th1DVect2.push_back(contents2.at(lI));
      }
      else if(isStrSame("TH2I", classes2.at(lI))){
	++th2ICounter2;
	th2IVect2.push_back(contents2.at(lI));
      }
      else if(isStrSame("TH2F", classes2.at(lI))){
	++th2FCounter2;
	th2FVect2.push_back(contents2.at(lI));
      }
      else if(isStrSame("TH2D", classes2.at(lI))){
	++th2DCounter2;
	th2DVect2.push_back(contents2.at(lI));
      }
      else if(isStrSame("RooUnfoldResponse", classes2.at(lI))){
	++rooUnfResCounter2;
	rooUnfResVect2.push_back(contents2.at(lI));
      }

      if(contents2.at(lI).find("generalHistDir") != std::string::npos){
	generalHistContents2.push_back(contents2.at(lI));
	generalHistClasses2.push_back(classes2.at(lI));
      }
    }

    std::cout << "Contents of \'" << inFileName2 << "\': " << std::endl;
    std::cout << " TH1I 2: " << th1ICounter2 << std::endl;
    std::cout << " TH1F 2: " << th1FCounter2 << std::endl;
    std::cout << " TH1D 2: " << th1DCounter2 << std::endl;
    std::cout << " TH2I 2: " << th2ICounter2 << std::endl;
    std::cout << " TH2F 2: " << th2FCounter2 << std::endl;
    std::cout << " TH2D 2: " << th2DCounter2 << std::endl;
    std::cout << " rooUnfRes 2: " << rooUnfResCounter2 << std::endl;

    std::cout << "GeneralHistDir in file1 \'" << inFileName << "\': " << std::endl;
    for(unsigned int gI = 0; gI < generalHistContents.size(); ++gI){
      std::cout << " " << gI << "/" << generalHistContents.size() << ": " << generalHistContents.at(gI) << ", \'" << generalHistClasses.at(gI) << "\'" << std::endl;
    }

    std::cout << "GeneralHistDir in file2 \'" << inFileName2 << "\': " << std::endl;
    for(unsigned int gI = 0; gI < generalHistContents2.size(); ++gI){
      std::cout << " " << gI << "/" << generalHistContents2.size() << ": " << generalHistContents2.at(gI) << ", \'" << generalHistClasses2.at(gI) << "\'" << std::endl;
    }

    std::vector<std::string> types = {"TH1I", "TH1F", "TH1D", "TH2I", "TH2F", "TH2D", "RooUnfRes"};
    std::vector<std::vector<std::string > > file1 = {th1IVect, th1FVect, th1DVect, th2IVect, th2FVect, th2DVect, rooUnfResVect};
    std::vector<std::vector<std::string > > file2 = {th1IVect2, th1FVect2, th1DVect2, th2IVect2, th2FVect2, th2DVect2, rooUnfResVect2};

    for(unsigned int tI = 0; tI < types.size(); ++tI){
      std::cout << "Checking type: " << types.at(tI) << std::endl;
      
      std::cout << " Checking all in file1 \'" << inFileName << "\'. " << std::endl;
      for(unsigned int fI1 = 0; fI1 < file1.at(tI).size(); ++fI1){
	bool isFound = false;

	for(unsigned int fI2 = 0; fI2 < file2.at(tI).size(); ++fI2){
	  if(isStrSame(file1.at(tI).at(fI1), file2.at(tI).at(fI2))){
	    isFound = true;
	    break;
	  }
	}
      
	if(!isFound) std::cout << "  Item " << file1.at(tI).at(fI1) << " not found." << std::endl;;
      }
  
      std::cout << " Checking all in file2 \'" << inFileName2 << "\'. " << std::endl;
      for(unsigned int fI2 = 0; fI2 < file2.at(tI).size(); ++fI2){
	bool isFound = false;

	for(unsigned int fI1 = 0; fI1 < file1.at(tI).size(); ++fI1){
	  if(isStrSame(file2.at(tI).at(fI2), file1.at(tI).at(fI1))){
	    isFound = true;
	    break;
	  }
	}
      
	if(!isFound) std::cout << "  Item " << file2.at(tI).at(fI2) << " not found." << std::endl;
      }
    }
  }

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc !=  2 && argc != 3){
    std::cout << "Usage: ./bin/checkFileNClassContents.exe <inFileName> <inFileName2-optional>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 2) retVal += checkFileNClassContents(argv[1]);
  else if(argc == 3) retVal += checkFileNClassContents(argv[1], argv[2]);
  return retVal;
}
