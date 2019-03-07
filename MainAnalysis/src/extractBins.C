//cpp dependencies
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TMath.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"

int extractBins(const std::string inFileName, double startVal, double epsilon)
{
  if(!checkFile(inFileName) || inFileName.find(".txt") == std::string::npos){
    std::cout << "inFileName \'" << inFileName << "\' is invalid. return 1." << std::endl;
    return 1;
  }
  
  std::vector<double> xVals;
  std::vector<double> yVals;
  std::vector<double> xBins = {100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 350, 400, 500, 700};

  std::ifstream inFile(inFileName.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    std::string firstStr = tempStr.substr(0, tempStr.find(","));
    tempStr.replace(0, tempStr.find(",")+1, "");
    xVals.push_back(std::stod(firstStr));
    yVals.push_back(std::stod(tempStr));
  }   
  inFile.close();


  std::vector<double> resBins = {startVal};
  while(resBins[resBins.size()-1] < 1000.){
    int binPos = -1;
    for(unsigned int bIX = 0; bIX < xVals.size(); ++bIX){
      if(resBins[resBins.size()-1] >= xBins[bIX] && resBins[resBins.size()-1] < xBins[bIX+1]){
	binPos = bIX;
	break;
      }
    }

    if(binPos < 0){
      if(resBins[resBins.size()-1] > xBins[xBins.size()-1]){
	binPos = xVals.size()-1;
      }
      else{
	std::cout << "Cannot find \'" << resBins[resBins.size()-1] << "\'. return 1" << std::endl;
	return 1;
      }
    }

    resBins.push_back(resBins[resBins.size()-1]*(1 + yVals[binPos]*epsilon));
  }

  while(resBins[resBins.size()-2] > 600){
    resBins.erase(resBins.begin() + (resBins.size()-2));
  }  
  resBins[resBins.size()-2] = 600;

  Int_t min = startVal;
  Int_t max = 600;
  Int_t nBins = (max - min)/5;
  Double_t bins[nBins+1];
  getLinBins(min, max, nBins, bins);
  
  for(unsigned int rI = 1; rI < resBins.size()-2; ++rI){

    Int_t pos = -1;
    Double_t maxDelta = 10000;
    for(int bI = 0; bI < nBins+1; ++bI){
      if(TMath::Abs(resBins[rI] - bins[bI]) < maxDelta){
	pos = bI;
	maxDelta = TMath::Abs(resBins[rI] - bins[bI]);
      }
    }

    resBins[rI] = bins[pos];
  }
  
  resBins[resBins.size()-1] = 1000;
  resBins.push_back(1500);
  
  for(unsigned int rI = 0; rI < resBins.size(); ++rI){
    std::cout << resBins[rI] << std::endl;
  }
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 4){
    std::cout << "Usage: ./bin/extractBins.exe <inFileName> <startVal> <epsilon>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += extractBins(argv[1], std::stod(argv[2]), std::stod(argv[3]));
  return retVal;
}
