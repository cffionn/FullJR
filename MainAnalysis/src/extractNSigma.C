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

int extractNSigma(const std::string inFileName, double startVal, double epsilon, std::string paramsStr = "")
{
  if(!checkFile(inFileName) || inFileName.find(".txt") == std::string::npos){
    std::cout << "inFileName \'" << inFileName << "\' is invalid. return 1." << std::endl;
    return 1;
  }

  std::vector<double> params;
  if(paramsStr.size() != 0){
    if(paramsStr.substr(paramsStr.size()-1,1).find(",") == std::string::npos) paramsStr = paramsStr + ",";
    while(paramsStr.find(",") != std::string::npos){
      params.push_back(std::stod(paramsStr.substr(0, paramsStr.find(","))));
      paramsStr.replace(0, paramsStr.find(",")+1, "");
    }
    std::cout << "USING PARAMS: " << std::endl;
    std::cout << " C=" << params[0] << std::endl;
    std::cout << " S=" << params[1] << std::endl;
    std::cout << " N=" << params[2] << std::endl;
  }
  
  std::vector<double> xVals;
  std::vector<double> yVals;
  std::vector<double> xBins = {100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 350, 400, 500, 700};

  std::string binsStr = "";
  
  std::ifstream inFile(inFileName.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    std::string firstStr = tempStr.substr(0, tempStr.find(","));
    tempStr.replace(0, tempStr.find(",")+1, "");
    xVals.push_back(std::stod(firstStr));
    yVals.push_back(std::stod(tempStr));
  }   
  inFile.close();


  Int_t min = 0;
  Int_t max = startVal;
  const Int_t nMaxBins = 500;
  Int_t nBins = (max - min)/10;
  Double_t bins[nMaxBins+1];
  getLinBins(min, max, nBins, bins);

  
  if(params.size() == 0){
    Int_t binPos = -1;
    for(Int_t bIX = 0; bIX < (Int_t)xVals.size(); ++bIX){
      if(startVal >= xBins[bIX] && startVal < xBins[bIX+1]){
	binPos = bIX;
      }
    }

    if(binPos > 1){
      Double_t res = (yVals[binPos] + yVals[binPos-1])/2.;
      std::cout << startVal - startVal*res*epsilon << std::endl;    
    }
    else{
      std::cout << "No bin pos found for " << startVal << " in " << xBins[0] << "-" << xBins[xBins.size()-1] << std::endl;
      return 1;
    }
  }
  else{
    double currentVal = startVal;
    binsStr = std::to_string((int)currentVal) + ",1000";
    
    while(currentVal > 100){
      double res = TMath::Sqrt(params[0]*params[0] + params[1]*params[1]/currentVal + params[2]*params[2]/(currentVal*currentVal));
      currentVal = currentVal - res*currentVal*epsilon;

      for(Int_t bI = 0; bI < nBins; ++bI){
	if(currentVal >= bins[bI] && currentVal < bins[bI+1]){
	  currentVal = bins[bI];
	  break;
	}
      }

      if(currentVal < 100) break;
      binsStr = std::to_string((int)currentVal) + "," + binsStr;
    }
    
    std::cout << binsStr << std::endl;
  }
  
  
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 4 && argc != 5){
    std::cout << "Usage: ./bin/extractNSigma.exe <inFileName> <startVal> <epsilon> <params-opt>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  if(argc == 4) retVal += extractNSigma(argv[1], std::stod(argv[2]), std::stod(argv[3]));
  else if(argc == 5) retVal += extractNSigma(argv[1], std::stod(argv[2]), std::stod(argv[3]), argv[4]);
  return retVal;
}
