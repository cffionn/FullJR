//cpp dependencies
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"

int extractSigmaFit(const std::string inFileName)
{
  if(!checkFile(inFileName) || inFileName.find(".txt") == std::string::npos){
    std::cout << "inFileName \'" << inFileName << "\' is invalid. return 1." << std::endl;
    return 1;
  }
  
  std::vector<double> xVals;
  std::vector<double> yVals;
  std::vector<double> xBins = {100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 350, 400, 500, 700};

  const int nMaxBins = 100;
  int nXBins = xBins.size()-1; 
  if(nMaxBins < (int)xBins.size()){
    std::cout << "XBins: " << nXBins << " is less than nMaxBins " << nMaxBins << ". return 1" << std::endl;
    return 1;
  }

  Double_t xBinsArr[nMaxBins];
  for(unsigned int xI = 0; xI < xBins.size(); ++xI){
    xBinsArr[xI] = xBins[xI];
  }
  
  std::ifstream inFile(inFileName.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    std::string firstStr = tempStr.substr(0, tempStr.find(","));
    tempStr.replace(0, tempStr.find(",")+1, "");
    xVals.push_back(std::stod(firstStr));
    yVals.push_back(std::stod(tempStr));
  }   
  inFile.close();

  TH1D* fitHist_p = new TH1D("fitHist_p", "", nXBins, xBinsArr);
  TF1* fit_p = new TF1("fit_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 100, 700);
  fit_p->SetParameter(0, 0.06);
  fit_p->SetParameter(1, 0.9);
  fit_p->SetParameter(2, 20);

  fit_p->SetParLimits(0, 0.04, 0.08);
  fit_p->SetParLimits(1, 0.7, 1.1);
  fit_p->SetParLimits(2, 0, 100);

  for(Int_t bIX = 0; bIX < nXBins; ++bIX){
    fitHist_p->SetBinContent(bIX+1, yVals[bIX]);
    fitHist_p->SetBinError(bIX+1, yVals[bIX]*.1);
  }
  
  fitHist_p->Print("ALL");  
  fitHist_p->Fit("fit_p", "Q E M N", "", 100, 700);

  std::cout << fit_p->GetParameter(0) << "," << fit_p->GetParameter(1) << "," << fit_p->GetParameter(2) << std::endl;
  
  delete fitHist_p;
    
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/extractSigmaFit.exe <inFileName>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += extractSigmaFit(argv[1]);
  return retVal;
}
