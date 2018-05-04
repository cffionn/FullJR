#ifndef SMEARINGFUNCS_H
#define SMEARINGFUNCS_H

#include <string>
#include "TF1.h"
#include "TMath.h"

#include "include/doGlobalDebug.h"

int getCentPos(int cent, std::vector<int> centBinsLow, std::vector<int> centBinsHi);

void getResForAKCs3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs3PU3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs3PU3PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs4PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs4PU4PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs4PU4PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs4PU3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs4PU3PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs6PU3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs6PU3PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs8PU3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs8PU3PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs10PU3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKCs10PU3PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKPu3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAKPu4PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAK1PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAK2PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAK3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAK4PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAK5PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAK6PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAK8PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
void getResForAK10PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn);
	      
double getResForPtEtaCentAlgo(double pt, double eta, int cent, const std::string algo, const bool doErr = false, const std::string relAlgo = "", const bool doErrRel = false)
{
  double smearingErrFactor = 1.15;

  //  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<int> centBinsLow;
  std::vector<int> centBinsHi;
  double cParam = 0;
  double sParam = 0;
  std::vector<double> nParam;

  double resVal = 0;

  //  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(algo.find("akCs3PU3PFFlow") != std::string::npos) getResForAKCs3PU3PFFlow(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs3PU3PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs3PU3PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs3PF") != std::string::npos) getResForAKCs3PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs4PU4PFFlow") != std::string::npos) getResForAKCs4PU4PFFlow(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs4PU4PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs4PU4PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs4PU3PFFlow") != std::string::npos) getResForAKCs4PU3PFFlow(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs4PU3PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs4PU3PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs6PU3PFFlow") != std::string::npos) getResForAKCs6PU3PFFlow(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs6PU3PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs6PU3PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs8PU3PFFlow") != std::string::npos) getResForAKCs8PU3PFFlow(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs8PU3PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs8PU3PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs10PU3PFFlow") != std::string::npos) getResForAKCs10PU3PFFlow(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs10PU3PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs10PU3PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akCs4PF") != std::string::npos) getResForAKCs4PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akPu3PF") != std::string::npos) getResForAKPu3PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("akPu4PF") != std::string::npos) getResForAKPu4PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("ak1PF") != std::string::npos) getResForAK1PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("ak2PF") != std::string::npos) getResForAK2PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("ak3PF") != std::string::npos) getResForAK3PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("ak4PF") != std::string::npos) getResForAK4PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("ak5PF") != std::string::npos) getResForAK5PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("ak6PF") != std::string::npos) getResForAK6PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("ak8PF") != std::string::npos) getResForAK8PF(centBinsLow, centBinsHi, cParam, sParam, nParam);
  else if(algo.find("ak10PF") != std::string::npos) getResForAK10PF(centBinsLow, centBinsHi, cParam, sParam, nParam);

  //  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  int centPos = getCentPos(cent, centBinsLow, centBinsHi);
  if(TMath::Abs(eta) < 2.0){
    
    //    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(relAlgo.size() == 0){
      //      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      //      if(doGlobalDebug) std::cout << "nparam size, centpos, algo: " << nParam.size() << ", " << centPos << ", " << algo << std::endl;

      if(!doErr) resVal = TMath::Sqrt(cParam*cParam + sParam*sParam/pt + nParam.at(centPos)*nParam.at(centPos)/(pt*pt));
      else resVal = TMath::Sqrt((smearingErrFactor*smearingErrFactor*cParam*cParam - cParam*cParam + (smearingErrFactor*smearingErrFactor*sParam*sParam - sParam*sParam)/pt + (smearingErrFactor*smearingErrFactor*nParam.at(centPos)*nParam.at(centPos) - nParam.at(centPos)*nParam.at(centPos))/(pt*pt)));
    }
    else{
      std::vector<int> centBinsLowRel;
      std::vector<int> centBinsHiRel;;
      double cParamRel = 0;
      double sParamRel = 0;
      std::vector<double> nParamRel;

      //      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      if(algo.find("akCs3PU3PFFlow") != std::string::npos) getResForAKCs3PU3PFFlow(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs3PU3PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs3PU3PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs3PF") != std::string::npos) getResForAKCs3PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs4PU4PFFlow") != std::string::npos) getResForAKCs4PU4PFFlow(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs4PU4PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs4PU4PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs4PU3PFFlow") != std::string::npos) getResForAKCs4PU3PFFlow(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs4PU3PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs4PU3PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs6PU3PFFlow") != std::string::npos) getResForAKCs6PU3PFFlow(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs6PU3PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs6PU3PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs8PU3PFFlow") != std::string::npos) getResForAKCs8PU3PFFlow(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs8PU3PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs8PU3PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs10PU3PFFlow") != std::string::npos) getResForAKCs10PU3PFFlow(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs10PU3PF") != std::string::npos && algo.find("Flow") == std::string::npos) getResForAKCs10PU3PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akCs4PF") != std::string::npos) getResForAKCs4PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akPu3PF") != std::string::npos) getResForAKPu3PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("akPu4PF") != std::string::npos) getResForAKPu4PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("ak1PF") != std::string::npos) getResForAK1PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("ak2PF") != std::string::npos) getResForAK2PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("ak3PF") != std::string::npos) getResForAK3PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);
      else if(algo.find("ak4PF") != std::string::npos) getResForAK4PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);      
      else if(algo.find("ak5PF") != std::string::npos) getResForAK5PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);      
      else if(algo.find("ak6PF") != std::string::npos) getResForAK6PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);      
      else if(algo.find("ak8PF") != std::string::npos) getResForAK8PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);      
      else if(algo.find("ak10PF") != std::string::npos) getResForAK10PF(centBinsLowRel, centBinsHiRel, cParamRel, sParamRel, nParamRel);      

      //      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      if(doErrRel) resVal = TMath::Sqrt(cParam*cParam - cParamRel*cParamRel + (sParam*sParam - sParamRel*sParamRel)/pt + (nParam.at(centPos)*nParam.at(centPos) - nParamRel.at(centPos)*nParamRel.at(centPos))/(pt*pt));
      else resVal = TMath::Sqrt(cParam*cParam - cParamRel*cParamRel + (sParam*sParam - sParamRel*sParamRel)/pt + (nParam.at(centPos)*nParam.at(centPos) - nParamRel.at(centPos)*nParamRel.at(centPos))/(pt*pt));
    }
  }

  //  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  return resVal;
}


void getResForAKCs3PU3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = 1.1;
  Double_t nParam[nCentBins] = {14.142, 10.023, 5.474, 0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}

void getResForAKCs3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = 1.1;
  Double_t nParam[nCentBins] = {14.142, 10.023, 5.474, 0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}


void getResForAKCs3PU3PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {14.3448, 10.3240, 6.7062, 0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}


void getResForAKCs4PU4PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = 1.1;
  Double_t nParam[nCentBins] = {14.142, 10.023, 5.474, 0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}

void getResForAKCs4PU3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = 1.1;
  Double_t nParam[nCentBins] = {14.142, 10.023, 5.474, 0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}


void getResForAKCs4PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {14.348, 10.324, 6.7062, 4.5972};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}


void getResForAKCs4PU4PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {17.5242, 12.2411, 7.2944, 3.5210};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}

void getResForAKCs4PU3PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {17.5242, 12.2411, 7.2944, 3.5210};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}


void getResForAKCs6PU3PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {25.4117, 17.7627, 10.2510, 4.6043};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}


void getResForAKCs6PU3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {25.4117, 17.7627, 10.2510, 4.6043};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}

void getResForAKCs8PU3PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {33.1414, 23.0186, 13.6921, 6.55};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}

void getResForAKCs8PU3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {33.1414, 23.0186, 13.6921, 6.55};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}


void getResForAKCs10PU3PFFlow(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {39.5597, 27.7144, 16.7752, 9.3324};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}


void getResForAKCs10PU3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {39.5597, 27.7144, 16.7752, 9.3324};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;
}



void getResForAKPu3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = 1.24;
  Double_t nParam[nCentBins] = {8.2, 5., 2.4, 0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;  
}

void getResForAKPu4PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 100};

  Double_t cParam = 0.06;
  Double_t sParam = 1.24;
  Double_t nParam[nCentBins] = {8.2, 5., 2.4, 0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;  
}

void getResForAK1PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 1;
  Int_t centBinsLow[nCentBins] = {0};
  Int_t centBinsHi[nCentBins] = {100};

  Double_t cParam = 0.06;
  Double_t sParam = .95;
  Double_t nParam[nCentBins] = {0.0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;  
}


void getResForAK2PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 1;
  Int_t centBinsLow[nCentBins] = {0};
  Int_t centBinsHi[nCentBins] = {100};

  Double_t cParam = 0.06;
  Double_t sParam = .95;
  Double_t nParam[nCentBins] = {0.0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;  
}


void getResForAK3PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 1;
  Int_t centBinsLow[nCentBins] = {0};
  Int_t centBinsHi[nCentBins] = {100};

  Double_t cParam = 0.06;
  Double_t sParam = .95;
  Double_t nParam[nCentBins] = {0.0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;  
}

void getResForAK4PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 1;
  Int_t centBinsLow[nCentBins] = {0};
  Int_t centBinsHi[nCentBins] = {100};

  Double_t cParam = 0.06;
  Double_t sParam = .95;
  Double_t nParam[nCentBins] = {0.0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;  
}


void getResForAK5PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 1;
  Int_t centBinsLow[nCentBins] = {0};
  Int_t centBinsHi[nCentBins] = {100};

  Double_t cParam = 0.06;
  Double_t sParam = .95;
  Double_t nParam[nCentBins] = {0.0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;  
}


void getResForAK6PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 1;
  Int_t centBinsLow[nCentBins] = {0};
  Int_t centBinsHi[nCentBins] = {100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {0.0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;  
}


void getResForAK8PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 1;
  Int_t centBinsLow[nCentBins] = {0};
  Int_t centBinsHi[nCentBins] = {100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {0.0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;  
}


void getResForAK10PF(std::vector<int>& centBinsLowIn, std::vector<int>& centBinsHiIn, double& cParamIn, double& sParamIn, std::vector<double>& nParamIn)
{
  const Int_t nCentBins = 1;
  Int_t centBinsLow[nCentBins] = {0};
  Int_t centBinsHi[nCentBins] = {100};

  Double_t cParam = 0.06;
  Double_t sParam = .93;
  Double_t nParam[nCentBins] = {0.0};

  cParamIn = cParam;
  sParamIn = sParam;
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    centBinsLowIn.push_back(centBinsLow[centIter]);
    centBinsHiIn.push_back(centBinsHi[centIter]);
    nParamIn.push_back(nParam[centIter]);
  }

  return;  
}


int getCentPos(int cent, std::vector<int> centBinsLow, std::vector<int> centBinsHi)
{
  int centPos = -1;
  int nCentBins = centBinsLow.size();
  for(int centIter = 0; centIter < nCentBins; ++centIter){
    if(cent >= centBinsLow.at(centIter) && cent < centBinsHi.at(centIter)){
      centPos = centIter;
      break;
    }
  }

  return centPos;
}



#endif
