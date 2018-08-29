#ifndef SMALLORLARGER_H
#define SMALLORLARGER_H

#include <iostream>
#include <vector>

class smallOrLargeR
{
 public:
  const double deltaBin = 0.1;
  //Chosen loose since bins will be effective integer but double because of TH1 initialization

  static const int nSmallR = 6;
  std::vector<int> smallRVals = {1, 2, 3, 4, 5, 6};
  
  static const int nLargeR = 4;
  std::vector<int> largeRVals = {7, 8, 9, 10};

  static const int nRecoJtPtBinsSmallR = 7;
  std::vector<double> recoJtPtBinsSmallR = {150., 200., 250., 300., 400., 600., 1000., 1500.};

  static const int nGenJtPtBinsSmallR = 9;
  std::vector<double> genJtPtBinsSmallR = {100., 150., 200., 250., 300., 400., 600., 1000., 1500., 2000.};

  static const int nRecoJtPtBinsLargeR = 5;
  std::vector<double> recoJtPtBinsLargeR = {250., 300., 400., 600., 1000., 1500.};

  static const int nGenJtPtBinsLargeR = 7;
  std::vector<double> genJtPtBinsLargeR = {200., 250., 300., 400., 600., 1000., 1500., 2000.};

  smallOrLargeR();

  int GetNSmallR(){return nSmallR;}
  int GetNLargeR(){return nLargeR;}

  std::vector<int> GetSmallRVals(){return smallRVals;}
  std::vector<int> GetLargeRVals(){return largeRVals;}

  int GetNRecoJtPtBinsSmallR(){return nRecoJtPtBinsSmallR;}
  int GetNRecoJtPtBinsLargeR(){return nRecoJtPtBinsLargeR;}
  int GetNGenJtPtBinsSmallR(){return nGenJtPtBinsSmallR;}
  int GetNGenJtPtBinsLargeR(){return nGenJtPtBinsLargeR;}

  //bins are are double for convenience; int would be fine and delta is loose to account
  std::vector<double> GetRecoJtPtBinsSmallR(){return recoJtPtBinsSmallR;}
  std::vector<double> GetRecoJtPtBinsLargeR(){return recoJtPtBinsLargeR;}
  std::vector<double> GetGenJtPtBinsSmallR(){return genJtPtBinsSmallR;}
  std::vector<double> GetGenJtPtBinsLargeR(){return genJtPtBinsLargeR;}

  bool CheckNSmallR(int inNSmallR);
  bool CheckNLargeR(int inNLargeR);

  bool CheckSmallRVals(std::vector<int> inSmallRVals);
  bool CheckLargeRVals(std::vector<int> inLargeRVals);

  bool CheckNRecoJtPtBinsSmallR(int inNRecoJtPtBinsSmallR);
  bool CheckNRecoJtPtBinsLargeR(int inNRecoJtPtBinsLargeR);
  bool CheckNGenJtPtBinsSmallR(int inNGenJtPtBinsSmallR);
  bool CheckNGenJtPtBinsLargeR(int inNGenJtPtBinsLargeR);

  bool CheckRecoJtPtBinsSmallR(std::vector<double> inRecoJtPtBinsSmallR);
  bool CheckRecoJtPtBinsLargeR(std::vector<double> inRecoJtPtBinsLargeR);
  bool CheckGenJtPtBinsSmallR(std::vector<double> inGenJtPtBinsSmallR);
  bool CheckGenJtPtBinsLargeR(std::vector<double> inGenJtPtBinsLargeR);

  bool GetIsSmallR(int rVal); 
  bool GetIsLargeR(int rVal);
  bool GetIsSmallOrLargeR(std::vector<int> vals, int rVal);
  bool GetIsSmallOrLargeR(int nVals, int vals[], int rVal);

  int GetSmallOrLargeRNBins(bool isSmallR, bool isGen);
  void GetSmallOrLargeRBins(bool isSmallR, bool isGen, int nBins, double bins[]);
};


smallOrLargeR::smallOrLargeR()
{
  if(nRecoJtPtBinsSmallR+1 != recoJtPtBinsSmallR.size()){
    std::cout << "Warning in smallOrLargeR: nRecoJtPtBinsSmallR+1, " << nRecoJtPtBinsSmallR+1 << ", is incompatible w/ recoJtPtBinsSmallR.size(), " << recoJtPtBinsSmallR.size() << "." << std::endl;
  }
  if(nRecoJtPtBinsLargeR+1 != recoJtPtBinsLargeR.size()){
    std::cout << "Warning in smallOrLargeR: nRecoJtPtBinsLargeR+1, " << nRecoJtPtBinsLargeR+1 << ", is incompatible w/ recoJtPtBinsLargeR.size(), " << recoJtPtBinsLargeR.size() << "." << std::endl;
  }
  if(nGenJtPtBinsSmallR+1 != genJtPtBinsSmallR.size()){
    std::cout << "Warning in smallOrLargeR: nGenJtPtBinsSmallR+1, " << nGenJtPtBinsSmallR+1 << ", is incompatible w/ genJtPtBinsSmallR.size(), " << genJtPtBinsSmallR.size() << "." << std::endl;
  }
  if(nGenJtPtBinsLargeR+1 != genJtPtBinsLargeR.size()){
    std::cout << "Warning in smallOrLargeR: nGenJtPtBinsLargeR+1, " << nGenJtPtBinsLargeR+1 << ", is incompatible w/ genJtPtBinsLargeR.size(), " << genJtPtBinsLargeR.size() << "." << std::endl;
  }
  
  return;
}


bool smallOrLargeR::CheckNSmallR(int inNSmallR){return inNSmallR == nSmallR;}
bool smallOrLargeR::CheckNLargeR(int inNLargeR){return inNLargeR == nLargeR;}

bool smallOrLargeR::CheckSmallRVals(std::vector<int> inSmallRVals)
{
  if(inSmallRVals.size() != smallRVals.size()) return false;
  for(unsigned int i = 0; i < smallRVals.size(); ++i){
    if(smallRVals.at(i) != inSmallRVals.at(i)) return false;
  }
  return true;
}
bool smallOrLargeR::CheckLargeRVals(std::vector<int> inLargeRVals)
{
  if(inLargeRVals.size() != largeRVals.size()) return false;
  for(unsigned int i = 0; i < largeRVals.size(); ++i){
    if(largeRVals.at(i) != inLargeRVals.at(i)) return false;
  }
  return true;
}

bool smallOrLargeR::CheckNRecoJtPtBinsSmallR(int inNRecoJtPtBinsSmallR){return inNRecoJtPtBinsSmallR == nRecoJtPtBinsSmallR;}
bool smallOrLargeR::CheckNRecoJtPtBinsLargeR(int inNRecoJtPtBinsLargeR){return inNRecoJtPtBinsLargeR == nRecoJtPtBinsLargeR;}
bool smallOrLargeR::CheckNGenJtPtBinsSmallR(int inNGenJtPtBinsSmallR){return inNGenJtPtBinsSmallR == nGenJtPtBinsSmallR;}
bool smallOrLargeR::CheckNGenJtPtBinsLargeR(int inNGenJtPtBinsLargeR){return inNGenJtPtBinsLargeR == nGenJtPtBinsLargeR;}

bool smallOrLargeR::CheckRecoJtPtBinsSmallR(std::vector<double> inRecoJtPtBinsSmallR)
{
  if(inRecoJtPtBinsSmallR.size() != recoJtPtBinsSmallR.size()) return false;
  for(unsigned int i = 0; i < recoJtPtBinsSmallR.size(); ++i){
    if(TMath::Abs(inRecoJtPtBinsSmallR.at(i) - recoJtPtBinsSmallR.at(i)) >= deltaBin) return false;
  }
  return true;
}
bool smallOrLargeR::CheckRecoJtPtBinsLargeR(std::vector<double> inRecoJtPtBinsLargeR)
{
  if(inRecoJtPtBinsLargeR.size() != recoJtPtBinsLargeR.size()) return false;
  for(unsigned int i = 0; i < recoJtPtBinsLargeR.size(); ++i){
    if(TMath::Abs(inRecoJtPtBinsLargeR.at(i) - recoJtPtBinsLargeR.at(i)) >= deltaBin) return false;
  }
  return true;
}
bool smallOrLargeR::CheckGenJtPtBinsSmallR(std::vector<double> inGenJtPtBinsSmallR)
{
  if(inGenJtPtBinsSmallR.size() != genJtPtBinsSmallR.size()) return false;
  for(unsigned int i = 0; i < genJtPtBinsSmallR.size(); ++i){
    if(TMath::Abs(inGenJtPtBinsSmallR.at(i) - genJtPtBinsSmallR.at(i)) >= deltaBin) return false;
  }
  return true;
}
bool smallOrLargeR::CheckGenJtPtBinsLargeR(std::vector<double> inGenJtPtBinsLargeR)
{
  if(inGenJtPtBinsLargeR.size() != genJtPtBinsLargeR.size()) return false;
  for(unsigned int i = 0; i < genJtPtBinsLargeR.size(); ++i){
    if(TMath::Abs(inGenJtPtBinsLargeR.at(i) - genJtPtBinsLargeR.at(i)) >= deltaBin) return false;
  }
  return true;
}


bool smallOrLargeR::GetIsSmallOrLargeR(std::vector<int> vals, int rVal)
{
  for(unsigned int bI = 0; bI < vals.size(); ++bI){
    if(vals.at(bI) == rVal) return true;
  }
  return false;
}
bool smallOrLargeR::GetIsSmallOrLargeR(int nVals, int vals[], int rVal)
{
  for(int bI = 0; bI < nVals; ++bI){
    if(vals[bI] == rVal) return true;
  }
  return false;
}
bool smallOrLargeR::GetIsSmallR(int rVal){return GetIsSmallOrLargeR(smallRVals, rVal);}
bool smallOrLargeR::GetIsLargeR(int rVal){return GetIsSmallOrLargeR(largeRVals, rVal);}

int smallOrLargeR::GetSmallOrLargeRNBins(bool isSmallR, bool isGen)
{
  if(isSmallR && isGen) return nGenJtPtBinsSmallR;
  else if(isSmallR && !isGen) return nRecoJtPtBinsSmallR;
  else if(!isSmallR && isGen) return nGenJtPtBinsLargeR;
  else if(!isSmallR && !isGen) return nRecoJtPtBinsLargeR;
  
  return -1;
}

void smallOrLargeR::GetSmallOrLargeRBins(bool isSmallR, bool isGen, int nBins, double bins[])
{ 
  for(int i = 0; i < nBins; ++i){
    if(isSmallR && isGen) bins[i] = genJtPtBinsSmallR[i];
    else if(isSmallR && !isGen) bins[i] = recoJtPtBinsSmallR[i];
    else if(!isSmallR && isGen) bins[i] = genJtPtBinsLargeR[i];
    else if(!isSmallR && !isGen) bins[i] = recoJtPtBinsLargeR[i];   
  }
  return;
}

#endif
