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

  bool CheckIntVect(std::vector<int> vect1, std::vector<int> vect2);
  bool CheckSmallRVals(std::vector<int> inSmallRVals);
  bool CheckLargeRVals(std::vector<int> inLargeRVals);

  bool CheckNRecoJtPtBinsSmallR(int inNRecoJtPtBinsSmallR);
  bool CheckNRecoJtPtBinsLargeR(int inNRecoJtPtBinsLargeR);
  bool CheckNGenJtPtBinsSmallR(int inNGenJtPtBinsSmallR);
  bool CheckNGenJtPtBinsLargeR(int inNGenJtPtBinsLargeR);

  bool CheckDoubleVect(std::vector<double> vect1, std::vector<double> vect2);
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


bool smallOrLargeR::CheckNSmallR(int inNSmallR)
{
  bool isGood = inNSmallR == nSmallR;
  if(!isGood) std::cout << "WARNING: nSmallR propagated \'" << inNSmallR << "\' doesn't match rReader \'" << nSmallR << "\'. return false" << std::endl;
  return isGood;
}

bool smallOrLargeR::CheckNLargeR(int inNLargeR)
{
  bool isGood = inNLargeR == nLargeR;
  if(!isGood) std::cout << "WARNING: nLargeR propagated \'" << inNLargeR << "\' doesn't match rReader \'" << nLargeR << "\'. return false" << std::endl;
  return isGood;
}

bool smallOrLargeR::CheckIntVect(std::vector<int> vect1, std::vector<int> vect2)
{
  bool isGood = vect1.size() == vect2.size();
  if(isGood){
    for(unsigned int i = 0; i < vect2.size(); ++i){
      if(vect2.at(i) != vect1.at(i)){
	isGood = false;
	break;
      }
    }
  }
  return isGood;
}
bool smallOrLargeR::CheckSmallRVals(std::vector<int> inSmallRVals)
{
  bool isGood = CheckIntVect(inSmallRVals, smallRVals);
  if(!isGood) std::cout << "WARNING: smallRVals propagated doesn't match rReader . return false" << std::endl;
  return isGood;
}
bool smallOrLargeR::CheckLargeRVals(std::vector<int> inLargeRVals)
{
  bool isGood = CheckIntVect(inLargeRVals, largeRVals);
  if(!isGood) std::cout << "WARNING: largeRVals propagated doesn't match rReader . return false" << std::endl;
  return isGood;
}

bool smallOrLargeR::CheckNRecoJtPtBinsSmallR(int inNRecoJtPtBinsSmallR)
{
  bool isGood = inNRecoJtPtBinsSmallR == nRecoJtPtBinsSmallR;
  if(!isGood) std::cout << "WARNING: nRecoJtPtBinsSmallR propagated \'" << inNRecoJtPtBinsSmallR << "\' doesn't match rReader \'" << nRecoJtPtBinsSmallR << "\'. return false" << std::endl;
  return isGood;
}
bool smallOrLargeR::CheckNRecoJtPtBinsLargeR(int inNRecoJtPtBinsLargeR)
{
  bool isGood = inNRecoJtPtBinsLargeR == nRecoJtPtBinsLargeR;
  if(!isGood) std::cout << "WARNING: nRecoJtPtBinsLargeR propagated \'" << inNRecoJtPtBinsLargeR << "\' doesn't match rReader \'" << nRecoJtPtBinsLargeR << "\'. return false" << std::endl;
  return isGood;
}
bool smallOrLargeR::CheckNGenJtPtBinsSmallR(int inNGenJtPtBinsSmallR)
{
  bool isGood = inNGenJtPtBinsSmallR == nGenJtPtBinsSmallR;
  if(!isGood) std::cout << "WARNING: nGenJtPtBinsSmallR propagated \'" << inNGenJtPtBinsSmallR << "\' doesn't match rReader \'" << nGenJtPtBinsSmallR << "\'. return false" << std::endl;
  return isGood;
}
bool smallOrLargeR::CheckNGenJtPtBinsLargeR(int inNGenJtPtBinsLargeR)
{
  bool isGood = inNGenJtPtBinsLargeR == nGenJtPtBinsLargeR;
  if(!isGood) std::cout << "WARNING: nGenJtPtBinsLargeR propagated \'" << inNGenJtPtBinsLargeR << "\' doesn't match rReader \'" << nGenJtPtBinsLargeR << "\'. return false" << std::endl;
  return isGood;
}

bool smallOrLargeR::CheckDoubleVect(std::vector<double> vect1, std::vector<double> vect2)
{
  bool isGood = vect1.size() == vect2.size();
  if(isGood){
    for(unsigned int i = 0; i < vect2.size(); ++i){
      if(TMath::Abs(vect2.at(i) - vect1.at(i)) >= deltaBin){
	isGood = false;
	break;
      }
    }
  }
  return isGood;
}
bool smallOrLargeR::CheckRecoJtPtBinsSmallR(std::vector<double> inRecoJtPtBinsSmallR)
{
  bool isGood = CheckDoubleVect(recoJtPtBinsSmallR, inRecoJtPtBinsSmallR);
  if(!isGood) std::cout << "WARNING: recoJtPtBinsSmallR propagated doesn't match rReader. return false" << std::endl;
  return isGood;
}
bool smallOrLargeR::CheckRecoJtPtBinsLargeR(std::vector<double> inRecoJtPtBinsLargeR)
{
  bool isGood = CheckDoubleVect(recoJtPtBinsLargeR, inRecoJtPtBinsLargeR);
  if(!isGood) std::cout << "WARNING: recoJtPtBinsLargeR propagated doesn't match rReader. return false" << std::endl;
  return isGood;
}
bool smallOrLargeR::CheckGenJtPtBinsSmallR(std::vector<double> inGenJtPtBinsSmallR)
{
  bool isGood = CheckDoubleVect(genJtPtBinsSmallR, inGenJtPtBinsSmallR);
  if(!isGood) std::cout << "WARNING: genJtPtBinsSmallR propagated doesn't match rReader. return false" << std::endl;
  return isGood;
}
bool smallOrLargeR::CheckGenJtPtBinsLargeR(std::vector<double> inGenJtPtBinsLargeR)
{
  bool isGood = CheckDoubleVect(genJtPtBinsLargeR, inGenJtPtBinsLargeR);
  if(!isGood) std::cout << "WARNING: genJtPtBinsLargeR propagated doesn't match rReader. return false" << std::endl;
  return isGood;
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
