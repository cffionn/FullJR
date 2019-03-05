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

  static const int nRecoJtPtBinsSmallRCent0to10 = 13;
  static const int nRecoJtPtBinsSmallRCent10to30 = 12;
  static const int nRecoJtPtBinsSmallRCent30to50 = 11;
  static const int nRecoJtPtBinsSmallRCent50to90 = 11;
  std::vector<double> recoJtPtBinsSmallR = {150., 175., 200., 225, 250., 275, 300., 350., 400., 450., 500., 600., 800., 1000.};

  static const int nGenJtPtSmallBinsSmallRCent0to10 = 16;
  static const int nGenJtPtSmallBinsSmallRCent10to30 = 15;
  static const int nGenJtPtSmallBinsSmallRCent30to50 = 14;
  static const int nGenJtPtSmallBinsSmallRCent50to90 = 14;
  std::vector<double> genJtPtSmallBinsSmallR = {100., 150., 175., 200., 225., 250., 275., 300., 350., 400., 450., 500., 600., 800., 1000., 1200., 1500.};

  static const int nGenJtPtLargeBinsSmallRCent0to10 = 6;
  static const int nGenJtPtLargeBinsSmallRCent10to30 = 5;
  static const int nGenJtPtLargeBinsSmallRCent30to50 = 4;
  static const int nGenJtPtLargeBinsSmallRCent50to90 = 4;
  std::vector<double> genJtPtLargeBinsSmallR = {100., 200., 300., 400, 600., 1000., 1500.};

  static const int nRecoJtPtBinsLargeRCent0to10 = 7;
  static const int nRecoJtPtBinsLargeRCent10to30 = 6;
  static const int nRecoJtPtBinsLargeRCent30to50 = 5;
  static const int nRecoJtPtBinsLargeRCent50to90 = 5;
  std::vector<double> recoJtPtBinsLargeR = {250., 300., 350., 400., 500., 600., 800, 1000.};

  static const int nGenJtPtSmallBinsLargeRCent0to10 = 10;
  static const int nGenJtPtSmallBinsLargeRCent10to30 = 9;
  static const int nGenJtPtSmallBinsLargeRCent30to50 = 8;
  static const int nGenJtPtSmallBinsLargeRCent50to90 = 8;
  std::vector<double> genJtPtSmallBinsLargeR = {200., 250., 300., 350., 400., 500., 600., 800., 1000., 1200., 1500.};

  static const int nGenJtPtLargeBinsLargeRCent0to10 = 4;
  static const int nGenJtPtLargeBinsLargeRCent10to30 = 3;
  static const int nGenJtPtLargeBinsLargeRCent30to50 = 2;
  static const int nGenJtPtLargeBinsLargeRCent50to90 = 2;
  std::vector<double> genJtPtLargeBinsLargeR = {200., 300., 600., 1000., 1500.};

  smallOrLargeR();

  int GetNSmallR(){return nSmallR;}
  int GetNLargeR(){return nLargeR;}

  std::vector<int> GetSmallRVals(){return smallRVals;}
  std::vector<int> GetLargeRVals(){return largeRVals;}

  int GetNRecoJtPtBinsSmallRCent0to10(){return nRecoJtPtBinsSmallRCent0to10;}
  int GetNRecoJtPtBinsLargeRCent0to10(){return nRecoJtPtBinsLargeRCent0to10;}
  int GetNGenJtPtSmallBinsSmallRCent0to10(){return nGenJtPtSmallBinsSmallRCent0to10;}
  int GetNGenJtPtSmallBinsLargeRCent0to10(){return nGenJtPtSmallBinsLargeRCent0to10;}
  int GetNGenJtPtLargeBinsSmallRCent0to10(){return nGenJtPtLargeBinsSmallRCent0to10;}
  int GetNGenJtPtLargeBinsLargeRCent0to10(){return nGenJtPtLargeBinsLargeRCent0to10;}

  int GetNRecoJtPtBinsSmallRCent10to30(){return nRecoJtPtBinsSmallRCent10to30;}
  int GetNRecoJtPtBinsLargeRCent10to30(){return nRecoJtPtBinsLargeRCent10to30;}
  int GetNGenJtPtSmallBinsSmallRCent10to30(){return nGenJtPtSmallBinsSmallRCent10to30;}
  int GetNGenJtPtSmallBinsLargeRCent10to30(){return nGenJtPtSmallBinsLargeRCent10to30;}
  int GetNGenJtPtLargeBinsSmallRCent10to30(){return nGenJtPtLargeBinsSmallRCent10to30;}
  int GetNGenJtPtLargeBinsLargeRCent10to30(){return nGenJtPtLargeBinsLargeRCent10to30;}

  int GetNRecoJtPtBinsSmallRCent30to50(){return nRecoJtPtBinsSmallRCent30to50;}
  int GetNRecoJtPtBinsLargeRCent30to50(){return nRecoJtPtBinsLargeRCent30to50;}
  int GetNGenJtPtSmallBinsSmallRCent30to50(){return nGenJtPtSmallBinsSmallRCent30to50;}
  int GetNGenJtPtSmallBinsLargeRCent30to50(){return nGenJtPtSmallBinsLargeRCent30to50;}
  int GetNGenJtPtLargeBinsSmallRCent30to50(){return nGenJtPtLargeBinsSmallRCent30to50;}
  int GetNGenJtPtLargeBinsLargeRCent30to50(){return nGenJtPtLargeBinsLargeRCent30to50;}

  int GetNRecoJtPtBinsSmallRCent50to90(){return nRecoJtPtBinsSmallRCent50to90;}
  int GetNRecoJtPtBinsLargeRCent50to90(){return nRecoJtPtBinsLargeRCent50to90;}
  int GetNGenJtPtSmallBinsSmallRCent50to90(){return nGenJtPtSmallBinsSmallRCent50to90;}
  int GetNGenJtPtSmallBinsLargeRCent50to90(){return nGenJtPtSmallBinsLargeRCent50to90;}
  int GetNGenJtPtLargeBinsSmallRCent50to90(){return nGenJtPtLargeBinsSmallRCent50to90;}
  int GetNGenJtPtLargeBinsLargeRCent50to90(){return nGenJtPtLargeBinsLargeRCent50to90;}

  //bins are are double for convenience; int would be fine and delta is loose to account
  std::vector<double> GetRecoJtPtBinsSmallR(){return recoJtPtBinsSmallR;}
  std::vector<double> GetRecoJtPtBinsLargeR(){return recoJtPtBinsLargeR;}
  std::vector<double> GetGenJtPtSmallBinsSmallR(){return genJtPtSmallBinsSmallR;}
  std::vector<double> GetGenJtPtSmallBinsLargeR(){return genJtPtSmallBinsLargeR;}
  std::vector<double> GetGenJtPtLargeBinsSmallR(){return genJtPtLargeBinsSmallR;}
  std::vector<double> GetGenJtPtLargeBinsLargeR(){return genJtPtLargeBinsLargeR;}

  bool CheckNSmallR(int inNSmallR);
  bool CheckNLargeR(int inNLargeR);

  bool CheckIntVect(std::vector<int> vect1, std::vector<int> vect2);
  bool CheckSmallRVals(std::vector<int> inSmallRVals);
  bool CheckLargeRVals(std::vector<int> inLargeRVals);

  bool CheckNBins(int inNBins, int compNBins, std::string nBinsStr);
  void SizeIsLTOrEQ(int inNBins, std::vector<double> inBins, std::string nBinsStr, std::string binsStr);
  
  bool CheckNRecoJtPtBinsSmallRCent0to10(int inNRecoJtPtBinsSmallRCent0to10);
  bool CheckNRecoJtPtBinsLargeRCent0to10(int inNRecoJtPtBinsLargeRCent0to10);
  bool CheckNGenJtPtSmallBinsSmallRCent0to10(int inNGenJtPtSmallBinsSmallRCent0to10);
  bool CheckNGenJtPtSmallBinsLargeRCent0to10(int inNGenJtPtSmallBinsLargeRCent0to10);
  bool CheckNGenJtPtLargeBinsSmallRCent0to10(int inNGenJtPtLargeBinsSmallRCent0to10);
  bool CheckNGenJtPtLargeBinsLargeRCent0to10(int inNGenJtPtLargeBinsLargeRCent0to10);

  bool CheckNRecoJtPtBinsSmallRCent10to30(int inNRecoJtPtBinsSmallRCent10to30);
  bool CheckNRecoJtPtBinsLargeRCent10to30(int inNRecoJtPtBinsLargeRCent10to30);
  bool CheckNGenJtPtSmallBinsSmallRCent10to30(int inNGenJtPtSmallBinsSmallRCent10to30);
  bool CheckNGenJtPtSmallBinsLargeRCent10to30(int inNGenJtPtSmallBinsLargeRCent10to30);
  bool CheckNGenJtPtLargeBinsSmallRCent10to30(int inNGenJtPtLargeBinsSmallRCent10to30);
  bool CheckNGenJtPtLargeBinsLargeRCent10to30(int inNGenJtPtLargeBinsLargeRCent10to30);

  bool CheckNRecoJtPtBinsSmallRCent30to50(int inNRecoJtPtBinsSmallRCent30to50);
  bool CheckNRecoJtPtBinsLargeRCent30to50(int inNRecoJtPtBinsLargeRCent30to50);
  bool CheckNGenJtPtSmallBinsSmallRCent30to50(int inNGenJtPtSmallBinsSmallRCent30to50);
  bool CheckNGenJtPtSmallBinsLargeRCent30to50(int inNGenJtPtSmallBinsLargeRCent30to50);
  bool CheckNGenJtPtLargeBinsSmallRCent30to50(int inNGenJtPtLargeBinsSmallRCent30to50);
  bool CheckNGenJtPtLargeBinsLargeRCent30to50(int inNGenJtPtLargeBinsLargeRCent30to50);

  bool CheckNRecoJtPtBinsSmallRCent50to90(int inNRecoJtPtBinsSmallRCent50to90);
  bool CheckNRecoJtPtBinsLargeRCent50to90(int inNRecoJtPtBinsLargeRCent50to90);
  bool CheckNGenJtPtSmallBinsSmallRCent50to90(int inNGenJtPtSmallBinsSmallRCent50to90);
  bool CheckNGenJtPtSmallBinsLargeRCent50to90(int inNGenJtPtSmallBinsLargeRCent50to90);
  bool CheckNGenJtPtLargeBinsSmallRCent50to90(int inNGenJtPtLargeBinsSmallRCent50to90);
  bool CheckNGenJtPtLargeBinsLargeRCent50to90(int inNGenJtPtLargeBinsLargeRCent50to90);

  bool CheckDoubleVect(std::vector<double> vect1, std::vector<double> vect2);
  bool CheckRecoJtPtBinsSmallR(std::vector<double> inRecoJtPtBinsSmallR);
  bool CheckRecoJtPtBinsLargeR(std::vector<double> inRecoJtPtBinsLargeR);
  bool CheckGenJtPtSmallBinsSmallR(std::vector<double> inGenJtPtSmallBinsSmallR);
  bool CheckGenJtPtSmallBinsLargeR(std::vector<double> inGenJtPtSmallBinsLargeR);
  bool CheckGenJtPtLargeBinsSmallR(std::vector<double> inGenJtPtLargeBinsSmallR);
  bool CheckGenJtPtLargeBinsLargeR(std::vector<double> inGenJtPtLargeBinsLargeR);

  bool GetIsSmallR(int rVal); 
  bool GetIsLargeR(int rVal);
  bool GetIsSmallOrLargeR(std::vector<int> vals, int rVal);
  bool GetIsSmallOrLargeR(int nVals, int vals[], int rVal);

  int GetSmallOrLargeRNBins(bool isSmallR, bool isGen, bool getSmallBins, std::string centStr);
  void GetSmallOrLargeRBins(bool isSmallR, bool isGen, int nBins, double bins[], bool getSmallBins);
};


smallOrLargeR::smallOrLargeR()
{
  SizeIsLTOrEQ(nRecoJtPtBinsSmallRCent0to10+1, recoJtPtBinsSmallR, "nRecoJtPtBinsSmallRCent0to10+1", "recoJtPtBinsSmallR.size()");
  SizeIsLTOrEQ(nRecoJtPtBinsLargeRCent0to10+1, recoJtPtBinsLargeR, "nRecoJtPtBinsLargeRCent0to10+1", "recoJtPtBinsLargeR.size()");
  SizeIsLTOrEQ(nGenJtPtSmallBinsSmallRCent0to10+1, genJtPtSmallBinsSmallR, "nGenJtPtSmallBinsSmallRCent0to10+1", "genJtPtSmallBinsSmallR.size()");
  SizeIsLTOrEQ(nGenJtPtSmallBinsLargeRCent0to10+1, genJtPtSmallBinsLargeR, "nGenJtPtSmallBinsLargeRCent0to10+1", "genJtPtSmallBinsLargeR.size()");
  SizeIsLTOrEQ(nGenJtPtLargeBinsSmallRCent0to10+1, genJtPtLargeBinsSmallR, "nGenJtPtLargeBinsSmallRCent0to10+1", "genJtPtLargeBinsSmallR.size()");
  SizeIsLTOrEQ(nGenJtPtLargeBinsLargeRCent0to10+1, genJtPtLargeBinsLargeR, "nGenJtPtLargeBinsLargeRCent0to10+1", "genJtPtLargeBinsLargeR.size()");

  SizeIsLTOrEQ(nRecoJtPtBinsSmallRCent10to30+1, recoJtPtBinsSmallR, "nRecoJtPtBinsSmallRCent10to30+1", "recoJtPtBinsSmallR.size()");
  SizeIsLTOrEQ(nRecoJtPtBinsLargeRCent10to30+1, recoJtPtBinsLargeR, "nRecoJtPtBinsLargeRCent10to30+1", "recoJtPtBinsLargeR.size()");
  SizeIsLTOrEQ(nGenJtPtSmallBinsSmallRCent10to30+1, genJtPtSmallBinsSmallR, "nGenJtPtSmallBinsSmallRCent10to30+1", "genJtPtSmallBinsSmallR.size()");
  SizeIsLTOrEQ(nGenJtPtSmallBinsLargeRCent10to30+1, genJtPtSmallBinsLargeR, "nGenJtPtSmallBinsLargeRCent10to30+1", "genJtPtSmallBinsLargeR.size()");
  SizeIsLTOrEQ(nGenJtPtLargeBinsSmallRCent10to30+1, genJtPtLargeBinsSmallR, "nGenJtPtLargeBinsSmallRCent10to30+1", "genJtPtLargeBinsSmallR.size()");
  SizeIsLTOrEQ(nGenJtPtLargeBinsLargeRCent10to30+1, genJtPtLargeBinsLargeR, "nGenJtPtLargeBinsLargeRCent10to30+1", "genJtPtLargeBinsLargeR.size()");
 
  SizeIsLTOrEQ(nRecoJtPtBinsSmallRCent30to50+1, recoJtPtBinsSmallR, "nRecoJtPtBinsSmallRCent30to50+1", "recoJtPtBinsSmallR.size()");
  SizeIsLTOrEQ(nRecoJtPtBinsLargeRCent30to50+1, recoJtPtBinsLargeR, "nRecoJtPtBinsLargeRCent30to50+1", "recoJtPtBinsLargeR.size()");
  SizeIsLTOrEQ(nGenJtPtSmallBinsSmallRCent30to50+1, genJtPtSmallBinsSmallR, "nGenJtPtSmallBinsSmallRCent30to50+1", "genJtPtSmallBinsSmallR.size()");
  SizeIsLTOrEQ(nGenJtPtSmallBinsLargeRCent30to50+1, genJtPtSmallBinsLargeR, "nGenJtPtSmallBinsLargeRCent30to50+1", "genJtPtSmallBinsLargeR.size()");
  SizeIsLTOrEQ(nGenJtPtLargeBinsSmallRCent30to50+1, genJtPtLargeBinsSmallR, "nGenJtPtLargeBinsSmallRCent30to50+1", "genJtPtLargeBinsSmallR.size()");
  SizeIsLTOrEQ(nGenJtPtLargeBinsLargeRCent30to50+1, genJtPtLargeBinsLargeR, "nGenJtPtLargeBinsLargeRCent30to50+1", "genJtPtLargeBinsLargeR.size()");

  SizeIsLTOrEQ(nRecoJtPtBinsSmallRCent50to90+1, recoJtPtBinsSmallR, "nRecoJtPtBinsSmallRCent50to90+1", "recoJtPtBinsSmallR.size()");
  SizeIsLTOrEQ(nRecoJtPtBinsLargeRCent50to90+1, recoJtPtBinsLargeR, "nRecoJtPtBinsLargeRCent50to90+1", "recoJtPtBinsLargeR.size()");
  SizeIsLTOrEQ(nGenJtPtSmallBinsSmallRCent50to90+1, genJtPtSmallBinsSmallR, "nGenJtPtSmallBinsSmallRCent50to90+1", "genJtPtSmallBinsSmallR.size()");
  SizeIsLTOrEQ(nGenJtPtSmallBinsLargeRCent50to90+1, genJtPtSmallBinsLargeR, "nGenJtPtSmallBinsLargeRCent50to90+1", "genJtPtSmallBinsLargeR.size()");
  SizeIsLTOrEQ(nGenJtPtLargeBinsSmallRCent50to90+1, genJtPtLargeBinsSmallR, "nGenJtPtLargeBinsSmallRCent50to90+1", "genJtPtLargeBinsSmallR.size()");
  SizeIsLTOrEQ(nGenJtPtLargeBinsLargeRCent50to90+1, genJtPtLargeBinsLargeR, "nGenJtPtLargeBinsLargeRCent50to90+1", "genJtPtLargeBinsLargeR.size()");
 
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

bool smallOrLargeR::CheckNBins(int inNBins, int compNBins, std::string nBinsStr)
{
  bool isGood = inNBins == compNBins;
  if(!isGood) std::cout << "WARNING: " << nBinsStr << " propagated \'" << inNBins << "\' doesn't match rReader \'" << compNBins << "\'. return false" << std::endl;
  return isGood;
}

void smallOrLargeR::SizeIsLTOrEQ(int inNBins, std::vector<double> inBins, std::string nBinsStr, std::string binsStr)
{
  if(((unsigned int)inNBins) > inBins.size()){
    std::cout << "Warning in smallOrLargeR: " << nBinsStr << ", " << inNBins << ", is incompatible w/ " << binsStr << ", " << inBins.size() << "." << std::endl;
  }
  return;
}

bool smallOrLargeR::CheckNRecoJtPtBinsSmallRCent0to10(int inNRecoJtPtBinsSmallRCent0to10){return CheckNBins(inNRecoJtPtBinsSmallRCent0to10, nRecoJtPtBinsSmallRCent0to10, "nRecoJtPtBinsSmallRCent0to10");}
bool smallOrLargeR::CheckNRecoJtPtBinsLargeRCent0to10(int inNRecoJtPtBinsLargeRCent0to10){return CheckNBins(inNRecoJtPtBinsLargeRCent0to10, nRecoJtPtBinsLargeRCent0to10, "nRecoJtPtBinsLargeRCent0to10");}
bool smallOrLargeR::CheckNGenJtPtSmallBinsSmallRCent0to10(int inNGenJtPtSmallBinsSmallRCent0to10){return CheckNBins(inNGenJtPtSmallBinsSmallRCent0to10, nGenJtPtSmallBinsSmallRCent0to10, "nGenJtPtSmallBinsSmallRCent0to10");}
bool smallOrLargeR::CheckNGenJtPtSmallBinsLargeRCent0to10(int inNGenJtPtSmallBinsLargeRCent0to10){return CheckNBins(inNGenJtPtSmallBinsLargeRCent0to10, nGenJtPtSmallBinsLargeRCent0to10, "nGenJtPtSmallBinsLargeRCent0to10");}
bool smallOrLargeR::CheckNGenJtPtLargeBinsSmallRCent0to10(int inNGenJtPtLargeBinsSmallRCent0to10){return CheckNBins(inNGenJtPtLargeBinsSmallRCent0to10, nGenJtPtLargeBinsSmallRCent0to10, "nGenJtPtLargeBinsSmallRCent0to10");}
bool smallOrLargeR::CheckNGenJtPtLargeBinsLargeRCent0to10(int inNGenJtPtLargeBinsLargeRCent0to10){return CheckNBins(inNGenJtPtLargeBinsLargeRCent0to10, nGenJtPtLargeBinsLargeRCent0to10, "nGenJtPtLargeBinsLargeRCent0to10");}
  

bool smallOrLargeR::CheckNRecoJtPtBinsSmallRCent10to30(int inNRecoJtPtBinsSmallRCent10to30){return CheckNBins(inNRecoJtPtBinsSmallRCent10to30, nRecoJtPtBinsSmallRCent10to30, "nRecoJtPtBinsSmallRCent10to30");}
bool smallOrLargeR::CheckNRecoJtPtBinsLargeRCent10to30(int inNRecoJtPtBinsLargeRCent10to30){return CheckNBins(inNRecoJtPtBinsLargeRCent10to30, nRecoJtPtBinsLargeRCent10to30, "nRecoJtPtBinsLargeRCent10to30");}
bool smallOrLargeR::CheckNGenJtPtSmallBinsSmallRCent10to30(int inNGenJtPtSmallBinsSmallRCent10to30){return CheckNBins(inNGenJtPtSmallBinsSmallRCent10to30, nGenJtPtSmallBinsSmallRCent10to30, "nGenJtPtSmallBinsSmallRCent10to30");}
bool smallOrLargeR::CheckNGenJtPtSmallBinsLargeRCent10to30(int inNGenJtPtSmallBinsLargeRCent10to30){return CheckNBins(inNGenJtPtSmallBinsLargeRCent10to30, nGenJtPtSmallBinsLargeRCent10to30, "nGenJtPtSmallBinsLargeRCent10to30");}
bool smallOrLargeR::CheckNGenJtPtLargeBinsSmallRCent10to30(int inNGenJtPtLargeBinsSmallRCent10to30){return CheckNBins(inNGenJtPtLargeBinsSmallRCent10to30, nGenJtPtLargeBinsSmallRCent10to30, "nGenJtPtLargeBinsSmallRCent10to30");}
bool smallOrLargeR::CheckNGenJtPtLargeBinsLargeRCent10to30(int inNGenJtPtLargeBinsLargeRCent10to30){return CheckNBins(inNGenJtPtLargeBinsLargeRCent10to30, nGenJtPtLargeBinsLargeRCent10to30, "nGenJtPtLargeBinsLargeRCent10to30");}
  

bool smallOrLargeR::CheckNRecoJtPtBinsSmallRCent30to50(int inNRecoJtPtBinsSmallRCent30to50){return CheckNBins(inNRecoJtPtBinsSmallRCent30to50, nRecoJtPtBinsSmallRCent30to50, "nRecoJtPtBinsSmallRCent30to50");}
bool smallOrLargeR::CheckNRecoJtPtBinsLargeRCent30to50(int inNRecoJtPtBinsLargeRCent30to50){return CheckNBins(inNRecoJtPtBinsLargeRCent30to50, nRecoJtPtBinsLargeRCent30to50, "nRecoJtPtBinsLargeRCent30to50");}
bool smallOrLargeR::CheckNGenJtPtSmallBinsSmallRCent30to50(int inNGenJtPtSmallBinsSmallRCent30to50){return CheckNBins(inNGenJtPtSmallBinsSmallRCent30to50, nGenJtPtSmallBinsSmallRCent30to50, "nGenJtPtSmallBinsSmallRCent30to50");}
bool smallOrLargeR::CheckNGenJtPtSmallBinsLargeRCent30to50(int inNGenJtPtSmallBinsLargeRCent30to50){return CheckNBins(inNGenJtPtSmallBinsLargeRCent30to50, nGenJtPtSmallBinsLargeRCent30to50, "nGenJtPtSmallBinsLargeRCent30to50");}
bool smallOrLargeR::CheckNGenJtPtLargeBinsSmallRCent30to50(int inNGenJtPtLargeBinsSmallRCent30to50){return CheckNBins(inNGenJtPtLargeBinsSmallRCent30to50, nGenJtPtLargeBinsSmallRCent30to50, "nGenJtPtLargeBinsSmallRCent30to50");}
bool smallOrLargeR::CheckNGenJtPtLargeBinsLargeRCent30to50(int inNGenJtPtLargeBinsLargeRCent30to50){return CheckNBins(inNGenJtPtLargeBinsLargeRCent30to50, nGenJtPtLargeBinsLargeRCent30to50, "nGenJtPtLargeBinsLargeRCent30to50");}
  

bool smallOrLargeR::CheckNRecoJtPtBinsSmallRCent50to90(int inNRecoJtPtBinsSmallRCent50to90){return CheckNBins(inNRecoJtPtBinsSmallRCent50to90, nRecoJtPtBinsSmallRCent50to90, "nRecoJtPtBinsSmallRCent50to90");}
bool smallOrLargeR::CheckNRecoJtPtBinsLargeRCent50to90(int inNRecoJtPtBinsLargeRCent50to90){return CheckNBins(inNRecoJtPtBinsLargeRCent50to90, nRecoJtPtBinsLargeRCent50to90, "nRecoJtPtBinsLargeRCent50to90");}
bool smallOrLargeR::CheckNGenJtPtSmallBinsSmallRCent50to90(int inNGenJtPtSmallBinsSmallRCent50to90){return CheckNBins(inNGenJtPtSmallBinsSmallRCent50to90, nGenJtPtSmallBinsSmallRCent50to90, "nGenJtPtSmallBinsSmallRCent50to90");}
bool smallOrLargeR::CheckNGenJtPtSmallBinsLargeRCent50to90(int inNGenJtPtSmallBinsLargeRCent50to90){return CheckNBins(inNGenJtPtSmallBinsLargeRCent50to90, nGenJtPtSmallBinsLargeRCent50to90, "nGenJtPtSmallBinsLargeRCent50to90");}
bool smallOrLargeR::CheckNGenJtPtLargeBinsSmallRCent50to90(int inNGenJtPtLargeBinsSmallRCent50to90){return CheckNBins(inNGenJtPtLargeBinsSmallRCent50to90, nGenJtPtLargeBinsSmallRCent50to90, "nGenJtPtLargeBinsSmallRCent50to90");}
bool smallOrLargeR::CheckNGenJtPtLargeBinsLargeRCent50to90(int inNGenJtPtLargeBinsLargeRCent50to90){return CheckNBins(inNGenJtPtLargeBinsLargeRCent50to90, nGenJtPtLargeBinsLargeRCent50to90, "nGenJtPtLargeBinsLargeRCent50to90");}  


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
bool smallOrLargeR::CheckGenJtPtSmallBinsSmallR(std::vector<double> inGenJtPtSmallBinsSmallR)
{
  bool isGood = CheckDoubleVect(genJtPtSmallBinsSmallR, inGenJtPtSmallBinsSmallR);
  if(!isGood) std::cout << "WARNING: genJtPtSmallBinsSmallR propagated doesn't match rReader. return false" << std::endl;
  return isGood;
}
bool smallOrLargeR::CheckGenJtPtSmallBinsLargeR(std::vector<double> inGenJtPtSmallBinsLargeR)
{
  bool isGood = CheckDoubleVect(genJtPtSmallBinsLargeR, inGenJtPtSmallBinsLargeR);
  if(!isGood) std::cout << "WARNING: genJtPtSmallBinsLargeR propagated doesn't match rReader. return false" << std::endl;
  return isGood;
}

bool smallOrLargeR::CheckGenJtPtLargeBinsSmallR(std::vector<double> inGenJtPtLargeBinsSmallR)
{
  bool isGood = CheckDoubleVect(genJtPtLargeBinsSmallR, inGenJtPtLargeBinsSmallR);
  if(!isGood) std::cout << "WARNING: genJtPtLargeBinsSmallR propagated doesn't match rReader. return false" << std::endl;
  return isGood;
}
bool smallOrLargeR::CheckGenJtPtLargeBinsLargeR(std::vector<double> inGenJtPtLargeBinsLargeR)
{
  bool isGood = CheckDoubleVect(genJtPtLargeBinsLargeR, inGenJtPtLargeBinsLargeR);
  if(!isGood) std::cout << "WARNING: genJtPtLargeBinsLargeR propagated doesn't match rReader. return false" << std::endl;
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

int smallOrLargeR::GetSmallOrLargeRNBins(bool isSmallR, bool isGen, bool getSmallBins, std::string centStr)
{
  if(centStr.find("Cent0to10") != std::string::npos){
    if(isSmallR && isGen && getSmallBins) return nGenJtPtSmallBinsSmallRCent0to10;
    else if(isSmallR && isGen && !getSmallBins) return nGenJtPtLargeBinsSmallRCent0to10;
    else if(isSmallR && !isGen) return nRecoJtPtBinsSmallRCent0to10;
    else if(!isSmallR && isGen && getSmallBins) return nGenJtPtSmallBinsLargeRCent0to10;
    else if(!isSmallR && isGen && !getSmallBins) return nGenJtPtLargeBinsLargeRCent0to10;
    else if(!isSmallR && !isGen) return nRecoJtPtBinsLargeRCent0to10;
  }
  else if(centStr.find("Cent10to30") != std::string::npos){
    if(isSmallR && isGen && getSmallBins) return nGenJtPtSmallBinsSmallRCent10to30;
    else if(isSmallR && isGen && !getSmallBins) return nGenJtPtLargeBinsSmallRCent10to30;
    else if(isSmallR && !isGen) return nRecoJtPtBinsSmallRCent10to30;
    else if(!isSmallR && isGen && getSmallBins) return nGenJtPtSmallBinsLargeRCent10to30;
    else if(!isSmallR && isGen && !getSmallBins) return nGenJtPtLargeBinsLargeRCent10to30;
    else if(!isSmallR && !isGen) return nRecoJtPtBinsLargeRCent10to30;
  }
  else if(centStr.find("Cent30to50") != std::string::npos){
    if(isSmallR && isGen && getSmallBins) return nGenJtPtSmallBinsSmallRCent30to50;
    else if(isSmallR && isGen && !getSmallBins) return nGenJtPtLargeBinsSmallRCent30to50;
    else if(isSmallR && !isGen) return nRecoJtPtBinsSmallRCent30to50;
    else if(!isSmallR && isGen && getSmallBins) return nGenJtPtSmallBinsLargeRCent30to50;
    else if(!isSmallR && isGen && !getSmallBins) return nGenJtPtLargeBinsLargeRCent30to50;
    else if(!isSmallR && !isGen) return nRecoJtPtBinsLargeRCent30to50;
  }
  else if(centStr.find("Cent50to90") != std::string::npos){
    if(isSmallR && isGen && getSmallBins) return nGenJtPtSmallBinsSmallRCent50to90;
    else if(isSmallR && isGen && !getSmallBins) return nGenJtPtLargeBinsSmallRCent50to90;
    else if(isSmallR && !isGen) return nRecoJtPtBinsSmallRCent50to90;
    else if(!isSmallR && isGen && getSmallBins) return nGenJtPtSmallBinsLargeRCent50to90;
    else if(!isSmallR && isGen && !getSmallBins) return nGenJtPtLargeBinsLargeRCent50to90;
    else if(!isSmallR && !isGen) return nRecoJtPtBinsLargeRCent50to90;
  }

  return -1;
}

void smallOrLargeR::GetSmallOrLargeRBins(bool isSmallR, bool isGen, int nBins, double bins[], bool getSmallBins)
{ 
  for(int i = 0; i < nBins; ++i){
    if(isSmallR && isGen && getSmallBins) bins[i] = genJtPtSmallBinsSmallR[i];
    else if(isSmallR && isGen && !getSmallBins) bins[i] = genJtPtLargeBinsSmallR[i];
    else if(isSmallR && !isGen) bins[i] = recoJtPtBinsSmallR[i];
    else if(!isSmallR && isGen && getSmallBins) bins[i] = genJtPtSmallBinsLargeR[i];
    else if(!isSmallR && isGen && !getSmallBins) bins[i] = genJtPtLargeBinsLargeR[i];
    else if(!isSmallR && !isGen) bins[i] = recoJtPtBinsLargeR[i];   
  }
  return;
}

#endif
