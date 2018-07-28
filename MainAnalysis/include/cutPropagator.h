//cpp dependencies
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>

//ROOT dependencies
#include "TFile.h"
#include "TDirectory.h"
#include "TNamed.h"

#include "Utility/include/returnRootFileContentsList.h"

class cutPropagator
{
 public:
  std::vector<std::string> inFileNames;
  std::vector<std::string> inFullFileNames;

  bool isPP;

  double jtAbsEtaMax;
  int nResponseMod;
  std::vector<double> responseMod;
  std::vector<double> scaleFactor;

  int nJtPtBins;
  std::vector<double> jtPtBins;
  int nJtAbsEtaBins;
  std::vector<double> jtAbsEtaBinsLow;
  std::vector<double> jtAbsEtaBinsHi;
  
  int nPthats;
  std::vector<double> pthats;
  std::vector<double> pthatWeights;

  int nCentBins;
  std::vector<int> centBinsLow;
  std::vector<int> centBinsHi;

  int nID;
  std::vector<std::string> idStr;
  std::vector<double> jtPfCHMFCutLow;
  std::vector<double> jtPfCHMFCutHi;
  std::vector<double> jtPfMUMFCutLow;
  std::vector<double> jtPfMUMFCutHi;
  std::vector<double> jtPfNHFCutLow;
  std::vector<double> jtPfNHFCutHi;
  std::vector<double> jtPfNEFCutLow;
  std::vector<double> jtPfNEFCutHi;
  std::vector<double> jtPfMUFCutLow;
  std::vector<double> jtPfMUFCutHi;
  std::vector<double> jtPfCHFCutLow;
  std::vector<double> jtPfCHFCutHi;
  std::vector<double> jtPfCEFCutLow;
  std::vector<double> jtPfCEFCutHi;
  std::vector<int> jtPfMinMult;
  std::vector<int> jtPfMinChgMult;

  int nSyst;
  std::vector<std::string> systStr;

  void Clean();
  bool GetAllVarFromFile(TFile* inFile_p);
  bool WriteAllVarToFile(TFile* inFile_p, TDirectory* inDir_p, TDirectory* inSubDir_p);
  bool CheckPropagatorsMatch(cutPropagator inCutProp, bool doBothMCOrBothData, bool doBothPPOrBothPbPb);

  std::vector<std::string> GetInFileNames(){return inFileNames;}
  std::vector<std::string> GetInFullFileNames(){return inFullFileNames;}

  bool GetIsPP(){return isPP;}

  double GetJtAbsEtaMax(){return jtAbsEtaMax;}

  int GetNResponseMod(){return nResponseMod;}
  std::vector<double> GetResponseMod(){return responseMod;}
  std::vector<double> GetScaleFactor(){return scaleFactor;}

  int GetNJtPtBins(){return nJtPtBins;}
  std::vector<double> GetJtPtBins(){return jtPtBins;}
  int GetNJtAbsEtaBins(){return nJtAbsEtaBins;}
  std::vector<double> GetJtAbsEtaBinsLow(){return jtAbsEtaBinsLow;}
  std::vector<double> GetJtAbsEtaBinsHi(){return jtAbsEtaBinsHi;}

  int GetNPthats(){return nPthats;}
  std::vector<double> GetPthats(){return pthats;}
  std::vector<double> GetPthatWeights(){return pthatWeights;}

  int GetNCentBins(){return nCentBins;}
  std::vector<int> GetCentBinsLow(){return centBinsLow;}
  std::vector<int> GetCentBinsHi(){return centBinsHi;}

  int GetNID(){return nID;}
  std::vector<std::string> GetIdStr(){return idStr;}
  std::vector<double> GetJtPfCHMFCutLow(){return jtPfCHMFCutLow;}
  std::vector<double> GetJtPfCHMFCutHi(){return jtPfCHMFCutHi;}
  std::vector<double> GetJtPfMUMFCutLow(){return jtPfMUMFCutLow;}
  std::vector<double> GetJtPfMUMFCutHi(){return jtPfMUMFCutHi;}
  std::vector<double> GetJtPfNHFCutLow(){return jtPfNHFCutLow;}
  std::vector<double> GetJtPfNHFCutHi(){return jtPfNHFCutHi;}
  std::vector<double> GetJtPfNEFCutLow(){return jtPfNEFCutLow;}
  std::vector<double> GetJtPfNEFCutHi(){return jtPfNEFCutHi;}
  std::vector<double> GetJtPfMUFCutLow(){return jtPfMUFCutLow;}
  std::vector<double> GetJtPfMUFCutHi(){return jtPfMUFCutHi;}
  std::vector<double> GetJtPfCHFCutLow(){return jtPfCHFCutLow;}
  std::vector<double> GetJtPfCHFCutHi(){return jtPfCHFCutHi;}
  std::vector<double> GetJtPfCEFCutLow(){return jtPfCEFCutLow;}
  std::vector<double> GetJtPfCEFCutHi(){return jtPfCEFCutHi;}
  std::vector<int> GetJtPfMinMult(){return jtPfMinMult;}
  std::vector<int> GetJtPfMinChgMult(){return jtPfMinChgMult;}

  int GetNSyst(){return nSyst;}
  std::vector<std::string> GetSystStr(){return systStr;}

  void SetInFileNames(std::vector<std::string> inInFileNames){inFileNames = inInFileNames; return;}
  void SetInFullFileNames(std::vector<std::string> inInFullFileNames){inFullFileNames = inInFullFileNames; return;}

  void SetIsPP(bool inIsPP){isPP = inIsPP; return;}

  void SetJtAbsEtaMax(double inJtAbsEtaMax){jtAbsEtaMax = inJtAbsEtaMax; return;}

  void SetNResponseMod(int inNResponseMod){nResponseMod = inNResponseMod; return;}
  void SetResponseMod(std::vector<double> inResponseMod){responseMod = inResponseMod; return;}
  void SetResponseMod(int inN, const Double_t inResponseMod[]);
  void SetScaleFactor(std::vector<double> inScaleFactor){scaleFactor = inScaleFactor; return;}
  void SetScaleFactor(int inN, const Double_t inScaleFactor[]);

  void SetNJtPtBins(int inNJtPtBins){nJtPtBins = inNJtPtBins; return;}
  void SetJtPtBins(std::vector<double> inJtPtBins){jtPtBins = inJtPtBins; return;}
  void SetJtPtBins(int inN, const Double_t inJtPtBins[]);
  void SetNJtAbsEtaBins(int inNJtAbsEtaBins){nJtAbsEtaBins = inNJtAbsEtaBins; return;}
  void SetJtAbsEtaBinsLow(std::vector<double> inJtAbsEtaBinsLow){jtAbsEtaBinsLow = inJtAbsEtaBinsLow; return;}
  void SetJtAbsEtaBinsLow(int inN, const Double_t inJtAbsEtaBinsLow[]);  
  void SetJtAbsEtaBinsHi(std::vector<double> inJtAbsEtaBinsHi){jtAbsEtaBinsHi = inJtAbsEtaBinsHi; return;}
  void SetJtAbsEtaBinsHi(int inN, const Double_t inJtAbsEtaBinsHi[]);  

  void SetNPthats(int inNPthats){nPthats = inNPthats; return;}
  void SetPthats(std::vector<double> inPthats){pthats = inPthats; return;}
  void SetPthats(int inN, const Double_t inPthats[]);
  void SetPthatWeights(std::vector<double> inPthatWeights){pthatWeights = inPthatWeights; return;}
  void SetPthatWeights(int inN, const Double_t inPthatWeights[]);

  void SetNCentBins(int inNCentBins){nCentBins = inNCentBins; return;}
  void SetCentBinsLow(std::vector<int> inCentBinsLow){centBinsLow = inCentBinsLow; return;}
  void SetCentBinsLow(int inN, const Int_t inCentBinsLow[]);
  void SetCentBinsHi(std::vector<int> inCentBinsHi){centBinsHi = inCentBinsHi; return;}
  void SetCentBinsHi(int inN, const Int_t inCentBinsHi[]);

  void SetNID(int inNID){nID = inNID; return;}
  void SetIdStr(std::vector<std::string> inIdStr){idStr = inIdStr; return;};
  void SetIdStr(int inN, const std::string inIdStr[]);
  void SetJtPfCHMFCutLow(std::vector<double> inJtPfCHMFCutLow){jtPfCHMFCutLow = inJtPfCHMFCutLow; return;}
  void SetJtPfCHMFCutLow(int inN, const Double_t inJtPfCHMFCutLow[]);
  void SetJtPfCHMFCutHi(std::vector<double> inJtPfCHMFCutHi){jtPfCHMFCutHi = inJtPfCHMFCutHi; return;}
  void SetJtPfCHMFCutHi(int inN, const Double_t inJtPfCHMFCutHi[]);
  void SetJtPfMUMFCutLow(std::vector<double> inJtPfMUMFCutLow){jtPfMUMFCutLow = inJtPfMUMFCutLow; return;}
  void SetJtPfMUMFCutLow(int inN, const Double_t inJtPfMUMFCutLow[]);
  void SetJtPfMUMFCutHi(std::vector<double> inJtPfMUMFCutHi){jtPfMUMFCutHi = inJtPfMUMFCutHi; return;}
  void SetJtPfMUMFCutHi(int inN, const Double_t inJtPfMUMFCutHi[]);
  void SetJtPfNHFCutLow(std::vector<double> inJtPfNHFCutLow){jtPfNHFCutLow = inJtPfNHFCutLow; return;}
  void SetJtPfNHFCutLow(int inN, const Double_t inJtPfNHFCutLow[]);
  void SetJtPfNHFCutHi(std::vector<double> inJtPfNHFCutHi){jtPfNHFCutHi = inJtPfNHFCutHi; return;}
  void SetJtPfNHFCutHi(int inN, const Double_t inJtPfNHFCutHi[]);
  void SetJtPfNEFCutLow(std::vector<double> inJtPfNEFCutLow){jtPfNEFCutLow = inJtPfNEFCutLow; return;}
  void SetJtPfNEFCutLow(int inN, const Double_t inJtPfNEFCutLow[]);
  void SetJtPfNEFCutHi(std::vector<double> inJtPfNEFCutHi){jtPfNEFCutHi = inJtPfNEFCutHi; return;}
  void SetJtPfNEFCutHi(int inN, const Double_t inJtPfNEFCutHi[]);
  void SetJtPfMUFCutLow(std::vector<double> inJtPfMUFCutLow){jtPfMUFCutLow = inJtPfMUFCutLow; return;}
  void SetJtPfMUFCutLow(int inN, const Double_t inJtPfMUFCutLow[]);
  void SetJtPfMUFCutHi(std::vector<double> inJtPfMUFCutHi){jtPfMUFCutHi = inJtPfMUFCutHi; return;}
  void SetJtPfMUFCutHi(int inN, const Double_t inJtPfMUFCutHi[]);
  void SetJtPfCHFCutLow(std::vector<double> inJtPfCHFCutLow){jtPfCHFCutLow = inJtPfCHFCutLow; return;}
  void SetJtPfCHFCutLow(int inN, const Double_t inJtPfCHFCutLow[]);
  void SetJtPfCHFCutHi(std::vector<double> inJtPfCHFCutHi){jtPfCHFCutHi = inJtPfCHFCutHi; return;}
  void SetJtPfCHFCutHi(int inN, const Double_t inJtPfCHFCutHi[]);
  void SetJtPfCEFCutLow(std::vector<double> inJtPfCEFCutLow){jtPfCEFCutLow = inJtPfCEFCutLow; return;}
  void SetJtPfCEFCutLow(int inN, const Double_t inJtPfCEFCutLow[]);
  void SetJtPfCEFCutHi(std::vector<double> inJtPfCEFCutHi){jtPfCEFCutHi = inJtPfCEFCutHi; return;}
  void SetJtPfCEFCutHi(int inN, const Double_t inJtPfCEFCutHi[]);
  void SetJtPfMinMult(std::vector<int> inJtPfMinMult){jtPfMinMult = inJtPfMinMult; return;}
  void SetJtPfMinMult(int inN, const Int_t inJtPfMinMult[]);
  void SetJtPfMinChgMult(std::vector<int> inJtPfMinChgMult){jtPfMinChgMult = inJtPfMinChgMult; return;}
  void SetJtPfMinChgMult(int inN, const Int_t inJtPfMinChgMult[]);

  void SetNSyst(int inNSyst){nSyst = inNSyst; return;}
  void SetSystStr(std::vector<std::string> inSystStr){systStr = inSystStr; return;};
  void SetSystStr(int inN, const std::string inSystStr[]);

  std::string to_string_with_precision(double a_value, const int n);
};


void cutPropagator::Clean()
{
  inFileNames.clear();
  inFullFileNames.clear();

  isPP = false;

  jtAbsEtaMax = -99;

  nResponseMod = -1;
  responseMod.clear();
  scaleFactor.clear();

  nJtPtBins = -1;
  jtPtBins.clear();
  nJtAbsEtaBins = -1;
  jtAbsEtaBinsLow.clear();
  jtAbsEtaBinsHi.clear();

  nPthats = -1;
  pthats.clear();
  pthatWeights.clear();

  nCentBins = -1;
  centBinsLow.clear();
  centBinsHi.clear();

  nID = -1;
  idStr.clear();
  jtPfCHMFCutLow.clear();
  jtPfCHMFCutHi.clear();
  jtPfMUMFCutLow.clear();
  jtPfMUMFCutHi.clear();
  jtPfNHFCutLow.clear();
  jtPfNHFCutHi.clear();
  jtPfNEFCutLow.clear();
  jtPfNEFCutHi.clear();
  jtPfMUFCutLow.clear();
  jtPfMUFCutHi.clear();
  jtPfCHFCutLow.clear();
  jtPfCHFCutHi.clear();
  jtPfCEFCutLow.clear();
  jtPfCEFCutHi.clear();
  jtPfMinMult.clear();
  jtPfMinChgMult.clear();

  nSyst = -1;
  systStr.clear();

  return;
}

bool cutPropagator::GetAllVarFromFile(TFile* inFile_p)
{
  inFile_p->cd();
  std::vector<std::string> cutDirList = returnRootFileContentsList(inFile_p, "TNamed", "");
  unsigned int pos = 0;
  while(pos < cutDirList.size()){
    if(cutDirList.at(pos).find("cutDir/") != std::string::npos) ++pos;
    else cutDirList.erase(cutDirList.begin()+pos);
  }

  if(cutDirList.size() == 0){
    std::cout << "cutPropagator::GetAllVarFromFile - given inFile_p \'" << inFile_p->GetName() << "\' contains no cutDir. return 1" << std::endl;
    return false;
  }

  for(unsigned int cI = 0; cI < cutDirList.size(); ++cI){
    std::string tempStr = cutDirList.at(cI);
    while(tempStr.find("/") != std::string::npos){tempStr.replace(0, tempStr.find("/")+1, "");}

    if(tempStr.find("nCentBins") != std::string::npos && tempStr.size() == std::string("nCentBins").size()) nCentBins = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("centBinsLow") != std::string::npos && tempStr.size() == std::string("centBinsLow").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        centBinsLow.push_back(std::stoi(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) centBinsLow.push_back(std::stoi(tempStr2));
    }
    else if(tempStr.find("centBinsHi") != std::string::npos && tempStr.size() == std::string("centBinsHi").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        centBinsHi.push_back(std::stoi(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) centBinsHi.push_back(std::stoi(tempStr2));
    }
    else if(tempStr.find("jtAbsEtaMax") != std::string::npos && tempStr.size() == std::string("jtAbsEtaMax").size()) jtAbsEtaMax = std::stof(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("nResponseMod") != std::string::npos && tempStr.size() == std::string("nResponseMod").size()) nResponseMod = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("nJtPtBins") != std::string::npos && tempStr.size() == std::string("nJtPtBins").size()) nJtPtBins = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("nJtAbsEtaBins") != std::string::npos && tempStr.size() == std::string("nJtAbsEtaBins").size()) nJtAbsEtaBins = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("inFileNames") != std::string::npos && tempStr.size() == std::string("inFileNames").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        inFileNames.push_back(tempStr2.substr(0, tempStr2.find(",")));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) inFileNames.push_back(tempStr2);
    }
    else if(tempStr.find("inFullFileNames_") != std::string::npos){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      inFullFileNames.push_back(tempStr2);
    }
    else if(tempStr.find("isPP") != std::string::npos && tempStr.size() == std::string("isPP").size()) isPP = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("nPthats") != std::string::npos && tempStr.size() == std::string("nPthats").size()) nPthats = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("pthats") != std::string::npos && tempStr.size() == std::string("pthats").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        pthats.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) pthats.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("pthatWeights") != std::string::npos && tempStr.size() == std::string("pthatWeights").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        pthatWeights.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) pthatWeights.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("responseMod") != std::string::npos && tempStr.size() == std::string("responseMod").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        responseMod.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) responseMod.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("scaleFactor") != std::string::npos && tempStr.size() == std::string("scaleFactor").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        scaleFactor.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) scaleFactor.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPtBins") != std::string::npos && tempStr.size() == std::string("jtPtBins").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPtBins.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPtBins.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtAbsEtaBinsLow") != std::string::npos && tempStr.size() == std::string("jtAbsEtaBinsLow").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtAbsEtaBinsLow.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtAbsEtaBinsLow.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtAbsEtaBinsHi") != std::string::npos && tempStr.size() == std::string("jtAbsEtaBinsHi").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtAbsEtaBinsHi.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtAbsEtaBinsHi.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("nID") != std::string::npos && tempStr.size() == std::string("nID").size()) nID = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("idStr") != std::string::npos && tempStr.size() == std::string("idStr").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        idStr.push_back(tempStr2.substr(0, tempStr2.find(",")));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) idStr.push_back(tempStr2);
    }
    else if(tempStr.find("jtPfCHMFCutLow") != std::string::npos && tempStr.size() == std::string("jtPfCHMFCutLow").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfCHMFCutLow.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfCHMFCutLow.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfCHMFCutHi") != std::string::npos && tempStr.size() == std::string("jtPfCHMFCutHi").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfCHMFCutHi.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfCHMFCutHi.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfMUMFCutLow") != std::string::npos && tempStr.size() == std::string("jtPfMUMFCutLow").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfMUMFCutLow.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfMUMFCutLow.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfMUMFCutHi") != std::string::npos && tempStr.size() == std::string("jtPfMUMFCutHi").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfMUMFCutHi.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfMUMFCutHi.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfNHFCutLow") != std::string::npos && tempStr.size() == std::string("jtPfNHFCutLow").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfNHFCutLow.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfNHFCutLow.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfNHFCutHi") != std::string::npos && tempStr.size() == std::string("jtPfNHFCutHi").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfNHFCutHi.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfNHFCutHi.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfNEFCutLow") != std::string::npos && tempStr.size() == std::string("jtPfNEFCutLow").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfNEFCutLow.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfNEFCutLow.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfNEFCutHi") != std::string::npos && tempStr.size() == std::string("jtPfNEFCutHi").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfNEFCutHi.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfNEFCutHi.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfMUFCutLow") != std::string::npos && tempStr.size() == std::string("jtPfMUFCutLow").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfMUFCutLow.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfMUFCutLow.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfMUFCutHi") != std::string::npos && tempStr.size() == std::string("jtPfMUFCutHi").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfMUFCutHi.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfMUFCutHi.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfCHFCutLow") != std::string::npos && tempStr.size() == std::string("jtPfCHFCutLow").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfCHFCutLow.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfCHFCutLow.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfCHFCutHi") != std::string::npos && tempStr.size() == std::string("jtPfCHFCutHi").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfCHFCutHi.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfCHFCutHi.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfCEFCutLow") != std::string::npos && tempStr.size() == std::string("jtPfCEFCutLow").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfCEFCutLow.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfCEFCutLow.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfCEFCutHi") != std::string::npos && tempStr.size() == std::string("jtPfCEFCutHi").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfCEFCutHi.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfCEFCutHi.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfMinMult") != std::string::npos && tempStr.size() == std::string("jtPfMinMult").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfMinMult.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfMinMult.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtPfMinChgMult") != std::string::npos && tempStr.size() == std::string("jtPfMinChgMult").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPfMinChgMult.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPfMinChgMult.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("nSyst") != std::string::npos && tempStr.size() == std::string("nSyst").size()) nSyst = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("systStr") != std::string::npos && tempStr.size() == std::string("systStr").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        systStr.push_back(tempStr2.substr(0, tempStr2.find(",")));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) systStr.push_back(tempStr2);
    }
  }

  return true;
}


bool cutPropagator::WriteAllVarToFile(TFile* inFile_p, TDirectory* inDir_p, TDirectory* inSubDir_p)
{
  if(inDir_p == NULL){
    std::cout << "cutPropafator::WriteAllVarToFile - Given indir is null. return false" << std::endl;
    return false;
  }
  else if(inSubDir_p == NULL){
    std::cout << "cutPropafator::WriteAllVarToFile - Given insubdir is null. return false" << std::endl;
    return false;
  }
  else if(std::string(inDir_p->GetName()).size() != std::string("cutDir").size()){
    std::cout << "cutPropafator::WriteAllVarToFile - Given indir is not \'cutDir\'. return false" << std::endl;
    return false;
  }
  else if(std::string(inDir_p->GetName()).find("cutDir") == std::string::npos){
    std::cout << "cutPropafator::WriteAllVarToFile - Given indir is not \'cutDir\'. return false" << std::endl;
    return false;
  }

  inFile_p->cd();
  inDir_p->cd();

  std::string inFileNames2 = "";
  for(unsigned int i = 0; i < inFileNames.size(); ++i){
    inFileNames2 = inFileNames2 + inFileNames.at(i) + ",";
  }

  std::string responseModStr = "";
  std::string scaleFactorStr = "";
  for(int jI = 0; jI < nResponseMod; ++jI){
    responseModStr = responseModStr + std::to_string(responseMod.at(jI)) + ",";
    scaleFactorStr = scaleFactorStr + std::to_string(scaleFactor.at(jI)) + ",";
  }

  std::string jtPtBinsStr = "";
  for(int jI = 0; jI < nJtPtBins+1; ++jI){
    jtPtBinsStr = jtPtBinsStr + std::to_string(jtPtBins.at(jI)) + ",";
  }

  std::string jtAbsEtaBinsLowStr = "";
  std::string jtAbsEtaBinsHiStr = "";
  for(int jI = 0; jI < nJtAbsEtaBins; ++jI){
    jtAbsEtaBinsLowStr = jtAbsEtaBinsLowStr + std::to_string(jtAbsEtaBinsLow.at(jI)) + ",";
    jtAbsEtaBinsHiStr = jtAbsEtaBinsHiStr + std::to_string(jtAbsEtaBinsHi.at(jI)) + ",";
  }

  std::string nPthatsStr = std::to_string(nPthats);
  std::string pthatsStr = "";
  for(unsigned int jI = 0; jI < pthats.size(); ++jI){
    pthatsStr = pthatsStr + std::to_string(pthats.at(jI)) + ",";
  }

  std::string pthatWeightsStr = "";
  for(unsigned int jI = 0; jI < pthatWeights.size(); ++jI){
    pthatWeightsStr = pthatWeightsStr + to_string_with_precision(pthatWeights.at(jI), 12) + ",";
    std::cout << "Weight orig, string: " << pthatWeights.at(jI) << ", " << to_string_with_precision(pthatWeights.at(jI), 12) << std::endl;
  }

  std::string centBinsLowStr = "";
  for(int jI = 0; jI < nCentBins; ++jI){
    centBinsLowStr = centBinsLowStr + std::to_string(centBinsLow.at(jI)) + ",";
  }

  std::string centBinsHiStr = "";
  for(int jI = 0; jI < nCentBins; ++jI){
    centBinsHiStr = centBinsHiStr + std::to_string(centBinsHi.at(jI)) + ",";
  }

  std::string nIDStr = std::to_string(nID);
  std::string idStr2 = "";
  std::string jtPfCHMFCutLowStr = "";
  std::string jtPfCHMFCutHiStr = "";
  std::string jtPfMUMFCutLowStr = "";
  std::string jtPfMUMFCutHiStr = "";
  std::string jtPfNHFCutLowStr = "";
  std::string jtPfNHFCutHiStr = "";
  std::string jtPfNEFCutLowStr = "";
  std::string jtPfNEFCutHiStr = "";
  std::string jtPfMUFCutLowStr = "";
  std::string jtPfMUFCutHiStr = "";
  std::string jtPfCHFCutLowStr = "";
  std::string jtPfCHFCutHiStr = "";
  std::string jtPfCEFCutLowStr = "";
  std::string jtPfCEFCutHiStr = "";
  std::string jtPfMinMultStr = "";
  std::string jtPfMinChgMultStr = "";

  for(int iI = 0; iI < nID; ++iI){
    idStr2 = idStr2 + idStr.at(iI) + ",";
    jtPfCHMFCutLowStr = jtPfCHMFCutLowStr + std::to_string(jtPfCHMFCutLow.at(iI)) + ",";
    jtPfCHMFCutHiStr = jtPfCHMFCutHiStr + std::to_string(jtPfCHMFCutHi.at(iI)) + ",";
    jtPfMUMFCutLowStr = jtPfMUMFCutLowStr + std::to_string(jtPfMUMFCutLow.at(iI)) + ",";
    jtPfMUMFCutHiStr = jtPfMUMFCutHiStr + std::to_string(jtPfMUMFCutHi.at(iI)) + ",";
    jtPfNHFCutLowStr = jtPfNHFCutLowStr + std::to_string(jtPfNHFCutLow.at(iI)) + ",";
    jtPfNHFCutHiStr = jtPfNHFCutHiStr + std::to_string(jtPfNHFCutHi.at(iI)) + ",";
    jtPfNEFCutLowStr = jtPfNEFCutLowStr + std::to_string(jtPfNEFCutLow.at(iI)) + ",";
    jtPfNEFCutHiStr = jtPfNEFCutHiStr + std::to_string(jtPfNEFCutHi.at(iI)) + ",";
    jtPfMUFCutLowStr = jtPfMUFCutLowStr + std::to_string(jtPfMUFCutLow.at(iI)) + ",";
    jtPfMUFCutHiStr = jtPfMUFCutHiStr + std::to_string(jtPfMUFCutHi.at(iI)) + ",";
    jtPfCHFCutLowStr = jtPfCHFCutLowStr + std::to_string(jtPfCHFCutLow.at(iI)) + ",";
    jtPfCHFCutHiStr = jtPfCHFCutHiStr + std::to_string(jtPfCHFCutHi.at(iI)) + ",";
    jtPfCEFCutLowStr = jtPfCEFCutLowStr + std::to_string(jtPfCEFCutLow.at(iI)) + ",";
    jtPfCEFCutHiStr = jtPfCEFCutHiStr + std::to_string(jtPfCEFCutHi.at(iI)) + ",";
    jtPfMinMultStr = jtPfMinMultStr + std::to_string(jtPfMinMult.at(iI)) + ",";
    jtPfMinChgMultStr = jtPfMinChgMultStr + std::to_string(jtPfMinChgMult.at(iI)) + ",";
  }

  std::string nSystStr = std::to_string(nSyst);
  std::string systStr2 = "";

  for(int sI = 0; sI < nSyst; ++sI){
    systStr2 = systStr2 + systStr.at(sI) + ",";
  }

  TNamed inFileNamesName("inFileNames", inFileNames2.c_str());
  TNamed isPPName("isPP", std::to_string(isPP));
  TNamed jtAbsEtaMaxName("jtAbsEtaMax", std::to_string(jtAbsEtaMax).c_str());
  TNamed nResponseModName("nResponseMod", std::to_string(nResponseMod).c_str());
  TNamed responseModName("responseMod", responseModStr.c_str());
  TNamed scaleFactorName("scaleFactor", scaleFactorStr.c_str());
  TNamed nJtPtBinsName("nJtPtBins", std::to_string(nJtPtBins).c_str());
  TNamed jtPtBinsName("jtPtBins", jtPtBinsStr.c_str());
  TNamed nJtAbsEtaBinsName("nJtAbsEtaBins", std::to_string(nJtAbsEtaBins).c_str());
  TNamed jtAbsEtaBinsLowName("jtAbsEtaBinsLow", jtAbsEtaBinsLowStr.c_str());
  TNamed jtAbsEtaBinsHiName("jtAbsEtaBinsHi", jtAbsEtaBinsHiStr.c_str());
  TNamed nPthatsName("nPthats", nPthatsStr.c_str());
  TNamed pthatsName("pthats", pthatsStr.c_str());
  TNamed pthatWeightsName("pthatWeights", pthatWeightsStr.c_str());
  TNamed nCentBinsName("nCentBins", std::to_string(nCentBins).c_str());
  TNamed centBinsLowName("centBinsLow", centBinsLowStr.c_str());
  TNamed centBinsHiName("centBinsHi", centBinsHiStr.c_str());
  TNamed nIDName("nID", nIDStr.c_str());
  TNamed idStrName("idStr", idStr2.c_str());
  TNamed jtPfCHMFCutLowName("jtPfCHMFCutLow", jtPfCHMFCutLowStr.c_str());
  TNamed jtPfCHMFCutHiName("jtPfCHMFCutHi", jtPfCHMFCutHiStr.c_str());
  TNamed jtPfMUMFCutLowName("jtPfMUMFCutLow", jtPfMUMFCutLowStr.c_str());
  TNamed jtPfMUMFCutHiName("jtPfMUMFCutHi", jtPfMUMFCutHiStr.c_str());
  TNamed jtPfNHFCutLowName("jtPfNHFCutLow", jtPfNHFCutLowStr.c_str());
  TNamed jtPfNHFCutHiName("jtPfNHFCutHi", jtPfNHFCutHiStr.c_str());
  TNamed jtPfNEFCutLowName("jtPfNEFCutLow", jtPfNEFCutLowStr.c_str());
  TNamed jtPfNEFCutHiName("jtPfNEFCutHi", jtPfNEFCutHiStr.c_str());
  TNamed jtPfMUFCutLowName("jtPfMUFCutLow", jtPfMUFCutLowStr.c_str());
  TNamed jtPfMUFCutHiName("jtPfMUFCutHi", jtPfMUFCutHiStr.c_str());
  TNamed jtPfCHFCutLowName("jtPfCHFCutLow", jtPfCHFCutLowStr.c_str());
  TNamed jtPfCHFCutHiName("jtPfCHFCutHi", jtPfCHFCutHiStr.c_str());
  TNamed jtPfCEFCutLowName("jtPfCEFCutLow", jtPfCEFCutLowStr.c_str());
  TNamed jtPfCEFCutHiName("jtPfCEFCutHi", jtPfCEFCutHiStr.c_str());
  TNamed jtPfMinMultName("jtPfMinMult", jtPfMinMultStr.c_str());
  TNamed jtPfMinChgMultName("jtPfMinChgMult", jtPfMinChgMultStr.c_str());
  TNamed nSystName("nSyst", nSystStr.c_str());
  TNamed systStrName("systStr", systStr2.c_str());

  inFileNamesName.Write("", TObject::kOverwrite);
  isPPName.Write("", TObject::kOverwrite);
  jtAbsEtaMaxName.Write("", TObject::kOverwrite);
  nResponseModName.Write("", TObject::kOverwrite);
  responseModName.Write("", TObject::kOverwrite);
  scaleFactorName.Write("", TObject::kOverwrite);
  nJtPtBinsName.Write("", TObject::kOverwrite);
  jtPtBinsName.Write("", TObject::kOverwrite);
  nJtAbsEtaBinsName.Write("", TObject::kOverwrite);
  jtAbsEtaBinsLowName.Write("", TObject::kOverwrite);
  jtAbsEtaBinsHiName.Write("", TObject::kOverwrite);
  nPthatsName.Write("", TObject::kOverwrite);
  pthatsName.Write("", TObject::kOverwrite);
  pthatWeightsName.Write("", TObject::kOverwrite);
  nCentBinsName.Write("", TObject::kOverwrite);
  centBinsLowName.Write("", TObject::kOverwrite);
  centBinsHiName.Write("", TObject::kOverwrite);
  nIDName.Write("", TObject::kOverwrite);
  idStrName.Write("", TObject::kOverwrite);
  jtPfCHMFCutLowName.Write("", TObject::kOverwrite);
  jtPfCHMFCutHiName.Write("", TObject::kOverwrite);
  jtPfMUMFCutLowName.Write("", TObject::kOverwrite);
  jtPfMUMFCutHiName.Write("", TObject::kOverwrite);
  jtPfNHFCutLowName.Write("", TObject::kOverwrite);
  jtPfNHFCutHiName.Write("", TObject::kOverwrite);
  jtPfNEFCutLowName.Write("", TObject::kOverwrite);
  jtPfNEFCutHiName.Write("", TObject::kOverwrite);
  jtPfMUFCutLowName.Write("", TObject::kOverwrite);
  jtPfMUFCutHiName.Write("", TObject::kOverwrite);
  jtPfCHFCutLowName.Write("", TObject::kOverwrite);
  jtPfCHFCutHiName.Write("", TObject::kOverwrite);
  jtPfCEFCutLowName.Write("", TObject::kOverwrite);
  jtPfCEFCutHiName.Write("", TObject::kOverwrite);
  jtPfCEFCutHiName.Write("", TObject::kOverwrite);
  jtPfMinMultName.Write("", TObject::kOverwrite);
  jtPfMinChgMultName.Write("", TObject::kOverwrite);
  nSystName.Write("", TObject::kOverwrite);
  systStrName.Write("", TObject::kOverwrite);


  inFile_p->cd();
  inDir_p->cd();
  inSubDir_p->cd();  

  std::vector<TNamed> inFullFileNames2;
  for(unsigned int i = 0; i < inFullFileNames.size(); ++i){
    inFullFileNames2.push_back(TNamed(("inFullFileNames_" + std::to_string(i)).c_str(), inFullFileNames.at(i).c_str()));
  }

  for(unsigned int i = 0; i < inFullFileNames2.size(); ++i){
    inFullFileNames2.at(i).Write("", TObject::kOverwrite);
  }


  return true;
}


bool cutPropagator::CheckPropagatorsMatch(cutPropagator inCutProp, bool doBothMCOrBothData, bool doBothPPOrBothPbPb)
{
  const double delta = 0.0001;
  if(jtAbsEtaMax - delta > inCutProp.GetJtAbsEtaMax()) return false;
  if(jtAbsEtaMax + delta < inCutProp.GetJtAbsEtaMax()) return false;
  if(nJtPtBins != inCutProp.GetNJtPtBins()) return false;
  if(nJtAbsEtaBins != inCutProp.GetNJtAbsEtaBins()) return false;
  if(nID != inCutProp.GetNID()) return false;
  if(jtPtBins.size() != inCutProp.GetJtPtBins().size()) return false;
  if(jtAbsEtaBinsLow.size() != inCutProp.GetJtAbsEtaBinsLow().size()) return false;
  if(jtAbsEtaBinsHi.size() != inCutProp.GetJtAbsEtaBinsHi().size()) return false;
  if(idStr.size() != inCutProp.GetIdStr().size()) return false;
  if(jtPfCHMFCutLow.size() != inCutProp.GetJtPfCHMFCutLow().size()) return false;
  if(jtPfCHMFCutHi.size() != inCutProp.GetJtPfCHMFCutHi().size()) return false;
  if(jtPfMUMFCutLow.size() != inCutProp.GetJtPfMUMFCutLow().size()) return false;
  if(jtPfMUMFCutHi.size() != inCutProp.GetJtPfMUMFCutHi().size()) return false;
  if(jtPfNHFCutLow.size() != inCutProp.GetJtPfNHFCutLow().size()) return false;
  if(jtPfNHFCutHi.size() != inCutProp.GetJtPfNHFCutHi().size()) return false;
  if(jtPfNEFCutLow.size() != inCutProp.GetJtPfNEFCutLow().size()) return false;
  if(jtPfNEFCutHi.size() != inCutProp.GetJtPfNEFCutHi().size()) return false;
  if(jtPfMUFCutLow.size() != inCutProp.GetJtPfMUFCutLow().size()) return false;
  if(jtPfMUFCutHi.size() != inCutProp.GetJtPfMUFCutHi().size()) return false;
  if(jtPfCHFCutLow.size() != inCutProp.GetJtPfCHFCutLow().size()) return false;
  if(jtPfCHFCutHi.size() != inCutProp.GetJtPfCHFCutHi().size()) return false;
  if(jtPfCEFCutLow.size() != inCutProp.GetJtPfCEFCutLow().size()) return false;
  if(jtPfCEFCutHi.size() != inCutProp.GetJtPfCEFCutHi().size()) return false;
  if(jtPfMinMult.size() != inCutProp.GetJtPfMinMult().size()) return false;
  if(jtPfMinChgMult.size() != inCutProp.GetJtPfMinChgMult().size()) return false;
  if(nSyst != inCutProp.GetNSyst()) return false;
  if(systStr.size() != inCutProp.GetSystStr().size()) return false;

  for(unsigned int i = 0; i < jtPtBins.size(); ++i){
    if(jtPtBins.at(i) - delta > inCutProp.GetJtPtBins().at(i)) return false;
    if(jtPtBins.at(i) + delta < inCutProp.GetJtPtBins().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtAbsEtaBinsLow.size(); ++i){
    if(jtAbsEtaBinsLow.at(i) - delta > inCutProp.GetJtAbsEtaBinsLow().at(i)) return false;
    if(jtAbsEtaBinsLow.at(i) + delta < inCutProp.GetJtAbsEtaBinsLow().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtAbsEtaBinsHi.size(); ++i){
    if(jtAbsEtaBinsHi.at(i) - delta > inCutProp.GetJtAbsEtaBinsHi().at(i)) return false;
    if(jtAbsEtaBinsHi.at(i) + delta < inCutProp.GetJtAbsEtaBinsHi().at(i)) return false;
  }

  for(unsigned int i = 0; i < idStr.size(); ++i){
    if(idStr.at(i).size() != inCutProp.GetIdStr().at(i).size()) return false;
    if(idStr.at(i).find(inCutProp.GetIdStr().at(i)) == std::string::npos) return false;
  }

  for(unsigned int i = 0; i < jtPfCHMFCutLow.size(); ++i){
    if(jtPfCHMFCutLow.at(i) - delta > inCutProp.GetJtPfCHMFCutLow().at(i)) return false;
    if(jtPfCHMFCutLow.at(i) + delta < inCutProp.GetJtPfCHMFCutLow().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfCHMFCutHi.size(); ++i){
    if(jtPfCHMFCutHi.at(i) - delta > inCutProp.GetJtPfCHMFCutHi().at(i)) return false;
    if(jtPfCHMFCutHi.at(i) + delta < inCutProp.GetJtPfCHMFCutHi().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfMUMFCutLow.size(); ++i){
    if(jtPfMUMFCutLow.at(i) - delta > inCutProp.GetJtPfMUMFCutLow().at(i)) return false;
    if(jtPfMUMFCutLow.at(i) + delta < inCutProp.GetJtPfMUMFCutLow().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfMUMFCutHi.size(); ++i){
    if(jtPfMUMFCutHi.at(i) - delta > inCutProp.GetJtPfMUMFCutHi().at(i)) return false;
    if(jtPfMUMFCutHi.at(i) + delta < inCutProp.GetJtPfMUMFCutHi().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfNHFCutLow.size(); ++i){
    if(jtPfNHFCutLow.at(i) - delta > inCutProp.GetJtPfNHFCutLow().at(i)) return false;
    if(jtPfNHFCutLow.at(i) + delta < inCutProp.GetJtPfNHFCutLow().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfNHFCutHi.size(); ++i){
    if(jtPfNHFCutHi.at(i) - delta > inCutProp.GetJtPfNHFCutHi().at(i)) return false;
    if(jtPfNHFCutHi.at(i) + delta < inCutProp.GetJtPfNHFCutHi().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfNEFCutLow.size(); ++i){
    if(jtPfNEFCutLow.at(i) - delta > inCutProp.GetJtPfNEFCutLow().at(i)) return false;
    if(jtPfNEFCutLow.at(i) + delta < inCutProp.GetJtPfNEFCutLow().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfNEFCutHi.size(); ++i){
    if(jtPfNEFCutHi.at(i) - delta > inCutProp.GetJtPfNEFCutHi().at(i)) return false;
    if(jtPfNEFCutHi.at(i) + delta < inCutProp.GetJtPfNEFCutHi().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfMUFCutLow.size(); ++i){
    if(jtPfMUFCutLow.at(i) - delta > inCutProp.GetJtPfMUFCutLow().at(i)) return false;
    if(jtPfMUFCutLow.at(i) + delta < inCutProp.GetJtPfMUFCutLow().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfMUFCutHi.size(); ++i){
    if(jtPfMUFCutHi.at(i) - delta > inCutProp.GetJtPfMUFCutHi().at(i)) return false;
    if(jtPfMUFCutHi.at(i) + delta < inCutProp.GetJtPfMUFCutHi().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfCHFCutLow.size(); ++i){
    if(jtPfCHFCutLow.at(i) - delta > inCutProp.GetJtPfCHFCutLow().at(i)) return false;
    if(jtPfCHFCutLow.at(i) + delta < inCutProp.GetJtPfCHFCutLow().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfCHFCutHi.size(); ++i){
    if(jtPfCHFCutHi.at(i) - delta > inCutProp.GetJtPfCHFCutHi().at(i)) return false;
    if(jtPfCHFCutHi.at(i) + delta < inCutProp.GetJtPfCHFCutHi().at(i)) return false;
  }


  for(unsigned int i = 0; i < jtPfCEFCutLow.size(); ++i){
    if(jtPfCEFCutLow.at(i) - delta > inCutProp.GetJtPfCEFCutLow().at(i)) return false;
    if(jtPfCEFCutLow.at(i) + delta < inCutProp.GetJtPfCEFCutLow().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfCEFCutHi.size(); ++i){
    if(jtPfCEFCutHi.at(i) - delta > inCutProp.GetJtPfCEFCutHi().at(i)) return false;
    if(jtPfCEFCutHi.at(i) + delta < inCutProp.GetJtPfCEFCutHi().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfMinMult.size(); ++i){
    if(jtPfMinMult.at(i) != inCutProp.GetJtPfMinMult().at(i)) return false;
  }

  for(unsigned int i = 0; i < jtPfMinChgMult.size(); ++i){
    if(jtPfMinChgMult.at(i) != inCutProp.GetJtPfMinChgMult().at(i)) return false;
  }

  for(unsigned int i = 0; i < systStr.size(); ++i){
    if(systStr.at(i).size() != inCutProp.GetSystStr().at(i).size()) return false;
    if(systStr.at(i).find(inCutProp.GetSystStr().at(i)) == std::string::npos) return false;
  }


  if(doBothMCOrBothData){
    if(nResponseMod != inCutProp.GetNResponseMod()) return false;
    if(responseMod.size() != inCutProp.GetResponseMod().size()) return false;
    if(scaleFactor.size() != inCutProp.GetScaleFactor().size()) return false;
    if(nPthats != inCutProp.GetNPthats()) return false;
    if(pthats.size() != inCutProp.GetPthats().size()) return false;
    if(pthatWeights.size() != inCutProp.GetPthatWeights().size()) return false;

    for(unsigned int i = 0; i < responseMod.size(); ++i){
      if(responseMod.at(i) - delta > inCutProp.GetResponseMod().at(i)) return false;
      if(responseMod.at(i) + delta < inCutProp.GetResponseMod().at(i)) return false;
    }

    for(unsigned int i = 0; i < scaleFactor.size(); ++i){
      if(scaleFactor.at(i) - delta > inCutProp.GetScaleFactor().at(i)) return false;
      if(scaleFactor.at(i) + delta < inCutProp.GetScaleFactor().at(i)) return false;
    }

    for(unsigned int i = 0; i < pthats.size(); ++i){
      if(pthats.at(i) - delta > inCutProp.GetPthats().at(i)) return false;
      if(pthats.at(i) + delta < inCutProp.GetPthats().at(i)) return false;
    }

    for(unsigned int i = 0; i < pthatWeights.size(); ++i){
      if(pthatWeights.at(i) - delta > inCutProp.GetPthatWeights().at(i)) return false;
      if(pthatWeights.at(i) + delta < inCutProp.GetPthatWeights().at(i)) return false;
    }
  }

  if(doBothPPOrBothPbPb){
    if(isPP != inCutProp.GetIsPP()) return false;
    if(nCentBins != inCutProp.GetNCentBins()) return false;
    if(centBinsLow.size() != inCutProp.GetCentBinsLow().size()) return false;
    if(centBinsHi.size() != inCutProp.GetCentBinsHi().size()) return false;

    for(unsigned int i = 0; i < centBinsLow.size(); ++i){
      if(centBinsLow.at(i) != inCutProp.GetCentBinsLow().at(i)) return false;
    }
    
    for(unsigned int i = 0; i < centBinsHi.size(); ++i){
      if(centBinsHi.at(i) != inCutProp.GetCentBinsHi().at(i)) return false;
    }    
  }

  return true;
}


void cutPropagator::SetResponseMod(int inN, const Double_t inResponseMod[])
{
  for(int i = 0; i < inN; ++i){
    responseMod.push_back(inResponseMod[i]);
  }

  return;
}

void cutPropagator::SetScaleFactor(int inN, const Double_t inScaleFactor[])
{
  for(int i = 0; i < inN; ++i){
    scaleFactor.push_back(inScaleFactor[i]);
  }

  return;
}

void cutPropagator::SetJtPtBins(int inN, const Double_t inJtPtBins[])
{
  for(int i = 0; i < inN; ++i){
    jtPtBins.push_back(inJtPtBins[i]);
  }

  return;
}

void cutPropagator::SetJtAbsEtaBinsLow(int inN, const Double_t inJtAbsEtaBinsLow[])
{
  for(int i = 0; i < inN; ++i){
    jtAbsEtaBinsLow.push_back(inJtAbsEtaBinsLow[i]);
  }

  return;
}


void cutPropagator::SetJtAbsEtaBinsHi(int inN, const Double_t inJtAbsEtaBinsHi[])
{
  for(int i = 0; i < inN; ++i){
    jtAbsEtaBinsHi.push_back(inJtAbsEtaBinsHi[i]);
  }

  return;
}


void cutPropagator::SetPthats(int inN, const Double_t inPthats[])
{
  for(int i = 0; i < inN; ++i){
    pthats.push_back(inPthats[i]);
  }

  return;
}


void cutPropagator::SetPthatWeights(int inN, const Double_t inPthatWeights[])
{
  for(int i = 0; i < inN; ++i){
    pthatWeights.push_back(inPthatWeights[i]);
  }

  return;
}


void cutPropagator::SetCentBinsLow(int inN, const Int_t inCentBinsLow[])
{
  for(int i = 0; i < inN; ++i){
    centBinsLow.push_back(inCentBinsLow[i]);
  }

  return;
}


void cutPropagator::SetCentBinsHi(int inN, const Int_t inCentBinsHi[])
{
  for(int i = 0; i < inN; ++i){
    centBinsHi.push_back(inCentBinsHi[i]);
  }

  return;
}


void cutPropagator::SetIdStr(int inN, const std::string inIdStr[])
{
  for(int i = 0; i < inN; ++i){
    idStr.push_back(inIdStr[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHMFCutLow(int inN, const Double_t inJtPfCHMFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHMFCutLow.push_back(inJtPfCHMFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHMFCutHi(int inN, const Double_t inJtPfCHMFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHMFCutHi.push_back(inJtPfCHMFCutHi[i]);
  }

  return;
}


void cutPropagator::SetJtPfMUMFCutLow(int inN, const Double_t inJtPfMUMFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUMFCutLow.push_back(inJtPfMUMFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfMUMFCutHi(int inN, const Double_t inJtPfMUMFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUMFCutHi.push_back(inJtPfMUMFCutHi[i]);
  }

  return;
}


void cutPropagator::SetJtPfNHFCutLow(int inN, const Double_t inJtPfNHFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNHFCutLow.push_back(inJtPfNHFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfNHFCutHi(int inN, const Double_t inJtPfNHFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNHFCutHi.push_back(inJtPfNHFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfNEFCutLow(int inN, const Double_t inJtPfNEFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNEFCutLow.push_back(inJtPfNEFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfNEFCutHi(int inN, const Double_t inJtPfNEFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNEFCutHi.push_back(inJtPfNEFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfMUFCutLow(int inN, const Double_t inJtPfMUFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUFCutLow.push_back(inJtPfMUFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfMUFCutHi(int inN, const Double_t inJtPfMUFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUFCutHi.push_back(inJtPfMUFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHFCutLow(int inN, const Double_t inJtPfCHFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHFCutLow.push_back(inJtPfCHFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHFCutHi(int inN, const Double_t inJtPfCHFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHFCutHi.push_back(inJtPfCHFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfCEFCutLow(int inN, const Double_t inJtPfCEFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCEFCutLow.push_back(inJtPfCEFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfCEFCutHi(int inN, const Double_t inJtPfCEFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCEFCutHi.push_back(inJtPfCEFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfMinMult(int inN, const Int_t inJtPfMinMult[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMinMult.push_back(inJtPfMinMult[i]);
  }

  return;
}

void cutPropagator::SetJtPfMinChgMult(int inN, const Int_t inJtPfMinChgMult[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMinChgMult.push_back(inJtPfMinChgMult[i]);
  }

  return;
}

void cutPropagator::SetSystStr(int inN, const std::string inSystStr[])
{
  for(int i = 0; i < inN; ++i){
    systStr.push_back(inSystStr[i]);
  }

  return;
}

std::string cutPropagator::to_string_with_precision(double a_value, const int n)
{
  std::ostringstream out;
  out << std::setprecision(n) << a_value;
  return out.str();
}
