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
  bool isPP;

  double jtAbsEtaMax;
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

  void Clean();
  bool GetAllVarFromFile(TFile* inFile_p);
  bool WriteAllVarToFile(TFile* inFile_p, TDirectory* inDir_p);
  bool CheckPropagatorMatch(cutPropagator inCutProp, bool doBothMCOrBothData, bool doBothPPOrBothPbPb);

  bool GetIsPP(){return isPP;}

  double GetJtAbsEtaMax(){return jtAbsEtaMax;}
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


  void SetIsPP(bool inIsPP){isPP = inIsPP; return;}

  void SetJtAbsEtaMax(double inJtAbsEtaMax){jtAbsEtaMax = inJtAbsEtaMax; return;}
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

  std::string to_string_with_precision(double a_value, const int n);
};


void cutPropagator::Clean()
{
  isPP = false;

  jtAbsEtaMax = -99;
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
    else if(tempStr.find("nJtPtBins") != std::string::npos && tempStr.size() == std::string("nJtPtBins").size()) nJtPtBins = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("nJtAbsEtaBins") != std::string::npos && tempStr.size() == std::string("nJtAbsEtaBins").size()) nJtAbsEtaBins = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("isPP") != std::string::npos && tempStr.size() == std::string("isPP").size()) isPP = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
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
  }

  return true;
}


bool cutPropagator::WriteAllVarToFile(TFile* inFile_p, TDirectory* inDir_p)
{
  if(inDir_p == NULL){
    std::cout << "cutPropafator::WriteAllVarToFile - Given indir is null. return false" << std::endl;
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

  for(int iI = 0; iI < nID; ++iI){
    idStr2 = idStr2 + idStr.at(iI) + ",";
    jtPfCHMFCutLowStr = jtPfCHMFCutLowStr + std::to_string(jtPfCHMFCutLow.at(iI)) + ",";
    jtPfCHMFCutHiStr = jtPfCHMFCutHiStr + std::to_string(jtPfCHMFCutHi.at(iI)) + ",";
    jtPfMUMFCutLowStr = jtPfMUMFCutLowStr + std::to_string(jtPfMUMFCutLow.at(iI)) + ",";
    jtPfMUMFCutHiStr = jtPfMUMFCutHiStr + std::to_string(jtPfMUMFCutHi.at(iI)) + ",";
  }

  TNamed isPPName("isPP", std::to_string(isPP));
  TNamed jtAbsEtaMaxName("jtAbsEtaMax", std::to_string(jtAbsEtaMax).c_str());
  TNamed nJtPtBinsName("nJtPtBins", std::to_string(nJtPtBins).c_str());
  TNamed jtPtBinsName("jtPtBins", jtPtBinsStr.c_str());
  TNamed nJtAbsEtaBinsName("nJtAbsEtaBins", std::to_string(nJtAbsEtaBins).c_str());
  TNamed jtAbsEtaBinsLowName("jtAbsEtaBinsLow", jtAbsEtaBinsLowStr.c_str());
  TNamed jtAbsEtaBinsHiName("jtAbsEtaBinsHi", jtAbsEtaBinsHiStr.c_str());
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

  isPPName.Write("", TObject::kOverwrite);
  jtAbsEtaMaxName.Write("", TObject::kOverwrite);
  nJtPtBinsName.Write("", TObject::kOverwrite);
  jtPtBinsName.Write("", TObject::kOverwrite);
  nJtAbsEtaBinsName.Write("", TObject::kOverwrite);
  jtAbsEtaBinsLowName.Write("", TObject::kOverwrite);
  jtAbsEtaBinsHiName.Write("", TObject::kOverwrite);
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

  if(doBothMCOrBothData){
    if(nPthats != inCutProp.GetNPthats()) return false;
    if(pthats.size() != inCutProp.GetPthats().size()) return false;
    if(pthatWeights.size() != inCutProp.GetPthatWeights().size()) return false;

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


std::string cutPropagator::to_string_with_precision(double a_value, const int n)
{
  std::ostringstream out;
  out << std::setprecision(n) << a_value;
  return out.str();
}
