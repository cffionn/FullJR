#ifndef JETDATA_H
#define JETDATA_H

#include <string>
#include <vector>

#include "TTree.h"

class jetData{
 public:
  static const int nMaxJets = 500;

  int hiBin_ = 0;
  float fullWeight_ = 0.;
  float pthatWeight_ = 0.;
  int nref_ = 0.;
  float jtpt_[nMaxJets] = {0};
  float rawpt_[nMaxJets] = {0};
  float refpt_[nMaxJets] = {0};
  float jteta_[nMaxJets] = {0};

  static const int nVar = 8;
  std::string varStr[nVar] = {"hiBin",
			      "fullWeight",
			      "pthatWeight",
			      "nref",
			      "jtpt",
			      "rawpt",
			      "refpt",
			      "jteta"};

  bool varIsGood[nVar] = {0};

  jetData();
  void SetStatusAndAddressRead(TTree* inTree_p,  std::vector<std::string> inList);
};

jetData::jetData()
{
  hiBin_ = -999;
  fullWeight_ = -999.;
  pthatWeight_ = -999.;
  nref_ = -999;
  
  for(Int_t i = 0; i < nMaxJets; ++i){
    jtpt_[i] = -999.;
    rawpt_[i] = -999.;
    refpt_[i] = -999.;
    jteta_[i] = -999.;
  }

  for(Int_t i = 0; i < nVar; ++i){varIsGood[i] = false;}
  return;
}

void jetData::SetStatusAndAddressRead(TTree* inTree_p,  std::vector<std::string> inList)
{
  if(inList.size() != 0){
    for(int i = 0; i < nVar; ++i){varIsGood[i] = false;}

    for(unsigned int i = 0; i < inList.size(); ++i){

      for(Int_t j = 0; j < nVar; ++j){
        if(inList.at(i).size() == varStr[j].size() && inList.at(i).find(varStr[j]) != std::string::npos){
          varIsGood[j] = true;
          break;
        }
      }
    }
  }

  for(Int_t i = 0; i < nVar; ++i){if(varIsGood[i]) inTree_p->SetBranchStatus(varStr[i].c_str(), 1);}

  if(varIsGood[0]) inTree_p->SetBranchAddress("hiBin", &hiBin_);
  if(varIsGood[1]) inTree_p->SetBranchAddress("fullWeight", &fullWeight_);
  if(varIsGood[2]) inTree_p->SetBranchAddress("pthatWeight", &pthatWeight_);
  if(varIsGood[3]) inTree_p->SetBranchAddress("nref", &nref_);
  if(varIsGood[4]) inTree_p->SetBranchAddress("jtpt", jtpt_);
  if(varIsGood[5]) inTree_p->SetBranchAddress("rawpt", rawpt_);
  if(varIsGood[6]) inTree_p->SetBranchAddress("refpt", refpt_);
  if(varIsGood[7]) inTree_p->SetBranchAddress("jteta", jteta_);

  return;
}


#endif
