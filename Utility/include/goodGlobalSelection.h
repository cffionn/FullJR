#ifndef GOODGLOBALSELECTION_H
#define GOODGLOBALSELECTION_H

#include "TMath.h"

class goodGlobalSelection{
 public:
  goodGlobalSelection();

  bool isPbPb_;
  
  float vz_;
  float hiHF_;

  int pprimaryVertexFilter_;
  int pBeamScrapingFilter_;
  int pcollisionEventSelection_;
  int HBHENoiseFilterResultRun2Loose_;
  int pclusterCompatibilityFilter_;
  int phfCoincFilter3_;

  void setIsPbPb(bool inIsPbPb){isPbPb_ = inIsPbPb; return;}
  void setVz(float inVz){vz_ = inVz; return;}
  void setHiHF(float inHiHF){hiHF_ = inHiHF; return;}
  void setPprimaryVertexFilter(int inPprimaryVertexFilter_){pprimaryVertexFilter_ = inPprimaryVertexFilter_; return;}
  void setPBeamScrapingFilter(int inPBeamScrapingFilter_){pBeamScrapingFilter_ = inPBeamScrapingFilter_; return;}
  void setPcollisionEventSelection(int inPcollisionEventSelection_){pcollisionEventSelection_ = inPcollisionEventSelection_; return;}
  void setHBHENoiseFilterResultRun2Loose(int inHBHENoiseFilterResultRun2Loose_){HBHENoiseFilterResultRun2Loose_ = inHBHENoiseFilterResultRun2Loose_; return;}
  void setPclusterCompatibilityFilter(int inPclusterCompatibilityFilter_){pclusterCompatibilityFilter_ = inPclusterCompatibilityFilter_; return;}
  void setPhfCoincFilter3(int inPhfCoincFilter3_){phfCoincFilter3_ = inPhfCoincFilter3_; return;}

  bool isGood();
};

goodGlobalSelection::goodGlobalSelection()
{
  isPbPb_ = false;
  vz_ = -999;
  hiHF_ = -999;

  pprimaryVertexFilter_ = -1;
  pBeamScrapingFilter_ = -1;
  pcollisionEventSelection_ = -1;
  HBHENoiseFilterResultRun2Loose_ = -1;
  pclusterCompatibilityFilter_ = -1;
  phfCoincFilter3_ = -1;

  return;
}

bool goodGlobalSelection::isGood()
{
  if(TMath::Abs(vz_) > 15.) return false;
  else if(hiHF_ > 5500. && isPbPb_) return false;
  else if(HBHENoiseFilterResultRun2Loose_ == 0) return false;
  else if(pprimaryVertexFilter_ == 0) return false;
  else if(!isPbPb_ && pBeamScrapingFilter_ == 0) return false;
  else if(isPbPb_ && phfCoincFilter3_ == 0) return false;
  else if(isPbPb_ && pclusterCompatibilityFilter_ == 0) return false;

  return true;
}

#endif
