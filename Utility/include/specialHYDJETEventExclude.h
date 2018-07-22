#ifndef SPECIALHYDJETEVENTEXCLUDE_H
#define SPECIALHYDJETEVENTEXCLUDE_H

//For excluding events that produced anamolous high pt hard scatterings from the underlying event, to simplify problems in evaluating performance
//Events are id'd by specific generator level dijet pairs, each with subid non-zero
//since multiple events are embedded with the same hydjet, it is better just to exclude these extreme anomolies altogether 
//Alternative handling would be looking for a better match in cases where high pt jets look fake
//but this comes with its own problems

class specialHYDJETEventExclude
{
 public:
  specialHYDJETEventExclude(){};

  //Following are defined with respect to ak4GenJets
  //First exclusion 600 gev dijet
  //Second exclusion 276 gev single jet, no compensating neutrino, no compensating large eta jet, complete mystery
  static const int nEvtToExclude = 2;
  Float_t leadingGenJtPt[nEvtToExclude] = {650.40368, 276.91653};
  Float_t subleadingGenJtPt[nEvtToExclude] = {616.23840, 30.184194};
  Float_t leadingGenJtPhi[nEvtToExclude] = {-0.294406, 1.9887496};
  Float_t subleadingGenJtPhi[nEvtToExclude] = {2.8826510, 0.3335092};
  Float_t leadingGenJtEta[nEvtToExclude] = {-0.222331, -0.214549};
  Float_t subleadingGenJtEta[nEvtToExclude] = {-0.346403, 0.0162385};

  double deltaVal = 0.01;

  bool CheckEventBadJet(const int ngen_, float genpt_[], float genphi_[], float geneta_[], int gensubid_[]);
};

bool specialHYDJETEventExclude::CheckEventBadJet(const int ngen_, float genpt_[], float genphi_[], float geneta_[], int gensubid_[])
{
  bool isLeadFound = false;
  bool isSubleadFound = false;

  for(int j = 0; j < nEvtToExclude; ++j){
    for(int i = 0; i < ngen_; ++i){
      if(gensubid_[i] == 0) continue;

      bool leadDeltaPt = genpt_[i] > leadingGenJtPt[j] - deltaVal && genpt_[i] < leadingGenJtPt[j] + deltaVal;
      bool subleadDeltaPt = genpt_[i] > subleadingGenJtPt[j] - deltaVal && genpt_[i] < subleadingGenJtPt[j] + deltaVal;
      
      bool leadDeltaPhi = genphi_[i] > leadingGenJtPhi[j] - deltaVal && genphi_[i] < leadingGenJtPhi[j] + deltaVal;
      bool subleadDeltaPhi = genphi_[i] > subleadingGenJtPhi[j] - deltaVal && genphi_[i] < subleadingGenJtPhi[j] + deltaVal;

      bool leadDeltaEta = geneta_[i] > leadingGenJtEta[j] - deltaVal && geneta_[i] < leadingGenJtEta[j] + deltaVal;
      bool subleadDeltaEta = geneta_[i] > subleadingGenJtEta[j] - deltaVal && geneta_[i] < subleadingGenJtEta[j] + deltaVal;

      if(leadDeltaPt && leadDeltaPhi && leadDeltaEta) isLeadFound = true;
      if(subleadDeltaPt && subleadDeltaPhi && subleadDeltaEta) isSubleadFound = true;
    }

    if(isLeadFound && isSubleadFound) break;
    else{
      isLeadFound = false;
      isSubleadFound = false;
    }
  }

  return isLeadFound && isSubleadFound;
}

#endif
