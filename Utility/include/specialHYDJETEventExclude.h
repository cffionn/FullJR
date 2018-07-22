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
  //Second exclusion 300 gev dijet
  //Third and fourth exclusion 317 and 276 gev single jet, no compensating neutrino, no compensating large eta jet, complete mystery
  static const int nEvtToExclude = 14;
  Float_t leadingGenJtPt[nEvtToExclude] = {650.40368, 365.46627, 319.86566, 317.03808, 305.15612, 276.91653, 237.58552, 221.78674, 216.37747, 212.6747, 203.78250, 198.67369, 196.15103, 193.30279};
  Float_t subleadingGenJtPt[nEvtToExclude] = {616.23840, 20.398469, 319.5466, 18.843473, 276.77761, 30.184194, 25.196847, 23.066457, 25.259368, 31.915657, 30.73803, 26.124143, 20.469167, 19.318740};
  Float_t leadingGenJtPhi[nEvtToExclude] = {-0.294406, -3.092175, 2.3265273, 1.1733226, 1.8576631, 1.9887496, 2.9177346, -0.913193, 2.1785547, -2.992421, -1.963628, 2.9601671, -1.495971, 0.6731430};
  Float_t subleadingGenJtPhi[nEvtToExclude] = {2.8826510, -0.782366, -0.575561, 2.5958199, -1.292150, 0.3335092, 1.6835835, 0.7550130, 2.2273616, 0.2651389, 0.9244814, -2.288487, 1.0721495, 1.6285703};
  Float_t leadingGenJtEta[nEvtToExclude] = {-0.222331, -0.697133, 0.7997897, -1.317555, 0.0763501, -0.214549, -0.992048, -0.620375, 1.3119092, -1.020294, -0.953381, -0.523874, 1.4934288, -0.163160};
  Float_t subleadingGenJtEta[nEvtToExclude] = {-0.346403, -1.855122, 0.8341307, -2.253091, 1.0136429, 0.0162385, -0.787025, -0.361022, -0.073839, 2.7714095, 1.4386611, -0.764248, -1.910048, -1.704711};
  Int_t excluded[nEvtToExclude] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  double deltaVal = 0.01;

  bool CheckEventBadJet(const int ngen_, float genpt_[], float genphi_[], float geneta_[], int gensubid_[]);
  void PrintExcludedNumbers();
  int TotalExcludedNumbers();
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

    if(isLeadFound && isSubleadFound){
      excluded[j] += 1;
      break;
    }
    else{
      isLeadFound = false;
      isSubleadFound = false;
    }
  }

  return isLeadFound && isSubleadFound;
}


void specialHYDJETEventExclude::PrintExcludedNumbers()
{
  Int_t total = 0;
  std::cout << "Events excluded per jet:" << std::endl;
  for(Int_t i = 0; i < nEvtToExclude; ++i){
    std::cout << " Leading Jet p_{T} " << leadingGenJtPt[i] << ": " << excluded[i] << std::endl;
    total += excluded[i];
  }
  std::cout << "Total: " << total << std::endl;

  return;
}

int specialHYDJETEventExclude::TotalExcludedNumbers()
{
  Int_t total = 0;
  for(Int_t i = 0; i < nEvtToExclude; ++i){
    total += excluded[i];
  }
  return total;
}

#endif
