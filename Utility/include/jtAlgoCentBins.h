#ifndef JTALGOCENTBINS_H
#define JTALGOCENTBINS_H

#include <vector>

class jtAlgoCentBins{
 public:
  static const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {50, 30, 10, 0};
  Int_t centBinsHi[nCentBins] = {90, 50, 30, 10};

  static const Int_t nAbsEtaBins = 1;
  Float_t absEtaBinsLow[nAbsEtaBins] = {0.0};
  Float_t absEtaBinsHi[nAbsEtaBins] = {2.0};

  static const Int_t nTruthBins[];
  static const Float_t truthBins5090[];
  static const Float_t truthBins3050[];
  static const Float_t truthBins1030[];
  static const Float_t truthBins010[];

  static const Int_t nRecoBins[];
  static const Float_t recoBins5090[];
  static const Float_t recoBins3050[];
  static const Float_t recoBins1030[];
  static const Float_t recoBins010[];

  static const Int_t nTruthBinsReduced[];
  static const Float_t truthBinsReduced5090[];
  static const Float_t truthBinsReduced3050[];
  static const Float_t truthBinsReduced1030[];
  static const Float_t truthBinsReduced010[];

  static const Int_t nRecoBinsReduced[];
  static const Float_t recoBinsReduced5090[];
  static const Float_t recoBinsReduced3050[];
  static const Float_t recoBinsReduced1030[];
  static const Float_t recoBinsReduced010[];

  static const Int_t rValBinSplit = 6;

  jtAlgoCentBins(){};
  
  Int_t getNBins(const Int_t centBinLow, const Int_t centBinHi, const Int_t rVal, const Bool_t isReco);
  void getBins(const Int_t centBinLow, const Int_t centBinHi, const Int_t rVal, const Bool_t isReco, const Int_t nBins, Float_t bins[]);

 protected:
  Int_t getCentPos(const Int_t centBinLow, const Int_t centBinHi);
};

const Int_t jtAlgoCentBins::nTruthBins[] = {12, 12, 12, 12};
const Float_t jtAlgoCentBins::truthBins5090[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700};
const Float_t jtAlgoCentBins::truthBins3050[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700};
const Float_t jtAlgoCentBins::truthBins1030[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700};
const Float_t jtAlgoCentBins::truthBins010[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700};

const Int_t jtAlgoCentBins::nRecoBins[] = {10, 10, 10, 10};
const Float_t jtAlgoCentBins::recoBins5090[] = {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650};
const Float_t jtAlgoCentBins::recoBins3050[] = {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650};
const Float_t jtAlgoCentBins::recoBins1030[] = {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650};
const Float_t jtAlgoCentBins::recoBins010[] = {150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650};

const Int_t jtAlgoCentBins::nTruthBinsReduced[] = {7, 7, 7, 7};
const Float_t jtAlgoCentBins::truthBinsReduced5090[] = {100, 200, 300, 400, 500, 600, 700, 800};
const Float_t jtAlgoCentBins::truthBinsReduced3050[] = {100, 200, 300, 400, 500, 600, 700, 800};
const Float_t jtAlgoCentBins::truthBinsReduced1030[] = {100, 200, 300, 400, 500, 600, 700, 800};
const Float_t jtAlgoCentBins::truthBinsReduced010[] = {100, 200, 300, 400, 500, 600, 700, 800};

const Int_t jtAlgoCentBins::nRecoBinsReduced[] = {5, 5, 5, 5};
const Float_t jtAlgoCentBins::recoBinsReduced5090[] = {200, 300, 400, 500, 600, 700};
const Float_t jtAlgoCentBins::recoBinsReduced3050[] = {200, 300, 400, 500, 600, 700};
const Float_t jtAlgoCentBins::recoBinsReduced1030[] = {200, 300, 400, 500, 600, 700};
const Float_t jtAlgoCentBins::recoBinsReduced010[] = {200, 300, 400, 500, 600, 700};


Int_t jtAlgoCentBins::getNBins(const Int_t centBinLow, const Int_t centBinHi, const Int_t rVal, const Bool_t isReco)
{
  Int_t nBins = -1;
  Int_t centPos = getCentPos(centBinLow, centBinHi);
  if(centPos == -1) return nBins;
  
  if(isReco){
    if(rVal <=  rValBinSplit) nBins = nRecoBins[centPos];
    else nBins = nRecoBinsReduced[centPos];
  }
  else{
    if(rVal <=  rValBinSplit) nBins = nTruthBins[centPos];
    else nBins = nTruthBinsReduced[centPos];
  }
  
  return nBins;
}

void jtAlgoCentBins::getBins(const Int_t centBinLow, const Int_t centBinHi, const Int_t rVal, const Bool_t isReco, const Int_t nBins, Float_t bins[])
{
  Int_t centPos = getCentPos(centBinLow, centBinHi);
  if(centPos == -1) return;

  if(isReco){
    if(rVal <=  rValBinSplit){
      for(Int_t bI = 0; bI < nBins+1; ++bI){
	if(centPos == 0) bins[bI] = recoBins5090[bI];
	else if(centPos == 1) bins[bI] = recoBins3050[bI];
	else if(centPos == 2) bins[bI] = recoBins1030[bI];
	else if(centPos == 3) bins[bI] = recoBins010[bI];
      }
    }
    else{
      for(Int_t bI = 0; bI < nBins+1; ++bI){
	if(centPos == 0) bins[bI] = recoBinsReduced5090[bI];
	else if(centPos == 1) bins[bI] = recoBinsReduced3050[bI];
	else if(centPos == 2) bins[bI] = recoBinsReduced1030[bI];
	else if(centPos == 3) bins[bI] = recoBinsReduced010[bI];
      }
    }
  }
  else{
    if(rVal <=  rValBinSplit){
      for(Int_t bI = 0; bI < nBins+1; ++bI){
	if(centPos == 0) bins[bI] = truthBins5090[bI];
	else if(centPos == 1) bins[bI] = truthBins3050[bI];
	else if(centPos == 2) bins[bI] = truthBins1030[bI];
	else if(centPos == 3) bins[bI] = truthBins010[bI];
      }
    }
    else{
      for(Int_t bI = 0; bI < nBins+1; ++bI){
        if(centPos == 0) bins[bI] = truthBinsReduced5090[bI];
        else if(centPos == 1) bins[bI] = truthBinsReduced3050[bI];
        else if(centPos == 2) bins[bI] = truthBinsReduced1030[bI];
        else if(centPos == 3) bins[bI] = truthBinsReduced010[bI];
      }
    }
  }

  return;
}


Int_t jtAlgoCentBins::getCentPos(const Int_t centBinLow, const Int_t centBinHi)
{
  Int_t centPos = -1;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    if(centBinsLow[cI] == centBinLow && centBinHi == centBinsHi[cI]){
      centPos = cI;
      break;
    }
  }

  if(centPos == -1){
    std::cout << "WARNING: given centLow, centHi \'" << centBinLow << ", " << centBinHi << "\' is not found in available bins: " << std::endl;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::cout << " \'" << centBinsLow[cI] << "-" << centBinsHi[cI] << "\'" << std::endl;
    }
    std::cout << "Returning -1" << std::endl;
  }
  
  return centPos;
}

#endif
