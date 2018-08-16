#include <iostream>
#include <string>

#include "TFile.h"
#include "TDatime.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TMath.h"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"

int altSmearingExample(std::string outFileName)
{
  TRandom3* randGen_p = new TRandom3(752604581);

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  if(outFileName.size() == 0) outFileName = "out";
  if(outFileName.find("output/") == std::string::npos) outFileName = "output/" + outFileName;
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
  outFileName = outFileName + "_AltSmearingExample_" + dateStr + ".root";

  const Int_t nBins = 100;
  const Double_t ptLowVal = 100;
  const Double_t ptHiVal = 1100;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* genSpect_h = new TH1F("genSpect_h", ";p_{T};Counts", nBins, ptLowVal, ptHiVal);
  TH1F* smearSpect_h = new TH1F("smearSpect_h", ";p_{T};Counts", nBins, ptLowVal, ptHiVal);
  
  const Double_t offVal = 200;
  const Double_t lowGenVal = ptLowVal + offVal;
  const Double_t hiGenVal = ptHiVal - offVal;

  const Int_t nResBins = (hiGenVal - lowGenVal)/50;
  Double_t resBins[nResBins+1];
  getLinBins(lowGenVal, hiGenVal, nResBins, resBins);

  TH1F* resHist_h[nResBins];
  for(Int_t rI = 0; rI < nResBins; ++rI){
    std::string name = "resHist_" + prettyString(resBins[rI], 1, true) + "to" + prettyString(resBins[rI+1], 1, true) + "_h";
    resHist_h[rI] = new TH1F(name.c_str(), ";Reco./Gen.;Counts", 101, 0, 2);
  }

  std::cout << "Res bins: " << nResBins << std::endl;
  for(Int_t rI = 0; rI < nResBins; ++rI){
    std::cout << " " << rI << "/" << nResBins << ": " << resBins[rI] << "-" << resBins[rI+1] << std::endl;
  }


  const Double_t cParam = 0.06;
  const Double_t sParam = 0.95;
  const Double_t nParam = 20.0;

  const Int_t nFills = 10000000;
  while(genSpect_h->GetEntries() < nFills){
    Double_t val = randGen_p->Uniform(0, 1);
    if(val == 0) genSpect_h->Fill(hiGenVal);
    else{
      val = TMath::Sqrt(val);
      val = TMath::Sqrt(val);
      val = lowGenVal/val;
      if(val > hiGenVal) continue;
      genSpect_h->Fill(val);

      Double_t sigma = TMath::Sqrt(cParam*cParam + sParam*sParam/val + nParam*nParam/(val*val));

      double newVal = val*randGen_p->Gaus(1., sigma);
      if(newVal < ptLowVal) newVal = ptLowVal;
      else if(newVal > ptHiVal) newVal = ptHiVal;

      smearSpect_h->Fill(newVal);

      Int_t resPos = -1;
      for(Int_t rI = 0; rI < nResBins; ++rI){
	if(val >= resBins[rI] && val < resBins[rI+1]){
	  resPos = rI;
	  break;
	}
      }

      resHist_h[resPos]->Fill(newVal/val);
    }
  }

  outFile_p->cd();

  genSpect_h->Write("", TObject::kOverwrite);
  delete genSpect_h;

  smearSpect_h->Write("", TObject::kOverwrite);
  delete smearSpect_h;

  for(Int_t rI = 0; rI < nResBins; ++rI){
    resHist_h[rI]->Write("", TObject::kOverwrite);
    delete resHist_h[rI];
  }

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/altSmearingExample.exe <outFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += altSmearingExample(argv[1]);
  return retVal;
}
