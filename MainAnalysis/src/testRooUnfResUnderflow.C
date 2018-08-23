//cpp dependencies
#include <iostream>

//ROOT dependencies
#include "TH1D.h"

//RooUnfold dependencies
#include "src/RooUnfoldResponse.h"

//Non-local (Utility) dependencies
#include "Utility/include/getLinBins.h"

int testRooUnfResUnderflow()
{
  const Double_t lowVal = 100;
  const Double_t hiVal = 2000;
  const Int_t nBins = 8;
  Double_t bins[nBins+1] = {lowVal, 150., 200., 250., 300., 400., 600., 1000., hiVal};

  TH1D* reco_p = new TH1D("reco_h", "", nBins, bins);
  TH1D* truth_p = new TH1D("truth_h", "", nBins, bins);

  RooUnfoldResponse* rooRes_p = new RooUnfoldResponse("rooRes_p","");
  rooRes_p->Setup(reco_p, truth_p);

  //Reco first, gen second
  rooRes_p->Fill(250, 300);
  rooRes_p->Fill(200, 175);
  rooRes_p->Fill(600, 555);
  rooRes_p->Fill(155, 125);
  rooRes_p->Fill(155, 175);
  
  std::cout << "Print 1" << std::endl;
  rooRes_p->Print("ALL");

  rooRes_p->Miss(125);
  std::cout << "Print 2" << std::endl;
  rooRes_p->Print("ALL");

  rooRes_p->Fill(125, 1500);
  rooRes_p->Fill(1500, 125);
  std::cout << "Print 3" << std::endl;
  rooRes_p->Print("ALL");

  TH1D* initTrue_p = (TH1D*)rooRes_p->Htruth()->Clone("initTrue_p");
  std::cout << "Print truth" << std::endl;
  initTrue_p->Print("ALL");

  delete initTrue_p;

  delete rooRes_p;

  delete truth_p;
  delete reco_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += testRooUnfResUnderflow();
  return retVal;
}
