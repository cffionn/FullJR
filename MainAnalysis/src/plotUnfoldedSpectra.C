//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TDatime.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"

int plotUnfoldedSpectra(const std::string inFileNamePP, const std::string inFileNamePbPb)
{
  if(!checkFile(inFileNamePP) || !checkFile(inFileNamePbPb)){
    if(!checkFile(inFileNamePP)) std::cout << "inFileNamePP, \'" << inFileNamePP << "\', is invalid. return 1" << std::endl;
    if(!checkFile(inFileNamePbPb)) std::cout << "inFileNamePbPb, \'" << inFileNamePbPb << "\', is invalid. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  const std::string dirStr = "pdfDir/" + dateStr;
  checkMakeDir(dirStr);

  TFile* inFilePP_p = new TFile(inFileNamePP.c_str(), "READ");
  cutPropagator cutPropPP;
  cutPropPP.Clean();
  cutPropPP.GetAllVarFromFile(inFilePP_p);

  TFile* inFilePbPb_p = new TFile(inFileNamePbPb.c_str(), "READ");
  cutPropagator cutPropPbPb;
  cutPropPbPb.Clean();
  cutPropPbPb.GetAllVarFromFile(inFilePbPb_p);

  if(!cutPropPP.GetIsPP() || cutPropPbPb.GetIsPP() || !cutPropPP.CheckPropagatorsMatch(cutPropPbPb, true, false)){
    if(!cutPropPP.GetIsPP()) std::cout << "inFileNamePP \'" << inFileNamePP << "\' is not pp, return 1" << std::endl;
    if(cutPropPbPb.GetIsPP()) std::cout << "inFileNamePbPb \'" << inFileNamePbPb << "\' is not PbPb, return 1" << std::endl;
    if(!cutPropPP.CheckPropagatorsMatch(cutPropPbPb, true, false)) std::cout << "Cut propagators of pp and pbpb do not match. inspect. return 1" << std::endl;
  
    inFilePP_p->Close();
    delete inFilePP_p;
    
    inFilePbPb_p->Close();
    delete inFilePbPb_p;
    
    return 1;
  }


  const Int_t nSyst = cutPropPbPb.GetNSyst();
  std::vector<std::string> systStr = cutPropPbPb.GetSystStr();

  const Int_t nResponseMod = cutPropPbPb.GetNResponseMod();
  std::vector<double> responseMod = cutPropPbPb.GetResponseMod();

  const Int_t nCentBins = cutPropPbPb.GetNCentBins();
  std::vector<Int_t> centBinsLow = cutPropPbPb.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutPropPbPb.GetCentBinsHi();

  const Int_t nJtPtBins = cutPropPbPb.GetNJtPtBins();
  std::vector<Double_t> jtPtBinsTemp = cutPropPbPb.GetJtPtBins();

  const Int_t nJtAbsEtaBins = cutPropPbPb.GetNJtAbsEtaBins();
  std::vector<Double_t> jtAbsEtaBinsLowTemp = cutPropPbPb.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHiTemp = cutPropPbPb.GetJtAbsEtaBinsHi();

  const Int_t nID = cutPropPbPb.GetNID();
  std::vector<std::string> idStr = cutPropPbPb.GetIdStr();
  std::vector<double> jtPfCHMFCutLow = cutPropPbPb.GetJtPfCHMFCutLow();
  std::vector<double> jtPfCHMFCutHi = cutPropPbPb.GetJtPfCHMFCutHi();
  std::vector<double> jtPfMUMFCutLow = cutPropPbPb.GetJtPfMUMFCutLow();
  std::vector<double> jtPfMUMFCutHi = cutPropPbPb.GetJtPfMUMFCutHi();


  



  inFilePP_p->Close();
  delete inFilePP_p;

  inFilePbPb_p->Close();
  delete inFilePbPb_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/plotUnfoldedSpectra.exe <inFileNamePP> <inFileNamePbPb>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotUnfoldedSpectra(argv[1], argv[2]);
  return retVal;
}
