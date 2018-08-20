#ifndef LUMIANDTAAUTIL_H
#define LUMIANDTAAUTIL_H

#include "TH1.h"
#include "TMath.h"

double getLumiFactor()
{
  return 27400.;
}

double getLumiPercentError()
{
  return .023;
}

double getLumiAbsError()
{
  return getLumiFactor()*getLumiPercentError();
}

double getNMBEvents()
{
  const Double_t nMBEvents_JSON_V2_DEFAULT = 55.1803882404*(49.454572);
  Double_t nMBEvents_JSON_V2_TEMP = nMBEvents_JSON_V2_DEFAULT;

  const Double_t nMBEvents_JSON_V2 = nMBEvents_JSON_V2_TEMP;
  return nMBEvents_JSON_V2;
}

double getTAAScaleFactor(const std::string inStr)
{
  double scaleFactor = -1.;
  if(inStr.find("Cent0to10") != std::string::npos) scaleFactor = 23.22;
  else if(inStr.find("Cent10to30") != std::string::npos) scaleFactor = 11.51;
  else if(inStr.find("Cent30to50") != std::string::npos) scaleFactor = 3.819;
  else if(inStr.find("Cent50to90") != std::string::npos) scaleFactor = 0.543 ;
  else{
    std::cout << "WARNING: \'" << inStr << "\' is not found, return -1" << std::endl;
  }
  return scaleFactor;
}

double getTAAScaleFactorUp(const std::string inStr)
{
  double scaleFactor = -1.;
  if(inStr.find("Cent0to10") != std::string::npos) scaleFactor = .019;
  else if(inStr.find("Cent10to30") != std::string::npos) scaleFactor = .026;
  else if(inStr.find("Cent30to50") != std::string::npos) scaleFactor = .054;
  else if(inStr.find("Cent50to90") != std::string::npos) scaleFactor = .112;
  else{
    std::cout << "WARNING: \'" << inStr << "\' is not found, return -1" << std::endl;
  }
  return scaleFactor;
}

double getTAAScaleFactorDown(const std::string inStr)
{
  double scaleFactor = -1.;
  if(inStr.find("Cent0to10") != std::string::npos) scaleFactor = .03;
  else if(inStr.find("Cent10to30") != std::string::npos) scaleFactor = .034;
  else if(inStr.find("Cent30to50") != std::string::npos) scaleFactor = .054;
  else if(inStr.find("Cent50to90") != std::string::npos) scaleFactor = .073;
  else{
    std::cout << "WARNING: \'" << inStr << "\' is not found, return -1" << std::endl;
  }
  return scaleFactor;
}

void divBinWidth(TH1* inHist_p)
{
  for(Int_t bIX = 0; bIX < inHist_p->GetNbinsX(); ++bIX){
    Double_t val = inHist_p->GetBinContent(bIX+1);
    Double_t err = inHist_p->GetBinError(bIX+1);
    val /= inHist_p->GetBinWidth(bIX+1);
    err /= inHist_p->GetBinWidth(bIX+1);
    inHist_p->SetBinContent(bIX+1, val);
    inHist_p->SetBinError(bIX+1, err);
  }

  return;
}


void scaleCentralValues(TH1* inHist_p, Double_t scaleVal)
{
  for(Int_t bIX = 0; bIX < inHist_p->GetNbinsX(); ++bIX){
    inHist_p->SetBinContent(bIX+1, inHist_p->GetBinContent(bIX+1)*scaleVal);
  }

  return;
}


void scaleErrorValues(TH1* inHist_p, Double_t scaleVal)
{
  for(Int_t bIX = 0; bIX < inHist_p->GetNbinsX(); ++bIX){
    inHist_p->SetBinError(bIX+1, inHist_p->GetBinError(bIX+1)*scaleVal);
  }

  return;
}

void scaleCentralAndErrorValues(TH1* inHist_p, Double_t scaleVal)
{
  scaleCentralValues(inHist_p, scaleVal); 
  scaleErrorValues(inHist_p, scaleVal);
  return;
}

double getMin(TH1* inHist_p)
{
  double minVal = inHist_p->GetMaximum();
  for(Int_t bIX = 0; bIX < inHist_p->GetNbinsX(); ++bIX){
    if(inHist_p->GetBinContent(bIX+1) < minVal) minVal = inHist_p->GetBinContent(bIX+1);
  }
  return minVal;
}

double getMinGTZero(TH1* inHist_p)
{
  double minVal = inHist_p->GetMaximum();
  for(Int_t bIX = 0; bIX < inHist_p->GetNbinsX(); ++bIX){
    if(inHist_p->GetBinContent(bIX+1) <= 0) continue;
    if(inHist_p->GetBinContent(bIX+1) < minVal) minVal = inHist_p->GetBinContent(bIX+1);
  }
  return minVal;
}

double getMax(TH1* inHist_p)
{
  double maxVal = inHist_p->GetMaximum();
  for(Int_t bIX = 0; bIX < inHist_p->GetNbinsX(); ++bIX){
    if(inHist_p->GetBinContent(bIX+1) > maxVal) maxVal = inHist_p->GetBinContent(bIX+1);
  }
  return maxVal;
}

int getColorPosFromAlgo(const std::string algoStr)
{
  Int_t colPos = 0;
  if(algoStr.find("akCs4P") != std::string::npos) colPos = 4;
  else if(algoStr.find("ak4PF") != std::string::npos) colPos = 4;
  else if(algoStr.find("akCs6P") != std::string::npos) colPos = 6;
  else if(algoStr.find("ak6PF") != std::string::npos) colPos = 6;
  else if(algoStr.find("akCs8P") != std::string::npos) colPos = 5;
  else if(algoStr.find("ak8PF") != std::string::npos) colPos = 5;
  else if(algoStr.find("akCs10P") != std::string::npos) colPos = 2;
  else if(algoStr.find("ak10PF") != std::string::npos) colPos = 2;

  return colPos;
}

int getColorPosFromCent(const std::string centStr, bool isPP)
{
  Int_t colPos = 0;
  if(isPP) colPos = 2;
  else if(centStr.find("Cent0to10") != std::string::npos) colPos = 0;
  else if(centStr.find("Cent10to30") != std::string::npos) colPos = 1;
  else if(centStr.find("Cent30to50") != std::string::npos) colPos = 3;
  else if(centStr.find("Cent50to90") != std::string::npos) colPos = 4;

  return colPos;
}

int getStyleFromAlgo(const std::string algoStr)
{
  Int_t colPos = 20;
  if(algoStr.find("akCs4P") != std::string::npos) colPos = 21;
  else if(algoStr.find("ak4PF") != std::string::npos) colPos = 21;
  else if(algoStr.find("akCs6P") != std::string::npos) colPos = 34;
  else if(algoStr.find("ak6PF") != std::string::npos) colPos = 34;
  else if(algoStr.find("akCs8P") != std::string::npos) colPos = 47;
  else if(algoStr.find("ak8PF") != std::string::npos) colPos = 47;
  else if(algoStr.find("akCs10P") != std::string::npos) colPos = 33;
  else if(algoStr.find("ak10PF") != std::string::npos) colPos = 33;

  return colPos;
}

int getAltStyleFromAlgo(const std::string algoStr)
{
  Int_t colPos = 24;
  if(algoStr.find("akCs4P") != std::string::npos) colPos = 25;
  else if(algoStr.find("ak4PF") != std::string::npos) colPos = 25;
  else if(algoStr.find("akCs6P") != std::string::npos) colPos = 28;
  else if(algoStr.find("ak6PF") != std::string::npos) colPos = 28;
  else if(algoStr.find("akCs8P") != std::string::npos) colPos = 46;
  else if(algoStr.find("ak8PF") != std::string::npos) colPos = 46;
  else if(algoStr.find("akCs10P") != std::string::npos) colPos = 27;
  else if(algoStr.find("ak10PF") != std::string::npos) colPos = 27;

  return colPos;
}

int getStyleFromCent(const std::string centStr, bool isPP)
{
  Int_t colPos = 0;
  if(isPP) colPos = 20;
  else if(centStr.find("Cent0to10") != std::string::npos) colPos = 34;
  else if(centStr.find("Cent10to30") != std::string::npos) colPos = 47;
  else if(centStr.find("Cent30to50") != std::string::npos) colPos = 33;
  else if(centStr.find("Cent50to90") != std::string::npos) colPos = 21;

  return colPos;
}


int getAltStyleFromCent(const std::string centStr, bool isPP)
{
  Int_t colPos = 0;
  if(isPP) colPos = 24;
  else if(centStr.find("Cent0to10") != std::string::npos) colPos = 28;
  else if(centStr.find("Cent10to30") != std::string::npos) colPos = 46;
  else if(centStr.find("Cent30to50") != std::string::npos) colPos = 27;
  else if(centStr.find("Cent50to90") != std::string::npos) colPos = 25;

  return colPos;
}

#endif
