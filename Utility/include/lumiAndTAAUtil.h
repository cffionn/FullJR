#ifndef LUMIANDTAAUTIL_H
#define LUMIANDTAAUTIL_H

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


#endif
