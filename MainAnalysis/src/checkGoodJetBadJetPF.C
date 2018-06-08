#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDatime.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/ncollFunctions_5TeV.h"
#include "Utility/include/specialHYDJETEventExclude.h"
#include "Utility/include/vanGoghPalette.h"


int checkGoodJetBadJetPF(const std::string inName, bool isPP = false)
{
  std::vector<std::string> fileList;
  std::vector<double> pthats;
  std::vector<double> pthatWeights;

  if(inName.find(".root") != std::string::npos){
    fileList.push_back(inName);
    pthats.push_back(1.);
    pthatWeights.push_back(1.);
  }
  else if(inName.find(".txt") != std::string::npos){
    std::ifstream file(inName.c_str());
    std::string tempStr;

    while(std::getline(file, tempStr)){
      while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
      else{
        if(tempStr.substr(0, std::string("PTHAT=").size()).find("PTHAT=") != std::string::npos){
          tempStr.replace(0, std::string("PTHAT=").size(), "");
          while(tempStr.find(",") != std::string::npos){
            pthats.push_back(std::stod(tempStr.substr(0, tempStr.find(","))));
            tempStr.replace(0, tempStr.find(",")+1, "");
          }
          if(tempStr.size() != 0) pthats.push_back(std::stod(tempStr));
        }
        else if(tempStr.substr(0, std::string("PTHATWEIGHTS=").size()).find("PTHATWEIGHTS=") != std::string::npos){
          tempStr.replace(0, std::string("PTHATWEIGHTS=").size(), "");
          while(tempStr.find(",") != std::string::npos){
            pthatWeights.push_back(std::stod(tempStr.substr(0, tempStr.find(","))));
            tempStr.replace(0, tempStr.find(",")+1, "");
          }
          if(tempStr.size() != 0) pthatWeights.push_back(std::stod(tempStr));
        }
        else if(tempStr.substr(0, std::string("ISPP=").size()).find("ISPP=") != std::string::npos){
          tempStr.replace(0,tempStr.find("=")+1, "");
          isPP = std::stoi(tempStr);
        }
        else std::cout << "WARNING: Line in \'" << inName << "\', \'" << tempStr << "\' is invalid. check input" << std::endl;
      }
    }

    file.close();
  }
  else{
    std::cout << "Given inName \'" << inName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  if(fileList.size() == 0){
    std::cout << "Given inName \'" << inName << "\' is gives no valid root files. return 1" << std::endl;
    return 1;
  }
  else if(pthats.size() == 0){
    std::cout << "Given inName \'" << inName << "\' contains no pthat list. return 1" << std::endl;
    return 1;
  }
  else if(pthatWeights.size() == 0){
    std::cout << "Given inName \'" << inName << "\' contains no pthatWeights list. return 1" << std::endl;
    return 1;
  }

  std::cout << "Pthats and weights: " << std::endl;
  for(unsigned int pI = 0; pI < pthats.size(); ++pI){
    std::cout << " " << pI << "/" << pthats.size() << ": " << pthats.at(pI) << ", " << pthatWeights.at(pI) << std::endl;
  }

  std::cout << "IsPP: " << isPP << std::endl;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  std::string outFileName = inName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
  else if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  outFileName = "output/" + outFileName + "_CheckGoodJetBadJetPF_" + dateStr + ".root";
  checkMakeDir("output");

  const Double_t jtRecoPtMin = 500.;
  const Double_t jtRefGoodPtMin = 350.;
  const Double_t jtRefBadPtMax = 200.;
  const Double_t jtAbsEtaMax = 2.;
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  TH1D* jtPt_Good_Unweighted_h = new TH1D("jtPt_Good_Unweighted_h", ";Jet p_{T} (Good);Counts", 10, 500, 1000);
  TH1D* jtPt_Bad_Unweighted_h = new TH1D("jtPt_Bad_Unweighted_h", ";Jet p_{T} (Bad);Counts", 10, 500, 1000);

  TH1D* jtPt_Good_Unweighted_NoSpecialCut_h = new TH1D("jtPt_Good_Unweighted_NoSpecialCut_h", ";Jet p_{T} (Good);Counts", 10, 500, 1000);
  TH1D* jtPt_Bad_Unweighted_NoSpecialCut_h = new TH1D("jtPt_Bad_Unweighted_NoSpecialCut_h", ";Jet p_{T} (Bad);Counts", 10, 500, 1000);

  TH1D* jtPt_Good_Unweighted_Cut_h = new TH1D("jtPt_Good_Unweighted_Cut_h", ";Jet p_{T} (Good);Counts", 10, 500, 1000);
  TH1D* jtPt_Bad_Unweighted_Cut_h = new TH1D("jtPt_Bad_Unweighted_Cut_h", ";Jet p_{T} (Bad);Counts", 10, 500, 1000);

  TH1D* jtPt_Good_Weighted_h = new TH1D("jtPt_Good_Weighted_h", ";Jet p_{T} (Good);Counts (Counts)", 10, 500, 1000);
  TH1D* jtPt_Bad_Weighted_h = new TH1D("jtPt_Bad_Weighted_h", ";Jet p_{T} (Bad);Counts (Counts)", 10, 500, 1000);

  TH1D* jtPt_Good_Weighted_NoSpecialCut_h = new TH1D("jtPt_Good_Weighted_NoSpecialCut_h", ";Jet p_{T} (Good);Counts (Counts)", 10, 500, 1000);
  TH1D* jtPt_Bad_Weighted_NoSpecialCut_h = new TH1D("jtPt_Bad_Weighted_NoSpecialCut_h", ";Jet p_{T} (Bad);Counts (Counts)", 10, 500, 1000);

  TH1D* jtPt_Good_Weighted_Cut_h = new TH1D("jtPt_Good_Weighted_Cut_h", ";Jet p_{T} (Good);Counts (Counts)", 10, 500, 1000);
  TH1D* jtPt_Bad_Weighted_Cut_h = new TH1D("jtPt_Bad_Weighted_Cut_h", ";Jet p_{T} (Bad);Counts (Counts)", 10, 500, 1000);

  TH1D* jtPfCHF_Good_Unweighted_h = new TH1D("jtPfCHF_Good_Unweighted_h", ";jtPfCHF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfCEF_Good_Unweighted_h = new TH1D("jtPfCEF_Good_Unweighted_h", ";jtPfCEF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfNHF_Good_Unweighted_h = new TH1D("jtPfNHF_Good_Unweighted_h", ";jtPfNHF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfNEF_Good_Unweighted_h = new TH1D("jtPfNEF_Good_Unweighted_h", ";jtPfNEF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfMUF_Good_Unweighted_h = new TH1D("jtPfMUF_Good_Unweighted_h", ";jtPfMUF;Counts", 10, 0.0, 1.0);

  TH2D* jtPfCHF_jtPfCEF_Good_Unweighted_h = new TH2D("jtPfCHF_jtPfCEF_Good_Unweighted_h", ";jtPfCHF;jtPfCEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHF_jtPfNHF_Good_Unweighted_h = new TH2D("jtPfCHF_jtPfNHF_Good_Unweighted_h", ";jtPfCHF;jtPfNHF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHF_jtPfNEF_Good_Unweighted_h = new TH2D("jtPfCHF_jtPfNEF_Good_Unweighted_h", ";jtPfCHF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHF_jtPfMUF_Good_Unweighted_h = new TH2D("jtPfCHF_jtPfMUF_Good_Unweighted_h", ";jtPfCHF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfCEF_jtPfNHF_Good_Unweighted_h = new TH2D("jtPfCEF_jtPfNHF_Good_Unweighted_h", ";jtPfCEF;jtPfNHF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCEF_jtPfNEF_Good_Unweighted_h = new TH2D("jtPfCEF_jtPfNEF_Good_Unweighted_h", ";jtPfCEF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCEF_jtPfMUF_Good_Unweighted_h = new TH2D("jtPfCEF_jtPfMUF_Good_Unweighted_h", ";jtPfCEF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfNHF_jtPfNEF_Good_Unweighted_h = new TH2D("jtPfNHF_jtPfNEF_Good_Unweighted_h", ";jtPfNHF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfNHF_jtPfMUF_Good_Unweighted_h = new TH2D("jtPfNHF_jtPfMUF_Good_Unweighted_h", ";jtPfNHF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfNEF_jtPfMUF_Good_Unweighted_h = new TH2D("jtPfNEF_jtPfMUF_Good_Unweighted_h", ";jtPfNEF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH1D* jtPfCHMF_Good_Unweighted_h = new TH1D("jtPfCHMF_Good_Unweighted_h", ";jtPfCHMF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfCEMF_Good_Unweighted_h = new TH1D("jtPfCEMF_Good_Unweighted_h", ";jtPfCEMF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfNHMF_Good_Unweighted_h = new TH1D("jtPfNHMF_Good_Unweighted_h", ";jtPfNHMF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfNEMF_Good_Unweighted_h = new TH1D("jtPfNEMF_Good_Unweighted_h", ";jtPfNEMF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfMUMF_Good_Unweighted_h = new TH1D("jtPfMUMF_Good_Unweighted_h", ";jtPfMUMF;Counts", 10, 0.0, 1.0);

  TH2D* jtPfCHMF_jtPfCEMF_Good_Unweighted_h = new TH2D("jtPfCHMF_jtPfCEMF_Good_Unweighted_h", ";jtPfCHMF;jtPfCEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHMF_jtPfNHMF_Good_Unweighted_h = new TH2D("jtPfCHMF_jtPfNHMF_Good_Unweighted_h", ";jtPfCHMF;jtPfNHMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHMF_jtPfNEMF_Good_Unweighted_h = new TH2D("jtPfCHMF_jtPfNEMF_Good_Unweighted_h", ";jtPfCHMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHMF_jtPfMUMF_Good_Unweighted_h = new TH2D("jtPfCHMF_jtPfMUMF_Good_Unweighted_h", ";jtPfCHMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfCEMF_jtPfNHMF_Good_Unweighted_h = new TH2D("jtPfCEMF_jtPfNHMF_Good_Unweighted_h", ";jtPfCEMF;jtPfNHMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCEMF_jtPfNEMF_Good_Unweighted_h = new TH2D("jtPfCEMF_jtPfNEMF_Good_Unweighted_h", ";jtPfCEMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCEMF_jtPfMUMF_Good_Unweighted_h = new TH2D("jtPfCEMF_jtPfMUMF_Good_Unweighted_h", ";jtPfCEMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfNHMF_jtPfNEMF_Good_Unweighted_h = new TH2D("jtPfNHMF_jtPfNEMF_Good_Unweighted_h", ";jtPfNHMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfNHMF_jtPfMUMF_Good_Unweighted_h = new TH2D("jtPfNHMF_jtPfMUMF_Good_Unweighted_h", ";jtPfNHMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfNEMF_jtPfMUMF_Good_Unweighted_h = new TH2D("jtPfNEMF_jtPfMUMF_Good_Unweighted_h", ";jtPfNEMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);


  TH1D* jtPfCHM_Good_Unweighted_h = new TH1D("jtPfCHM_Good_Unweighted_h", ";jtPfCHM;Counts", 20, 0.0, 40.0);
  TH1D* jtPfCEM_Good_Unweighted_h = new TH1D("jtPfCEM_Good_Unweighted_h", ";jtPfCEM;Counts", 20, 0.0, 40.0);
  TH1D* jtPfNHM_Good_Unweighted_h = new TH1D("jtPfNHM_Good_Unweighted_h", ";jtPfNHM;Counts", 20, 0.0, 40.0);
  TH1D* jtPfNEM_Good_Unweighted_h = new TH1D("jtPfNEM_Good_Unweighted_h", ";jtPfNEM;Counts", 20, 0.0, 40.0);
  TH1D* jtPfMUM_Good_Unweighted_h = new TH1D("jtPfMUM_Good_Unweighted_h", ";jtPfMUM;Counts", 20, 0.0, 40.0);

  TH1D* jtPfCHF_Bad_Unweighted_h = new TH1D("jtPfCHF_Bad_Unweighted_h", ";jtPfCHF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfCEF_Bad_Unweighted_h = new TH1D("jtPfCEF_Bad_Unweighted_h", ";jtPfCEF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfNHF_Bad_Unweighted_h = new TH1D("jtPfNHF_Bad_Unweighted_h", ";jtPfNHF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfNEF_Bad_Unweighted_h = new TH1D("jtPfNEF_Bad_Unweighted_h", ";jtPfNEF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfMUF_Bad_Unweighted_h = new TH1D("jtPfMUF_Bad_Unweighted_h", ";jtPfMUF;Counts", 10, 0.0, 1.0);

  TH2D* jtPfCHF_jtPfCEF_Bad_Unweighted_h = new TH2D("jtPfCHF_jtPfCEF_Bad_Unweighted_h", ";jtPfCHF;jtPfCEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHF_jtPfNHF_Bad_Unweighted_h = new TH2D("jtPfCHF_jtPfNHF_Bad_Unweighted_h", ";jtPfCHF;jtPfNHF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHF_jtPfNEF_Bad_Unweighted_h = new TH2D("jtPfCHF_jtPfNEF_Bad_Unweighted_h", ";jtPfCHF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHF_jtPfMUF_Bad_Unweighted_h = new TH2D("jtPfCHF_jtPfMUF_Bad_Unweighted_h", ";jtPfCHF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfCEF_jtPfNHF_Bad_Unweighted_h = new TH2D("jtPfCEF_jtPfNHF_Bad_Unweighted_h", ";jtPfCEF;jtPfNHF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCEF_jtPfNEF_Bad_Unweighted_h = new TH2D("jtPfCEF_jtPfNEF_Bad_Unweighted_h", ";jtPfCEF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCEF_jtPfMUF_Bad_Unweighted_h = new TH2D("jtPfCEF_jtPfMUF_Bad_Unweighted_h", ";jtPfCEF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfNHF_jtPfNEF_Bad_Unweighted_h = new TH2D("jtPfNHF_jtPfNEF_Bad_Unweighted_h", ";jtPfNHF;jtPfNEF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfNHF_jtPfMUF_Bad_Unweighted_h = new TH2D("jtPfNHF_jtPfMUF_Bad_Unweighted_h", ";jtPfNHF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfNEF_jtPfMUF_Bad_Unweighted_h = new TH2D("jtPfNEF_jtPfMUF_Bad_Unweighted_h", ";jtPfNEF;jtPfMUF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH1D* jtPfCHMF_Bad_Unweighted_h = new TH1D("jtPfCHMF_Bad_Unweighted_h", ";jtPfCHMF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfCEMF_Bad_Unweighted_h = new TH1D("jtPfCEMF_Bad_Unweighted_h", ";jtPfCEMF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfNHMF_Bad_Unweighted_h = new TH1D("jtPfNHMF_Bad_Unweighted_h", ";jtPfNHMF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfNEMF_Bad_Unweighted_h = new TH1D("jtPfNEMF_Bad_Unweighted_h", ";jtPfNEMF;Counts", 10, 0.0, 1.0);
  TH1D* jtPfMUMF_Bad_Unweighted_h = new TH1D("jtPfMUMF_Bad_Unweighted_h", ";jtPfMUMF;Counts", 10, 0.0, 1.0);

  TH2D* jtPfCHMF_jtPfCEMF_Bad_Unweighted_h = new TH2D("jtPfCHMF_jtPfCEMF_Bad_Unweighted_h", ";jtPfCHMF;jtPfCEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHMF_jtPfNHMF_Bad_Unweighted_h = new TH2D("jtPfCHMF_jtPfNHMF_Bad_Unweighted_h", ";jtPfCHMF;jtPfNHMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHMF_jtPfNEMF_Bad_Unweighted_h = new TH2D("jtPfCHMF_jtPfNEMF_Bad_Unweighted_h", ";jtPfCHMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCHMF_jtPfMUMF_Bad_Unweighted_h = new TH2D("jtPfCHMF_jtPfMUMF_Bad_Unweighted_h", ";jtPfCHMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfCEMF_jtPfNHMF_Bad_Unweighted_h = new TH2D("jtPfCEMF_jtPfNHMF_Bad_Unweighted_h", ";jtPfCEMF;jtPfNHMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCEMF_jtPfNEMF_Bad_Unweighted_h = new TH2D("jtPfCEMF_jtPfNEMF_Bad_Unweighted_h", ";jtPfCEMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfCEMF_jtPfMUMF_Bad_Unweighted_h = new TH2D("jtPfCEMF_jtPfMUMF_Bad_Unweighted_h", ";jtPfCEMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfNHMF_jtPfNEMF_Bad_Unweighted_h = new TH2D("jtPfNHMF_jtPfNEMF_Bad_Unweighted_h", ";jtPfNHMF;jtPfNEMF", 10, 0.0, 1.0, 10, 0.0, 1.0);
  TH2D* jtPfNHMF_jtPfMUMF_Bad_Unweighted_h = new TH2D("jtPfNHMF_jtPfMUMF_Bad_Unweighted_h", ";jtPfNHMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);

  TH2D* jtPfNEMF_jtPfMUMF_Bad_Unweighted_h = new TH2D("jtPfNEMF_jtPfMUMF_Bad_Unweighted_h", ";jtPfNEMF;jtPfMUMF", 10, 0.0, 1.0, 10, 0.0, 1.0);


  TH1D* jtPfCHM_Bad_Unweighted_h = new TH1D("jtPfCHM_Bad_Unweighted_h", ";jtPfCHM;Counts", 20, 0.0, 40.0);
  TH1D* jtPfCEM_Bad_Unweighted_h = new TH1D("jtPfCEM_Bad_Unweighted_h", ";jtPfCEM;Counts", 20, 0.0, 40.0);
  TH1D* jtPfNHM_Bad_Unweighted_h = new TH1D("jtPfNHM_Bad_Unweighted_h", ";jtPfNHM;Counts", 20, 0.0, 40.0);
  TH1D* jtPfNEM_Bad_Unweighted_h = new TH1D("jtPfNEM_Bad_Unweighted_h", ";jtPfNEM;Counts", 20, 0.0, 40.0);
  TH1D* jtPfMUM_Bad_Unweighted_h = new TH1D("jtPfMUM_Bad_Unweighted_h", ";jtPfMUM;Counts", 20, 0.0, 40.0);


  centerTitles({jtPt_Good_Unweighted_h, jtPt_Good_Unweighted_Cut_h, jtPt_Good_Unweighted_NoSpecialCut_h, jtPt_Good_Weighted_h, jtPt_Good_Weighted_Cut_h, jtPt_Good_Weighted_NoSpecialCut_h, jtPfCHF_Good_Unweighted_h, jtPfCEF_Good_Unweighted_h, jtPfNHF_Good_Unweighted_h, jtPfNEF_Good_Unweighted_h, jtPfMUF_Good_Unweighted_h, jtPt_Bad_Unweighted_h, jtPt_Bad_Unweighted_Cut_h, jtPt_Bad_Unweighted_NoSpecialCut_h, jtPt_Bad_Weighted_h, jtPt_Bad_Weighted_Cut_h, jtPt_Bad_Weighted_NoSpecialCut_h, jtPfCHF_Bad_Unweighted_h, jtPfCEF_Bad_Unweighted_h, jtPfNHF_Bad_Unweighted_h, jtPfNEF_Bad_Unweighted_h, jtPfMUF_Bad_Unweighted_h, jtPfCHF_jtPfCEF_Good_Unweighted_h, jtPfCHF_jtPfNHF_Good_Unweighted_h, jtPfCHF_jtPfNEF_Good_Unweighted_h, jtPfCHF_jtPfMUF_Good_Unweighted_h, jtPfCEF_jtPfNHF_Good_Unweighted_h, jtPfCEF_jtPfNEF_Good_Unweighted_h, jtPfCEF_jtPfMUF_Good_Unweighted_h, jtPfNHF_jtPfNEF_Good_Unweighted_h, jtPfNHF_jtPfMUF_Good_Unweighted_h, jtPfNEF_jtPfMUF_Good_Unweighted_h, jtPfCHF_jtPfCEF_Bad_Unweighted_h, jtPfCHF_jtPfNHF_Bad_Unweighted_h, jtPfCHF_jtPfNEF_Bad_Unweighted_h, jtPfCHF_jtPfMUF_Bad_Unweighted_h, jtPfCEF_jtPfNHF_Bad_Unweighted_h, jtPfCEF_jtPfNEF_Bad_Unweighted_h, jtPfCEF_jtPfMUF_Bad_Unweighted_h, jtPfNHF_jtPfNEF_Bad_Unweighted_h, jtPfNHF_jtPfMUF_Bad_Unweighted_h, jtPfNEF_jtPfMUF_Bad_Unweighted_h, jtPfCHMF_Good_Unweighted_h, jtPfCEMF_Good_Unweighted_h, jtPfNHMF_Good_Unweighted_h, jtPfNEMF_Good_Unweighted_h, jtPfMUMF_Good_Unweighted_h, jtPfCHM_Good_Unweighted_h, jtPfCEM_Good_Unweighted_h, jtPfNHM_Good_Unweighted_h, jtPfNEM_Good_Unweighted_h, jtPfMUM_Good_Unweighted_h, jtPt_Bad_Unweighted_h, jtPfCHMF_Bad_Unweighted_h, jtPfCEMF_Bad_Unweighted_h, jtPfNHMF_Bad_Unweighted_h, jtPfNEMF_Bad_Unweighted_h, jtPfMUMF_Bad_Unweighted_h, jtPfCHM_Bad_Unweighted_h, jtPfCEM_Bad_Unweighted_h, jtPfNHM_Bad_Unweighted_h, jtPfNEM_Bad_Unweighted_h, jtPfMUM_Bad_Unweighted_h, jtPfCHMF_jtPfCEMF_Good_Unweighted_h, jtPfCHMF_jtPfNHMF_Good_Unweighted_h, jtPfCHMF_jtPfNEMF_Good_Unweighted_h, jtPfCHMF_jtPfMUMF_Good_Unweighted_h, jtPfCEMF_jtPfNHMF_Good_Unweighted_h, jtPfCEMF_jtPfNEMF_Good_Unweighted_h, jtPfCEMF_jtPfMUMF_Good_Unweighted_h, jtPfNHMF_jtPfNEMF_Good_Unweighted_h, jtPfNHMF_jtPfMUMF_Good_Unweighted_h, jtPfNEMF_jtPfMUMF_Good_Unweighted_h, jtPfCHMF_jtPfCEMF_Bad_Unweighted_h, jtPfCHMF_jtPfNHMF_Bad_Unweighted_h, jtPfCHMF_jtPfNEMF_Bad_Unweighted_h, jtPfCHMF_jtPfMUMF_Bad_Unweighted_h, jtPfCEMF_jtPfNHMF_Bad_Unweighted_h, jtPfCEMF_jtPfNEMF_Bad_Unweighted_h, jtPfCEMF_jtPfMUMF_Bad_Unweighted_h, jtPfNHMF_jtPfNEMF_Bad_Unweighted_h, jtPfNHMF_jtPfMUMF_Bad_Unweighted_h, jtPfNEMF_jtPfMUMF_Bad_Unweighted_h});

  setSumW2({jtPt_Good_Unweighted_h, jtPt_Good_Unweighted_Cut_h, jtPt_Good_Unweighted_NoSpecialCut_h, jtPt_Good_Weighted_h, jtPt_Good_Weighted_Cut_h, jtPt_Good_Weighted_NoSpecialCut_h, jtPfCHF_Good_Unweighted_h, jtPfCEF_Good_Unweighted_h, jtPfNHF_Good_Unweighted_h, jtPfNEF_Good_Unweighted_h, jtPfMUF_Good_Unweighted_h, jtPt_Bad_Unweighted_h, jtPt_Bad_Unweighted_Cut_h, jtPt_Bad_Unweighted_NoSpecialCut_h, jtPt_Bad_Weighted_h, jtPt_Bad_Weighted_Cut_h, jtPt_Bad_Weighted_NoSpecialCut_h, jtPfCHF_Bad_Unweighted_h, jtPfCEF_Bad_Unweighted_h, jtPfNHF_Bad_Unweighted_h, jtPfNEF_Bad_Unweighted_h, jtPfMUF_Bad_Unweighted_h, jtPfCHF_jtPfCEF_Good_Unweighted_h, jtPfCHF_jtPfNHF_Good_Unweighted_h, jtPfCHF_jtPfNEF_Good_Unweighted_h, jtPfCHF_jtPfMUF_Good_Unweighted_h, jtPfCEF_jtPfNHF_Good_Unweighted_h, jtPfCEF_jtPfNEF_Good_Unweighted_h, jtPfCEF_jtPfMUF_Good_Unweighted_h, jtPfNHF_jtPfNEF_Good_Unweighted_h, jtPfNHF_jtPfMUF_Good_Unweighted_h, jtPfNEF_jtPfMUF_Good_Unweighted_h, jtPfCHF_jtPfCEF_Bad_Unweighted_h, jtPfCHF_jtPfNHF_Bad_Unweighted_h, jtPfCHF_jtPfNEF_Bad_Unweighted_h, jtPfCHF_jtPfMUF_Bad_Unweighted_h, jtPfCEF_jtPfNHF_Bad_Unweighted_h, jtPfCEF_jtPfNEF_Bad_Unweighted_h, jtPfCEF_jtPfMUF_Bad_Unweighted_h, jtPfNHF_jtPfNEF_Bad_Unweighted_h, jtPfNHF_jtPfMUF_Bad_Unweighted_h, jtPfNEF_jtPfMUF_Bad_Unweighted_h, jtPfCHMF_Good_Unweighted_h, jtPfCEMF_Good_Unweighted_h, jtPfNHMF_Good_Unweighted_h, jtPfNEMF_Good_Unweighted_h, jtPfMUMF_Good_Unweighted_h, jtPfCHM_Good_Unweighted_h, jtPfCEM_Good_Unweighted_h, jtPfNHM_Good_Unweighted_h, jtPfNEM_Good_Unweighted_h, jtPfMUM_Good_Unweighted_h, jtPt_Bad_Unweighted_h, jtPfCHMF_Bad_Unweighted_h, jtPfCEMF_Bad_Unweighted_h, jtPfNHMF_Bad_Unweighted_h, jtPfNEMF_Bad_Unweighted_h, jtPfMUMF_Bad_Unweighted_h, jtPfCHM_Bad_Unweighted_h, jtPfCEM_Bad_Unweighted_h, jtPfNHM_Bad_Unweighted_h, jtPfNEM_Bad_Unweighted_h, jtPfMUM_Bad_Unweighted_h, jtPfCHMF_jtPfCEMF_Good_Unweighted_h, jtPfCHMF_jtPfNHMF_Good_Unweighted_h, jtPfCHMF_jtPfNEMF_Good_Unweighted_h, jtPfCHMF_jtPfMUMF_Good_Unweighted_h, jtPfCEMF_jtPfNHMF_Good_Unweighted_h, jtPfCEMF_jtPfNEMF_Good_Unweighted_h, jtPfCEMF_jtPfMUMF_Good_Unweighted_h, jtPfNHMF_jtPfNEMF_Good_Unweighted_h, jtPfNHMF_jtPfMUMF_Good_Unweighted_h, jtPfNEMF_jtPfMUMF_Good_Unweighted_h, jtPfCHMF_jtPfCEMF_Bad_Unweighted_h, jtPfCHMF_jtPfNHMF_Bad_Unweighted_h, jtPfCHMF_jtPfNEMF_Bad_Unweighted_h, jtPfCHMF_jtPfMUMF_Bad_Unweighted_h, jtPfCEMF_jtPfNHMF_Bad_Unweighted_h, jtPfCEMF_jtPfNEMF_Bad_Unweighted_h, jtPfCEMF_jtPfMUMF_Bad_Unweighted_h, jtPfNHMF_jtPfNEMF_Bad_Unweighted_h, jtPfNHMF_jtPfMUMF_Bad_Unweighted_h, jtPfNEMF_jtPfMUMF_Bad_Unweighted_h});

  const Int_t nMaxJet_ = 500;
  Float_t pthat_;
  Int_t nref_;
  Float_t jtpt_[nMaxJet_];
  Float_t refpt_[nMaxJet_];
  Float_t jteta_[nMaxJet_];
  Float_t jtphi_[nMaxJet_];
  Float_t jtPfCHF_[nMaxJet_];
  Float_t jtPfCEF_[nMaxJet_];
  Float_t jtPfNHF_[nMaxJet_];
  Float_t jtPfNEF_[nMaxJet_];
  Float_t jtPfMUF_[nMaxJet_];

  Float_t jtPfCHMF_[nMaxJet_];
  Float_t jtPfCEMF_[nMaxJet_];
  Float_t jtPfNHMF_[nMaxJet_];
  Float_t jtPfNEMF_[nMaxJet_];
  Float_t jtPfMUMF_[nMaxJet_];

  Int_t jtPfCHM_[nMaxJet_];
  Int_t jtPfCEM_[nMaxJet_];
  Int_t jtPfNHM_[nMaxJet_];
  Int_t jtPfNEM_[nMaxJet_];
  Int_t jtPfMUM_[nMaxJet_];

  Int_t ngen_;
  Float_t genpt_[nMaxJet_];
  Float_t geneta_[nMaxJet_];
  Float_t genphi_[nMaxJet_];
  Int_t gensubid_[nMaxJet_];

  unsigned int run_, lumi_;
  unsigned long long evt_;
  Int_t hiBin_;
  Float_t vz_;
  Float_t hiHF_;

  Int_t HBHENoiseFilterResultRun2Loose_ = -1;
  Int_t pprimaryVertexFilter_ = -1;
  Int_t pBeamScrapingFilter_ = -1;
  Int_t phfCoincFilter3_ = -1;
  Int_t pclusterCompatibilityFilter_ = -1;

  goodGlobalSelection globalSel;
  globalSel.setIsPbPb(!isPP);

  specialHYDJETEventExclude specialSel;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "Processing file " << fI << "/" << fileList.size() << ": \'" << fileList.at(fI) << "\'" << std::endl;

    TFile* inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(fI)).c_str(), "READ");
    TTree* jetTree_p = (TTree*)inFile_p->Get("akCs4PU3PFFlowJetAnalyzer/t");
    jetTree_p->SetBranchStatus("*", 0);
    jetTree_p->SetBranchStatus("pthat", 1);
    jetTree_p->SetBranchStatus("nref", 1);
    jetTree_p->SetBranchStatus("jtpt", 1);
    jetTree_p->SetBranchStatus("refpt", 1);
    jetTree_p->SetBranchStatus("jteta", 1);
    jetTree_p->SetBranchStatus("jtphi", 1);
    jetTree_p->SetBranchStatus("jtPfCHF", 1);
    jetTree_p->SetBranchStatus("jtPfCEF", 1);
    jetTree_p->SetBranchStatus("jtPfNHF", 1);
    jetTree_p->SetBranchStatus("jtPfNEF", 1);
    jetTree_p->SetBranchStatus("jtPfMUF", 1);
    jetTree_p->SetBranchStatus("jtPfCHMF", 1);
    jetTree_p->SetBranchStatus("jtPfCEMF", 1);
    jetTree_p->SetBranchStatus("jtPfNHMF", 1);
    jetTree_p->SetBranchStatus("jtPfNEMF", 1);
    jetTree_p->SetBranchStatus("jtPfMUMF", 1);
    jetTree_p->SetBranchStatus("jtPfCHM", 1);
    jetTree_p->SetBranchStatus("jtPfCEM", 1);
    jetTree_p->SetBranchStatus("jtPfNHM", 1);
    jetTree_p->SetBranchStatus("jtPfNEM", 1);
    jetTree_p->SetBranchStatus("jtPfMUM", 1);
    jetTree_p->SetBranchStatus("ngen", 1);
    jetTree_p->SetBranchStatus("genpt", 1);
    jetTree_p->SetBranchStatus("geneta", 1);
    jetTree_p->SetBranchStatus("genphi", 1);
    jetTree_p->SetBranchStatus("gensubid", 1);

    jetTree_p->SetBranchAddress("pthat", &pthat_);
    jetTree_p->SetBranchAddress("nref", &nref_);
    jetTree_p->SetBranchAddress("jtpt", jtpt_);
    jetTree_p->SetBranchAddress("refpt", refpt_);
    jetTree_p->SetBranchAddress("jteta", jteta_);
    jetTree_p->SetBranchAddress("jtphi", jtphi_);
    jetTree_p->SetBranchAddress("jtPfCHF", jtPfCHF_);
    jetTree_p->SetBranchAddress("jtPfCEF", jtPfCEF_);
    jetTree_p->SetBranchAddress("jtPfNHF", jtPfNHF_);
    jetTree_p->SetBranchAddress("jtPfNEF", jtPfNEF_);
    jetTree_p->SetBranchAddress("jtPfMUF", jtPfMUF_);
    jetTree_p->SetBranchAddress("jtPfCHMF", jtPfCHMF_);
    jetTree_p->SetBranchAddress("jtPfCEMF", jtPfCEMF_);
    jetTree_p->SetBranchAddress("jtPfNHMF", jtPfNHMF_);
    jetTree_p->SetBranchAddress("jtPfNEMF", jtPfNEMF_);
    jetTree_p->SetBranchAddress("jtPfMUMF", jtPfMUMF_);
    jetTree_p->SetBranchAddress("jtPfCHM", jtPfCHM_);
    jetTree_p->SetBranchAddress("jtPfCEM", jtPfCEM_);
    jetTree_p->SetBranchAddress("jtPfNHM", jtPfNHM_);
    jetTree_p->SetBranchAddress("jtPfNEM", jtPfNEM_);
    jetTree_p->SetBranchAddress("jtPfMUM", jtPfMUM_);
    jetTree_p->SetBranchAddress("ngen", &ngen_);
    jetTree_p->SetBranchAddress("genpt", genpt_);
    jetTree_p->SetBranchAddress("geneta", geneta_);
    jetTree_p->SetBranchAddress("genphi", genphi_);
    jetTree_p->SetBranchAddress("gensubid", gensubid_);

    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");

    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("hiHF", 1);
    hiTree_p->SetBranchStatus("vz", 1);

    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("evt", &evt_);
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("hiHF", &hiHF_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    
    TTree* skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");

    skimTree_p->SetBranchStatus("*", 0);
    skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);

    if(!isPP){
      skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
      skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
      skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);

      skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
      skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
      skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);
    }
    else{
      skimTree_p->SetBranchStatus("pBeamScrapingFilter", 1);
      skimTree_p->SetBranchStatus("pPAprimaryVertexFilter", 1);

      skimTree_p->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter_);
      skimTree_p->SetBranchAddress("pPAprimaryVertexFilter", &pprimaryVertexFilter_);
    }

    const Int_t nEntries = jetTree_p->GetEntries();
    const Int_t printInterval = TMath::Max(1, nEntries/20);

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(nEntries >= 50000 && entry%printInterval == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;

      jetTree_p->GetEntry(entry);
      hiTree_p->GetEntry(entry);
      skimTree_p->GetEntry(entry);

      globalSel.setVz(vz_);
      globalSel.setHiHF(hiHF_);
      globalSel.setPprimaryVertexFilter(pprimaryVertexFilter_);
      globalSel.setPBeamScrapingFilter(pBeamScrapingFilter_);
      globalSel.setPhfCoincFilter3(phfCoincFilter3_);
      globalSel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
      globalSel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);

      if(!globalSel.isGood()) continue;
      if(hiBin_ > 60) continue;

      bool badJetSpecialSel = specialSel.CheckEventBadJet(ngen_, genpt_, genphi_, geneta_, gensubid_);

      Double_t pthatWeight_ = -1;
      for(unsigned int pI = 0; pI < pthats.size()-1; ++pI){
        if(pthats.at(pI) <= pthat_ && pthat_ < pthats.at(pI+1)){
          pthatWeight_ = pthatWeights.at(pI);
          break;
        }
      }
      if(pthat_ > pthats.at(pthats.size()-1)) pthatWeight_ = pthatWeights.at(pthatWeights.size()-1);

      if(pthatWeight_ < 0){
	std::cout << "WARNING - NO WEIGHT FOR pthat \'" << pthat_ << "\'. Set to 1" << std::endl;
        pthatWeight_ = 1.;
      }

      //      Double_t ncollWeight_ = findNcoll_Renorm(hiBin_);
      //      Double_t fullWeight_ = ncollWeight_*pthatWeight_;

      for(Int_t jI = 0; jI < nref_; ++jI){
	if(TMath::Abs(jteta_[jI]) > jtAbsEtaMax) continue;
	if(jtpt_[jI] < jtRecoPtMin) continue;

	bool isGood = refpt_[jI] > jtRefGoodPtMin;
	bool isBad = refpt_[jI] < jtRefBadPtMax;
	bool isNotCut = jtPfMUMF_[jI] < 0.6 && jtPfCHMF_[jI] < 0.9;

	if(isBad && jtpt_[jI] >= 600 && jtpt_[jI] < 750 && isNotCut){
	  std::cout << "Bad jet at entry: " << entry << std::endl;
	  std::cout << " pt, refpt, eta, phi: " << jtpt_[jI] << ", " << refpt_[jI] << ", " << jteta_[jI] << ", " << jtphi_[jI] << std::endl;
	}

	if(isGood){
	  jtPt_Good_Unweighted_NoSpecialCut_h->Fill(jtpt_[jI]);
	  jtPt_Good_Weighted_NoSpecialCut_h->Fill(jtpt_[jI], pthatWeight_);	  
	}
	else{
	  jtPt_Bad_Unweighted_NoSpecialCut_h->Fill(jtpt_[jI]);
	  jtPt_Bad_Weighted_NoSpecialCut_h->Fill(jtpt_[jI], pthatWeight_);	  
	}

	if(!isPP && badJetSpecialSel) continue;

	if(isGood){
	  jtPt_Good_Unweighted_h->Fill(jtpt_[jI]);
	  jtPt_Good_Weighted_h->Fill(jtpt_[jI], pthatWeight_);

	  if(isNotCut){
	    jtPt_Good_Unweighted_Cut_h->Fill(jtpt_[jI]);
	    jtPt_Good_Weighted_Cut_h->Fill(jtpt_[jI], pthatWeight_);
	  }

	  jtPfCHF_Good_Unweighted_h->Fill(jtPfCHF_[jI]);
	  jtPfCEF_Good_Unweighted_h->Fill(jtPfCEF_[jI]);
	  jtPfNHF_Good_Unweighted_h->Fill(jtPfNHF_[jI]);
	  jtPfNEF_Good_Unweighted_h->Fill(jtPfNEF_[jI]);
	  jtPfMUF_Good_Unweighted_h->Fill(jtPfMUF_[jI]);

	  jtPfCHM_Good_Unweighted_h->Fill(jtPfCHM_[jI]);
	  jtPfCEM_Good_Unweighted_h->Fill(jtPfCEM_[jI]);
	  jtPfNHM_Good_Unweighted_h->Fill(jtPfNHM_[jI]);
	  jtPfNEM_Good_Unweighted_h->Fill(jtPfNEM_[jI]);
	  jtPfMUM_Good_Unweighted_h->Fill(jtPfMUM_[jI]);

	  jtPfCHF_jtPfCEF_Good_Unweighted_h->Fill(jtPfCHF_[jI], jtPfCEF_[jI]);
	  jtPfCHF_jtPfNHF_Good_Unweighted_h->Fill(jtPfCHF_[jI], jtPfNHF_[jI]);
	  jtPfCHF_jtPfNEF_Good_Unweighted_h->Fill(jtPfCHF_[jI], jtPfNEF_[jI]);
	  jtPfCHF_jtPfMUF_Good_Unweighted_h->Fill(jtPfCHF_[jI], jtPfMUF_[jI]);

	  jtPfCEF_jtPfNHF_Good_Unweighted_h->Fill(jtPfCEF_[jI], jtPfNHF_[jI]);
	  jtPfCEF_jtPfNEF_Good_Unweighted_h->Fill(jtPfCEF_[jI], jtPfNEF_[jI]);
	  jtPfCEF_jtPfMUF_Good_Unweighted_h->Fill(jtPfCEF_[jI], jtPfMUF_[jI]);

	  jtPfNHF_jtPfNEF_Good_Unweighted_h->Fill(jtPfNHF_[jI], jtPfNEF_[jI]);
	  jtPfNHF_jtPfMUF_Good_Unweighted_h->Fill(jtPfNHF_[jI], jtPfMUF_[jI]);

	  jtPfNEF_jtPfMUF_Good_Unweighted_h->Fill(jtPfNEF_[jI], jtPfMUF_[jI]);

	  jtPfCHMF_Good_Unweighted_h->Fill(jtPfCHMF_[jI]);
	  jtPfCEMF_Good_Unweighted_h->Fill(jtPfCEMF_[jI]);
	  jtPfNHMF_Good_Unweighted_h->Fill(jtPfNHMF_[jI]);
	  jtPfNEMF_Good_Unweighted_h->Fill(jtPfNEMF_[jI]);
	  jtPfMUMF_Good_Unweighted_h->Fill(jtPfMUMF_[jI]);

	  jtPfCHMF_jtPfCEMF_Good_Unweighted_h->Fill(jtPfCHMF_[jI], jtPfCEMF_[jI]);
	  jtPfCHMF_jtPfNHMF_Good_Unweighted_h->Fill(jtPfCHMF_[jI], jtPfNHMF_[jI]);
	  jtPfCHMF_jtPfNEMF_Good_Unweighted_h->Fill(jtPfCHMF_[jI], jtPfNEMF_[jI]);
	  jtPfCHMF_jtPfMUMF_Good_Unweighted_h->Fill(jtPfCHMF_[jI], jtPfMUMF_[jI]);

	  jtPfCEMF_jtPfNHMF_Good_Unweighted_h->Fill(jtPfCEMF_[jI], jtPfNHMF_[jI]);
	  jtPfCEMF_jtPfNEMF_Good_Unweighted_h->Fill(jtPfCEMF_[jI], jtPfNEMF_[jI]);
	  jtPfCEMF_jtPfMUMF_Good_Unweighted_h->Fill(jtPfCEMF_[jI], jtPfMUMF_[jI]);

	  jtPfNHMF_jtPfNEMF_Good_Unweighted_h->Fill(jtPfNHMF_[jI], jtPfNEMF_[jI]);
	  jtPfNHMF_jtPfMUMF_Good_Unweighted_h->Fill(jtPfNHMF_[jI], jtPfMUMF_[jI]);

	  jtPfNEMF_jtPfMUMF_Good_Unweighted_h->Fill(jtPfNEMF_[jI], jtPfMUMF_[jI]);
	}
	else if(isBad){
	  jtPt_Bad_Unweighted_h->Fill(jtpt_[jI]);
	  jtPt_Bad_Weighted_h->Fill(jtpt_[jI], pthatWeight_);

	  if(isNotCut){
	    jtPt_Bad_Unweighted_Cut_h->Fill(jtpt_[jI]);
	    jtPt_Bad_Weighted_Cut_h->Fill(jtpt_[jI], pthatWeight_);
	  }


	  jtPfCHF_Bad_Unweighted_h->Fill(jtPfCHF_[jI]);
	  jtPfCEF_Bad_Unweighted_h->Fill(jtPfCEF_[jI]);
	  jtPfNHF_Bad_Unweighted_h->Fill(jtPfNHF_[jI]);
	  jtPfNEF_Bad_Unweighted_h->Fill(jtPfNEF_[jI]);
	  jtPfMUF_Bad_Unweighted_h->Fill(jtPfMUF_[jI]);

	  jtPfCHM_Bad_Unweighted_h->Fill(jtPfCHM_[jI]);
	  jtPfCEM_Bad_Unweighted_h->Fill(jtPfCEM_[jI]);
	  jtPfNHM_Bad_Unweighted_h->Fill(jtPfNHM_[jI]);
	  jtPfNEM_Bad_Unweighted_h->Fill(jtPfNEM_[jI]);
	  jtPfMUM_Bad_Unweighted_h->Fill(jtPfMUM_[jI]);


	  jtPfCHF_jtPfCEF_Bad_Unweighted_h->Fill(jtPfCHF_[jI], jtPfCEF_[jI]);
	  jtPfCHF_jtPfNHF_Bad_Unweighted_h->Fill(jtPfCHF_[jI], jtPfNHF_[jI]);
	  jtPfCHF_jtPfNEF_Bad_Unweighted_h->Fill(jtPfCHF_[jI], jtPfNEF_[jI]);
	  jtPfCHF_jtPfMUF_Bad_Unweighted_h->Fill(jtPfCHF_[jI], jtPfMUF_[jI]);

	  jtPfCEF_jtPfNHF_Bad_Unweighted_h->Fill(jtPfCEF_[jI], jtPfNHF_[jI]);
	  jtPfCEF_jtPfNEF_Bad_Unweighted_h->Fill(jtPfCEF_[jI], jtPfNEF_[jI]);
	  jtPfCEF_jtPfMUF_Bad_Unweighted_h->Fill(jtPfCEF_[jI], jtPfMUF_[jI]);

	  jtPfNHF_jtPfNEF_Bad_Unweighted_h->Fill(jtPfNHF_[jI], jtPfNEF_[jI]);
	  jtPfNHF_jtPfMUF_Bad_Unweighted_h->Fill(jtPfNHF_[jI], jtPfMUF_[jI]);

	  jtPfNEF_jtPfMUF_Bad_Unweighted_h->Fill(jtPfNEF_[jI], jtPfMUF_[jI]);	  

	  jtPfCHMF_Bad_Unweighted_h->Fill(jtPfCHMF_[jI]);
	  jtPfCEMF_Bad_Unweighted_h->Fill(jtPfCEMF_[jI]);
	  jtPfNHMF_Bad_Unweighted_h->Fill(jtPfNHMF_[jI]);
	  jtPfNEMF_Bad_Unweighted_h->Fill(jtPfNEMF_[jI]);
	  jtPfMUMF_Bad_Unweighted_h->Fill(jtPfMUMF_[jI]);

	  jtPfCHMF_jtPfCEMF_Bad_Unweighted_h->Fill(jtPfCHMF_[jI], jtPfCEMF_[jI]);
	  jtPfCHMF_jtPfNHMF_Bad_Unweighted_h->Fill(jtPfCHMF_[jI], jtPfNHMF_[jI]);
	  jtPfCHMF_jtPfNEMF_Bad_Unweighted_h->Fill(jtPfCHMF_[jI], jtPfNEMF_[jI]);
	  jtPfCHMF_jtPfMUMF_Bad_Unweighted_h->Fill(jtPfCHMF_[jI], jtPfMUMF_[jI]);

	  jtPfCEMF_jtPfNHMF_Bad_Unweighted_h->Fill(jtPfCEMF_[jI], jtPfNHMF_[jI]);
	  jtPfCEMF_jtPfNEMF_Bad_Unweighted_h->Fill(jtPfCEMF_[jI], jtPfNEMF_[jI]);
	  jtPfCEMF_jtPfMUMF_Bad_Unweighted_h->Fill(jtPfCEMF_[jI], jtPfMUMF_[jI]);

	  jtPfNHMF_jtPfNEMF_Bad_Unweighted_h->Fill(jtPfNHMF_[jI], jtPfNEMF_[jI]);
	  jtPfNHMF_jtPfMUMF_Bad_Unweighted_h->Fill(jtPfNHMF_[jI], jtPfMUMF_[jI]);

	  jtPfNEMF_jtPfMUMF_Bad_Unweighted_h->Fill(jtPfNEMF_[jI], jtPfMUMF_[jI]);	  
	}
      }
    }

    inFile_p->Close();
    delete inFile_p;
  }

  std::vector<TH1*> goodHists;
  std::vector<TH1*> badHists;

  std::vector<TH2*> goodHistsTH2;
  std::vector<TH2*> badHistsTH2;

  goodHists.push_back(jtPt_Good_Unweighted_h);
  goodHists.push_back(jtPt_Good_Unweighted_Cut_h);
  goodHists.push_back(jtPt_Good_Unweighted_NoSpecialCut_h);
  goodHists.push_back(jtPt_Good_Weighted_h);
  goodHists.push_back(jtPt_Good_Weighted_Cut_h);
  goodHists.push_back(jtPt_Good_Weighted_NoSpecialCut_h);

  goodHists.push_back(jtPfCHF_Good_Unweighted_h);
  goodHists.push_back(jtPfCEF_Good_Unweighted_h);
  goodHists.push_back(jtPfNHF_Good_Unweighted_h);
  goodHists.push_back(jtPfNEF_Good_Unweighted_h);
  goodHists.push_back(jtPfMUF_Good_Unweighted_h);

  goodHists.push_back(jtPfCHM_Good_Unweighted_h);
  goodHists.push_back(jtPfCEM_Good_Unweighted_h);
  goodHists.push_back(jtPfNHM_Good_Unweighted_h);
  goodHists.push_back(jtPfNEM_Good_Unweighted_h);
  goodHists.push_back(jtPfMUM_Good_Unweighted_h);

  goodHistsTH2.push_back(jtPfCHF_jtPfCEF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfCHF_jtPfNHF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfCHF_jtPfNEF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfCHF_jtPfMUF_Good_Unweighted_h);

  goodHistsTH2.push_back(jtPfCEF_jtPfNHF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfCEF_jtPfNEF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfCEF_jtPfMUF_Good_Unweighted_h);

  goodHistsTH2.push_back(jtPfNHF_jtPfNEF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfNHF_jtPfMUF_Good_Unweighted_h);

  goodHistsTH2.push_back(jtPfNEF_jtPfMUF_Good_Unweighted_h);

  goodHists.push_back(jtPfCHMF_Good_Unweighted_h);
  goodHists.push_back(jtPfCEMF_Good_Unweighted_h);
  goodHists.push_back(jtPfNHMF_Good_Unweighted_h);
  goodHists.push_back(jtPfNEMF_Good_Unweighted_h);
  goodHists.push_back(jtPfMUMF_Good_Unweighted_h);

  goodHistsTH2.push_back(jtPfCHMF_jtPfCEMF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfCHMF_jtPfNHMF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfCHMF_jtPfNEMF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfCHMF_jtPfMUMF_Good_Unweighted_h);

  goodHistsTH2.push_back(jtPfCEMF_jtPfNHMF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfCEMF_jtPfNEMF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfCEMF_jtPfMUMF_Good_Unweighted_h);

  goodHistsTH2.push_back(jtPfNHMF_jtPfNEMF_Good_Unweighted_h);
  goodHistsTH2.push_back(jtPfNHMF_jtPfMUMF_Good_Unweighted_h);

  goodHistsTH2.push_back(jtPfNEMF_jtPfMUMF_Good_Unweighted_h);

  badHists.push_back(jtPt_Bad_Unweighted_h);
  badHists.push_back(jtPt_Bad_Unweighted_Cut_h);
  badHists.push_back(jtPt_Bad_Unweighted_NoSpecialCut_h);
  badHists.push_back(jtPt_Bad_Weighted_h);
  badHists.push_back(jtPt_Bad_Weighted_Cut_h);
  badHists.push_back(jtPt_Bad_Weighted_NoSpecialCut_h);

  badHists.push_back(jtPfCHF_Bad_Unweighted_h);
  badHists.push_back(jtPfCEF_Bad_Unweighted_h);
  badHists.push_back(jtPfNHF_Bad_Unweighted_h);
  badHists.push_back(jtPfNEF_Bad_Unweighted_h);
  badHists.push_back(jtPfMUF_Bad_Unweighted_h);

  badHists.push_back(jtPfCHM_Bad_Unweighted_h);
  badHists.push_back(jtPfCEM_Bad_Unweighted_h);
  badHists.push_back(jtPfNHM_Bad_Unweighted_h);
  badHists.push_back(jtPfNEM_Bad_Unweighted_h);
  badHists.push_back(jtPfMUM_Bad_Unweighted_h);

  badHistsTH2.push_back(jtPfCHF_jtPfCEF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfCHF_jtPfNHF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfCHF_jtPfNEF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfCHF_jtPfMUF_Bad_Unweighted_h);

  badHistsTH2.push_back(jtPfCEF_jtPfNHF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfCEF_jtPfNEF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfCEF_jtPfMUF_Bad_Unweighted_h);

  badHistsTH2.push_back(jtPfNHF_jtPfNEF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfNHF_jtPfMUF_Bad_Unweighted_h);

  badHistsTH2.push_back(jtPfNEF_jtPfMUF_Bad_Unweighted_h);

  badHists.push_back(jtPfCHMF_Bad_Unweighted_h);
  badHists.push_back(jtPfCEMF_Bad_Unweighted_h);
  badHists.push_back(jtPfNHMF_Bad_Unweighted_h);
  badHists.push_back(jtPfNEMF_Bad_Unweighted_h);
  badHists.push_back(jtPfMUMF_Bad_Unweighted_h);

  badHistsTH2.push_back(jtPfCHMF_jtPfCEMF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfCHMF_jtPfNHMF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfCHMF_jtPfNEMF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfCHMF_jtPfMUMF_Bad_Unweighted_h);

  badHistsTH2.push_back(jtPfCEMF_jtPfNHMF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfCEMF_jtPfNEMF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfCEMF_jtPfMUMF_Bad_Unweighted_h);

  badHistsTH2.push_back(jtPfNHMF_jtPfNEMF_Bad_Unweighted_h);
  badHistsTH2.push_back(jtPfNHMF_jtPfMUMF_Bad_Unweighted_h);

  badHistsTH2.push_back(jtPfNEMF_jtPfMUMF_Bad_Unweighted_h);

  outFile_p->cd();

  checkMakeDir("pdfDir");

  vanGoghPalette vg;

  const int nPlots = 2;
  const int styles[nPlots] = {20, 21};
  const int colors[nPlots] = {vg.getColor(0), vg.getColor(1)};

  double maxValPt = -1;
  double minValPt = 100000;

  for(unsigned int i = 0; i < goodHists.size(); ++i){
    std::string canvStr = goodHists.at(i)->GetName();
    if(canvStr.find("jtPt_") == std::string::npos) continue;

    for(Int_t bIX = 0; bIX < goodHists.at(i)->GetNbinsX(); ++bIX){
      if(maxValPt < goodHists.at(i)->GetBinContent(bIX+1)) maxValPt = goodHists.at(i)->GetBinContent(bIX+1);
      if(minValPt > goodHists.at(i)->GetBinContent(bIX+1) && goodHists.at(i)->GetBinContent(bIX+1) > 0) minValPt = goodHists.at(i)->GetBinContent(bIX+1);

      if(maxValPt < badHists.at(i)->GetBinContent(bIX+1)) maxValPt = badHists.at(i)->GetBinContent(bIX+1);
      if(minValPt > badHists.at(i)->GetBinContent(bIX+1) && badHists.at(i)->GetBinContent(bIX+1) > 0) minValPt = badHists.at(i)->GetBinContent(bIX+1);
    }
  }

  maxValPt *= 5.;
  minValPt /= 5.;

  for(unsigned int i = 0; i < goodHists.size(); ++i){
    std::string canvStr = goodHists.at(i)->GetName();
    canvStr.replace(canvStr.find("_Good"), std::string("_Good").size(), "");

    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.01);
    
    if(canvStr.find("jtPt_") == std::string::npos){
      goodHists.at(i)->Scale(1./goodHists.at(i)->Integral());
      badHists.at(i)->Scale(1./badHists.at(i)->Integral());

      goodHists.at(i)->GetYaxis()->SetTitle("Counts (Norm. To 1)");
      badHists.at(i)->GetYaxis()->SetTitle("Counts (Norm. To 1)");
    }
    else{
      goodHists.at(i)->SetMaximum(maxValPt);
      goodHists.at(i)->SetMinimum(minValPt);

      badHists.at(i)->SetMaximum(maxValPt);
      badHists.at(i)->SetMinimum(minValPt);
    }


    if(canvStr.find("jtPt_") == std::string::npos){
      double maxVal = -1;
      double minVal = 9999;
      
      for(Int_t bIX = 0; bIX < goodHists.at(i)->GetNbinsX(); ++bIX){
	if(maxVal < goodHists.at(i)->GetBinContent(bIX+1)) maxVal = goodHists.at(i)->GetBinContent(bIX+1);
	if(maxVal < badHists.at(i)->GetBinContent(bIX+1)) maxVal = badHists.at(i)->GetBinContent(bIX+1);
	
	if(minVal > goodHists.at(i)->GetBinContent(bIX+1) && goodHists.at(i)->GetBinContent(bIX+1) > 0) minVal = goodHists.at(i)->GetBinContent(bIX+1);
	if(minVal > badHists.at(i)->GetBinContent(bIX+1) && badHists.at(i)->GetBinContent(bIX+1) > 0) minVal = badHists.at(i)->GetBinContent(bIX+1);
      }
      
      goodHists.at(i)->SetMaximum(maxVal*5.);
      goodHists.at(i)->SetMinimum(minVal/5.);
    }

    goodHists.at(i)->SetMarkerColor(colors[0]);
    goodHists.at(i)->SetMarkerStyle(styles[0]);
    goodHists.at(i)->SetMarkerSize(0.8);
    goodHists.at(i)->SetLineColor(colors[0]);

    badHists.at(i)->SetMarkerColor(colors[1]);
    badHists.at(i)->SetMarkerStyle(styles[1]);
    badHists.at(i)->SetMarkerSize(0.6);
    badHists.at(i)->SetLineColor(colors[1]);
    
    canv_p->cd();

    goodHists.at(i)->DrawCopy("HIST E1 P");
    badHists.at(i)->DrawCopy("HIST E1 P SAME");
    gStyle->SetOptStat(0);

    gPad->SetLogy();

    canv_p->SaveAs(("pdfDir/" + canvStr + "_" + dateStr + ".pdf").c_str());
    delete canv_p;

    goodHists.at(i)->Write("", TObject::kOverwrite);
    badHists.at(i)->Write("", TObject::kOverwrite);

    delete goodHists.at(i);
    delete badHists.at(i);
  }


  for(unsigned int i = 0; i < goodHistsTH2.size(); ++i){
    std::string canvStr = goodHistsTH2.at(i)->GetName();
    canvStr.replace(canvStr.find("_Good"), std::string("_Good").size(), "");

    TCanvas* canv_p = new TCanvas("canv_p", "", 450*2, 450);
    
    goodHistsTH2.at(i)->Scale(1./goodHistsTH2.at(i)->Integral());
    badHistsTH2.at(i)->Scale(1./badHistsTH2.at(i)->Integral());

    double maxVal = -1;
    double minVal = 9999;

    for(Int_t bIX = 0; bIX < goodHistsTH2.at(i)->GetNbinsX(); ++bIX){
      for(Int_t bIY = 0; bIY < goodHistsTH2.at(i)->GetNbinsY(); ++bIY){
	if(maxVal < goodHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1)) maxVal = goodHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1);
	if(maxVal < badHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1)) maxVal = badHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1);
	
	if(minVal > goodHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1) && goodHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1) > 0) minVal = goodHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1);
	if(minVal > badHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1) && badHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1) > 0) minVal = badHistsTH2.at(i)->GetBinContent(bIX+1, bIY+1);
      }
    }

    goodHistsTH2.at(i)->SetMaximum(maxVal*5.);
    goodHistsTH2.at(i)->SetMinimum(minVal/5.);

    canv_p->Divide(2, 1);
    canv_p->cd();
    canv_p->cd(1);

    goodHistsTH2.at(i)->SetTitle("GOOD");
    goodHistsTH2.at(i)->DrawCopy("COLZ");
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    canv_p->cd();
    canv_p->cd(2);

    badHistsTH2.at(i)->SetTitle("BAD");
    badHistsTH2.at(i)->DrawCopy("COLZ");
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    canv_p->SaveAs(("pdfDir/" + canvStr + "_" + dateStr + ".pdf").c_str());
    delete canv_p;

    goodHistsTH2.at(i)->Write("", TObject::kOverwrite);
    badHistsTH2.at(i)->Write("", TObject::kOverwrite);

    delete goodHistsTH2.at(i);
    delete badHistsTH2.at(i);
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage: ./bin/checkGoodJetBadJetPF.exe <inFileName> <isPP-Opt>" << std::endl;
    return 1;
  }

  int retVal = 0; 
  if(argc == 2) retVal += checkGoodJetBadJetPF(argv[1]);
  else if(argc == 3) retVal += checkGoodJetBadJetPF(argv[1], std::stoi(argv[2]));
  return retVal;
}
