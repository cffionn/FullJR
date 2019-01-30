//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TBox.h"
#include "TCanvas.h"
#include "TCollection.h"
#include "TDatime.h"
#include "TDirectory.h"
#include "TError.h"
#include "TFile.h"
#include "TH1D.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/doLocalDebug.h"
#include "MainAnalysis/include/smallOrLargeR.h"
#include "MainAnalysis/include/systFunctions.h"
#include "MainAnalysis/include/texSlideCreator.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/doGlobalDebug.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/getLogBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vanGoghPalette.h"

int plotUnfoldedRAA(std::vector<std::string> inFileList)
{
  //Double check files
  std::cout << "Checking fileList..." << std::endl;
  unsigned int posFI = 0;

  std::vector<std::string> listOfPbPbAlgo;
  std::vector<std::string> fileOfPbPbAlgo;
  std::vector<int> ppMatchedToPbPb;
  std::vector<std::string> listOfPPAlgo;
  std::vector<std::string> fileOfPPAlgo;
  std::vector<int> isPPAlgoUsed;

  const Int_t nJtMax = 10;
  const Int_t nMaxCentBins = 4;
  const Int_t nMaxID = 6;
  const Int_t nMaxSyst = 20;

  while(posFI < inFileList.size()){
    bool keep = true;

    if(inFileList.at(posFI).find(".root") == std::string::npos) keep = false;
    else if(!checkFile(inFileList.at(posFI))) keep = false;
    else{
      TFile* inFile_p = new TFile(inFileList.at(posFI).c_str(), "READ");
      std::vector<std::string> dirList = returnRootFileContentsList(inFile_p, "TDirectory", "JetAnalyzer");
      std::vector<std::string> dirList2 = returnRootFileContentsList(inFile_p, "TDirectoryFile", "JetAnalyzer");
      dirList.insert(dirList.end(), dirList2.begin(), dirList2.end());

      if(dirList.size() == 0){
	keep = false;	
      }
      else{
	for(unsigned int dI = 0; dI < dirList.size(); ++dI){
	  bool isPbPb = dirList.at(dI).find("akCs") != std::string::npos;
	  bool isUnique = true;
	  int pos = -1;

	  if(isPbPb){
	    for(unsigned int dI2 = 0; dI2 < listOfPbPbAlgo.size(); ++dI2){
	      if(isStrSame(listOfPbPbAlgo.at(dI2), dirList.at(dI))){
		pos = dI2;
		isUnique = false;
		break;
	      }
	    }
	  }
	  else{
	    for(unsigned int dI2 = 0; dI2 < listOfPPAlgo.size(); ++dI2){
	      if(isStrSame(listOfPPAlgo.at(dI2), dirList.at(dI))){
		pos = dI2;
		isUnique = false;
		break;
	      }
	    }
	  }
	
	  
	  if(!isUnique){
	    std::cout << "WARNING: ALGO \'" << dirList.at(dI) << "\' is not unique!" << std::endl;
	    if(isPbPb) std::cout << " First found in file \'" << fileOfPbPbAlgo.at(pos) << "\'" << std::endl;
	    else std::cout << " First found in file \'" << fileOfPPAlgo.at(pos) << "\'" << std::endl;
	    std::cout << " Second found in file \'" << inFileList.at(posFI) << "\'" << std::endl;
	    std::cout << " return 1" << std::endl;

	    inFile_p->Close();
	    delete inFile_p;
	    return 1;
	  }
	  else if(isPbPb){
	    listOfPbPbAlgo.push_back(dirList.at(dI));
	    fileOfPbPbAlgo.push_back(inFileList.at(posFI));
	  }
	  else{
	    listOfPPAlgo.push_back(dirList.at(dI));
	    fileOfPPAlgo.push_back(inFileList.at(posFI));
	    isPPAlgoUsed.push_back(0);
	  }
	}
      }
    
      inFile_p->Close();
      delete inFile_p;
    }

    if(keep) ++posFI;
    else{
      std::cout << " File \'" << inFileList.at(posFI) << "\' invalid. erasing" << std::endl;
      inFileList.erase(inFileList.begin()+posFI);      
    }
  }

  posFI = 0;
  while(posFI < listOfPbPbAlgo.size()){
    int rValPbPb = getRVal(listOfPbPbAlgo.at(posFI));

    int ppMatchPos = -1;
    
    for(unsigned int tI2 = 0; tI2 < listOfPPAlgo.size(); ++tI2){
      int rValPP = getRVal(listOfPPAlgo.at(tI2));
      
      if(rValPbPb == rValPP){
	ppMatchPos = tI2;
	break;
      }
    }

    if(ppMatchPos >= 0){
      ppMatchedToPbPb.push_back(ppMatchPos);
      isPPAlgoUsed.at(ppMatchPos) = 1;
      ++posFI;
    }
    else{
      std::cout << "No pp match found for \'" << listOfPbPbAlgo.at(posFI) << "\', removing" << std::endl;
      listOfPbPbAlgo.erase(listOfPbPbAlgo.begin() + posFI);
      fileOfPbPbAlgo.erase(listOfPbPbAlgo.begin() + posFI);
    }
  }

  posFI = 0;
  while(posFI < listOfPPAlgo.size()){
    if(isPPAlgoUsed.at(posFI) == 0){
      std::cout << "pp algo \'" << listOfPPAlgo.at(posFI) << "\' is not used in matching to PbPb. removing" << std::endl;
      listOfPPAlgo.erase(listOfPPAlgo.begin() + posFI);
      fileOfPPAlgo.erase(fileOfPPAlgo.begin() + posFI);
      isPPAlgoUsed.erase(isPPAlgoUsed.begin() + posFI);
    }
    else ++posFI;
  }

  inFileList.clear();
  inFileList.insert(inFileList.end(), fileOfPbPbAlgo.begin(), fileOfPbPbAlgo.end());

  for(unsigned int fI = 0; fI < fileOfPPAlgo.size(); ++fI){
    bool isUnique = true;
    for(unsigned int fI2 = 0; fI2 < fileOfPbPbAlgo.size(); ++fI2){
      if(isStrSame(fileOfPPAlgo.at(fI), fileOfPbPbAlgo.at(fI2))){
	isUnique = false;
	break;
      }
    }

    if(isUnique) inFileList.push_back(fileOfPPAlgo.at(fI));
  }

  std::cout << "Final files to process: " << std::endl;
  for(unsigned int fI = 0; fI < inFileList.size(); ++fI){
    std::cout << " " << fI << "/" << inFileList.size() << ": " << inFileList.at(fI) << std::endl;    
  }

  std::cout << "PbPb algos to process" << std::endl;
  for(unsigned int tI = 0; tI < listOfPbPbAlgo.size(); ++tI){
    std::cout << " " << tI << "/" << listOfPbPbAlgo.size() << ": " <<  listOfPbPbAlgo.at(tI) << std::endl;
  }

  std::cout << "PP algos to process" << std::endl;
  for(unsigned int tI = 0; tI < listOfPPAlgo.size(); ++tI){
    std::cout << " " << tI << "/" << listOfPPAlgo.size() << ": " <<  listOfPPAlgo.at(tI) << std::endl;
  }
  
  std::vector<std::string> slideTitles;
  std::vector<std::vector<std::string > > pdfPerSlide;
  std::vector<std::string> slideTitlesMain;
  std::vector<std::vector<std::string > > pdfPerSlideMain;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  std::string totRStr = "";
  for(unsigned int rI = 0; rI < listOfPbPbAlgo.size(); ++rI){
    totRStr = totRStr + "R" + std::to_string(getRVal(listOfPbPbAlgo.at(rI)));
  }

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  checkMakeDir("pdfDir/" + dateStr + "/PlotUnfolded/");

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);

  std::string outFileName = "output/" + dateStr + "/plotUnfoldedRAA_" + totRStr + "_" + dateStr + ".root";
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  outFile_p->cd();
  

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/plotUnfoldedRAA.exe <commaSeparatedInputs>" << std::endl;
    return 1;
  }

  std::string commaSeperatedInputs = std::string(argv[1]) + ",";
  while(commaSeperatedInputs.find(",,") != std::string::npos){commaSeperatedInputs.replace(commaSeperatedInputs.find(",,"), 2, ",");}

  std::vector<std::string> fileList;
  while(commaSeperatedInputs.find(",") != std::string::npos){
    fileList.push_back(commaSeperatedInputs.substr(0, commaSeperatedInputs.find(",")));
    commaSeperatedInputs.replace(0, commaSeperatedInputs.find(",")+1, "");
  }
  
  if(fileList.size() == 0){
    std::cout << "Input \'" << argv[1] << "\' produced no valid files. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotUnfoldedRAA(fileList);
  return retVal;
}
