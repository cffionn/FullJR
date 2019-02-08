//cpp depdencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TCollection.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TKey.h"
#include "TNamed.h"
#include "TTree.h"

//RooUnfold dependencies
#include "src/RooUnfoldResponse.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/doLocalDebug.h"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/stringUtil.h"

int recursiveCombineDir(TFile* inFile_p, const std::string dirName, TFile* outFile_p, const Int_t nDir, TDirectory* dirs_p[])
{
  int retVal = 0;

  inFile_p->cd();
  TDirectoryFile* dir_p = (TDirectoryFile*)inFile_p->Get(dirName.c_str());

  Int_t dirPos = -1;
  Int_t maxVal = -1;
  for(Int_t dI = 0; dI < nDir; ++dI){
    std::string tempName = dirs_p[dI]->GetName();

    if(tempName.find(dirName) != std::string::npos && maxVal < (Int_t )tempName.size()){
      dirPos = dI;
      maxVal = tempName.size();
    }
  }

  inFile_p->cd();
  dir_p->cd();
  TIter next(dir_p->GetListOfKeys());
  TKey* key = NULL;

  while( (key = (TKey*)next()) ){
    const std::string name = key->GetName();
    const std::string className = key->GetClassName();

    if(className.find("TDirectory") != std::string::npos){retVal += recursiveCombineDir(inFile_p, dirName + "/" + name, outFile_p, nDir, dirs_p);}

    if(className.find("TH1I") != std::string::npos){
      TH1I* hist_p = (TH1I*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("TH1F") != std::string::npos){
      TH1F* hist_p = (TH1F*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("TH1D") != std::string::npos){
      TH1D* hist_p = (TH1D*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("TH2I") != std::string::npos){
      TH2I* hist_p = (TH2I*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("TH2F") != std::string::npos){
      TH2F* hist_p = (TH2F*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("TH2D") != std::string::npos){
      TH2D* hist_p = (TH2D*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("RooUnfoldResponse") != std::string::npos){
      RooUnfoldResponse* rooRes_p = (RooUnfoldResponse*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      rooRes_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete rooRes_p;
    }
  }

  return retVal;
}

int combineFiles(const std::string outFileName, std::vector<std::string> fileList)
{
  cppWatch globalWatch;
  globalWatch.start();

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    if(!checkFile(fileList.at(fI)) || fileList.at(fI).find(".root") == std::string::npos){
      std::cout << "File \'" << fileList.at(fI) << "\' is invalid. return 1" << std::endl;
      return 1;
    }
  }

  cutPropagator prevProp;
  prevProp.Clean();

  std::vector<std::string> combinedHistTag;
  std::vector<int> combinedBestBayes;

  std::vector<std::string> tdirNameVect;
  std::vector<std::vector<std::string> > tdirNamePerFile;

  std::cout << "Combining " << fileList.size() << " files..." << std::endl;
  
  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << " " << fI << "/" << fileList.size() << ": " << fileList.at(fI) << std::endl;
    TFile* inFile_p = new TFile(fileList.at(fI).c_str(), "READ");

    cutPropagator cutPropCurrent;
    cutPropCurrent.Clean();    
    cutPropCurrent.GetAllVarFromFile(inFile_p);
    
    std::vector<std::string> tempHistTag = cutPropCurrent.GetHistTag();
    std::vector<int> tempHistBestBayes = cutPropCurrent.GetHistBestBayes();

    if(tempHistTag.size() != 0){
      combinedHistTag.insert(combinedHistTag.end(), tempHistTag.begin(), tempHistTag.end());
      combinedBestBayes.insert(combinedBestBayes.end(), tempHistBestBayes.begin(), tempHistBestBayes.end());
    }

    std::vector<std::string> classList;
    std::vector<std::string> contentsList = returnRootFileContentsList(inFile_p, "", "", -1, &(classList));

    tdirNamePerFile.push_back({});

    for(unsigned int i = 0; i < contentsList.size(); ++i){
      if(classList.at(i).find("TDirectory") == std::string::npos) continue;
      
      bool isSame = false;
      for(unsigned int j = 0; j < tdirNameVect.size(); ++j){
	if(isStrSame(contentsList.at(i), tdirNameVect.at(j))){
	  isSame = true;
	  break;
	}
      }

      if(!isSame) tdirNameVect.push_back(contentsList.at(i));

      tdirNamePerFile.at(tdirNamePerFile.size()-1).push_back(contentsList.at(i));
    }

    if(fI != 0){
      bool testProp = cutPropCurrent.CheckPropagatorsMatch(prevProp, true, true);
      if(!testProp){
	std::cout << "Warning: \'" << fileList.at(fI-1) << "\', and \'" << fileList.at(fI) << "\' doesn't match." << std::endl;
      }
    }
    else{
      prevProp.Clean();
      prevProp.GetAllVarFromFile(inFile_p);
    }

    inFile_p->Close();
    delete inFile_p;
  }

  unsigned int pos = 0;
  while(pos < tdirNameVect.size()){
    bool isSwap = false;
    for(unsigned int dI = pos+1; dI < tdirNameVect.size(); ++dI){
      if(tdirNameVect.at(pos).find("/") != std::string::npos && tdirNameVect.at(dI).find("/") == std::string::npos){
	std::string tempStr = tdirNameVect.at(pos);
        tdirNameVect.at(pos) = tdirNameVect.at(dI);
        tdirNameVect.at(dI) = tempStr;
        isSwap = true;
      }
    }

    for(unsigned int dI = pos+1; dI < tdirNameVect.size(); ++dI){
      if(tdirNameVect.at(dI).find("/") != std::string::npos) continue;

      if(getRVal(tdirNameVect.at(dI)) < getRVal(tdirNameVect.at(pos))){
	std::string tempStr = tdirNameVect.at(pos);
        tdirNameVect.at(pos) = tdirNameVect.at(dI);
        tdirNameVect.at(dI) = tempStr;
        isSwap = true;
      }
    }
    

    if(!isSwap) pos++;
  }

  pos = 0;
  while(pos < tdirNameVect.size()){
    bool isFound = false;
    for(unsigned int dI = pos+1; dI < tdirNameVect.size(); ++dI){
      if(tdirNameVect.at(pos).find(tdirNameVect.at(dI)) != std::string::npos){
	std::string tempStr = tdirNameVect.at(pos);
	tdirNameVect.at(pos) = tdirNameVect.at(dI);
	tdirNameVect.at(dI) = tempStr;
	isFound = true;
      }
    }

    if(!isFound) pos++;
  }

  std::cout << "Reproducing " << tdirNameVect.size() << " directories..." << std::endl;
  for(unsigned int i = 0; i < tdirNameVect.size(); ++i){
    std::cout << " " << i << "/" << tdirNameVect.size() << ": " << tdirNameVect.at(i) << std::endl;
  }

  const Int_t nDir = tdirNameVect.size();
  TDirectory* dirs_p[nDir];

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);

  
  int cutDirPos = -1;
  int subDirPos = -1;
  int unfoldDirPos = -1;
  
  for(Int_t dI = 0; dI < nDir; ++dI){
    outFile_p->cd();
    dirs_p[dI] = NULL;
    if(tdirNameVect.at(dI).find("/") == std::string::npos) dirs_p[dI] = (TDirectory*)outFile_p->mkdir(tdirNameVect.at(dI).c_str());
    else{
      Int_t prePos = -1;
      Int_t tempMaxVal = -1;
      for(Int_t dI2 = 0; dI2 < dI; ++dI2){
	if(tdirNameVect.at(dI).find(tdirNameVect.at(dI2)) != std::string::npos && tempMaxVal < (int)tdirNameVect.at(dI2).size()){
	  prePos = dI2;
	  tempMaxVal = tdirNameVect.at(dI2).size();
	}
      }
    
      std::string tempDirName = tdirNameVect.at(dI);
      tempDirName.replace(0, tempDirName.rfind("/")+1, "");
      std::cout << "Prepos: " << prePos << ", " << dirs_p[prePos]->GetName() << std::endl;
      dirs_p[prePos]->cd();
      dirs_p[dI] = (TDirectory*)dirs_p[prePos]->mkdir(tempDirName.c_str());
    }
    dirs_p[dI]->cd();

    if(isStrSame("cutDir", tdirNameVect.at(dI))) cutDirPos = dI;
    else if(tdirNameVect.at(dI).find("subDir") != std::string::npos) subDirPos = dI;
    else if(tdirNameVect.at(dI).find("unfoldDir") != std::string::npos) unfoldDirPos = dI;
  }

  if(unfoldDirPos >= 0){
    std::cout << "UnfoldDirPos: " << unfoldDirPos << std::endl;
    std::cout << " " << dirs_p[unfoldDirPos]->GetName() << std::endl;
  }
  else{
    std::cout << "No unfolding dir found warning" << std::endl;
  }

  outFile_p->cd();

  cppWatch writingWatch;
  writingWatch.start();

  int retVal = 0;

  std::cout << "Writing..." << std::endl;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    TFile* inFile_p = new TFile(fileList.at(fI).c_str(), "READ");

    TIter next(inFile_p->GetListOfKeys());
    TKey* key = NULL;
    
    while( (key = (TKey*)next() ) ){
      const std::string name = key->GetName();
      const std::string className = key->GetClassName();
      
      if(className.find("TDirectory") != std::string::npos){retVal += recursiveCombineDir(inFile_p, name, outFile_p, nDir, dirs_p);}
      
      if(className.find("TH1I") != std::string::npos){      
	TH1I* hist_p = (TH1I*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("TH1F") != std::string::npos){      
	TH1F* hist_p = (TH1F*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("TH1D") != std::string::npos){      
	TH1D* hist_p = (TH1D*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("TH2I") != std::string::npos){      
	TH2I* hist_p = (TH2I*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("TH2F") != std::string::npos){      
	TH2F* hist_p = (TH2F*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("TH2D") != std::string::npos){      
	TH2D* hist_p = (TH2D*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("RooUnfoldResponse") != std::string::npos){      
	RooUnfoldResponse* rooRes_p = (RooUnfoldResponse*)key->ReadObj();
      
	outFile_p->cd();
	rooRes_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete rooRes_p;
      }
    }

    inFile_p->Close();
    delete inFile_p;
  }

  if(retVal != 0) std::cout << "Uhoh rec nonzero something bad" << std::endl;

  writingWatch.stop();
  
  std::cout << "Write time: " << writingWatch.totalWall() << std::endl;

  prevProp.SetNHistDim(combinedHistTag.size());
  prevProp.SetHistTag(combinedHistTag);
  prevProp.SetHistBestBayes(combinedBestBayes);

  std::cout << " Best bayes: " << combinedHistTag.size() << ", " << combinedBestBayes.size() << std::endl;
  std::cout << "  unfoldDirPos: " << unfoldDirPos << std::endl;
  
  if(unfoldDirPos < 0){
    if(!prevProp.WriteAllVarToFile(outFile_p, dirs_p[cutDirPos], dirs_p[subDirPos], NULL)) std::cout << "Warning: Cut writing has failed" << std::endl;
  }
  else{
    if(!prevProp.WriteAllVarToFile(outFile_p, dirs_p[cutDirPos], dirs_p[subDirPos], dirs_p[unfoldDirPos])) std::cout << "Warning: Cut writing has failed" << std::endl;
  }

  outFile_p->Close();
  delete outFile_p;

  globalWatch.stop();

  std::cout << "Job complete." << std::endl;
  std::cout << "Total time: " << globalWatch.totalWall() << std::endl;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc < 3){
    std::cout << "Usage: ./bin/combineFiles.exe <outFileName> <Long> <List> <Of> <Input> <Files> ..." << std::endl;
    return 1;
  }

  std::vector<std::string> fileList;
  for(int i = 2; i < argc; ++i){
    fileList.push_back(argv[i]);
  }

  int retVal = 0;
  retVal += combineFiles(argv[1], fileList);
  return retVal;
}
