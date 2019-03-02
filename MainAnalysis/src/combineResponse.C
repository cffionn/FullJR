//cpp dependencies
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT dependencies
#include "TDatime.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TKey.h"
#include "TMath.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/stringUtil.h"

bool checkAllNames(std::vector<std::string> names1, std::vector<std::string> names2, std::string file1, std::string file2)
{
  for(auto const & name1: names1){
    bool name1Found = false;
    for(auto const & name2: names2){
      if(isStrSame(name1, name2)){
	name1Found = true;
	break;
      }
    }
    if(!name1Found){
      std::cout << "Hist \'" << name1 << "\' from file \'" << file1 << "\' is not found in file \'" << file2 << "\'" << std::endl;
      return false;
    }
  }
  return true;
}

void constructHist(std::vector<TH2D*>* hists_p, std::map<std::string, unsigned int>* histNameToPos, TFile* file_p, std::string dirName, bool doNew = false)
{
  TDirectory* dir_p = (TDirectory*)file_p->Get(dirName.c_str());
  TIter next(dir_p->GetListOfKeys());
  TKey* key=NULL;

  const Int_t nMaxBins = 100;
  Double_t binsX[nMaxBins];
  Double_t binsY[nMaxBins];
  

  while((key = (TKey*)next())){
    const std::string name = key->GetName();
    const std::string name2 = dirName + "/" + name;
    const std::string className = key->GetClassName();

    if(className.find("TH2D") == std::string::npos) continue;
    if(histNameToPos->count(name2) == 0) continue;

    unsigned int pos = (*histNameToPos)[name2];
    TH2D* tempHist_p = (TH2D*)file_p->Get((dirName + "/" + name).c_str());

    if(doNew){
      std::string newName = name;
      while(newName.find("/") != std::string::npos){newName.replace(0, newName.find("/")+1, "");}

      const Int_t nBinsX = tempHist_p->GetNbinsX();
      const Int_t nBinsY = tempHist_p->GetNbinsY();

      for(Int_t bIX = 0; bIX < nBinsX+1; ++bIX){
	binsX[bIX] = tempHist_p->GetXaxis()->GetBinLowEdge(bIX+1);
      }

      for(Int_t bIY = 0; bIY < nBinsY+1; ++bIY){
	binsY[bIY] = tempHist_p->GetYaxis()->GetBinLowEdge(bIY+1);
      }
      std::string title = ";" + std::string(tempHist_p->GetXaxis()->GetTitle()) + ";" + std::string(tempHist_p->GetYaxis()->GetTitle());
      (*hists_p)[pos] = new TH2D(newName.c_str(), title.c_str(), nBinsX, binsX, nBinsY, binsY);
    }

    for(Int_t bIX = 0; bIX < tempHist_p->GetNbinsX(); ++bIX){
      for(Int_t bIY = 0; bIY < tempHist_p->GetNbinsY(); ++bIY){
	Double_t val = (*hists_p)[pos]->GetBinContent(bIX+1, bIY+1);
	Double_t err = (*hists_p)[pos]->GetBinError(bIX+1, bIY+1);

	val += tempHist_p->GetBinContent(bIX+1, bIY+1);
	err = TMath::Sqrt(err*err + tempHist_p->GetBinError(bIX+1, bIY+1)*tempHist_p->GetBinError(bIX+1, bIY+1));

	(*hists_p)[pos]->SetBinContent(bIX+1, bIY+1, val);
	(*hists_p)[pos]->SetBinError(bIX+1, bIY+1, err);
      }
    } 
  }

  return;
}

int combineResponse(std::string inFileNames, std::string tagStr)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  const std::string permaFileNames = inFileNames;

  if(inFileNames.size() != 0){
    if(inFileNames.substr(inFileNames.size()-1, 1).find(",") == std::string::npos) inFileNames = inFileNames + ",";
  }

  std::vector<std::string> fileNames;
  while(inFileNames.find(",") != std::string::npos){
    std::string tempFileName = inFileNames.substr(0, inFileNames.find(","));
    if(!checkFile(tempFileName)) continue;
    if(tempFileName.find(".root") == std::string::npos) continue;
    fileNames.push_back(tempFileName);
    inFileNames.replace(0, inFileNames.find(",")+1, "");
  }

  if(fileNames.size() == 0){
    std::cout << "Given inputs give no valid files. \'" << permaFileNames << "\'. return 1." << std::endl;
    return 1;
  }

  std::string outFileName = "output/" + dateStr + "/combinedResponse_" + tagStr + "_" + dateStr + ".root";
  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);  

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  std::vector<TH2D*> hists_p;
  std::map<std::string, unsigned int> histNameToPos;

  ULong64_t totalNResponseTH2 = 0;
  TFile* inFile_p = new TFile(fileNames[0].c_str(), "READ");
  std::vector<std::string> responseNames = returnRootFileContentsList(inFile_p, "TH2D", "");
  
  inFile_p->Close();
  delete inFile_p;

  for(unsigned int fI = 1; fI < fileNames.size(); ++fI){
    inFile_p = new TFile(fileNames[fI].c_str(), "READ");

    std::vector<std::string> responseNames2 = returnRootFileContentsList(inFile_p, "TH2D", "");

    bool test1 = checkAllNames(responseNames, responseNames2, fileNames[0], fileNames[fI]);
    bool test2 = checkAllNames(responseNames2, responseNames, fileNames[fI], fileNames[0]);

    if(!test1 || !test2){
      std::cout << "Name in file \'" << fileNames[fI] << "\', not found in \'" << fileNames[0] << "\'. return 1" << std::endl;	
      inFile_p->Close();
      delete inFile_p;
      
      outFile_p->Close();
      delete outFile_p;
      
      return 1;
    }

    inFile_p->Close();
    delete inFile_p;
  }

  std::cout << "N responseNames: " << responseNames.size() << std::endl;
  for(unsigned int nI = 0; nI < responseNames.size(); ++nI){
    ++totalNResponseTH2;
    histNameToPos[responseNames[nI]] = nI;
  }
  hists_p.reserve(totalNResponseTH2);
  for(unsigned int nI = 0; nI < totalNResponseTH2; ++nI){
    hists_p.push_back(NULL);
  }

  for(unsigned int fI = 0; fI < fileNames.size(); ++fI){
    std::cout << "Processing file \'" << fileNames[fI] << "\'.." << std::endl;
    inFile_p = new TFile(fileNames[fI].c_str(), "READ");
    std::vector<std::string> dirNames = returnRootFileContentsList(inFile_p, "TDirectory", "");
    std::vector<std::string> dirNames2 = returnRootFileContentsList(inFile_p, "TDirectoryFile", "");
    dirNames.insert(std::end(dirNames), std::begin(dirNames2), std::end(dirNames2));

    for(auto const& dir : dirNames){
      std::cout << "Dir: " << dir << std::endl;
      if(fI == 0) constructHist(&hists_p, &histNameToPos, inFile_p, dir, true);
      else constructHist(&hists_p, &histNameToPos, inFile_p, dir, false);
    }
  }

  outFile_p->cd();
  std::cout << "HISTS: " << hists_p.size() << std::endl;

  for(unsigned int hI = 0; hI < hists_p.size(); ++hI){
    std::cout << "WRITING: " << hists_p[hI]->GetName() << std::endl;
    hists_p[hI]->Write("", TObject::kOverwrite);
    delete hists_p[hI];
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/combineResponse.exe <commaSeparatedInFileNames> <tagStr>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += combineResponse(argv[1], argv[2]);
  return retVal;
}
