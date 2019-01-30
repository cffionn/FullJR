//cpp dependencies
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT dependencies
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"

//Non-Local dependencies
#include "MainAnalysis/include/cutPropagator.h"

template <class T>
bool compWithWarning(const std::string typeStr, T a, T b)
{
  bool val = a == b;
  if(!val) std::cout << typeStr << " \'" << a << "\' from one file does not match \'" << b << "\' from other. fail";
  return val;
}

ULong64_t getKey(ULong64_t jI, ULong64_t cI, ULong64_t idI, ULong64_t rI, ULong64_t aI, ULong64_t sI){return jI + 100*cI + 10000*idI + 1000000*rI + 100000000*aI + 10000000000*sI;}

ULong64_t getKey(ULong64_t jI, ULong64_t cI, ULong64_t idI, ULong64_t rI, ULong64_t aI, ULong64_t sI, ULong64_t bI){return jI + 100*cI + 10000*idI + 1000000*rI + 100000000*aI + 10000000000*sI + 10000000000000*bI;}

int posInStrVect(std::string strToCheck, std::string front, std::vector<std::string> vect, std::string back)
{
  int pos = -1;
  for(unsigned i = 0; i < vect.size(); ++i){
    if(strToCheck.find(front + vect[i] + back) != std::string::npos){
      pos = i;
      break;
    }
  }

  return pos;
}

int plotUnfoldedAll(const std::string inFileNamePP, const std::string inFileNamePbPb)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  const Int_t nFiles = 2;
  TFile* inFile_p[nFiles];
  std::vector<bool> isPbPb = {true, false};
  std::vector<std::string> fileNames = {inFileNamePbPb, inFileNamePP};
  std::vector<cutPropagator*> cutProps;

  std::vector<std::string> jtAlgosPbPb;
  std::vector<std::string> jtAlgosPP;
  std::vector<std::string> dirListPbPb;
  std::vector<std::string> dirListPP;

  std::vector<std::string> histTagsPbPb;
  std::vector<std::string> histTagsPP;
  std::vector<int> histBestBayesPbPb;
  std::vector<int> histBestBayesPP;
  std::map<std::string, int> histTagsToBestBayesPbPb;
  std::map<std::string, int> histTagsToBestBayesPP;
  std::map<ULong64_t, ULong64_t> histKeyToBestBayesKeyPbPb;
  std::map<ULong64_t, ULong64_t> histKeyToBestBayesKeyPP;

  Int_t nCentBins = -1;
  std::vector<Int_t> centBinsLow;
  std::vector<Int_t> centBinsHi;
  std::vector<std::string> centBinsStr;

  Int_t nID = -1;
  std::vector<std::string> idStr;
  std::vector<bool> goodID;

  Int_t nResponseMod = -1;
  std::vector<double> responseMod;
  std::vector<std::string> responseModStr;
  std::vector<bool> goodResponseMod;

  Int_t nJtAbsEtaBins = -1;  
  std::vector<double> jtAbsEtaBinsLow;
  std::vector<double> jtAbsEtaBinsHi;
  std::vector<std::string> jtAbsEtaBinsStr;
  std::vector<bool> goodJtAbsEtaBins;

  Int_t nSyst = -1;
  std::vector<std::string> systStr;

  Int_t nBayes = -1;
  std::vector<int> bayesVal;
  std::vector<std::string> bayesValStr;

  std::cout << "Processing inputs..." << std::endl;
  unsigned int pos = 0;
  for(auto const & fName : fileNames){
    inFile_p[pos] = new TFile(fName.c_str(), "READ");

    cutProps.push_back(new cutPropagator());
    cutProps[pos]->Clean();
    cutProps[pos]->GetAllVarFromFile(inFile_p[pos]);        
    
    if(isPbPb[pos]){
      dirListPbPb = returnRootFileContentsList(inFile_p[pos], "TDirectoryFile", "JetAnalyzer", 1);
      histTagsPbPb = cutProps[pos]->GetHistTag();
      histBestBayesPbPb = cutProps[pos]->GetHistBestBayes();

      nCentBins = cutProps[pos]->GetNCentBins();
      centBinsLow = cutProps[pos]->GetCentBinsLow();
      centBinsHi = cutProps[pos]->GetCentBinsHi();
      jtAlgosPbPb = cutProps[pos]->GetJtAlgos();
    }      
    else{
      dirListPP = returnRootFileContentsList(inFile_p[pos], "TDirectoryFile", "JetAnalyzer", 1);
      histTagsPP = cutProps[pos]->GetHistTag();
      histBestBayesPP = cutProps[pos]->GetHistBestBayes();

      jtAlgosPP = cutProps[pos]->GetJtAlgos();      
    }
    
    if(pos == 0){
      nID = cutProps[pos]->GetNID();
      nResponseMod = cutProps[pos]->GetNResponseMod();
      nSyst = cutProps[pos]->GetNSyst();
      nJtAbsEtaBins = cutProps[pos]->GetNJtAbsEtaBins();
      nBayes = cutProps[pos]->GetNBayes();

      idStr = cutProps[pos]->GetIdStr();
      responseMod = cutProps[pos]->GetResponseMod();
      jtAbsEtaBinsLow = cutProps[pos]->GetJtAbsEtaBinsLow();
      jtAbsEtaBinsHi = cutProps[pos]->GetJtAbsEtaBinsHi();
      systStr = cutProps[pos]->GetSystStr();
      bayesVal = cutProps[pos]->GetBayesVal();
    }
    else{
      if(!compWithWarning("nID", nID, cutProps[pos]->GetNID())) return 1;
      if(!compWithWarning("nResponseMod", nResponseMod, cutProps[pos]->GetNResponseMod())) return 1;
      if(!compWithWarning("nSyst", nSyst, cutProps[pos]->GetNSyst())) return 1;
      if(!compWithWarning("nJtAbsEtaBins", nJtAbsEtaBins, cutProps[pos]->GetNJtAbsEtaBins())) return 1;

      bool propMatch = cutProps[0]->CheckPropagatorsMatch(*(cutProps[pos]), true, false, false);

      if(!propMatch) std::cout << " Mismatch in propagator 0 to " << pos << std::endl;
      else std::cout << " Good propagator " << pos << std::endl;
    }

    inFile_p[pos]->Close();
    delete inFile_p[pos];
    ++pos;
  }
  
  for(unsigned int hI = 0; hI < histTagsPbPb.size(); ++hI){
    if(histTagsPbPb[hI].find("nHistDim") != std::string::npos) continue;
    histTagsToBestBayesPbPb[histTagsPbPb[hI]] = histBestBayesPbPb[hI];
  }
  for(unsigned int hI = 0; hI < histTagsPP.size(); ++hI){
    if(histTagsPP[hI].find("nHistDim") != std::string::npos) continue;
    histTagsToBestBayesPP[histTagsPP[hI]] = histBestBayesPP[hI];
  }

  for(auto const & res : responseMod){responseModStr.push_back("ResponseMod" + prettyString(res, 2, true));}

  for(unsigned int aI = 0; aI < jtAbsEtaBinsLow.size(); ++aI){
    jtAbsEtaBinsStr.push_back("AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true));
  }

  for(auto const & b : bayesVal){bayesValStr.push_back("Bayes" + std::to_string(b));}

  for(unsigned int cI = 0; cI < centBinsLow.size(); ++cI){
    centBinsStr.push_back("Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]));
  }

  std::cout << "nCentBins: " << nCentBins << std::endl;
  std::cout << "nBayes: " << nBayes << std::endl;

  std::cout << "PbPb Algos..." << std::endl;
  for(auto const & algo : jtAlgosPbPb){
    std::cout << " " << algo << std::endl;
  }
  std::cout << "PbPb Dir..." << std::endl;
  for(auto const & dir : dirListPbPb){
    std::cout << " " << dir << std::endl;
  }

  std::cout << "PP Algos..." << std::endl;
  for(auto const & algo : jtAlgosPP){
    std::cout << " " << algo << std::endl;
  }
  std::cout << "PP Dir..." << std::endl;
  for(auto const & dir : dirListPP){
    std::cout << " " << dir << std::endl;
  }

  const UInt_t nJtAlgos = jtAlgosPbPb.size();
  //nJtAlgos + nCentBins*100 + nMaxID*10000 + nMaxResponseMod*1000000 + nMaxJtAbsEtaBins*100000000 + nMaxSyst*10000000000

  std::map<ULong64_t, ULong64_t> keyToVectPosPbPb;
  std::map<ULong64_t, ULong64_t> keyPbPbToKeyPP;
  UInt_t nKeyPbPb = 0;
  std::map<ULong64_t, ULong64_t> keyToVectPosPP;
  UInt_t nKeyPP = 0;
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
      goodID.push_back(false);

      for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	goodResponseMod.push_back(false);

	for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	  goodJtAbsEtaBins.push_back(false);

	  for(ULong64_t sI = 0; sI < (ULong64_t)nSyst; ++sI){
	    for(ULong64_t bI = 0; bI < (ULong64_t)nBayes; ++bI){

	      ULong64_t keyPP = getKey(jI, 0, idI, rI, aI, sI, bI);
	      keyToVectPosPP[keyPP] = nKeyPP;
	      ++nKeyPP;

	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		ULong64_t keyPbPb = getKey(jI, cI, idI, rI, aI, sI, bI);
		keyToVectPosPbPb[keyPbPb] = nKeyPbPb;
		keyPbPbToKeyPP[keyPbPb] = keyPP;
		++nKeyPbPb;
	      }	      
	    }
	  }
	}
      }
    }
  }
  std::cout << "nKeys: " << nKeyPbPb << std::endl;
  std::vector<TH1D*> histVectPbPb;
  std::vector<TH1D*> histVectPP;

  histVectPbPb.reserve(nKeyPbPb);
  for(ULong64_t kI = 0; kI < nKeyPbPb; ++kI){histVectPbPb.push_back(NULL);}

  histVectPP.reserve(nKeyPP);
  for(ULong64_t kI = 0; kI < nKeyPP; ++kI){histVectPP.push_back(NULL);}

  std::cout << "Getting all hists...." << std::endl;
  for(Int_t fI = 0; fI < nFiles; ++fI){
    std::vector<std::string> dirs;
    if(isPbPb[fI]) dirs = dirListPbPb;
    else dirs = dirListPP;

    inFile_p[fI] = new TFile(fileNames[fI].c_str(), "READ");

    for(auto const & dir : dirs){
      TDirectory* dir_p = (TDirectory*)inFile_p[fI]->Get(dir.c_str());
      TIter next(dir_p->GetListOfKeys());
      TKey* key = NULL;
    
      int algoPos = -1;
      std::vector<std::string> algos;
      if(isPbPb[fI]) algos = jtAlgosPbPb;
      else algos = jtAlgosPP;

      for(UInt_t algoI = 0; algoI < algos.size(); ++algoI){
	if(algos[algoI].find(dir) != std::string::npos){
	  algoPos = algoI;
	  break;
	}
      }

      while((key = (TKey*)next())){
	const std::string name = key->GetName();
	const std::string className = key->GetClassName();
	if(className.find("TH1") == std::string::npos) continue;
	
	int idPos = posInStrVect(name, "_", idStr, "_");
	int modPos = posInStrVect(name, "_", responseModStr, "_");
	int absEtaPos = posInStrVect(name, "_", jtAbsEtaBinsStr, "_");
	int systPos = posInStrVect(name, "_", systStr, "_");
	if(systPos < 0 && name.find("FlatPrior") == std::string::npos) systPos = 0;
	int bayesPos = posInStrVect(name, "_", bayesValStr, "_");
	int centPos = 0;
	if(isPbPb[fI]) centPos = posInStrVect(name, "_", centBinsStr, "_");

	if(algoPos < 0) std::cout << "Missing algoPos in name \'" << name << "\'" << std::endl;
	if(idPos < 0) std::cout << "Missing idPos in name \'" << name << "\'" << std::endl;
	if(modPos < 0) std::cout << "Missing modPos in name \'" << name << "\'" << std::endl;
	if(absEtaPos < 0) std::cout << "Missing absEtaPos in name \'" << name << "\'" << std::endl;
	if(systPos < 0) std::cout << "Missing systPos in name \'" << name << "\'" << std::endl;
	if(bayesPos < 0) std::cout << "Missing bayesPos in name \'" << name << "\'" << std::endl;
	if(centPos < 0) std::cout << "Missing centPos in name \'" << name << "\'" << std::endl;	

	goodID[idPos] = true;
	goodResponseMod[modPos] = true;
	goodJtAbsEtaBins[absEtaPos] = true;

	ULong64_t idKey = getKey(algoPos, centPos, idPos, modPos, absEtaPos, systPos, bayesPos);
	if(isPbPb[fI]){
	  ULong64_t vectPos = keyToVectPosPbPb[idKey];	  
	  histVectPbPb[vectPos] = (TH1D*)key->ReadObj();
	}
	else{
	  ULong64_t vectPos = keyToVectPosPP[idKey];
	  histVectPP[vectPos] = (TH1D*)key->ReadObj();
	}
      }
    }
  }

  std::cout << "Processing best bayes PbPb..." << std::endl;
  for(auto const & tagToBB : histTagsToBestBayesPbPb){
    int algoPos = posInStrVect(tagToBB.first, "", dirListPbPb, "_");
    int idPos = posInStrVect(tagToBB.first, "_", idStr, "_");
    int modPos = posInStrVect(tagToBB.first, "_", responseModStr, "_");
    int absEtaPos = posInStrVect(tagToBB.first, "_", jtAbsEtaBinsStr, "_");
    int systPos = posInStrVect(tagToBB.first, "_", systStr, "_");
    if(systPos < 0 && tagToBB.first.find("FlatPrior") == std::string::npos) systPos = 0;
    int centPos = posInStrVect(tagToBB.first, "_", centBinsStr, "_");
    int bayesPos = tagToBB.second;
    
    if(algoPos < 0) std::cout << "Missing algoPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(idPos < 0) std::cout << "Missing idPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(modPos < 0) std::cout << "Missing modPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(absEtaPos < 0) std::cout << "Missing absEtaPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(systPos < 0) std::cout << "Missing systPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(centPos < 0) std::cout << "Missing centPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;	
    if(bayesPos < 0){
      std::cout << "Missing bayesPos in \'" << tagToBB.first << "\'. set to default 10" << std::endl;
      bayesPos = 10;
    }

    ULong64_t idKey = getKey(algoPos, centPos, idPos, modPos, absEtaPos, systPos);
    ULong64_t idKeyBB = getKey(algoPos, centPos, idPos, modPos, absEtaPos, systPos, (ULong64_t)bayesPos);
    histKeyToBestBayesKeyPbPb[idKey] = idKeyBB;
  }

  std::cout << "Processing best bayes PP..." << std::endl;
  for(auto const & tagToBB : histTagsToBestBayesPP){
    int algoPos = posInStrVect(tagToBB.first, "", dirListPP, "_");
    int idPos = posInStrVect(tagToBB.first, "_", idStr, "_");
    int modPos = posInStrVect(tagToBB.first, "_", responseModStr, "_");
    int absEtaPos = posInStrVect(tagToBB.first, "_", jtAbsEtaBinsStr, "_");
    int systPos = posInStrVect(tagToBB.first, "_", systStr, "_");
    if(systPos < 0 && tagToBB.first.find("FlatPrior") == std::string::npos) systPos = 0;
    int centPos = 0;
    int bayesPos = tagToBB.second;
    
    if(algoPos < 0) std::cout << "Missing algoPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(idPos < 0) std::cout << "Missing idPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(modPos < 0) std::cout << "Missing modPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(absEtaPos < 0) std::cout << "Missing absEtaPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(systPos < 0) std::cout << "Missing systPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(centPos < 0) std::cout << "Missing centPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;	
    if(bayesPos < 0){
      std::cout << "Missing bayesPos in \'" << tagToBB.first << "\'. set to default 10" << std::endl;
      bayesPos = 10;
    }

    ULong64_t idKey = getKey(algoPos, centPos, idPos, modPos, absEtaPos, systPos);
    ULong64_t idKeyBB = getKey(algoPos, centPos, idPos, modPos, absEtaPos, systPos, (ULong64_t)bayesPos);
    histKeyToBestBayesKeyPP[idKey] = idKeyBB;
  }
  
  const UInt_t nDummyStr = 256;
  char dummyStr[nDummyStr];
  std::cin.get(dummyStr, nDummyStr);

  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
      if(!goodID[idI]) continue;

      for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	if(!goodResponseMod[rI]) continue;

	for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	  if(!goodJtAbsEtaBins[aI]) continue;	 

	  std::cout << jtAlgosPbPb[jI] << ", " << jtAlgosPP[jI] << ", " << idStr[idI] << ", " << responseModStr[rI] << ", " << jtAbsEtaBinsStr[aI] << std::endl;

	  std::vector<TH1D*> histVectPbPbClones;
	  histVectPbPbClones.reserve(nSyst*nCentBins);
	  std::vector<TH1D*> histVectPPClones;	  
	  histVectPPClones.reserve(nSyst);

	  std::cout << __LINE__ << std::endl;
	  for(ULong64_t sI = 0; sI < (ULong64_t)nSyst; ++sI){
	    std::cout << __LINE__ << ", sI " << sI << std::endl;

	    ULong64_t idKey = getKey(jI, 0, idI, rI, aI, sI);
	    ULong64_t idKeyBB = histKeyToBestBayesKeyPP[idKey];
	    std::cout << __LINE__ << ", sI " << sI << std::endl;
	    histVectPPClones.push_back(NULL);
	    ULong64_t vectPos = keyToVectPosPP[idKeyBB];
	    std::string tempName = std::string(histVectPP[vectPos]->GetName()) + "_Clone";
	    std::cout << __LINE__ << ", sI " << sI << std::endl;
	    histVectPPClones[sI] = (TH1D*)histVectPP[vectPos]->Clone(tempName.c_str());
	    std::cout << __LINE__ << ", sI " << sI << std::endl;

	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      std::cout << __LINE__ << ", sI " << sI << ", cI " << cI << std::endl;
	      idKey = getKey(jI, cI, idI, rI, aI, sI);
	      idKeyBB = histKeyToBestBayesKeyPbPb[idKey];
	      histVectPbPbClones.push_back(NULL);
	      vectPos = keyToVectPosPbPb[idKeyBB];
	      tempName = std::string(histVectPbPb[vectPos]->GetName()) + "_Clone";
	      histVectPbPbClones[sI + cI*nSyst] = (TH1D*)histVectPbPb[vectPos]->Clone(tempName.c_str());
	    }
	  }


	  /*
	  for(ULong64_t sI = 0; sI < (ULong64_t)nSyst; ++sI){
	    delete histVectPPClones[sI];
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      delete histVectPbPbClones[sI + cI*nSyst];
	    }
	  }
	  histVectPPClones.clear();
	  histVectPbPbClones.clear();
	  */
	}
      }
    }
  }  

  for(Int_t fI = 0; fI < nFiles; ++fI){
    inFile_p[fI]->Close();
    delete inFile_p[fI];
  }

  for(unsigned int i = 0; i < cutProps.size(); ++i){
    delete (cutProps[i]);
  }

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/plotUnfoldedAll.exe <inFileNamePP> <inFileNamePbPb>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotUnfoldedAll(argv[1], argv[2]);
  return retVal;
}
