#ifndef CUTPROPAGATOR_H
#define CUTPROPAGATOR_H

//cpp dependencies
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>

//ROOT dependencies
#include "TFile.h"
#include "TDirectory.h"
#include "TNamed.h"
#include "TCollection.h"
#include "TList.h"
#include "TKey.h"
#include "TString.h"

//Local FullJR (MainAnalysis) dependencies
#include "MainAnalysis/include/doLocalDebug.h"

//Non-local FullJR (Utility, etc.) dependencies
#include "Utility/include/doGlobalDebug.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/stringUtil.h"

class cutPropagator
{
 public:
  std::vector<std::string> inFileNames;
  std::vector<std::string> inFullFileNames;

  bool isPP;

  double jtAbsEtaMax;

  std::string rcDiffFileName;
  std::string flatPriorFileName;
  double jecVarMC;
  double jerVarMC;

  double jecVarData;

  double delta = 0.0001;

  int nResponseMod;
  std::vector<double> responseMod;
  std::vector<double> jerVarData;

  int nJtAlgos;
  std::vector<std::string> jtAlgos;
  std::vector<double> minJtPtCut;
  std::vector<double> multiJtPtCut;
  std::vector<int> recoTruncPos;

  int nJtPtBins;
  std::vector<double> jtPtBins;
  int nJtAbsEtaBins;
  std::vector<double> jtAbsEtaBinsLow;
  std::vector<double> jtAbsEtaBinsHi;
  
  int nPthats;
  std::vector<double> pthats;
  std::vector<double> pthatWeights;

  int nCentBins;
  std::vector<int> centBinsLow;
  std::vector<int> centBinsHi;

  int nID;
  std::vector<std::string> idStr;
  std::vector<double> jtPfCHMFCutLow;
  std::vector<double> jtPfCHMFCutHi;
  std::vector<double> jtPfMUMFCutLow;
  std::vector<double> jtPfMUMFCutHi;
  std::vector<double> jtPfNHFCutLow;
  std::vector<double> jtPfNHFCutHi;
  std::vector<double> jtPfNEFCutLow;
  std::vector<double> jtPfNEFCutHi;
  std::vector<double> jtPfMUFCutLow;
  std::vector<double> jtPfMUFCutHi;
  std::vector<double> jtPfCHFCutLow;
  std::vector<double> jtPfCHFCutHi;
  std::vector<double> jtPfCEFCutLow;
  std::vector<double> jtPfCEFCutHi;
  std::vector<int> jtPfMinMult;
  std::vector<int> jtPfMinChgMult;

  int nSyst;
  std::vector<std::string> systStr;

  int nBayes;
  int nBigBayesSymm;
  std::vector<int> bayesVal;

  int nSuperBayes;

  int nHistDim;
  std::vector<std::string> histTag;
  std::vector<int> histBestBayes;

  void Clean();
  bool GetAllVarFromFile(TFile* inFile_p);
  bool WriteAllVarToFile(TFile* inFile_p, TDirectory* inDir_p, TDirectory* inSubDir_p, TDirectory* unfoldDir_p);
  bool CheckPropagatorsMatch(cutPropagator inCutProp, bool doBothMCOrBothData, bool doBothPPOrBothPbPb, bool skipJtAlgo);

  bool CheckDouble(double in1, double in2);
  bool CheckInt(int in1, int in2);
  bool CheckBool(bool in1, bool in2);
  bool CheckString(std::string in1, std::string in2);
  bool CheckVectDouble(std::vector<double> in1, std::vector<double> in2);
  bool CheckVectInt(std::vector<int> in1, std::vector<int> in2);
  bool CheckVectString(std::vector<std::string> in1, std::vector<std::string> in2);

  bool CheckJtAbsEtaMax(double inJtAbsEtaMax);
  bool CheckJtAbsEtaMax(cutPropagator inCutProp);
  bool CheckNJtPtBins(int inNJtPtBins);
  bool CheckNJtPtBins(cutPropagator inCutProp);
  bool CheckNJtAbsEtaBins(int inNJtAbsEtaBins);
  bool CheckNJtAbsEtaBins(cutPropagator inCutProp);
  bool CheckNID(int inNID);
  bool CheckNID(cutPropagator inCutProp);
  bool CheckJtPtBins(std::vector<double> inJtPtBins);
  bool CheckJtPtBins(cutPropagator inCutProp);
  bool CheckJtAbsEtaBinsLow(std::vector<double> inJtAbsEtaBinsLow);
  bool CheckJtAbsEtaBinsLow(cutPropagator inCutProp);
  bool CheckJtAbsEtaBinsHi(std::vector<double> inJtAbsEtaBinsHi);
  bool CheckJtAbsEtaBinsHi(cutPropagator inCutProp);
  bool CheckIdStr(std::vector<std::string> inIDStr);
  bool CheckIdStr(cutPropagator inCutProp);
  bool CheckJtPfCHMFCutLow(std::vector<double> inJtPfCHMFCutLow);
  bool CheckJtPfCHMFCutLow(cutPropagator inCutProp);
  bool CheckJtPfCHMFCutHi(std::vector<double> inJtPfCHMFCutHi);
  bool CheckJtPfCHMFCutHi(cutPropagator inCutProp);
  bool CheckJtPfMUMFCutLow(std::vector<double> inJtPfMUMFCutLow);
  bool CheckJtPfMUMFCutLow(cutPropagator inCutProp);
  bool CheckJtPfMUMFCutHi(std::vector<double> inJtPfMUMFCutHi);
  bool CheckJtPfMUMFCutHi(cutPropagator inCutProp);
  bool CheckJtPfNHFCutLow(std::vector<double> inJtPfNHFCutLow);
  bool CheckJtPfNHFCutLow(cutPropagator inCutProp);
  bool CheckJtPfNHFCutHi(std::vector<double> inJtPfNHFCutHi);
  bool CheckJtPfNHFCutHi(cutPropagator inCutProp);
  bool CheckJtPfNEFCutLow(std::vector<double> inJtPfNEFCutLow);
  bool CheckJtPfNEFCutLow(cutPropagator inCutProp);
  bool CheckJtPfNEFCutHi(std::vector<double> inJtPfNEFCutHi);
  bool CheckJtPfNEFCutHi(cutPropagator inCutProp);
  bool CheckJtPfMUFCutLow(std::vector<double> inJtPfMUFCutLow);
  bool CheckJtPfMUFCutLow(cutPropagator inCutProp);
  bool CheckJtPfMUFCutHi(std::vector<double> inJtPfMUFCutHi);
  bool CheckJtPfMUFCutHi(cutPropagator inCutProp);
  bool CheckJtPfCHFCutLow(std::vector<double> inJtPfCHFCutLow);
  bool CheckJtPfCHFCutLow(cutPropagator inCutProp);
  bool CheckJtPfCHFCutHi(std::vector<double> inJtPfCHFCutHi);
  bool CheckJtPfCHFCutHi(cutPropagator inCutProp);
  bool CheckJtPfCEFCutLow(std::vector<double> inJtPfCEFCutLow);
  bool CheckJtPfCEFCutLow(cutPropagator inCutProp);
  bool CheckJtPfCEFCutHi(std::vector<double> inJtPfCEFCutHi);
  bool CheckJtPfCEFCutHi(cutPropagator inCutProp);
  bool CheckJtPfMinMult(std::vector<int> inJtPfMinMult);
  bool CheckJtPfMinMult(cutPropagator inCutProp);
  bool CheckJtPfMinChgMult(std::vector<int> inJtPfMinChgMult);
  bool CheckJtPfMinChgMult(cutPropagator inCutProp);
  bool CheckNSyst(int inNSyst);
  bool CheckNSyst(cutPropagator inCutProp);
  bool CheckSystStr(std::vector<std::string> inSystStr);
  bool CheckSystStr(cutPropagator inCutProp);
  bool CheckNBayes(int inNBayes);
  bool CheckNBayes(cutPropagator inCutProp);
  bool CheckNBigBayesSymm(int inNBigBayesSymm);
  bool CheckNBigBayesSymm(cutPropagator inCutProp);
  bool CheckBayesVal(std::vector<int> inBayesVal);
  bool CheckBayesVal(cutPropagator inCutProp);
  bool CheckNSuperBayes(int inNSuperBayes);
  bool CheckNSuperBayes(cutPropagator inCutProp);
  bool CheckNHistDim(int inNHistDim);
  bool CheckNHistDim(cutPropagator inCutProp);
  bool CheckHistTag(std::vector<int> inHistTag);
  bool CheckHistTag(cutPropagator inCutProp);
  bool CheckHistBestBayes(std::vector<int> inHistBestBayes);
  bool CheckHistBestBayes(cutPropagator inCutProp);
  bool CheckRCDiffFileName(std::string inRCDiffFileName);
  bool CheckRCDiffFileName(cutPropagator inCutProp);
  bool CheckFlatPriorFileName(std::string inFlatPriorFileName);
  bool CheckFlatPriorFileName(cutPropagator inCutProp);
  bool CheckJECVarMC(double inJECVarMC);
  bool CheckJECVarMC(cutPropagator inCutProp);
  bool CheckJERVarMC(double inJERVarMC);
  bool CheckJERVarMC(cutPropagator inCutProp);
  bool CheckJECVarData(double inJECVarData);
  bool CheckJECVarData(cutPropagator inCutProp);
  bool CheckNResponseMod(int inNResponseMod);
  bool CheckNResponseMod(cutPropagator inCutProp);
  bool CheckResponseMod(std::vector<double> inResponseMod);
  bool CheckResponseMod(cutPropagator inCutProp);
  bool CheckJERVarData(std::vector<double> inJERVarData);
  bool CheckJERVarData(cutPropagator inCutProp);
  bool CheckNPthats(int inNPthats);
  bool CheckNPthats(cutPropagator inCutProp);
  bool CheckPthats(std::vector<double> inPthats);
  bool CheckPthats(cutPropagator inCutProp);
  bool CheckPthatWeights(std::vector<double> inPthatWeights);
  bool CheckPthatWeights(cutPropagator inCutProp);
  bool CheckNJtAlgos(int inNJtAlgos);
  bool CheckNJtAlgos(cutPropagator inCutProp);
  bool CheckJtAlgos(std::vector<std::string> inJtAlgos);
  bool CheckJtAlgos(cutPropagator inCutProp);
  bool CheckMinJtPtCut(std::vector<double> inMinJtPtCut);
  bool CheckMinJtPtCut(cutPropagator inCutProp);
  bool CheckMultiJtPtCut(std::vector<double> inMultiJtPtCut);
  bool CheckMultiJtPtCut(cutPropagator inCutProp);
  bool CheckRecoTruncPos(std::vector<int> inRecoTruncPos);
  bool CheckRecoTruncPos(cutPropagator inCutProp);
  bool CheckIsPP(bool inIsPP);
  bool CheckIsPP(cutPropagator inCutProp);
  bool CheckNCentBins(int inNCentBins);
  bool CheckNCentBins(cutPropagator inCutProp);
  bool CheckCentBinsLow(std::vector<int> inCentBinsLow);
  bool CheckCentBinsLow(cutPropagator inCutProp);
  bool CheckCentBinsHi(std::vector<int> inCentBinsHi);
  bool CheckCentBinsHi(cutPropagator inCutProp);

  std::vector<std::string> GetInFileNames(){return inFileNames;}
  std::vector<std::string> GetInFullFileNames(){return inFullFileNames;}

  bool GetIsPP(){return isPP;}

  double GetJtAbsEtaMax(){return jtAbsEtaMax;}

  std::string GetRCDiffFileName(){return rcDiffFileName;}
  std::string GetFlatPriorFileName(){return flatPriorFileName;}

  double GetJECVarMC(){return jecVarMC;}
  double GetJERVarMC(){return jerVarMC;}

  double GetJECVarData(){return jecVarData;}

  int GetNResponseMod(){return nResponseMod;}
  std::vector<double> GetResponseMod(){return responseMod;}
  std::vector<double> GetJERVarData(){return jerVarData;}

  int GetNJtAlgos(){return nJtAlgos;}
  std::vector<std::string> GetJtAlgos(){return jtAlgos;}
  std::vector<double> GetMinJtPtCut(){return minJtPtCut;}
  std::vector<double> GetMultiJtPtCut(){return multiJtPtCut;}
  std::vector<int> GetRecoTruncPos(){return recoTruncPos;}


  int GetNJtPtBins(){return nJtPtBins;}
  std::vector<double> GetJtPtBins(){return jtPtBins;}
  int GetNJtAbsEtaBins(){return nJtAbsEtaBins;}
  std::vector<double> GetJtAbsEtaBinsLow(){return jtAbsEtaBinsLow;}
  std::vector<double> GetJtAbsEtaBinsHi(){return jtAbsEtaBinsHi;}

  int GetNPthats(){return nPthats;}
  std::vector<double> GetPthats(){return pthats;}
  std::vector<double> GetPthatWeights(){return pthatWeights;}

  int GetNCentBins(){return nCentBins;}
  std::vector<int> GetCentBinsLow(){return centBinsLow;}
  std::vector<int> GetCentBinsHi(){return centBinsHi;}

  int GetNID(){return nID;}
  std::vector<std::string> GetIdStr(){return idStr;}
  std::vector<double> GetJtPfCHMFCutLow(){return jtPfCHMFCutLow;}
  std::vector<double> GetJtPfCHMFCutHi(){return jtPfCHMFCutHi;}
  std::vector<double> GetJtPfMUMFCutLow(){return jtPfMUMFCutLow;}
  std::vector<double> GetJtPfMUMFCutHi(){return jtPfMUMFCutHi;}
  std::vector<double> GetJtPfNHFCutLow(){return jtPfNHFCutLow;}
  std::vector<double> GetJtPfNHFCutHi(){return jtPfNHFCutHi;}
  std::vector<double> GetJtPfNEFCutLow(){return jtPfNEFCutLow;}
  std::vector<double> GetJtPfNEFCutHi(){return jtPfNEFCutHi;}
  std::vector<double> GetJtPfMUFCutLow(){return jtPfMUFCutLow;}
  std::vector<double> GetJtPfMUFCutHi(){return jtPfMUFCutHi;}
  std::vector<double> GetJtPfCHFCutLow(){return jtPfCHFCutLow;}
  std::vector<double> GetJtPfCHFCutHi(){return jtPfCHFCutHi;}
  std::vector<double> GetJtPfCEFCutLow(){return jtPfCEFCutLow;}
  std::vector<double> GetJtPfCEFCutHi(){return jtPfCEFCutHi;}
  std::vector<int> GetJtPfMinMult(){return jtPfMinMult;}
  std::vector<int> GetJtPfMinChgMult(){return jtPfMinChgMult;}

  int GetNSyst(){return nSyst;}
  std::vector<std::string> GetSystStr(){return systStr;}

  int GetNBayes(){return nBayes;}
  int GetNBigBayesSymm(){return nBigBayesSymm;}
  std::vector<int> GetBayesVal(){return bayesVal;}

  int GetNSuperBayes(){return nSuperBayes;}

  int GetNHistDim(){return nHistDim;}
  std::vector<std::string> GetHistTag(){return histTag;}
  std::vector<int> GetHistBestBayes(){return histBestBayes;}

  void SetInFileNames(std::vector<std::string> inInFileNames){inFileNames = inInFileNames; return;}
  void SetInFullFileNames(std::vector<std::string> inInFullFileNames){inFullFileNames = inInFullFileNames; return;}

  void SetIsPP(bool inIsPP){isPP = inIsPP; return;}

  void SetJtAbsEtaMax(double inJtAbsEtaMax){jtAbsEtaMax = inJtAbsEtaMax; return;}

  void SetRCDiffFileName(std::string inRCDiffFileName){rcDiffFileName = inRCDiffFileName; return;}
  void SetFlatPriorFileName(std::string inFlatPriorFileName){flatPriorFileName = inFlatPriorFileName; return;}
  
  void SetJECVarMC(double inJECVarMC){jecVarMC = inJECVarMC; return;}
  void SetJERVarMC(double inJERVarMC){jerVarMC = inJERVarMC; return;}

  void SetJECVarData(double inJECVarData){jecVarData = inJECVarData; return;}

  void SetNResponseMod(int inNResponseMod){nResponseMod = inNResponseMod; return;}
  void SetResponseMod(std::vector<double> inResponseMod){responseMod = inResponseMod; return;}
  void SetResponseMod(int inN, const Double_t inResponseMod[]);
  void SetJERVarData(std::vector<double> inJERVarData){jerVarData = inJERVarData; return;}
  void SetJERVarData(int inN, const Double_t inJERVarData[]);

  void SetNJtAlgos(int inNJtAlgos){nJtAlgos = inNJtAlgos; return;}
  void SetJtAlgos(std::vector<std::string> inJtAlgos){jtAlgos = inJtAlgos; return;}
  void SetJtAlgos(int inN, std::string inJtAlgos[]);
  void SetMinJtPtCut(std::vector<double> inMinJtPtCut){minJtPtCut = inMinJtPtCut; return;}
  void SetMinJtPtCut(int inN, double inMinJtPtCut[]);
  void SetMultiJtPtCut(std::vector<double> inMultiJtPtCut){multiJtPtCut = inMultiJtPtCut; return;}
  void SetMultiJtPtCut(int inN, double inMultiJtPtCut[]);
  void SetRecoTruncPos(std::vector<int> inRecoTruncPos){recoTruncPos = inRecoTruncPos; return;}
  void SetRecoTruncPos(int inN, int inRecoTruncPos[]);

  void SetNJtPtBins(int inNJtPtBins){nJtPtBins = inNJtPtBins; return;}
  void SetJtPtBins(std::vector<double> inJtPtBins){jtPtBins = inJtPtBins; return;}
  void SetJtPtBins(int inN, const Double_t inJtPtBins[]);
  void SetNJtAbsEtaBins(int inNJtAbsEtaBins){nJtAbsEtaBins = inNJtAbsEtaBins; return;}
  void SetJtAbsEtaBinsLow(std::vector<double> inJtAbsEtaBinsLow){jtAbsEtaBinsLow = inJtAbsEtaBinsLow; return;}
  void SetJtAbsEtaBinsLow(int inN, const Double_t inJtAbsEtaBinsLow[]);  
  void SetJtAbsEtaBinsHi(std::vector<double> inJtAbsEtaBinsHi){jtAbsEtaBinsHi = inJtAbsEtaBinsHi; return;}
  void SetJtAbsEtaBinsHi(int inN, const Double_t inJtAbsEtaBinsHi[]);  

  void SetNPthats(int inNPthats){nPthats = inNPthats; return;}
  void SetPthats(std::vector<double> inPthats){pthats = inPthats; return;}
  void SetPthats(int inN, const Double_t inPthats[]);
  void SetPthatWeights(std::vector<double> inPthatWeights){pthatWeights = inPthatWeights; return;}
  void SetPthatWeights(int inN, const Double_t inPthatWeights[]);

  void SetNCentBins(int inNCentBins){nCentBins = inNCentBins; return;}
  void SetCentBinsLow(std::vector<int> inCentBinsLow){centBinsLow = inCentBinsLow; return;}
  void SetCentBinsLow(int inN, const Int_t inCentBinsLow[]);
  void SetCentBinsHi(std::vector<int> inCentBinsHi){centBinsHi = inCentBinsHi; return;}
  void SetCentBinsHi(int inN, const Int_t inCentBinsHi[]);

  void SetNID(int inNID){nID = inNID; return;}
  void SetIdStr(std::vector<std::string> inIdStr){idStr = inIdStr; return;};
  void SetIdStr(int inN, const std::string inIdStr[]);
  void SetJtPfCHMFCutLow(std::vector<double> inJtPfCHMFCutLow){jtPfCHMFCutLow = inJtPfCHMFCutLow; return;}
  void SetJtPfCHMFCutLow(int inN, const Double_t inJtPfCHMFCutLow[]);
  void SetJtPfCHMFCutHi(std::vector<double> inJtPfCHMFCutHi){jtPfCHMFCutHi = inJtPfCHMFCutHi; return;}
  void SetJtPfCHMFCutHi(int inN, const Double_t inJtPfCHMFCutHi[]);
  void SetJtPfMUMFCutLow(std::vector<double> inJtPfMUMFCutLow){jtPfMUMFCutLow = inJtPfMUMFCutLow; return;}
  void SetJtPfMUMFCutLow(int inN, const Double_t inJtPfMUMFCutLow[]);
  void SetJtPfMUMFCutHi(std::vector<double> inJtPfMUMFCutHi){jtPfMUMFCutHi = inJtPfMUMFCutHi; return;}
  void SetJtPfMUMFCutHi(int inN, const Double_t inJtPfMUMFCutHi[]);
  void SetJtPfNHFCutLow(std::vector<double> inJtPfNHFCutLow){jtPfNHFCutLow = inJtPfNHFCutLow; return;}
  void SetJtPfNHFCutLow(int inN, const Double_t inJtPfNHFCutLow[]);
  void SetJtPfNHFCutHi(std::vector<double> inJtPfNHFCutHi){jtPfNHFCutHi = inJtPfNHFCutHi; return;}
  void SetJtPfNHFCutHi(int inN, const Double_t inJtPfNHFCutHi[]);
  void SetJtPfNEFCutLow(std::vector<double> inJtPfNEFCutLow){jtPfNEFCutLow = inJtPfNEFCutLow; return;}
  void SetJtPfNEFCutLow(int inN, const Double_t inJtPfNEFCutLow[]);
  void SetJtPfNEFCutHi(std::vector<double> inJtPfNEFCutHi){jtPfNEFCutHi = inJtPfNEFCutHi; return;}
  void SetJtPfNEFCutHi(int inN, const Double_t inJtPfNEFCutHi[]);
  void SetJtPfMUFCutLow(std::vector<double> inJtPfMUFCutLow){jtPfMUFCutLow = inJtPfMUFCutLow; return;}
  void SetJtPfMUFCutLow(int inN, const Double_t inJtPfMUFCutLow[]);
  void SetJtPfMUFCutHi(std::vector<double> inJtPfMUFCutHi){jtPfMUFCutHi = inJtPfMUFCutHi; return;}
  void SetJtPfMUFCutHi(int inN, const Double_t inJtPfMUFCutHi[]);
  void SetJtPfCHFCutLow(std::vector<double> inJtPfCHFCutLow){jtPfCHFCutLow = inJtPfCHFCutLow; return;}
  void SetJtPfCHFCutLow(int inN, const Double_t inJtPfCHFCutLow[]);
  void SetJtPfCHFCutHi(std::vector<double> inJtPfCHFCutHi){jtPfCHFCutHi = inJtPfCHFCutHi; return;}
  void SetJtPfCHFCutHi(int inN, const Double_t inJtPfCHFCutHi[]);
  void SetJtPfCEFCutLow(std::vector<double> inJtPfCEFCutLow){jtPfCEFCutLow = inJtPfCEFCutLow; return;}
  void SetJtPfCEFCutLow(int inN, const Double_t inJtPfCEFCutLow[]);
  void SetJtPfCEFCutHi(std::vector<double> inJtPfCEFCutHi){jtPfCEFCutHi = inJtPfCEFCutHi; return;}
  void SetJtPfCEFCutHi(int inN, const Double_t inJtPfCEFCutHi[]);
  void SetJtPfMinMult(std::vector<int> inJtPfMinMult){jtPfMinMult = inJtPfMinMult; return;}
  void SetJtPfMinMult(int inN, const Int_t inJtPfMinMult[]);
  void SetJtPfMinChgMult(std::vector<int> inJtPfMinChgMult){jtPfMinChgMult = inJtPfMinChgMult; return;}
  void SetJtPfMinChgMult(int inN, const Int_t inJtPfMinChgMult[]);

  void SetNSyst(int inNSyst){nSyst = inNSyst; return;}
  void SetSystStr(std::vector<std::string> inSystStr){systStr = inSystStr; return;};
  void SetSystStr(int inN, const std::string inSystStr[]);

  void SetNBayes(int inNBayes){nBayes = inNBayes; return;}
  void SetNBigBayesSymm(int inNBigBayesSymm){nBigBayesSymm = inNBigBayesSymm; return;}
  void SetBayesVal(std::vector<int> inBayesVal){bayesVal = inBayesVal; return;};
  void SetBayesVal(int inN, const int inBayesVal[]);
  void SetNSuperBayes(int inNSuperBayes){nSuperBayes = inNSuperBayes; return;}

  void SetNHistDim(int inNHistDim){nHistDim = inNHistDim; return;}
  void SetHistTag(std::vector<std::string> inHistTag){histTag = inHistTag; return;};
  void SetHistTag(int inN, const std::string inHistTag[]);
  void SetHistBestBayes(std::vector<int> inHistBestBayes){histBestBayes = inHistBestBayes; return;};
  void SetHistBestBayes(int inN, const int inHistBestBayes[]);

  std::vector<int> StringToIntVect(std::string inStr);
  std::vector<double> StringToDoubleVect(std::string inStr);
  std::vector<std::string> StringToStringVect(std::string inStr);

  std::string to_string_with_precision(double a_value, const int n);
};


void cutPropagator::Clean()
{
  inFileNames.clear();
  inFullFileNames.clear();

  isPP = false;

  jtAbsEtaMax = -99;

  rcDiffFileName = "";
  flatPriorFileName = "";

  jecVarMC = -99;
  jerVarMC = -99;

  jecVarData = -99;

  nResponseMod = -1;
  responseMod.clear();
  jerVarData.clear();

  nJtAlgos = -1;
  jtAlgos.clear();
  minJtPtCut.clear();
  multiJtPtCut.clear();
  recoTruncPos.clear();

  nJtPtBins = -1;
  jtPtBins.clear();
  nJtAbsEtaBins = -1;
  jtAbsEtaBinsLow.clear();
  jtAbsEtaBinsHi.clear();

  nPthats = -1;
  pthats.clear();
  pthatWeights.clear();

  nCentBins = -1;
  centBinsLow.clear();
  centBinsHi.clear();

  nID = -1;
  idStr.clear();
  jtPfCHMFCutLow.clear();
  jtPfCHMFCutHi.clear();
  jtPfMUMFCutLow.clear();
  jtPfMUMFCutHi.clear();
  jtPfNHFCutLow.clear();
  jtPfNHFCutHi.clear();
  jtPfNEFCutLow.clear();
  jtPfNEFCutHi.clear();
  jtPfMUFCutLow.clear();
  jtPfMUFCutHi.clear();
  jtPfCHFCutLow.clear();
  jtPfCHFCutHi.clear();
  jtPfCEFCutLow.clear();
  jtPfCEFCutHi.clear();
  jtPfMinMult.clear();
  jtPfMinChgMult.clear();

  nSyst = -1;
  systStr.clear();

  nBayes = -1;
  nBigBayesSymm = -1;
  bayesVal.clear();

  nSuperBayes = -1;

  nHistDim = -1;
  histTag.clear();
  histBestBayes.clear();

  return;
}

bool cutPropagator::GetAllVarFromFile(TFile* inFile_p)
{
  inFile_p->cd();

  std::vector<std::string> cutDirCheck = returnRootFileContentsList(inFile_p, "TDirectoryFile", "cutDir", 1);
  if(cutDirCheck.size() == 0){
    std::cout << "cutPropagator::GetAllVarFromFile - given inFile_p \'" << inFile_p->GetName() << "\' contains no cutDir. return 1" << std::endl;
    return false;
  }

  std::vector<std::string> cutDirTNameds = returnTDirContentsList(inFile_p, "cutDir", "TNamed", "cutDir", 0, -1);
  std::vector<std::string> cutDirTDir = returnTDirContentsList(inFile_p, "cutDir", "TDirectoryFile", "", 0, -1);
  //returnRootFileContentsList(inFile_p, "TNamed", "cutDir");
  unsigned int pos = 0;
  while(pos < cutDirTNameds.size()){
    if(cutDirTNameds.at(pos).find("cutDir/") != std::string::npos) ++pos;
    else cutDirTNameds.erase(cutDirTNameds.begin()+pos);
  }

  if(cutDirTNameds.size() == 0){
    std::cout << "cutPropagator::GetAllVarFromFile - given inFile_p \'" << inFile_p->GetName() << "\' contains empty cutDir. return 1" << std::endl;
    return false;
  }

  //  std::cout << "Printing " << cutDirTNameds.size() << " TNameds..." << std::endl;
  for(unsigned int cI = 0; cI < cutDirTNameds.size(); ++cI){
    std::string tempStr = cutDirTNameds.at(cI);
    while(tempStr.find("/") != std::string::npos){tempStr.replace(0, tempStr.find("/")+1, "");}

    //    std::cout << " " << cI << "/" << cutDirTNameds.size() << ": " << cutDirTNameds.at(cI) << std::endl;
    bool isHistBayesBest = cutDirTNameds.at(cI).find("/unfoldDir/") != std::string::npos;

    if(isHistBayesBest){
      
      histTag.push_back(tempStr);
      histBestBayes.push_back(std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle()));
      
      continue;
    }


    if(isStrSame("nCentBins", tempStr)) nCentBins = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("centBinsLow", tempStr)) centBinsLow = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("centBinsHi", tempStr)) centBinsHi = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtAbsEtaMax", tempStr)) jtAbsEtaMax = std::stof(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("rcDiffFileName", tempStr)) rcDiffFileName = ((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle();
    else if(isStrSame("flatPriorFileName", tempStr)) flatPriorFileName = ((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle();
    else if(isStrSame("jecVarMC", tempStr)) jecVarMC = std::stof(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jerVarMC", tempStr)) jerVarMC = std::stof(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jecVarData", tempStr)) jecVarData = std::stof(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nResponseMod", tempStr)) nResponseMod = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nJtAlgos", tempStr)) nJtAlgos = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtAlgos", tempStr)) jtAlgos = StringToStringVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("minJtPtCut", tempStr)) minJtPtCut = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("multiJtPtCut", tempStr)) multiJtPtCut = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoTruncPos", tempStr)) recoTruncPos = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nJtPtBins", tempStr)) nJtPtBins = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nJtAbsEtaBins", tempStr)) nJtAbsEtaBins = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("inFileNames", tempStr)) inFileNames = StringToStringVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("inFullFileNames_") != std::string::npos){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle();
      inFullFileNames.push_back(tempStr2);
    }
    else if(isStrSame("isPP", tempStr)) isPP = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nPthats", tempStr)) nPthats = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("pthats", tempStr)) pthats = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("pthatWeights", tempStr)) pthatWeights = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("responseMod", tempStr)) responseMod = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jerVarData", tempStr)) jerVarData = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPtBins", tempStr)) jtPtBins = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtAbsEtaBinsLow", tempStr)) jtAbsEtaBinsLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtAbsEtaBinsHi", tempStr)) jtAbsEtaBinsHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nID", tempStr)) nID = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("idStr", tempStr)) idStr = StringToStringVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCHMFCutLow", tempStr)) jtPfCHMFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCHMFCutHi", tempStr)) jtPfCHMFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMUMFCutLow", tempStr)) jtPfMUMFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMUMFCutHi", tempStr)) jtPfMUMFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfNHFCutLow", tempStr)) jtPfNHFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfNHFCutHi", tempStr)) jtPfNHFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfNEFCutLow", tempStr)) jtPfNEFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfNEFCutHi", tempStr)) jtPfNEFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMUFCutLow", tempStr)) jtPfMUFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMUFCutHi", tempStr)) jtPfMUFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCHFCutLow", tempStr)) jtPfCHFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCHFCutHi", tempStr)) jtPfCHFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCEFCutLow", tempStr)) jtPfCEFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCEFCutHi", tempStr)) jtPfCEFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMinMult", tempStr)) jtPfMinMult = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMinChgMult", tempStr)) jtPfMinChgMult = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nSyst", tempStr)) nSyst = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("systStr", tempStr)) systStr = StringToStringVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nBayes", tempStr)) nBayes = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nBigBayesSymm", tempStr)) nBigBayesSymm = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("bayesVal", tempStr)) bayesVal = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nSuperBayes", tempStr)) nSuperBayes = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nHistDim", tempStr)) nHistDim = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("histTag", tempStr)) histTag = StringToStringVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("histBestBayes", tempStr)) histBestBayes = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else{
      std::cout << "WARNING: TNAMED \'" << tempStr << "\' is unaccounted for in cutPropagator. Consider fixing" << std::endl;
    }
  }
  
  //  std::cout << "Handling cutPropDirs: " << std::endl;
  for(unsigned int dI = 0; dI < cutDirTDir.size(); ++dI){
    //    std::cout << " " << dI << "/" << cutDirTDir.size() << ": " << cutDirTDir.at(dI) << std::endl;
    TDirectory* tempDir_p = (TDirectory*)inFile_p->Get(cutDirTDir.at(dI).c_str());
    TIter next(tempDir_p->GetListOfKeys());
    TKey* tempKey = NULL;

    while( (tempKey = (TKey*)next()) ){
      const std::string name = tempKey->GetName();
      const std::string className = tempKey->GetClassName();

      if(className.find("TNamed") == std::string::npos) continue;

      TNamed* tempName = (TNamed*)tempKey->ReadObj();
      if(cutDirTDir.at(dI).find("unfoldDir") != std::string::npos){
	histTag.push_back(std::string(tempName->GetName()));
	histBestBayes.push_back(std::stoi(std::string(tempName->GetTitle())));
      }
      else if(cutDirTDir.at(dI).find("subDir") != std::string::npos){
	inFullFileNames.push_back(std::string(tempName->GetTitle()));
      }

    }

  }

  return true;
}


bool cutPropagator::WriteAllVarToFile(TFile* inFile_p, TDirectory* inDir_p, TDirectory* inSubDir_p, TDirectory* unfoldDir_p = NULL)
{
  if(inDir_p == NULL){
    std::cout << "cutPropafator::WriteAllVarToFile - Given indir is null. return false" << std::endl;
    return false;
  }
  else if(inSubDir_p == NULL){
    std::cout << "cutPropafator::WriteAllVarToFile - Given insubdir is null. return false" << std::endl;
    return false;
  }
  else if(std::string(inDir_p->GetName()).size() != std::string("cutDir").size()){
    std::cout << "cutPropafator::WriteAllVarToFile - Given indir is not \'cutDir\'. return false" << std::endl;
    return false;
  }
  else if(std::string(inDir_p->GetName()).find("cutDir") == std::string::npos){
    std::cout << "cutPropafator::WriteAllVarToFile - Given indir is not \'cutDir\'. return false" << std::endl;
    return false;
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inFile_p->cd();
  inDir_p->cd();

  std::string inFileNames2 = "";
  for(unsigned int i = 0; i < inFileNames.size(); ++i){
    inFileNames2 = inFileNames2 + inFileNames.at(i) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string responseModStr = "";
  std::string jerVarDataStr = "";
  for(int jI = 0; jI < nResponseMod; ++jI){
    responseModStr = responseModStr + std::to_string(responseMod.at(jI)) + ",";
    jerVarDataStr = jerVarDataStr + std::to_string(jerVarData.at(jI)) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string jtPtBinsStr = "";
  for(int jI = 0; jI < nJtPtBins+1; ++jI){
    jtPtBinsStr = jtPtBinsStr + std::to_string(jtPtBins.at(jI)) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string jtAlgosStr = "";
  std::string minJtPtCutStr = "";
  std::string multiJtPtCutStr = "";
  std::string recoTruncPosStr = "";

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(int jI = 0; jI < nJtAlgos; ++jI){
    if(jtAlgos.size() != 0) jtAlgosStr = jtAlgosStr + jtAlgos.at(jI) + ",";
    if(minJtPtCut.size() != 0) minJtPtCutStr = minJtPtCutStr + prettyString(minJtPtCut.at(jI), 1, false) + ",";
    if(multiJtPtCut.size() != 0) multiJtPtCutStr = multiJtPtCutStr + prettyString(multiJtPtCut.at(jI), 1, false) + ",";
    if(recoTruncPos.size() != 0) recoTruncPosStr = recoTruncPosStr + std::to_string((int)(recoTruncPos.at(jI))) + ",";
  }

  std::string jtAbsEtaBinsLowStr = "";
  std::string jtAbsEtaBinsHiStr = "";
  for(int jI = 0; jI < nJtAbsEtaBins; ++jI){
    jtAbsEtaBinsLowStr = jtAbsEtaBinsLowStr + std::to_string(jtAbsEtaBinsLow.at(jI)) + ",";
    jtAbsEtaBinsHiStr = jtAbsEtaBinsHiStr + std::to_string(jtAbsEtaBinsHi.at(jI)) + ",";
  }

  std::string nPthatsStr = std::to_string(nPthats);
  std::string pthatsStr = "";
  for(unsigned int jI = 0; jI < pthats.size(); ++jI){
    pthatsStr = pthatsStr + std::to_string(pthats.at(jI)) + ",";
  }

  std::string pthatWeightsStr = "";
  for(unsigned int jI = 0; jI < pthatWeights.size(); ++jI){
    pthatWeightsStr = pthatWeightsStr + to_string_with_precision(pthatWeights.at(jI), 12) + ",";
    std::cout << "Weight orig, string: " << pthatWeights.at(jI) << ", " << to_string_with_precision(pthatWeights.at(jI), 12) << std::endl;
  }

  std::string centBinsLowStr = "";
  for(int jI = 0; jI < nCentBins; ++jI){
    centBinsLowStr = centBinsLowStr + std::to_string(centBinsLow.at(jI)) + ",";
  }

  std::string centBinsHiStr = "";
  for(int jI = 0; jI < nCentBins; ++jI){
    centBinsHiStr = centBinsHiStr + std::to_string(centBinsHi.at(jI)) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string nIDStr = std::to_string(nID);
  std::string idStr2 = "";
  std::string jtPfCHMFCutLowStr = "";
  std::string jtPfCHMFCutHiStr = "";
  std::string jtPfMUMFCutLowStr = "";
  std::string jtPfMUMFCutHiStr = "";
  std::string jtPfNHFCutLowStr = "";
  std::string jtPfNHFCutHiStr = "";
  std::string jtPfNEFCutLowStr = "";
  std::string jtPfNEFCutHiStr = "";
  std::string jtPfMUFCutLowStr = "";
  std::string jtPfMUFCutHiStr = "";
  std::string jtPfCHFCutLowStr = "";
  std::string jtPfCHFCutHiStr = "";
  std::string jtPfCEFCutLowStr = "";
  std::string jtPfCEFCutHiStr = "";
  std::string jtPfMinMultStr = "";
  std::string jtPfMinChgMultStr = "";

  for(int iI = 0; iI < nID; ++iI){
    idStr2 = idStr2 + idStr.at(iI) + ",";
    jtPfCHMFCutLowStr = jtPfCHMFCutLowStr + std::to_string(jtPfCHMFCutLow.at(iI)) + ",";
    jtPfCHMFCutHiStr = jtPfCHMFCutHiStr + std::to_string(jtPfCHMFCutHi.at(iI)) + ",";
    jtPfMUMFCutLowStr = jtPfMUMFCutLowStr + std::to_string(jtPfMUMFCutLow.at(iI)) + ",";
    jtPfMUMFCutHiStr = jtPfMUMFCutHiStr + std::to_string(jtPfMUMFCutHi.at(iI)) + ",";
    jtPfNHFCutLowStr = jtPfNHFCutLowStr + std::to_string(jtPfNHFCutLow.at(iI)) + ",";
    jtPfNHFCutHiStr = jtPfNHFCutHiStr + std::to_string(jtPfNHFCutHi.at(iI)) + ",";
    jtPfNEFCutLowStr = jtPfNEFCutLowStr + std::to_string(jtPfNEFCutLow.at(iI)) + ",";
    jtPfNEFCutHiStr = jtPfNEFCutHiStr + std::to_string(jtPfNEFCutHi.at(iI)) + ",";
    jtPfMUFCutLowStr = jtPfMUFCutLowStr + std::to_string(jtPfMUFCutLow.at(iI)) + ",";
    jtPfMUFCutHiStr = jtPfMUFCutHiStr + std::to_string(jtPfMUFCutHi.at(iI)) + ",";
    jtPfCHFCutLowStr = jtPfCHFCutLowStr + std::to_string(jtPfCHFCutLow.at(iI)) + ",";
    jtPfCHFCutHiStr = jtPfCHFCutHiStr + std::to_string(jtPfCHFCutHi.at(iI)) + ",";
    jtPfCEFCutLowStr = jtPfCEFCutLowStr + std::to_string(jtPfCEFCutLow.at(iI)) + ",";
    jtPfCEFCutHiStr = jtPfCEFCutHiStr + std::to_string(jtPfCEFCutHi.at(iI)) + ",";
    jtPfMinMultStr = jtPfMinMultStr + std::to_string(jtPfMinMult.at(iI)) + ",";
    jtPfMinChgMultStr = jtPfMinChgMultStr + std::to_string(jtPfMinChgMult.at(iI)) + ",";
  }

  std::string nSystStr = std::to_string(nSyst);
  std::string systStr2 = "";

  for(int sI = 0; sI < nSyst; ++sI){
    systStr2 = systStr2 + systStr.at(sI) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string nBayesStr = std::to_string(nBayes);
  std::string nBigBayesSymmStr = std::to_string(nBigBayesSymm);
  std::string nSuperBayesStr = std::to_string(nSuperBayes);
  std::string bayesVal2 = "";

  for(int sI = 0; sI < nBayes; ++sI){
    bayesVal2 = bayesVal2 + std::to_string(bayesVal.at(sI)) + ",";
  }

  std::string nHistDimStr = std::to_string(nHistDim);
  

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TNamed inFileNamesName("inFileNames", inFileNames2.c_str());
  TNamed isPPName("isPP", std::to_string(isPP).c_str());
  TNamed jtAbsEtaMaxName("jtAbsEtaMax", std::to_string(jtAbsEtaMax).c_str());
  TNamed rcDiffFileNameName("rcDiffFileName", rcDiffFileName.c_str());
  TNamed flatPriorFileNameName("flatPriorFileName", flatPriorFileName.c_str());
  TNamed jecVarMCName("jecVarMC", std::to_string(jecVarMC).c_str());
  TNamed jerVarMCName("jerVarMC", std::to_string(jerVarMC).c_str());
  TNamed jecVarDataName("jecVarData", std::to_string(jecVarData).c_str());
  TNamed nResponseModName("nResponseMod", std::to_string(nResponseMod).c_str());
  TNamed responseModName("responseMod", responseModStr.c_str());
  TNamed jerVarDataName("jerVarData", jerVarDataStr.c_str());
  TNamed nJtAlgosName("nJtAlgos", std::to_string(nJtAlgos).c_str());
  TNamed jtAlgosName("jtAlgos", jtAlgosStr.c_str());
  TNamed minJtPtCutName("minJtPtCut", minJtPtCutStr.c_str());
  TNamed multiJtPtCutName("multiJtPtCut", multiJtPtCutStr.c_str());
  TNamed recoTruncPosName("recoTruncPos", recoTruncPosStr.c_str());
  TNamed nJtPtBinsName("nJtPtBins", std::to_string(nJtPtBins).c_str());
  TNamed jtPtBinsName("jtPtBins", jtPtBinsStr.c_str());
  TNamed nJtAbsEtaBinsName("nJtAbsEtaBins", std::to_string(nJtAbsEtaBins).c_str());
  TNamed jtAbsEtaBinsLowName("jtAbsEtaBinsLow", jtAbsEtaBinsLowStr.c_str());
  TNamed jtAbsEtaBinsHiName("jtAbsEtaBinsHi", jtAbsEtaBinsHiStr.c_str());
  TNamed nPthatsName("nPthats", nPthatsStr.c_str());
  TNamed pthatsName("pthats", pthatsStr.c_str());
  TNamed pthatWeightsName("pthatWeights", pthatWeightsStr.c_str());
  TNamed nCentBinsName("nCentBins", std::to_string(nCentBins).c_str());
  TNamed centBinsLowName("centBinsLow", centBinsLowStr.c_str());
  TNamed centBinsHiName("centBinsHi", centBinsHiStr.c_str());
  TNamed nIDName("nID", nIDStr.c_str());
  TNamed idStrName("idStr", idStr2.c_str());
  TNamed jtPfCHMFCutLowName("jtPfCHMFCutLow", jtPfCHMFCutLowStr.c_str());
  TNamed jtPfCHMFCutHiName("jtPfCHMFCutHi", jtPfCHMFCutHiStr.c_str());
  TNamed jtPfMUMFCutLowName("jtPfMUMFCutLow", jtPfMUMFCutLowStr.c_str());
  TNamed jtPfMUMFCutHiName("jtPfMUMFCutHi", jtPfMUMFCutHiStr.c_str());
  TNamed jtPfNHFCutLowName("jtPfNHFCutLow", jtPfNHFCutLowStr.c_str());
  TNamed jtPfNHFCutHiName("jtPfNHFCutHi", jtPfNHFCutHiStr.c_str());
  TNamed jtPfNEFCutLowName("jtPfNEFCutLow", jtPfNEFCutLowStr.c_str());
  TNamed jtPfNEFCutHiName("jtPfNEFCutHi", jtPfNEFCutHiStr.c_str());
  TNamed jtPfMUFCutLowName("jtPfMUFCutLow", jtPfMUFCutLowStr.c_str());
  TNamed jtPfMUFCutHiName("jtPfMUFCutHi", jtPfMUFCutHiStr.c_str());
  TNamed jtPfCHFCutLowName("jtPfCHFCutLow", jtPfCHFCutLowStr.c_str());
  TNamed jtPfCHFCutHiName("jtPfCHFCutHi", jtPfCHFCutHiStr.c_str());
  TNamed jtPfCEFCutLowName("jtPfCEFCutLow", jtPfCEFCutLowStr.c_str());
  TNamed jtPfCEFCutHiName("jtPfCEFCutHi", jtPfCEFCutHiStr.c_str());
  TNamed jtPfMinMultName("jtPfMinMult", jtPfMinMultStr.c_str());
  TNamed jtPfMinChgMultName("jtPfMinChgMult", jtPfMinChgMultStr.c_str());
  TNamed nSystName("nSyst", nSystStr.c_str());
  TNamed systStrName("systStr", systStr2.c_str());
  TNamed nBayesName("nBayes", nBayesStr.c_str());
  TNamed nBigBayesSymmName("nBigBayesSymm", nBigBayesSymmStr.c_str());
  TNamed nSuperBayesName("nSuperBayes", nSuperBayesStr.c_str());
  TNamed bayesValName("bayesVal", bayesVal2.c_str());
  TNamed nHistDimName("nHistDim", nHistDimStr.c_str());
  
  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inFileNamesName.Write("", TObject::kOverwrite);
  isPPName.Write("", TObject::kOverwrite);
  jtAbsEtaMaxName.Write("", TObject::kOverwrite);
  rcDiffFileNameName.Write("", TObject::kOverwrite);
  flatPriorFileNameName.Write("", TObject::kOverwrite);
  jecVarMCName.Write("", TObject::kOverwrite);
  jerVarMCName.Write("", TObject::kOverwrite);
  jecVarDataName.Write("", TObject::kOverwrite);
  nResponseModName.Write("", TObject::kOverwrite);
  responseModName.Write("", TObject::kOverwrite);
  jerVarDataName.Write("", TObject::kOverwrite);
  nJtAlgosName.Write("", TObject::kOverwrite);
  jtAlgosName.Write("", TObject::kOverwrite);
  minJtPtCutName.Write("", TObject::kOverwrite);
  multiJtPtCutName.Write("", TObject::kOverwrite);
  recoTruncPosName.Write("", TObject::kOverwrite);
  nJtPtBinsName.Write("", TObject::kOverwrite);
  jtPtBinsName.Write("", TObject::kOverwrite);
  nJtAbsEtaBinsName.Write("", TObject::kOverwrite);
  jtAbsEtaBinsLowName.Write("", TObject::kOverwrite);
  jtAbsEtaBinsHiName.Write("", TObject::kOverwrite);
  nPthatsName.Write("", TObject::kOverwrite);
  pthatsName.Write("", TObject::kOverwrite);
  pthatWeightsName.Write("", TObject::kOverwrite);
  nCentBinsName.Write("", TObject::kOverwrite);
  centBinsLowName.Write("", TObject::kOverwrite);
  centBinsHiName.Write("", TObject::kOverwrite);
  nIDName.Write("", TObject::kOverwrite);
  idStrName.Write("", TObject::kOverwrite);
  jtPfCHMFCutLowName.Write("", TObject::kOverwrite);
  jtPfCHMFCutHiName.Write("", TObject::kOverwrite);
  jtPfMUMFCutLowName.Write("", TObject::kOverwrite);
  jtPfMUMFCutHiName.Write("", TObject::kOverwrite);
  jtPfNHFCutLowName.Write("", TObject::kOverwrite);
  jtPfNHFCutHiName.Write("", TObject::kOverwrite);
  jtPfNEFCutLowName.Write("", TObject::kOverwrite);
  jtPfNEFCutHiName.Write("", TObject::kOverwrite);
  jtPfMUFCutLowName.Write("", TObject::kOverwrite);
  jtPfMUFCutHiName.Write("", TObject::kOverwrite);
  jtPfCHFCutLowName.Write("", TObject::kOverwrite);
  jtPfCHFCutHiName.Write("", TObject::kOverwrite);
  jtPfCEFCutLowName.Write("", TObject::kOverwrite);
  jtPfCEFCutHiName.Write("", TObject::kOverwrite);
  jtPfCEFCutHiName.Write("", TObject::kOverwrite);
  jtPfMinMultName.Write("", TObject::kOverwrite);
  jtPfMinChgMultName.Write("", TObject::kOverwrite);
  nSystName.Write("", TObject::kOverwrite);
  systStrName.Write("", TObject::kOverwrite);
  nBayesName.Write("", TObject::kOverwrite);
  nSuperBayesName.Write("", TObject::kOverwrite);
  nBigBayesSymmName.Write("", TObject::kOverwrite);
  bayesValName.Write("", TObject::kOverwrite);
  nHistDimName.Write("", TObject::kOverwrite);

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inFile_p->cd();
  inDir_p->cd();
  inSubDir_p->cd();  


  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<TNamed> inFullFileNames2;
  for(unsigned int i = 0; i < inFullFileNames.size(); ++i){
    inFullFileNames2.push_back(TNamed(("inFullFileNames_" + std::to_string(i)).c_str(), inFullFileNames.at(i).c_str()));
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int i = 0; i < inFullFileNames2.size(); ++i){
    inFullFileNames2.at(i).Write("", TObject::kOverwrite);
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(unfoldDir_p != NULL){
    std::cout << "Attempting a write of all unfold obj" << std::endl;

    inFile_p->cd();
    inDir_p->cd();

    std::cout << " " << unfoldDir_p->GetName() << std::endl;

    unfoldDir_p->cd();
    nHistDimName.Write("", TObject::kOverwrite);
    
    for(Int_t hI = 0; hI < nHistDim; ++hI){
      unfoldDir_p->cd();
      TNamed tempTName(histTag.at(hI).c_str(), std::to_string(histBestBayes.at(hI)).c_str());
      tempTName.Write("", TObject::kOverwrite);
    }
  }

  return true;
}


bool cutPropagator::CheckPropagatorsMatch(cutPropagator inCutProp, bool doBothMCOrBothData, bool doBothPPOrBothPbPb, bool skipJtAlgo = false)
{
  if(!CheckJtAbsEtaMax(inCutProp)) return false;
  if(!CheckNJtPtBins(inCutProp)) return false;
  if(!CheckNJtAbsEtaBins(inCutProp)) return false;
  if(!CheckNID(inCutProp)) return false;
  if(!CheckJtPtBins(inCutProp)) return false;
  if(!CheckJtAbsEtaBinsLow(inCutProp)) return false;
  if(!CheckJtAbsEtaBinsHi(inCutProp)) return false;
  if(!CheckIdStr(inCutProp)) return false;
  if(!CheckJtPfCHMFCutLow(inCutProp)) return false;
  if(!CheckJtPfCHMFCutHi(inCutProp)) return false;
  if(!CheckJtPfMUMFCutLow(inCutProp)) return false;
  if(!CheckJtPfMUMFCutHi(inCutProp)) return false;
  if(!CheckJtPfNHFCutLow(inCutProp)) return false;
  if(!CheckJtPfNHFCutHi(inCutProp)) return false;
  if(!CheckJtPfNEFCutLow(inCutProp)) return false;
  if(!CheckJtPfNEFCutHi(inCutProp)) return false;
  if(!CheckJtPfMUFCutLow(inCutProp)) return false;
  if(!CheckJtPfMUFCutHi(inCutProp)) return false;
  if(!CheckJtPfCHFCutLow(inCutProp)) return false;
  if(!CheckJtPfCHFCutHi(inCutProp)) return false;
  if(!CheckJtPfCEFCutLow(inCutProp)) return false;
  if(!CheckJtPfCEFCutHi(inCutProp)) return false;
  if(!CheckJtPfMinMult(inCutProp)) return false;
  if(!CheckJtPfMinChgMult(inCutProp)) return false;

  if(!CheckNBayes(inCutProp)) return false;
  if(!CheckNBigBayesSymm(inCutProp)) return false;
  if(!CheckBayesVal(inCutProp)) return false;
  if(!CheckNSuperBayes(inCutProp)) return false;

  if(doBothMCOrBothData){
    if(!CheckNSyst(inCutProp)) return false;
    if(!CheckSystStr(inCutProp)) return false;
    if(!CheckRCDiffFileName(inCutProp)) return false;
    if(!CheckFlatPriorFileName(inCutProp)) return false;
    if(!CheckJECVarMC(inCutProp)) return false;
    if(!CheckJERVarMC(inCutProp)) return false;
    if(!CheckJECVarData(inCutProp)) return false;
    if(!CheckNResponseMod(inCutProp)) return false;
    if(!CheckResponseMod(inCutProp)) return false;
    if(!CheckJERVarData(inCutProp)) return false;
    if(!CheckNPthats(inCutProp)) return false;
    if(!CheckPthats(inCutProp)) return false;
    if(!CheckPthatWeights(inCutProp)) return false;
  }

  if(doBothPPOrBothPbPb){
    if(skipJtAlgo){
      if(!CheckNJtAlgos(inCutProp)) return false;
      if(!CheckJtAlgos(inCutProp)) return false;
      if(!CheckMinJtPtCut(inCutProp)) return false;
      if(!CheckMultiJtPtCut(inCutProp)) return false;
      if(!CheckRecoTruncPos(inCutProp)) return false;
    }
    if(!CheckIsPP(inCutProp)) return false;
    if(!CheckNCentBins(inCutProp)) return false;
    if(!CheckCentBinsLow(inCutProp)) return false;
    if(!CheckCentBinsHi(inCutProp)) return false;
  }

  return true;
}

bool cutPropagator::CheckDouble(double in1, double in2)
{
  if(in1 + delta < in2) return false;
  else if(in1 - delta > in2) return false;
  return true;
}

bool cutPropagator::CheckInt(int in1, int in2){return in1 == in2;}

bool cutPropagator::CheckBool(bool in1, bool in2){return in1 == in2;}

bool cutPropagator::CheckString(std::string in1, std::string in2){return isStrSame(in1, in2);}

bool cutPropagator::CheckVectDouble(std::vector<double> in1, std::vector<double> in2)
{
  if(in1.size() != in2.size()) return false;
  for(unsigned int i = 0; i < in1.size(); ++i){
    if(!CheckDouble(in1.at(i), in2.at(i))) return false;
  }
  return true;
}

bool cutPropagator::CheckVectInt(std::vector<int> in1, std::vector<int> in2)
{
  if(in1.size() != in2.size()) return false;
  for(unsigned int i = 0; i < in1.size(); ++i){
    if(!CheckInt(in1.at(i), in2.at(i))) return false;
  }
  return true;
}

bool cutPropagator::CheckVectString(std::vector<std::string> in1, std::vector<std::string> in2)
{
  if(in1.size() != in2.size()) return false;
  for(unsigned int i = 0; i < in1.size(); ++i){
    if(!CheckString(in1.at(i), in2.at(i))) return false;
  }
  return true;
}

bool cutPropagator::CheckJtAbsEtaMax(double inJtAbsEtaMax){return CheckDouble(jtAbsEtaMax, inJtAbsEtaMax);}
bool cutPropagator::CheckJtAbsEtaMax(cutPropagator inCutProp){return CheckJtAbsEtaMax(inCutProp.GetJtAbsEtaMax());}
bool cutPropagator::CheckNJtPtBins(int inNJtPtBins){return CheckInt(nJtPtBins, inNJtPtBins);}
bool cutPropagator::CheckNJtPtBins(cutPropagator inCutProp){return CheckNJtPtBins(inCutProp.GetNJtPtBins());}
bool cutPropagator::CheckNJtAbsEtaBins(int inNJtAbsEtaBins){return CheckInt(nJtAbsEtaBins, inNJtAbsEtaBins);}
bool cutPropagator::CheckNJtAbsEtaBins(cutPropagator inCutProp){return CheckNJtAbsEtaBins(inCutProp.GetNJtAbsEtaBins());}
bool cutPropagator::CheckNID(int inNID){return CheckInt(nID, inNID);}
bool cutPropagator::CheckNID(cutPropagator inCutProp){return CheckNID(inCutProp.GetNID());}
bool cutPropagator::CheckJtPtBins(std::vector<double> inJtPtBins){return CheckVectDouble(jtPtBins, inJtPtBins);}
bool cutPropagator::CheckJtPtBins(cutPropagator inCutProp){return CheckJtPtBins(inCutProp.GetJtPtBins());}
bool cutPropagator::CheckJtAbsEtaBinsLow(std::vector<double> inJtAbsEtaBinsLow){return CheckVectDouble(jtAbsEtaBinsLow, inJtAbsEtaBinsLow);}
bool cutPropagator::CheckJtAbsEtaBinsLow(cutPropagator inCutProp){return CheckJtAbsEtaBinsLow(inCutProp.GetJtAbsEtaBinsLow());}
bool cutPropagator::CheckJtAbsEtaBinsHi(std::vector<double> inJtAbsEtaBinsHi){return CheckVectDouble(jtAbsEtaBinsHi, inJtAbsEtaBinsHi);}
bool cutPropagator::CheckJtAbsEtaBinsHi(cutPropagator inCutProp){return CheckJtAbsEtaBinsHi(inCutProp.GetJtAbsEtaBinsHi());}
bool cutPropagator::CheckIdStr(std::vector<std::string> inIDStr){return CheckVectString(idStr, inIDStr);}
bool cutPropagator::CheckIdStr(cutPropagator inCutProp){return CheckIdStr(inCutProp.GetIdStr());}
bool cutPropagator::CheckJtPfCHMFCutLow(std::vector<double> inJtPfCHMFCutLow){return CheckVectDouble(jtPfCHMFCutLow, inJtPfCHMFCutLow);}
bool cutPropagator::CheckJtPfCHMFCutLow(cutPropagator inCutProp){return CheckJtPfCHMFCutLow(inCutProp.GetJtPfCHMFCutLow());}
bool cutPropagator::CheckJtPfCHMFCutHi(std::vector<double> inJtPfCHMFCutHi){return CheckVectDouble(jtPfCHMFCutHi, inJtPfCHMFCutHi);}
bool cutPropagator::CheckJtPfCHMFCutHi(cutPropagator inCutProp){return CheckJtPfCHMFCutHi(inCutProp.GetJtPfCHMFCutHi());}
bool cutPropagator::CheckJtPfMUMFCutLow(std::vector<double> inJtPfMUMFCutLow){return CheckVectDouble(jtPfMUMFCutLow, inJtPfMUMFCutLow);}
bool cutPropagator::CheckJtPfMUMFCutLow(cutPropagator inCutProp){return CheckJtPfMUMFCutLow(inCutProp.GetJtPfMUMFCutLow());}
bool cutPropagator::CheckJtPfMUMFCutHi(std::vector<double> inJtPfMUMFCutHi){return CheckVectDouble(jtPfMUMFCutHi, inJtPfMUMFCutHi);}
bool cutPropagator::CheckJtPfMUMFCutHi(cutPropagator inCutProp){return CheckJtPfMUMFCutHi(inCutProp.GetJtPfMUMFCutHi());}
bool cutPropagator::CheckJtPfNHFCutLow(std::vector<double> inJtPfNHFCutLow){return CheckVectDouble(jtPfNHFCutLow, inJtPfNHFCutLow);}
bool cutPropagator::CheckJtPfNHFCutLow(cutPropagator inCutProp){return CheckJtPfNHFCutLow(inCutProp.GetJtPfNHFCutLow());}
bool cutPropagator::CheckJtPfNHFCutHi(std::vector<double> inJtPfNHFCutHi){return CheckVectDouble(jtPfNHFCutHi, inJtPfNHFCutHi);}
bool cutPropagator::CheckJtPfNHFCutHi(cutPropagator inCutProp){return CheckJtPfNHFCutHi(inCutProp.GetJtPfNHFCutHi());}
bool cutPropagator::CheckJtPfNEFCutLow(std::vector<double> inJtPfNEFCutLow){return CheckVectDouble(jtPfNEFCutLow, inJtPfNEFCutLow);}
bool cutPropagator::CheckJtPfNEFCutLow(cutPropagator inCutProp){return CheckJtPfNEFCutLow(inCutProp.GetJtPfNEFCutLow());}
bool cutPropagator::CheckJtPfNEFCutHi(std::vector<double> inJtPfNEFCutHi){return CheckVectDouble(jtPfNEFCutHi, inJtPfNEFCutHi);}
bool cutPropagator::CheckJtPfNEFCutHi(cutPropagator inCutProp){return CheckJtPfNEFCutHi(inCutProp.GetJtPfNEFCutHi());}
bool cutPropagator::CheckJtPfMUFCutLow(std::vector<double> inJtPfMUFCutLow){return CheckVectDouble(jtPfMUFCutLow, inJtPfMUFCutLow);}
bool cutPropagator::CheckJtPfMUFCutLow(cutPropagator inCutProp){return CheckJtPfMUFCutLow(inCutProp.GetJtPfMUFCutLow());}
bool cutPropagator::CheckJtPfMUFCutHi(std::vector<double> inJtPfMUFCutHi){return CheckVectDouble(jtPfMUFCutHi, inJtPfMUFCutHi);}
bool cutPropagator::CheckJtPfMUFCutHi(cutPropagator inCutProp){return CheckJtPfMUFCutHi(inCutProp.GetJtPfMUFCutHi());}
bool cutPropagator::CheckJtPfCHFCutLow(std::vector<double> inJtPfCHFCutLow){return CheckVectDouble(jtPfCHFCutLow, inJtPfCHFCutLow);}
bool cutPropagator::CheckJtPfCHFCutLow(cutPropagator inCutProp){return CheckJtPfCHFCutLow(inCutProp.GetJtPfCHFCutLow());}
bool cutPropagator::CheckJtPfCHFCutHi(std::vector<double> inJtPfCHFCutHi){return CheckVectDouble(jtPfCHFCutHi, inJtPfCHFCutHi);}
bool cutPropagator::CheckJtPfCHFCutHi(cutPropagator inCutProp){return CheckJtPfCHFCutHi(inCutProp.GetJtPfCHFCutHi());}
bool cutPropagator::CheckJtPfCEFCutLow(std::vector<double> inJtPfCEFCutLow){return CheckVectDouble(jtPfCEFCutLow, inJtPfCEFCutLow);}
bool cutPropagator::CheckJtPfCEFCutLow(cutPropagator inCutProp){return CheckJtPfCEFCutLow(inCutProp.GetJtPfCEFCutLow());}
bool cutPropagator::CheckJtPfCEFCutHi(std::vector<double> inJtPfCEFCutHi){return CheckVectDouble(jtPfCEFCutHi, inJtPfCEFCutHi);}
bool cutPropagator::CheckJtPfCEFCutHi(cutPropagator inCutProp){return CheckJtPfCEFCutHi(inCutProp.GetJtPfCEFCutHi());}
bool cutPropagator::CheckJtPfMinMult(std::vector<int> inJtPfMinMult){return CheckVectInt(jtPfMinMult, inJtPfMinMult);}
bool cutPropagator::CheckJtPfMinMult(cutPropagator inCutProp){return CheckJtPfMinMult(inCutProp.GetJtPfMinMult());}
bool cutPropagator::CheckJtPfMinChgMult(std::vector<int> inJtPfMinChgMult){return CheckVectInt(jtPfMinChgMult, inJtPfMinChgMult);}
bool cutPropagator::CheckJtPfMinChgMult(cutPropagator inCutProp){return CheckJtPfMinChgMult(inCutProp.GetJtPfMinChgMult());}
bool cutPropagator::CheckNSyst(int inNSyst){return CheckInt(inNSyst, nSyst);}
bool cutPropagator::CheckNSyst(cutPropagator inCutProp){return CheckNSyst(inCutProp.GetNSyst());}
bool cutPropagator::CheckSystStr(std::vector<std::string> inSystStr){return CheckVectString(inSystStr, systStr);}
bool cutPropagator::CheckSystStr(cutPropagator inCutProp){return CheckSystStr(inCutProp.GetSystStr());}
bool cutPropagator::CheckNBayes(int inNBayes){return CheckInt(inNBayes, nBayes);}
bool cutPropagator::CheckNBayes(cutPropagator inCutProp){return CheckNBayes(inCutProp.GetNBayes());}
bool cutPropagator::CheckNBigBayesSymm(int inNBigBayesSymm){return CheckInt(inNBigBayesSymm, nBigBayesSymm);}
bool cutPropagator::CheckNBigBayesSymm(cutPropagator inCutProp){return CheckNBigBayesSymm(inCutProp.GetNBigBayesSymm());}
bool cutPropagator::CheckBayesVal(std::vector<int> inBayesVal){return CheckVectInt(inBayesVal, bayesVal);}
bool cutPropagator::CheckBayesVal(cutPropagator inCutProp){return CheckBayesVal(inCutProp.GetBayesVal());}
bool cutPropagator::CheckNSuperBayes(int inNSuperBayes){return CheckInt(inNSuperBayes, nSuperBayes);}
bool cutPropagator::CheckNSuperBayes(cutPropagator inCutProp){return CheckNSuperBayes(inCutProp.GetNSuperBayes());}
bool cutPropagator::CheckRCDiffFileName(std::string inRCDiffFileName){return CheckString(inRCDiffFileName, rcDiffFileName);}
bool cutPropagator::CheckRCDiffFileName(cutPropagator inCutProp){return CheckRCDiffFileName(inCutProp.GetRCDiffFileName());}
bool cutPropagator::CheckFlatPriorFileName(std::string inFlatPriorFileName){return CheckString(inFlatPriorFileName, flatPriorFileName);}
bool cutPropagator::CheckFlatPriorFileName(cutPropagator inCutProp){return CheckFlatPriorFileName(inCutProp.GetFlatPriorFileName());}
bool cutPropagator::CheckJECVarMC(double inJECVarMC){return CheckDouble(inJECVarMC, jecVarMC);}
bool cutPropagator::CheckJECVarMC(cutPropagator inCutProp){return CheckJECVarMC(inCutProp.GetJECVarMC());}
bool cutPropagator::CheckJERVarMC(double inJERVarMC){return CheckDouble(inJERVarMC, jerVarMC);}
bool cutPropagator::CheckJERVarMC(cutPropagator inCutProp){return CheckJERVarMC(inCutProp.GetJERVarMC());}
bool cutPropagator::CheckJECVarData(double inJECVarData){return CheckDouble(inJECVarData, jecVarData);}
bool cutPropagator::CheckJECVarData(cutPropagator inCutProp){return CheckJECVarData(inCutProp.GetJECVarData());}
bool cutPropagator::CheckNResponseMod(int inNResponseMod){return CheckInt(inNResponseMod, nResponseMod);}
bool cutPropagator::CheckNResponseMod(cutPropagator inCutProp){return CheckNResponseMod(inCutProp.GetNResponseMod());}
bool cutPropagator::CheckResponseMod(std::vector<double> inResponseMod){return CheckVectDouble(inResponseMod, responseMod);}
bool cutPropagator::CheckResponseMod(cutPropagator inCutProp){return CheckResponseMod(inCutProp.GetResponseMod());}
bool cutPropagator::CheckJERVarData(std::vector<double> inJERVarData){return CheckVectDouble(inJERVarData, jerVarData);}
bool cutPropagator::CheckJERVarData(cutPropagator inCutProp){return CheckJERVarData(inCutProp.GetJERVarData());}
bool cutPropagator::CheckNPthats(int inNPthats){return CheckInt(inNPthats, nPthats);}
bool cutPropagator::CheckNPthats(cutPropagator inCutProp){return CheckNPthats(inCutProp.GetNPthats());}
bool cutPropagator::CheckPthats(std::vector<double> inPthats){return CheckVectDouble(inPthats, pthats);}
bool cutPropagator::CheckPthats(cutPropagator inCutProp){return CheckPthats(inCutProp.GetPthats());}
bool cutPropagator::CheckPthatWeights(std::vector<double> inPthatWeights){return CheckVectDouble(inPthatWeights, pthatWeights);}
bool cutPropagator::CheckPthatWeights(cutPropagator inCutProp){return CheckPthatWeights(inCutProp.GetPthatWeights());}
bool cutPropagator::CheckNJtAlgos(int inNJtAlgos){return CheckInt(inNJtAlgos, nJtAlgos);}
bool cutPropagator::CheckNJtAlgos(cutPropagator inCutProp){return CheckNJtAlgos(inCutProp.GetNJtAlgos());}
bool cutPropagator::CheckJtAlgos(std::vector<std::string> inJtAlgos){return CheckVectString(inJtAlgos, jtAlgos);}
bool cutPropagator::CheckJtAlgos(cutPropagator inCutProp){return CheckJtAlgos(inCutProp.GetJtAlgos());}
bool cutPropagator::CheckMinJtPtCut(std::vector<double> inMinJtPtCut){return CheckVectDouble(inMinJtPtCut, minJtPtCut);}
bool cutPropagator::CheckMinJtPtCut(cutPropagator inCutProp){return CheckMinJtPtCut(inCutProp.GetMinJtPtCut());}
bool cutPropagator::CheckMultiJtPtCut(std::vector<double> inMultiJtPtCut){return CheckVectDouble(inMultiJtPtCut, multiJtPtCut);}
bool cutPropagator::CheckMultiJtPtCut(cutPropagator inCutProp){return CheckMultiJtPtCut(inCutProp.GetMultiJtPtCut());}
bool cutPropagator::CheckRecoTruncPos(std::vector<int> inRecoTruncPos){return CheckVectInt(inRecoTruncPos, recoTruncPos);}
bool cutPropagator::CheckRecoTruncPos(cutPropagator inCutProp){return CheckRecoTruncPos(inCutProp.GetRecoTruncPos());}
bool cutPropagator::CheckIsPP(bool inIsPP){return CheckBool(inIsPP, isPP);}
bool cutPropagator::CheckIsPP(cutPropagator inCutProp){return CheckIsPP(inCutProp.GetIsPP());}
bool cutPropagator::CheckNCentBins(int inNCentBins){return CheckInt(inNCentBins, nCentBins);}
bool cutPropagator::CheckNCentBins(cutPropagator inCutProp){return CheckNCentBins(inCutProp.GetNCentBins());}
bool cutPropagator::CheckCentBinsLow(std::vector<int> inCentBinsLow){return CheckVectInt(inCentBinsLow, centBinsLow);}
bool cutPropagator::CheckCentBinsLow(cutPropagator inCutProp){return CheckCentBinsLow(inCutProp.GetCentBinsLow());}
bool cutPropagator::CheckCentBinsHi(std::vector<int> inCentBinsHi){return CheckVectInt(inCentBinsHi, centBinsHi);}
bool cutPropagator::CheckCentBinsHi(cutPropagator inCutProp){return CheckCentBinsHi(inCutProp.GetCentBinsHi());}

void cutPropagator::SetResponseMod(int inN, const Double_t inResponseMod[])
{
  for(int i = 0; i < inN; ++i){
    responseMod.push_back(inResponseMod[i]);
  }

  return;
}

void cutPropagator::SetJERVarData(int inN, const Double_t inJERVarData[])
{
  for(int i = 0; i < inN; ++i){
    jerVarData.push_back(inJERVarData[i]);
  }

  return;
}

void cutPropagator::SetJtAlgos(int inN, std::string inJtAlgos[])
{
  for(int i = 0; i < inN; ++i){
    jtAlgos.push_back(inJtAlgos[i]);
  }

  return;
}

void cutPropagator::SetMinJtPtCut(int inN, Double_t inMinJtPtCut[])
{
  for(int i = 0; i < inN; ++i){
    minJtPtCut.push_back(inMinJtPtCut[i]);
  }

  return;
}

void cutPropagator::SetMultiJtPtCut(int inN, Double_t inMultiJtPtCut[])
{
  for(int i = 0; i < inN; ++i){
    multiJtPtCut.push_back(inMultiJtPtCut[i]);
  }

  return;
}

void cutPropagator::SetRecoTruncPos(int inN, Int_t inRecoTruncPos[])
{
  for(int i = 0; i < inN; ++i){
    recoTruncPos.push_back(inRecoTruncPos[i]);
  }

  return;
}


void cutPropagator::SetJtPtBins(int inN, const Double_t inJtPtBins[])
{
  for(int i = 0; i < inN; ++i){
    jtPtBins.push_back(inJtPtBins[i]);
  }

  return;
}

void cutPropagator::SetJtAbsEtaBinsLow(int inN, const Double_t inJtAbsEtaBinsLow[])
{
  for(int i = 0; i < inN; ++i){
    jtAbsEtaBinsLow.push_back(inJtAbsEtaBinsLow[i]);
  }

  return;
}


void cutPropagator::SetJtAbsEtaBinsHi(int inN, const Double_t inJtAbsEtaBinsHi[])
{
  for(int i = 0; i < inN; ++i){
    jtAbsEtaBinsHi.push_back(inJtAbsEtaBinsHi[i]);
  }

  return;
}


void cutPropagator::SetPthats(int inN, const Double_t inPthats[])
{
  for(int i = 0; i < inN; ++i){
    pthats.push_back(inPthats[i]);
  }

  return;
}


void cutPropagator::SetPthatWeights(int inN, const Double_t inPthatWeights[])
{
  for(int i = 0; i < inN; ++i){
    pthatWeights.push_back(inPthatWeights[i]);
  }

  return;
}


void cutPropagator::SetCentBinsLow(int inN, const Int_t inCentBinsLow[])
{
  for(int i = 0; i < inN; ++i){
    centBinsLow.push_back(inCentBinsLow[i]);
  }

  return;
}


void cutPropagator::SetCentBinsHi(int inN, const Int_t inCentBinsHi[])
{
  for(int i = 0; i < inN; ++i){
    centBinsHi.push_back(inCentBinsHi[i]);
  }

  return;
}


void cutPropagator::SetIdStr(int inN, const std::string inIdStr[])
{
  for(int i = 0; i < inN; ++i){
    idStr.push_back(inIdStr[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHMFCutLow(int inN, const Double_t inJtPfCHMFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHMFCutLow.push_back(inJtPfCHMFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHMFCutHi(int inN, const Double_t inJtPfCHMFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHMFCutHi.push_back(inJtPfCHMFCutHi[i]);
  }

  return;
}


void cutPropagator::SetJtPfMUMFCutLow(int inN, const Double_t inJtPfMUMFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUMFCutLow.push_back(inJtPfMUMFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfMUMFCutHi(int inN, const Double_t inJtPfMUMFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUMFCutHi.push_back(inJtPfMUMFCutHi[i]);
  }

  return;
}


void cutPropagator::SetJtPfNHFCutLow(int inN, const Double_t inJtPfNHFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNHFCutLow.push_back(inJtPfNHFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfNHFCutHi(int inN, const Double_t inJtPfNHFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNHFCutHi.push_back(inJtPfNHFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfNEFCutLow(int inN, const Double_t inJtPfNEFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNEFCutLow.push_back(inJtPfNEFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfNEFCutHi(int inN, const Double_t inJtPfNEFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNEFCutHi.push_back(inJtPfNEFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfMUFCutLow(int inN, const Double_t inJtPfMUFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUFCutLow.push_back(inJtPfMUFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfMUFCutHi(int inN, const Double_t inJtPfMUFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUFCutHi.push_back(inJtPfMUFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHFCutLow(int inN, const Double_t inJtPfCHFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHFCutLow.push_back(inJtPfCHFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHFCutHi(int inN, const Double_t inJtPfCHFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHFCutHi.push_back(inJtPfCHFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfCEFCutLow(int inN, const Double_t inJtPfCEFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCEFCutLow.push_back(inJtPfCEFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfCEFCutHi(int inN, const Double_t inJtPfCEFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCEFCutHi.push_back(inJtPfCEFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfMinMult(int inN, const Int_t inJtPfMinMult[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMinMult.push_back(inJtPfMinMult[i]);
  }

  return;
}

void cutPropagator::SetJtPfMinChgMult(int inN, const Int_t inJtPfMinChgMult[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMinChgMult.push_back(inJtPfMinChgMult[i]);
  }

  return;
}

void cutPropagator::SetSystStr(int inN, const std::string inSystStr[])
{
  for(int i = 0; i < inN; ++i){
    systStr.push_back(inSystStr[i]);
  }

  return;
}

void cutPropagator::SetBayesVal(int inN, const int inBayesVal[])
{
  for(int i = 0; i < inN; ++i){
    bayesVal.push_back(inBayesVal[i]);
  }

  return;
}


void cutPropagator::SetHistTag(int inN, const std::string inHistTag[])
{
  for(int i = 0; i < inN; ++i){
    histTag.push_back(inHistTag[i]);
  }

  return;
}


void cutPropagator::SetHistBestBayes(int inN, const int inHistBestBayes[])
{
  for(int i = 0; i < inN; ++i){
    histBestBayes.push_back(inHistBestBayes[i]);
  }

  return;
}

std::vector<int> cutPropagator::StringToIntVect(std::string inStr)
{
  std::vector<int> intVect;
  while(inStr.find(",") != std::string::npos){
    intVect.push_back(std::stoi(inStr.substr(0, inStr.find(","))));
    inStr.replace(0, inStr.find(",")+1, "");
  }
  if(inStr.size() != 0) intVect.push_back(std::stoi(inStr));
  return intVect;
}

std::vector<double> cutPropagator::StringToDoubleVect(std::string inStr)
{
  std::vector<double> doubleVect;
  while(inStr.find(",") != std::string::npos){
    doubleVect.push_back(std::stod(inStr.substr(0, inStr.find(","))));
    inStr.replace(0, inStr.find(",")+1, "");
  }
  if(inStr.size() != 0) doubleVect.push_back(std::stod(inStr));
  return doubleVect;
}

std::vector<std::string> cutPropagator::StringToStringVect(std::string inStr)
{
  std::vector<std::string> stringVect;
  while(inStr.find(",") != std::string::npos){
    stringVect.push_back(inStr.substr(0, inStr.find(",")));
    inStr.replace(0, inStr.find(",")+1, "");
  }
  if(inStr.size() != 0) stringVect.push_back(inStr);
  return stringVect;
}


std::string cutPropagator::to_string_with_precision(double a_value, const int n)
{
  std::ostringstream out;
  out << std::setprecision(n) << a_value;
  return out.str();
}

#endif
