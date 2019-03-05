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

  int nSmallR;
  std::vector<int> smallRVals;
  int nLargeR;
  std::vector<int> largeRVals;

  int nRecoJtPtBinsSmallRCent0to10;
  int nGenJtPtSmallBinsSmallRCent0to10;
  int nGenJtPtLargeBinsSmallRCent0to10;
  int nRecoJtPtBinsSmallRCent10to30;
  int nGenJtPtSmallBinsSmallRCent10to30;
  int nGenJtPtLargeBinsSmallRCent10to30;
  int nRecoJtPtBinsSmallRCent30to50;
  int nGenJtPtSmallBinsSmallRCent30to50;
  int nGenJtPtLargeBinsSmallRCent30to50;
  int nRecoJtPtBinsSmallRCent50to90;
  int nGenJtPtSmallBinsSmallRCent50to90;
  int nGenJtPtLargeBinsSmallRCent50to90;

  int nRecoJtPtBinsLargeRCent0to10;
  int nGenJtPtSmallBinsLargeRCent0to10;
  int nGenJtPtLargeBinsLargeRCent0to10;
  int nRecoJtPtBinsLargeRCent10to30;
  int nGenJtPtSmallBinsLargeRCent10to30;
  int nGenJtPtLargeBinsLargeRCent10to30;
  int nRecoJtPtBinsLargeRCent30to50;
  int nGenJtPtSmallBinsLargeRCent30to50;
  int nGenJtPtLargeBinsLargeRCent30to50;
  int nRecoJtPtBinsLargeRCent50to90;
  int nGenJtPtSmallBinsLargeRCent50to90;
  int nGenJtPtLargeBinsLargeRCent50to90;

  std::vector<double> recoJtPtBinsSmallR;
  std::vector<double> genJtPtSmallBinsSmallR;
  std::vector<double> genJtPtLargeBinsSmallR;
  std::vector<double> recoJtPtBinsLargeR;
  std::vector<double> genJtPtSmallBinsLargeR;
  std::vector<double> genJtPtLargeBinsLargeR;

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

  bool CheckNSmallR(int inNSmallR);
  bool CheckNSmallR(cutPropagator cutProp);
  bool CheckSmallRVals(std::vector<int> inSmallRVals);
  bool CheckSmallRVals(cutPropagator cutProp);
  bool CheckNLargeR(int inNLargeR);
  bool CheckNLargeR(cutPropagator cutProp);
  bool CheckLargeRVals(std::vector<int> inLargeRVals);
  bool CheckLargeRVals(cutPropagator cutProp);

  bool CheckNBins(int inBins, int compBins, std::string binStr);

  bool CheckNRecoJtPtBinsSmallRCent0to10(int inNJtPtBinsSmallR);
  bool CheckNRecoJtPtBinsSmallRCent0to10(cutPropagator inCutProp);
  bool CheckNGenJtPtSmallBinsSmallRCent0to10(int inNJtPtSmallBinsSmallR);
  bool CheckNGenJtPtSmallBinsSmallRCent0to10(cutPropagator inCutProp);
  bool CheckNGenJtPtLargeBinsSmallRCent0to10(int inNJtPtLargeBinsSmallR);
  bool CheckNGenJtPtLargeBinsSmallRCent0to10(cutPropagator inCutProp);
  bool CheckNRecoJtPtBinsLargeRCent0to10(int inNJtPtBinsLargeR);
  bool CheckNRecoJtPtBinsLargeRCent0to10(cutPropagator inCutProp);
  bool CheckNGenJtPtSmallBinsLargeRCent0to10(int inNJtPtSmallBinsLargeR);
  bool CheckNGenJtPtSmallBinsLargeRCent0to10(cutPropagator inCutProp);
  bool CheckNGenJtPtLargeBinsLargeRCent0to10(int inNJtPtLargeBinsLargeR);
  bool CheckNGenJtPtLargeBinsLargeRCent0to10(cutPropagator inCutProp);

  bool CheckNRecoJtPtBinsSmallRCent10to30(int inNJtPtBinsSmallR);
  bool CheckNRecoJtPtBinsSmallRCent10to30(cutPropagator inCutProp);
  bool CheckNGenJtPtSmallBinsSmallRCent10to30(int inNJtPtSmallBinsSmallR);
  bool CheckNGenJtPtSmallBinsSmallRCent10to30(cutPropagator inCutProp);
  bool CheckNGenJtPtLargeBinsSmallRCent10to30(int inNJtPtLargeBinsSmallR);
  bool CheckNGenJtPtLargeBinsSmallRCent10to30(cutPropagator inCutProp);
  bool CheckNRecoJtPtBinsLargeRCent10to30(int inNJtPtBinsLargeR);
  bool CheckNRecoJtPtBinsLargeRCent10to30(cutPropagator inCutProp);
  bool CheckNGenJtPtSmallBinsLargeRCent10to30(int inNJtPtSmallBinsLargeR);
  bool CheckNGenJtPtSmallBinsLargeRCent10to30(cutPropagator inCutProp);
  bool CheckNGenJtPtLargeBinsLargeRCent10to30(int inNJtPtLargeBinsLargeR);
  bool CheckNGenJtPtLargeBinsLargeRCent10to30(cutPropagator inCutProp);

  bool CheckNRecoJtPtBinsSmallRCent30to50(int inNJtPtBinsSmallR);
  bool CheckNRecoJtPtBinsSmallRCent30to50(cutPropagator inCutProp);
  bool CheckNGenJtPtSmallBinsSmallRCent30to50(int inNJtPtSmallBinsSmallR);
  bool CheckNGenJtPtSmallBinsSmallRCent30to50(cutPropagator inCutProp);
  bool CheckNGenJtPtLargeBinsSmallRCent30to50(int inNJtPtLargeBinsSmallR);
  bool CheckNGenJtPtLargeBinsSmallRCent30to50(cutPropagator inCutProp);
  bool CheckNRecoJtPtBinsLargeRCent30to50(int inNJtPtBinsLargeR);
  bool CheckNRecoJtPtBinsLargeRCent30to50(cutPropagator inCutProp);
  bool CheckNGenJtPtSmallBinsLargeRCent30to50(int inNJtPtSmallBinsLargeR);
  bool CheckNGenJtPtSmallBinsLargeRCent30to50(cutPropagator inCutProp);
  bool CheckNGenJtPtLargeBinsLargeRCent30to50(int inNJtPtLargeBinsLargeR);
  bool CheckNGenJtPtLargeBinsLargeRCent30to50(cutPropagator inCutProp);

  bool CheckNRecoJtPtBinsSmallRCent50to90(int inNJtPtBinsSmallR);
  bool CheckNRecoJtPtBinsSmallRCent50to90(cutPropagator inCutProp);
  bool CheckNGenJtPtSmallBinsSmallRCent50to90(int inNJtPtSmallBinsSmallR);
  bool CheckNGenJtPtSmallBinsSmallRCent50to90(cutPropagator inCutProp);
  bool CheckNGenJtPtLargeBinsSmallRCent50to90(int inNJtPtLargeBinsSmallR);
  bool CheckNGenJtPtLargeBinsSmallRCent50to90(cutPropagator inCutProp);
  bool CheckNRecoJtPtBinsLargeRCent50to90(int inNJtPtBinsLargeR);
  bool CheckNRecoJtPtBinsLargeRCent50to90(cutPropagator inCutProp);
  bool CheckNGenJtPtSmallBinsLargeRCent50to90(int inNJtPtSmallBinsLargeR);
  bool CheckNGenJtPtSmallBinsLargeRCent50to90(cutPropagator inCutProp);
  bool CheckNGenJtPtLargeBinsLargeRCent50to90(int inNJtPtLargeBinsLargeR);
  bool CheckNGenJtPtLargeBinsLargeRCent50to90(cutPropagator inCutProp);

  bool CheckNJtAbsEtaBins(int inNJtAbsEtaBins);
  bool CheckNJtAbsEtaBins(cutPropagator inCutProp);
  bool CheckNID(int inNID);
  bool CheckNID(cutPropagator inCutProp);

  bool CheckRecoJtPtBinsSmallR(std::vector<double> inRecoJtPtBinsSmallR);
  bool CheckRecoJtPtBinsSmallR(cutPropagator inCutProp);
  bool CheckGenJtPtSmallBinsSmallR(std::vector<double> inGenJtPtSmallBinsSmallR);
  bool CheckGenJtPtSmallBinsSmallR(cutPropagator inCutProp);
  bool CheckGenJtPtLargeBinsSmallR(std::vector<double> inGenJtPtLargeBinsSmallR);
  bool CheckGenJtPtLargeBinsSmallR(cutPropagator inCutProp);
  bool CheckRecoJtPtBinsLargeR(std::vector<double> inRecoJtPtBinsLargeR);
  bool CheckRecoJtPtBinsLargeR(cutPropagator inCutProp);
  bool CheckGenJtPtSmallBinsLargeR(std::vector<double> inGenJtPtSmallBinsLargeR);
  bool CheckGenJtPtSmallBinsLargeR(cutPropagator inCutProp);
  bool CheckGenJtPtLargeBinsLargeR(std::vector<double> inGenJtPtLargeBinsLargeR);
  bool CheckGenJtPtLargeBinsLargeR(cutPropagator inCutProp);

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
  bool CheckPriorFlatFileName(std::string inPriorFlatFileName);
  bool CheckPriorFlatFileName(cutPropagator inCutProp);
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
  std::string GetPriorFlatFileName(){return flatPriorFileName;}

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

  int GetNSmallR(){return nSmallR;}
  std::vector<int> GetSmallRVals(){return smallRVals;}
  int GetNLargeR(){return nLargeR;}
  std::vector<int> GetLargeRVals(){return largeRVals;}

  int GetNRecoJtPtBinsSmallRCent0to10(){return nRecoJtPtBinsSmallRCent0to10;}
  int GetNGenJtPtSmallBinsSmallRCent0to10(){return nGenJtPtSmallBinsSmallRCent0to10;}
  int GetNGenJtPtLargeBinsSmallRCent0to10(){return nGenJtPtLargeBinsSmallRCent0to10;}
  int GetNRecoJtPtBinsLargeRCent0to10(){return nRecoJtPtBinsLargeRCent0to10;}
  int GetNGenJtPtSmallBinsLargeRCent0to10(){return nGenJtPtSmallBinsLargeRCent0to10;}
  int GetNGenJtPtLargeBinsLargeRCent0to10(){return nGenJtPtLargeBinsLargeRCent0to10;}

  int GetNRecoJtPtBinsSmallRCent10to30(){return nRecoJtPtBinsSmallRCent10to30;}
  int GetNGenJtPtSmallBinsSmallRCent10to30(){return nGenJtPtSmallBinsSmallRCent10to30;}
  int GetNGenJtPtLargeBinsSmallRCent10to30(){return nGenJtPtLargeBinsSmallRCent10to30;}
  int GetNRecoJtPtBinsLargeRCent10to30(){return nRecoJtPtBinsLargeRCent10to30;}
  int GetNGenJtPtSmallBinsLargeRCent10to30(){return nGenJtPtSmallBinsLargeRCent10to30;}
  int GetNGenJtPtLargeBinsLargeRCent10to30(){return nGenJtPtLargeBinsLargeRCent10to30;}

  int GetNRecoJtPtBinsSmallRCent30to50(){return nRecoJtPtBinsSmallRCent30to50;}
  int GetNGenJtPtSmallBinsSmallRCent30to50(){return nGenJtPtSmallBinsSmallRCent30to50;}
  int GetNGenJtPtLargeBinsSmallRCent30to50(){return nGenJtPtLargeBinsSmallRCent30to50;}
  int GetNRecoJtPtBinsLargeRCent30to50(){return nRecoJtPtBinsLargeRCent30to50;}
  int GetNGenJtPtSmallBinsLargeRCent30to50(){return nGenJtPtSmallBinsLargeRCent30to50;}
  int GetNGenJtPtLargeBinsLargeRCent30to50(){return nGenJtPtLargeBinsLargeRCent30to50;}

  int GetNRecoJtPtBinsSmallRCent50to90(){return nRecoJtPtBinsSmallRCent50to90;}
  int GetNGenJtPtSmallBinsSmallRCent50to90(){return nGenJtPtSmallBinsSmallRCent50to90;}
  int GetNGenJtPtLargeBinsSmallRCent50to90(){return nGenJtPtLargeBinsSmallRCent50to90;}
  int GetNRecoJtPtBinsLargeRCent50to90(){return nRecoJtPtBinsLargeRCent50to90;}
  int GetNGenJtPtSmallBinsLargeRCent50to90(){return nGenJtPtSmallBinsLargeRCent50to90;}
  int GetNGenJtPtLargeBinsLargeRCent50to90(){return nGenJtPtLargeBinsLargeRCent50to90;}

  std::vector<double> GetRecoJtPtBinsSmallR(){return recoJtPtBinsSmallR;}
  std::vector<double> GetGenJtPtSmallBinsSmallR(){return genJtPtSmallBinsSmallR;}
  std::vector<double> GetGenJtPtLargeBinsSmallR(){return genJtPtLargeBinsSmallR;}
  std::vector<double> GetRecoJtPtBinsLargeR(){return recoJtPtBinsLargeR;}
  std::vector<double> GetGenJtPtSmallBinsLargeR(){return genJtPtSmallBinsLargeR;}
  std::vector<double> GetGenJtPtLargeBinsLargeR(){return genJtPtLargeBinsLargeR;}

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
  void SetPriorFlatFileName(std::string inPriorFlatFileName){flatPriorFileName = inPriorFlatFileName; return;}
  
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

  void SetNSmallR(int inNSmallR){nSmallR = inNSmallR; return;}
  void SetSmallRVals(std::vector<int> inSmallRVals){smallRVals = inSmallRVals; return;}
  void SetSmallRVals(int inN, const Int_t inSmallRVals[]);

  void SetNLargeR(int inNLargeR){nLargeR = inNLargeR; return;}
  void SetLargeRVals(std::vector<int> inLargeRVals){largeRVals = inLargeRVals; return;}
  void SetLargeRVals(int inN, const Int_t inLargeRVals[]);

  void SetNRecoJtPtBinsSmallRCent0to10(int inNRecoJtPtBinsSmallRCent0to10){nRecoJtPtBinsSmallRCent0to10 = inNRecoJtPtBinsSmallRCent0to10; return;}
  void SetNGenJtPtSmallBinsSmallRCent0to10(int inNGenJtPtSmallBinsSmallRCent0to10){nGenJtPtSmallBinsSmallRCent0to10 = inNGenJtPtSmallBinsSmallRCent0to10; return;}
  void SetNGenJtPtLargeBinsSmallRCent0to10(int inNGenJtPtLargeBinsSmallRCent0to10){nGenJtPtLargeBinsSmallRCent0to10 = inNGenJtPtLargeBinsSmallRCent0to10; return;}
  void SetNRecoJtPtBinsLargeRCent0to10(int inNRecoJtPtBinsLargeRCent0to10){nRecoJtPtBinsLargeRCent0to10 = inNRecoJtPtBinsLargeRCent0to10; return;}
  void SetNGenJtPtSmallBinsLargeRCent0to10(int inNGenJtPtSmallBinsLargeRCent0to10){nGenJtPtSmallBinsLargeRCent0to10 = inNGenJtPtSmallBinsLargeRCent0to10; return;}
  void SetNGenJtPtLargeBinsLargeRCent0to10(int inNGenJtPtLargeBinsLargeRCent0to10){nGenJtPtLargeBinsLargeRCent0to10 = inNGenJtPtLargeBinsLargeRCent0to10; return;}

  void SetNRecoJtPtBinsSmallRCent10to30(int inNRecoJtPtBinsSmallRCent10to30){nRecoJtPtBinsSmallRCent10to30 = inNRecoJtPtBinsSmallRCent10to30; return;}
  void SetNGenJtPtSmallBinsSmallRCent10to30(int inNGenJtPtSmallBinsSmallRCent10to30){nGenJtPtSmallBinsSmallRCent10to30 = inNGenJtPtSmallBinsSmallRCent10to30; return;}
  void SetNGenJtPtLargeBinsSmallRCent10to30(int inNGenJtPtLargeBinsSmallRCent10to30){nGenJtPtLargeBinsSmallRCent10to30 = inNGenJtPtLargeBinsSmallRCent10to30; return;}
  void SetNRecoJtPtBinsLargeRCent10to30(int inNRecoJtPtBinsLargeRCent10to30){nRecoJtPtBinsLargeRCent10to30 = inNRecoJtPtBinsLargeRCent10to30; return;}
  void SetNGenJtPtSmallBinsLargeRCent10to30(int inNGenJtPtSmallBinsLargeRCent10to30){nGenJtPtSmallBinsLargeRCent10to30 = inNGenJtPtSmallBinsLargeRCent10to30; return;}
  void SetNGenJtPtLargeBinsLargeRCent10to30(int inNGenJtPtLargeBinsLargeRCent10to30){nGenJtPtLargeBinsLargeRCent10to30 = inNGenJtPtLargeBinsLargeRCent10to30; return;}

  void SetNRecoJtPtBinsSmallRCent30to50(int inNRecoJtPtBinsSmallRCent30to50){nRecoJtPtBinsSmallRCent30to50 = inNRecoJtPtBinsSmallRCent30to50; return;}
  void SetNGenJtPtSmallBinsSmallRCent30to50(int inNGenJtPtSmallBinsSmallRCent30to50){nGenJtPtSmallBinsSmallRCent30to50 = inNGenJtPtSmallBinsSmallRCent30to50; return;}
  void SetNGenJtPtLargeBinsSmallRCent30to50(int inNGenJtPtLargeBinsSmallRCent30to50){nGenJtPtLargeBinsSmallRCent30to50 = inNGenJtPtLargeBinsSmallRCent30to50; return;}
  void SetNRecoJtPtBinsLargeRCent30to50(int inNRecoJtPtBinsLargeRCent30to50){nRecoJtPtBinsLargeRCent30to50 = inNRecoJtPtBinsLargeRCent30to50; return;}
  void SetNGenJtPtSmallBinsLargeRCent30to50(int inNGenJtPtSmallBinsLargeRCent30to50){nGenJtPtSmallBinsLargeRCent30to50 = inNGenJtPtSmallBinsLargeRCent30to50; return;}
  void SetNGenJtPtLargeBinsLargeRCent30to50(int inNGenJtPtLargeBinsLargeRCent30to50){nGenJtPtLargeBinsLargeRCent30to50 = inNGenJtPtLargeBinsLargeRCent30to50; return;}

  void SetNRecoJtPtBinsSmallRCent50to90(int inNRecoJtPtBinsSmallRCent50to90){nRecoJtPtBinsSmallRCent50to90 = inNRecoJtPtBinsSmallRCent50to90; return;}
  void SetNGenJtPtSmallBinsSmallRCent50to90(int inNGenJtPtSmallBinsSmallRCent50to90){nGenJtPtSmallBinsSmallRCent50to90 = inNGenJtPtSmallBinsSmallRCent50to90; return;}
  void SetNGenJtPtLargeBinsSmallRCent50to90(int inNGenJtPtLargeBinsSmallRCent50to90){nGenJtPtLargeBinsSmallRCent50to90 = inNGenJtPtLargeBinsSmallRCent50to90; return;}
  void SetNRecoJtPtBinsLargeRCent50to90(int inNRecoJtPtBinsLargeRCent50to90){nRecoJtPtBinsLargeRCent50to90 = inNRecoJtPtBinsLargeRCent50to90; return;}
  void SetNGenJtPtSmallBinsLargeRCent50to90(int inNGenJtPtSmallBinsLargeRCent50to90){nGenJtPtSmallBinsLargeRCent50to90 = inNGenJtPtSmallBinsLargeRCent50to90; return;}
  void SetNGenJtPtLargeBinsLargeRCent50to90(int inNGenJtPtLargeBinsLargeRCent50to90){nGenJtPtLargeBinsLargeRCent50to90 = inNGenJtPtLargeBinsLargeRCent50to90; return;}

  void SetRecoJtPtBinsSmallR(std::vector<double> inRecoJtPtBinsSmallR){recoJtPtBinsSmallR = inRecoJtPtBinsSmallR; return;}
  void SetRecoJtPtBinsSmallR(int inN, const Double_t inRecoJtPtBinsSmallR[]);
  void SetGenJtPtSmallBinsSmallR(std::vector<double> inGenJtPtSmallBinsSmallR){genJtPtSmallBinsSmallR = inGenJtPtSmallBinsSmallR; return;}
  void SetGenJtPtSmallBinsSmallR(int inN, const Double_t inGenJtPtSmallBinsSmallR[]);
  void SetGenJtPtLargeBinsSmallR(std::vector<double> inGenJtPtLargeBinsSmallR){genJtPtLargeBinsSmallR = inGenJtPtLargeBinsSmallR; return;}
  void SetGenJtPtLargeBinsSmallR(int inN, const Double_t inGenJtPtLargeBinsSmallR[]);

  void SetRecoJtPtBinsLargeR(std::vector<double> inRecoJtPtBinsLargeR){recoJtPtBinsLargeR = inRecoJtPtBinsLargeR; return;}
  void SetRecoJtPtBinsLargeR(int inN, const Double_t inRecoJtPtBinsLargeR[]);
  void SetGenJtPtSmallBinsLargeR(std::vector<double> inGenJtPtSmallBinsLargeR){genJtPtSmallBinsLargeR = inGenJtPtSmallBinsLargeR; return;}
  void SetGenJtPtSmallBinsLargeR(int inN, const Double_t inGenJtPtSmallBinsLargeR[]);
  void SetGenJtPtLargeBinsLargeR(std::vector<double> inGenJtPtLargeBinsLargeR){genJtPtLargeBinsLargeR = inGenJtPtLargeBinsLargeR; return;}
  void SetGenJtPtLargeBinsLargeR(int inN, const Double_t inGenJtPtLargeBinsLargeR[]);

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

  nSmallR = -1;
  smallRVals.clear();

  nLargeR = -1;
  largeRVals.clear();

  nRecoJtPtBinsSmallRCent0to10 = -1;
  nGenJtPtSmallBinsSmallRCent0to10 = -1;
  nGenJtPtLargeBinsSmallRCent0to10 = -1;
  nRecoJtPtBinsLargeRCent0to10 = -1;
  nGenJtPtSmallBinsLargeRCent0to10 = -1;
  nGenJtPtLargeBinsLargeRCent0to10 = -1;

  nRecoJtPtBinsSmallRCent10to30 = -1;
  nGenJtPtSmallBinsSmallRCent10to30 = -1;
  nGenJtPtLargeBinsSmallRCent10to30 = -1;
  nRecoJtPtBinsLargeRCent10to30 = -1;
  nGenJtPtSmallBinsLargeRCent10to30 = -1;
  nGenJtPtLargeBinsLargeRCent10to30 = -1;

  nRecoJtPtBinsSmallRCent30to50 = -1;
  nGenJtPtSmallBinsSmallRCent30to50 = -1;
  nGenJtPtLargeBinsSmallRCent30to50 = -1;
  nRecoJtPtBinsLargeRCent30to50 = -1;
  nGenJtPtSmallBinsLargeRCent30to50 = -1;
  nGenJtPtLargeBinsLargeRCent30to50 = -1;

  nRecoJtPtBinsSmallRCent50to90 = -1;
  nGenJtPtSmallBinsSmallRCent50to90 = -1;
  nGenJtPtLargeBinsSmallRCent50to90 = -1;
  nRecoJtPtBinsLargeRCent50to90 = -1;
  nGenJtPtSmallBinsLargeRCent50to90 = -1;
  nGenJtPtLargeBinsLargeRCent50to90 = -1;

  recoJtPtBinsSmallR.clear();
  genJtPtSmallBinsSmallR.clear();
  genJtPtLargeBinsSmallR.clear();
  recoJtPtBinsLargeR.clear();
  genJtPtSmallBinsLargeR.clear();
  genJtPtLargeBinsLargeR.clear();

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
    else if(isStrSame("nSmallR", tempStr)) nSmallR = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("smallRVals", tempStr)) smallRVals = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nLargeR", tempStr)) nLargeR = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("largeRVals", tempStr)) largeRVals = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nRecoJtPtBinsSmallRCent0to10", tempStr)) nRecoJtPtBinsSmallRCent0to10 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtSmallBinsSmallRCent0to10", tempStr)) nGenJtPtSmallBinsSmallRCent0to10 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtLargeBinsSmallRCent0to10", tempStr)) nGenJtPtLargeBinsSmallRCent0to10 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nRecoJtPtBinsLargeRCent0to10", tempStr)) nRecoJtPtBinsLargeRCent0to10 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtSmallBinsLargeRCent0to10", tempStr)) nGenJtPtSmallBinsLargeRCent0to10 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtLargeBinsLargeRCent0to10", tempStr)) nGenJtPtLargeBinsLargeRCent0to10 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nRecoJtPtBinsSmallRCent10to30", tempStr)) nRecoJtPtBinsSmallRCent10to30 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtSmallBinsSmallRCent10to30", tempStr)) nGenJtPtSmallBinsSmallRCent10to30 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtLargeBinsSmallRCent10to30", tempStr)) nGenJtPtLargeBinsSmallRCent10to30 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nRecoJtPtBinsLargeRCent10to30", tempStr)) nRecoJtPtBinsLargeRCent10to30 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtSmallBinsLargeRCent10to30", tempStr)) nGenJtPtSmallBinsLargeRCent10to30 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtLargeBinsLargeRCent10to30", tempStr)) nGenJtPtLargeBinsLargeRCent10to30 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nRecoJtPtBinsSmallRCent30to50", tempStr)) nRecoJtPtBinsSmallRCent30to50 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtSmallBinsSmallRCent30to50", tempStr)) nGenJtPtSmallBinsSmallRCent30to50 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtLargeBinsSmallRCent30to50", tempStr)) nGenJtPtLargeBinsSmallRCent30to50 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nRecoJtPtBinsLargeRCent30to50", tempStr)) nRecoJtPtBinsLargeRCent30to50 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtSmallBinsLargeRCent30to50", tempStr)) nGenJtPtSmallBinsLargeRCent30to50 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtLargeBinsLargeRCent30to50", tempStr)) nGenJtPtLargeBinsLargeRCent30to50 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nRecoJtPtBinsSmallRCent50to90", tempStr)) nRecoJtPtBinsSmallRCent50to90 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtSmallBinsSmallRCent50to90", tempStr)) nGenJtPtSmallBinsSmallRCent50to90 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtLargeBinsSmallRCent50to90", tempStr)) nGenJtPtLargeBinsSmallRCent50to90 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nRecoJtPtBinsLargeRCent50to90", tempStr)) nRecoJtPtBinsLargeRCent50to90 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtSmallBinsLargeRCent50to90", tempStr)) nGenJtPtSmallBinsLargeRCent50to90 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGenJtPtLargeBinsLargeRCent50to90", tempStr)) nGenJtPtLargeBinsLargeRCent50to90 = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
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
    else if(isStrSame("recoJtPtBinsSmallR", tempStr)) recoJtPtBinsSmallR = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtSmallBinsSmallR", tempStr)) genJtPtSmallBinsSmallR = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtLargeBinsSmallR", tempStr)) genJtPtLargeBinsSmallR = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsLargeR", tempStr)) recoJtPtBinsLargeR = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtSmallBinsLargeR", tempStr)) genJtPtSmallBinsLargeR = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtLargeBinsLargeR", tempStr)) genJtPtLargeBinsLargeR = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
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

  std::string smallRValsStr = "";
  for(Int_t rI = 0; rI < nSmallR; ++rI){
    smallRValsStr = smallRValsStr + std::to_string(smallRVals.at(rI)) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string largeRValsStr = "";
  for(Int_t rI = 0; rI < nLargeR; ++rI){
    largeRValsStr = largeRValsStr + std::to_string(largeRVals.at(rI)) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string recoJtPtBinsSmallRStr = "";
  for(unsigned int jI = 0; jI < recoJtPtBinsSmallR.size(); ++jI){
    recoJtPtBinsSmallRStr = recoJtPtBinsSmallRStr + std::to_string(recoJtPtBinsSmallR.at(jI)) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string recoJtPtBinsLargeRStr = "";
  for(unsigned int jI = 0; jI < recoJtPtBinsLargeR.size(); ++jI){
    recoJtPtBinsLargeRStr = recoJtPtBinsLargeRStr + std::to_string(recoJtPtBinsLargeR.at(jI)) + ",";
  }

  std::string genJtPtSmallBinsSmallRStr = "";
  for(unsigned int jI = 0; jI < genJtPtSmallBinsSmallR.size(); ++jI){
    genJtPtSmallBinsSmallRStr = genJtPtSmallBinsSmallRStr + std::to_string(genJtPtSmallBinsSmallR.at(jI)) + ",";
  }
  std::string genJtPtLargeBinsSmallRStr = "";
  for(unsigned int jI = 0; jI < genJtPtLargeBinsSmallR.size(); ++jI){
    genJtPtLargeBinsSmallRStr = genJtPtLargeBinsSmallRStr + std::to_string(genJtPtLargeBinsSmallR.at(jI)) + ",";
  }

  std::string genJtPtSmallBinsLargeRStr = "";
  for(unsigned int jI = 0; jI < genJtPtSmallBinsLargeR.size(); ++jI){
    genJtPtSmallBinsLargeRStr = genJtPtSmallBinsLargeRStr + std::to_string(genJtPtSmallBinsLargeR.at(jI)) + ",";
  }
  std::string genJtPtLargeBinsLargeRStr = "";
  for(unsigned int jI = 0; jI < genJtPtLargeBinsLargeR.size(); ++jI){
    genJtPtLargeBinsLargeRStr = genJtPtLargeBinsLargeRStr + std::to_string(genJtPtLargeBinsLargeR.at(jI)) + ",";
  }

  std::string jtAlgosStr = "";
  std::string minJtPtCutStr = "";
  std::string multiJtPtCutStr = "";
  std::string recoTruncPosStr = "";

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

  std::vector<TNamed*> nameList;  
  nameList.push_back(new TNamed("inFileNames", inFileNames2.c_str()));
  nameList.push_back(new TNamed("isPP", std::to_string(isPP).c_str()));
  nameList.push_back(new TNamed("jtAbsEtaMax", std::to_string(jtAbsEtaMax).c_str()));
  nameList.push_back(new TNamed("rcDiffFileName", rcDiffFileName.c_str()));
  nameList.push_back(new TNamed("flatPriorFileName", flatPriorFileName.c_str()));
  nameList.push_back(new TNamed("jecVarMC", std::to_string(jecVarMC).c_str()));
  nameList.push_back(new TNamed("jerVarMC", std::to_string(jerVarMC).c_str()));
  nameList.push_back(new TNamed("jecVarData", std::to_string(jecVarData).c_str()));
  nameList.push_back(new TNamed("nResponseMod", std::to_string(nResponseMod).c_str()));
  nameList.push_back(new TNamed("responseMod", responseModStr.c_str()));
  nameList.push_back(new TNamed("jerVarData", jerVarDataStr.c_str()));
  nameList.push_back(new TNamed("nJtAlgos", std::to_string(nJtAlgos).c_str()));
  nameList.push_back(new TNamed("jtAlgos", jtAlgosStr.c_str()));
  nameList.push_back(new TNamed("minJtPtCut", minJtPtCutStr.c_str()));
  nameList.push_back(new TNamed("multiJtPtCut", multiJtPtCutStr.c_str()));
  nameList.push_back(new TNamed("recoTruncPos", recoTruncPosStr.c_str()));
  nameList.push_back(new TNamed("nSmallR", std::to_string(nSmallR).c_str()));
  nameList.push_back(new TNamed("smallRVals", smallRValsStr.c_str()));
  nameList.push_back(new TNamed("nLargeR", std::to_string(nLargeR).c_str()));
  nameList.push_back(new TNamed("largeRVals", largeRValsStr.c_str()));

  nameList.push_back(new TNamed("nRecoJtPtBinsSmallRCent0to10", std::to_string(nRecoJtPtBinsSmallRCent0to10).c_str()));
  nameList.push_back(new TNamed("nGenJtPtSmallBinsSmallRCent0to10", std::to_string(nGenJtPtSmallBinsSmallRCent0to10).c_str()));
  nameList.push_back(new TNamed("nGenJtPtLargeBinsSmallRCent0to10", std::to_string(nGenJtPtLargeBinsSmallRCent0to10).c_str()));
  nameList.push_back(new TNamed("nRecoJtPtBinsLargeRCent0to10", std::to_string(nRecoJtPtBinsLargeRCent0to10).c_str()));
  nameList.push_back(new TNamed("nGenJtPtSmallBinsLargeRCent0to10", std::to_string(nGenJtPtSmallBinsLargeRCent0to10).c_str()));
  nameList.push_back(new TNamed("nGenJtPtLargeBinsLargeRCent0to10", std::to_string(nGenJtPtLargeBinsLargeRCent0to10).c_str()));

  nameList.push_back(new TNamed("nRecoJtPtBinsSmallRCent10to30", std::to_string(nRecoJtPtBinsSmallRCent10to30).c_str()));
  nameList.push_back(new TNamed("nGenJtPtSmallBinsSmallRCent10to30", std::to_string(nGenJtPtSmallBinsSmallRCent10to30).c_str()));
  nameList.push_back(new TNamed("nGenJtPtLargeBinsSmallRCent10to30", std::to_string(nGenJtPtLargeBinsSmallRCent10to30).c_str()));
  nameList.push_back(new TNamed("nRecoJtPtBinsLargeRCent10to30", std::to_string(nRecoJtPtBinsLargeRCent10to30).c_str()));
  nameList.push_back(new TNamed("nGenJtPtSmallBinsLargeRCent10to30", std::to_string(nGenJtPtSmallBinsLargeRCent10to30).c_str()));
  nameList.push_back(new TNamed("nGenJtPtLargeBinsLargeRCent10to30", std::to_string(nGenJtPtLargeBinsLargeRCent10to30).c_str()));

  nameList.push_back(new TNamed("nRecoJtPtBinsSmallRCent30to50", std::to_string(nRecoJtPtBinsSmallRCent30to50).c_str()));
  nameList.push_back(new TNamed("nGenJtPtSmallBinsSmallRCent30to50", std::to_string(nGenJtPtSmallBinsSmallRCent30to50).c_str()));
  nameList.push_back(new TNamed("nGenJtPtLargeBinsSmallRCent30to50", std::to_string(nGenJtPtLargeBinsSmallRCent30to50).c_str()));
  nameList.push_back(new TNamed("nRecoJtPtBinsLargeRCent30to50", std::to_string(nRecoJtPtBinsLargeRCent30to50).c_str()));
  nameList.push_back(new TNamed("nGenJtPtSmallBinsLargeRCent30to50", std::to_string(nGenJtPtSmallBinsLargeRCent30to50).c_str()));
  nameList.push_back(new TNamed("nGenJtPtLargeBinsLargeRCent30to50", std::to_string(nGenJtPtLargeBinsLargeRCent30to50).c_str()));

  nameList.push_back(new TNamed("nRecoJtPtBinsSmallRCent50to90", std::to_string(nRecoJtPtBinsSmallRCent50to90).c_str()));
  nameList.push_back(new TNamed("nGenJtPtSmallBinsSmallRCent50to90", std::to_string(nGenJtPtSmallBinsSmallRCent50to90).c_str()));
  nameList.push_back(new TNamed("nGenJtPtLargeBinsSmallRCent50to90", std::to_string(nGenJtPtLargeBinsSmallRCent50to90).c_str()));
  nameList.push_back(new TNamed("nRecoJtPtBinsLargeRCent50to90", std::to_string(nRecoJtPtBinsLargeRCent50to90).c_str()));
  nameList.push_back(new TNamed("nGenJtPtSmallBinsLargeRCent50to90", std::to_string(nGenJtPtSmallBinsLargeRCent50to90).c_str()));
  nameList.push_back(new TNamed("nGenJtPtLargeBinsLargeRCent50to90", std::to_string(nGenJtPtLargeBinsLargeRCent50to90).c_str()));

  nameList.push_back(new TNamed("recoJtPtBinsSmallR", recoJtPtBinsSmallRStr.c_str()));
  nameList.push_back(new TNamed("genJtPtSmallBinsSmallR", genJtPtSmallBinsSmallRStr.c_str()));
  nameList.push_back(new TNamed("genJtPtLargeBinsSmallR", genJtPtLargeBinsSmallRStr.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsLargeR", recoJtPtBinsLargeRStr.c_str()));
  nameList.push_back(new TNamed("genJtPtSmallBinsLargeR", genJtPtSmallBinsLargeRStr.c_str()));
  nameList.push_back(new TNamed("genJtPtLargeBinsLargeR", genJtPtLargeBinsLargeRStr.c_str()));

  nameList.push_back(new TNamed("nJtAbsEtaBins", std::to_string(nJtAbsEtaBins).c_str()));
  nameList.push_back(new TNamed("jtAbsEtaBinsLow", jtAbsEtaBinsLowStr.c_str()));
  nameList.push_back(new TNamed("jtAbsEtaBinsHi", jtAbsEtaBinsHiStr.c_str()));
  nameList.push_back(new TNamed("nPthats", nPthatsStr.c_str()));
  nameList.push_back(new TNamed("pthats", pthatsStr.c_str()));
  nameList.push_back(new TNamed("pthatWeights", pthatWeightsStr.c_str()));
  nameList.push_back(new TNamed("nCentBins", std::to_string(nCentBins).c_str()));
  nameList.push_back(new TNamed("centBinsLow", centBinsLowStr.c_str()));
  nameList.push_back(new TNamed("centBinsHi", centBinsHiStr.c_str()));
  nameList.push_back(new TNamed("nID", nIDStr.c_str()));
  nameList.push_back(new TNamed("idStr", idStr2.c_str()));
  nameList.push_back(new TNamed("jtPfCHMFCutLow", jtPfCHMFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfCHMFCutHi", jtPfCHMFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfMUMFCutLow", jtPfMUMFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfMUMFCutHi", jtPfMUMFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfNHFCutLow", jtPfNHFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfNHFCutHi", jtPfNHFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfNEFCutLow", jtPfNEFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfNEFCutHi", jtPfNEFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfMUFCutLow", jtPfMUFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfMUFCutHi", jtPfMUFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfCHFCutLow", jtPfCHFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfCHFCutHi", jtPfCHFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfCEFCutLow", jtPfCEFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfCEFCutHi", jtPfCEFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfMinMult", jtPfMinMultStr.c_str()));
  nameList.push_back(new TNamed("jtPfMinChgMult", jtPfMinChgMultStr.c_str()));
  nameList.push_back(new TNamed("nSyst", nSystStr.c_str()));
  nameList.push_back(new TNamed("systStr", systStr2.c_str()));
  nameList.push_back(new TNamed("nBayes", nBayesStr.c_str()));
  nameList.push_back(new TNamed("nBigBayesSymm", nBigBayesSymmStr.c_str()));
  nameList.push_back(new TNamed("nSuperBayes", nSuperBayesStr.c_str()));
  nameList.push_back(new TNamed("bayesVal", bayesVal2.c_str()));
  nameList.push_back(new TNamed("nHistDim", nHistDimStr.c_str()));
  
  for(auto const & name : nameList){
    name->Write("", TObject::kOverwrite);
  }

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
    //    nHistDimName.Write("", TObject::kOverwrite);
    
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
  if(!CheckNSmallR(inCutProp)) return false;
  if(!CheckSmallRVals(inCutProp)) return false;
  if(!CheckNLargeR(inCutProp)) return false;
  if(!CheckLargeRVals(inCutProp)) return false;

  if(!CheckNRecoJtPtBinsSmallRCent0to10(inCutProp)) return false;
  if(!CheckNGenJtPtSmallBinsSmallRCent0to10(inCutProp)) return false;
  if(!CheckNGenJtPtLargeBinsSmallRCent0to10(inCutProp)) return false;
  if(!CheckNRecoJtPtBinsLargeRCent0to10(inCutProp)) return false;
  if(!CheckNGenJtPtSmallBinsLargeRCent0to10(inCutProp)) return false;
  if(!CheckNGenJtPtLargeBinsLargeRCent0to10(inCutProp)) return false;

  if(!CheckNRecoJtPtBinsSmallRCent10to30(inCutProp)) return false;
  if(!CheckNGenJtPtSmallBinsSmallRCent10to30(inCutProp)) return false;
  if(!CheckNGenJtPtLargeBinsSmallRCent10to30(inCutProp)) return false;
  if(!CheckNRecoJtPtBinsLargeRCent10to30(inCutProp)) return false;
  if(!CheckNGenJtPtSmallBinsLargeRCent10to30(inCutProp)) return false;
  if(!CheckNGenJtPtLargeBinsLargeRCent10to30(inCutProp)) return false;

  if(!CheckNRecoJtPtBinsSmallRCent30to50(inCutProp)) return false;
  if(!CheckNGenJtPtSmallBinsSmallRCent30to50(inCutProp)) return false;
  if(!CheckNGenJtPtLargeBinsSmallRCent30to50(inCutProp)) return false;
  if(!CheckNRecoJtPtBinsLargeRCent30to50(inCutProp)) return false;
  if(!CheckNGenJtPtSmallBinsLargeRCent30to50(inCutProp)) return false;
  if(!CheckNGenJtPtLargeBinsLargeRCent30to50(inCutProp)) return false;

  if(!CheckNRecoJtPtBinsSmallRCent50to90(inCutProp)) return false;
  if(!CheckNGenJtPtSmallBinsSmallRCent50to90(inCutProp)) return false;
  if(!CheckNGenJtPtLargeBinsSmallRCent50to90(inCutProp)) return false;
  if(!CheckNRecoJtPtBinsLargeRCent50to90(inCutProp)) return false;
  if(!CheckNGenJtPtSmallBinsLargeRCent50to90(inCutProp)) return false;
  if(!CheckNGenJtPtLargeBinsLargeRCent50to90(inCutProp)) return false;

  if(!CheckNJtAbsEtaBins(inCutProp)) return false;
  if(!CheckNID(inCutProp)) return false;

  if(!CheckRecoJtPtBinsSmallR(inCutProp)) return false;
  if(!CheckGenJtPtSmallBinsSmallR(inCutProp)) return false;
  if(!CheckGenJtPtLargeBinsSmallR(inCutProp)) return false;
  if(!CheckRecoJtPtBinsLargeR(inCutProp)) return false;
  if(!CheckGenJtPtSmallBinsLargeR(inCutProp)) return false;
  if(!CheckGenJtPtLargeBinsLargeR(inCutProp)) return false;

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
    if(!CheckPriorFlatFileName(inCutProp)) return false;
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

//Series of checks
bool cutPropagator::CheckJtAbsEtaMax(double inJtAbsEtaMax)
{
  bool checkVal = CheckDouble(jtAbsEtaMax, inJtAbsEtaMax);
  if(!checkVal) std::cout << "cutPropagator check failed on jtAbsEtaMax" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtAbsEtaMax(cutPropagator inCutProp){return CheckJtAbsEtaMax(inCutProp.GetJtAbsEtaMax());}

bool cutPropagator::CheckNSmallR(int inNSmallR)
{
  bool checkVal = CheckInt(nSmallR, inNSmallR);
  if(!checkVal) std::cout << "cutPropagator check failed on nSmallR" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNSmallR(cutPropagator inCutProp){return CheckNSmallR(inCutProp.GetNSmallR());}

bool cutPropagator::CheckNLargeR(int inNLargeR)
{
  bool checkVal = CheckInt(nLargeR, inNLargeR);
  if(!checkVal) std::cout << "cutPropagator check failed on nLargeR" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNLargeR(cutPropagator inCutProp){return CheckNLargeR(inCutProp.GetNLargeR());}

bool cutPropagator::CheckSmallRVals(std::vector<int> inSmallRVals)
{
  bool checkVal = CheckVectInt(smallRVals, inSmallRVals);
  if(!checkVal) std::cout << "cutPropagator check failed on smallRVals" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckSmallRVals(cutPropagator inCutProp){return CheckSmallRVals(inCutProp.GetSmallRVals());}

bool cutPropagator::CheckLargeRVals(std::vector<int> inLargeRVals)
{
  bool checkVal = CheckVectInt(largeRVals, inLargeRVals);
  if(!checkVal) std::cout << "cutPropagator check failed on largeRVals" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckLargeRVals(cutPropagator inCutProp){return CheckLargeRVals(inCutProp.GetLargeRVals());}

bool cutPropagator::CheckNBins(int inBins, int compBins, std::string binStr)
{
  bool checkVal = CheckInt(compBins, inBins);
  if(!checkVal) std::cout << "cutPropagator check failed on " << binStr << std::endl;
  return checkVal;
}

bool cutPropagator::CheckNRecoJtPtBinsSmallRCent0to10(int inNRecoJtPtBinsSmallRCent0to10){return CheckNBins(inNRecoJtPtBinsSmallRCent0to10, nRecoJtPtBinsSmallRCent0to10, "nRecoJtPtBinsSmallRCent0to10");}
bool cutPropagator::CheckNRecoJtPtBinsSmallRCent0to10(cutPropagator inCutProp){return CheckNRecoJtPtBinsSmallRCent0to10(inCutProp.GetNRecoJtPtBinsSmallRCent0to10());}
bool cutPropagator::CheckNRecoJtPtBinsSmallRCent10to30(int inNRecoJtPtBinsSmallRCent10to30){return CheckNBins(inNRecoJtPtBinsSmallRCent10to30, nRecoJtPtBinsSmallRCent10to30, "nRecoJtPtBinsSmallRCent10to30");}
bool cutPropagator::CheckNRecoJtPtBinsSmallRCent10to30(cutPropagator inCutProp){return CheckNRecoJtPtBinsSmallRCent10to30(inCutProp.GetNRecoJtPtBinsSmallRCent10to30());}
bool cutPropagator::CheckNRecoJtPtBinsSmallRCent30to50(int inNRecoJtPtBinsSmallRCent30to50){return CheckNBins(inNRecoJtPtBinsSmallRCent30to50, nRecoJtPtBinsSmallRCent30to50, "nRecoJtPtBinsSmallRCent30to50");}
bool cutPropagator::CheckNRecoJtPtBinsSmallRCent30to50(cutPropagator inCutProp){return CheckNRecoJtPtBinsSmallRCent30to50(inCutProp.GetNRecoJtPtBinsSmallRCent30to50());}
bool cutPropagator::CheckNRecoJtPtBinsSmallRCent50to90(int inNRecoJtPtBinsSmallRCent50to90){return CheckNBins(inNRecoJtPtBinsSmallRCent50to90, nRecoJtPtBinsSmallRCent50to90, "nRecoJtPtBinsSmallRCent50to90");}
bool cutPropagator::CheckNRecoJtPtBinsSmallRCent50to90(cutPropagator inCutProp){return CheckNRecoJtPtBinsSmallRCent50to90(inCutProp.GetNRecoJtPtBinsSmallRCent50to90());}

bool cutPropagator::CheckNGenJtPtSmallBinsSmallRCent0to10(int inNGenJtPtSmallBinsSmallRCent0to10){return CheckNBins(inNGenJtPtSmallBinsSmallRCent0to10, nGenJtPtSmallBinsSmallRCent0to10, "nGenJtPtSmallBinsSmallRCent0to10");}
bool cutPropagator::CheckNGenJtPtSmallBinsSmallRCent0to10(cutPropagator inCutProp){return CheckNGenJtPtSmallBinsSmallRCent0to10(inCutProp.GetNGenJtPtSmallBinsSmallRCent0to10());}
bool cutPropagator::CheckNGenJtPtSmallBinsSmallRCent10to30(int inNGenJtPtSmallBinsSmallRCent10to30){return CheckNBins(inNGenJtPtSmallBinsSmallRCent10to30, nGenJtPtSmallBinsSmallRCent10to30, "nGenJtPtSmallBinsSmallRCent10to30");}
bool cutPropagator::CheckNGenJtPtSmallBinsSmallRCent10to30(cutPropagator inCutProp){return CheckNGenJtPtSmallBinsSmallRCent10to30(inCutProp.GetNGenJtPtSmallBinsSmallRCent10to30());}
bool cutPropagator::CheckNGenJtPtSmallBinsSmallRCent30to50(int inNGenJtPtSmallBinsSmallRCent30to50){return CheckNBins(inNGenJtPtSmallBinsSmallRCent30to50, nGenJtPtSmallBinsSmallRCent30to50, "nGenJtPtSmallBinsSmallRCent30to50");}
bool cutPropagator::CheckNGenJtPtSmallBinsSmallRCent30to50(cutPropagator inCutProp){return CheckNGenJtPtSmallBinsSmallRCent30to50(inCutProp.GetNGenJtPtSmallBinsSmallRCent30to50());}
bool cutPropagator::CheckNGenJtPtSmallBinsSmallRCent50to90(int inNGenJtPtSmallBinsSmallRCent50to90){return CheckNBins(inNGenJtPtSmallBinsSmallRCent50to90, nGenJtPtSmallBinsSmallRCent50to90, "nGenJtPtSmallBinsSmallRCent50to90");}
bool cutPropagator::CheckNGenJtPtSmallBinsSmallRCent50to90(cutPropagator inCutProp){return CheckNGenJtPtSmallBinsSmallRCent50to90(inCutProp.GetNGenJtPtSmallBinsSmallRCent50to90());}

bool cutPropagator::CheckNGenJtPtLargeBinsSmallRCent0to10(int inNGenJtPtLargeBinsSmallRCent0to10){return CheckNBins(inNGenJtPtLargeBinsSmallRCent0to10, nGenJtPtLargeBinsSmallRCent0to10, "nGenJtPtLargeBinsSmallRCent0to10");}
bool cutPropagator::CheckNGenJtPtLargeBinsSmallRCent0to10(cutPropagator inCutProp){return CheckNGenJtPtLargeBinsSmallRCent0to10(inCutProp.GetNGenJtPtLargeBinsSmallRCent0to10());}
bool cutPropagator::CheckNGenJtPtLargeBinsSmallRCent10to30(int inNGenJtPtLargeBinsSmallRCent10to30){return CheckNBins(inNGenJtPtLargeBinsSmallRCent10to30, nGenJtPtLargeBinsSmallRCent10to30, "nGenJtPtLargeBinsSmallRCent10to30");}
bool cutPropagator::CheckNGenJtPtLargeBinsSmallRCent10to30(cutPropagator inCutProp){return CheckNGenJtPtLargeBinsSmallRCent10to30(inCutProp.GetNGenJtPtLargeBinsSmallRCent10to30());}
bool cutPropagator::CheckNGenJtPtLargeBinsSmallRCent30to50(int inNGenJtPtLargeBinsSmallRCent30to50){return CheckNBins(inNGenJtPtLargeBinsSmallRCent30to50, nGenJtPtLargeBinsSmallRCent30to50, "nGenJtPtLargeBinsSmallRCent30to50");}
bool cutPropagator::CheckNGenJtPtLargeBinsSmallRCent30to50(cutPropagator inCutProp){return CheckNGenJtPtLargeBinsSmallRCent30to50(inCutProp.GetNGenJtPtLargeBinsSmallRCent30to50());}
bool cutPropagator::CheckNGenJtPtLargeBinsSmallRCent50to90(int inNGenJtPtLargeBinsSmallRCent50to90){return CheckNBins(inNGenJtPtLargeBinsSmallRCent50to90, nGenJtPtLargeBinsSmallRCent50to90, "nGenJtPtLargeBinsSmallRCent50to90");}
bool cutPropagator::CheckNGenJtPtLargeBinsSmallRCent50to90(cutPropagator inCutProp){return CheckNGenJtPtLargeBinsSmallRCent50to90(inCutProp.GetNGenJtPtLargeBinsSmallRCent50to90());}


bool cutPropagator::CheckNRecoJtPtBinsLargeRCent0to10(int inNRecoJtPtBinsLargeRCent0to10){return CheckNBins(inNRecoJtPtBinsLargeRCent0to10, nRecoJtPtBinsLargeRCent0to10, "nRecoJtPtBinsLargeRCent0to10");}
bool cutPropagator::CheckNRecoJtPtBinsLargeRCent0to10(cutPropagator inCutProp){return CheckNRecoJtPtBinsLargeRCent0to10(inCutProp.GetNRecoJtPtBinsLargeRCent0to10());}
bool cutPropagator::CheckNRecoJtPtBinsLargeRCent10to30(int inNRecoJtPtBinsLargeRCent10to30){return CheckNBins(inNRecoJtPtBinsLargeRCent10to30, nRecoJtPtBinsLargeRCent10to30, "nRecoJtPtBinsLargeRCent10to30");}
bool cutPropagator::CheckNRecoJtPtBinsLargeRCent10to30(cutPropagator inCutProp){return CheckNRecoJtPtBinsLargeRCent10to30(inCutProp.GetNRecoJtPtBinsLargeRCent10to30());}
bool cutPropagator::CheckNRecoJtPtBinsLargeRCent30to50(int inNRecoJtPtBinsLargeRCent30to50){return CheckNBins(inNRecoJtPtBinsLargeRCent30to50, nRecoJtPtBinsLargeRCent30to50, "nRecoJtPtBinsLargeRCent30to50");}
bool cutPropagator::CheckNRecoJtPtBinsLargeRCent30to50(cutPropagator inCutProp){return CheckNRecoJtPtBinsLargeRCent30to50(inCutProp.GetNRecoJtPtBinsLargeRCent30to50());}
bool cutPropagator::CheckNRecoJtPtBinsLargeRCent50to90(int inNRecoJtPtBinsLargeRCent50to90){return CheckNBins(inNRecoJtPtBinsLargeRCent50to90, nRecoJtPtBinsLargeRCent50to90, "nRecoJtPtBinsLargeRCent50to90");}
bool cutPropagator::CheckNRecoJtPtBinsLargeRCent50to90(cutPropagator inCutProp){return CheckNRecoJtPtBinsLargeRCent50to90(inCutProp.GetNRecoJtPtBinsLargeRCent50to90());}

bool cutPropagator::CheckNGenJtPtSmallBinsLargeRCent0to10(int inNGenJtPtSmallBinsLargeRCent0to10){return CheckNBins(inNGenJtPtSmallBinsLargeRCent0to10, nGenJtPtSmallBinsLargeRCent0to10, "nGenJtPtSmallBinsLargeRCent0to10");}
bool cutPropagator::CheckNGenJtPtSmallBinsLargeRCent0to10(cutPropagator inCutProp){return CheckNGenJtPtSmallBinsLargeRCent0to10(inCutProp.GetNGenJtPtSmallBinsLargeRCent0to10());}
bool cutPropagator::CheckNGenJtPtSmallBinsLargeRCent10to30(int inNGenJtPtSmallBinsLargeRCent10to30){return CheckNBins(inNGenJtPtSmallBinsLargeRCent10to30, nGenJtPtSmallBinsLargeRCent10to30, "nGenJtPtSmallBinsLargeRCent10to30");}
bool cutPropagator::CheckNGenJtPtSmallBinsLargeRCent10to30(cutPropagator inCutProp){return CheckNGenJtPtSmallBinsLargeRCent10to30(inCutProp.GetNGenJtPtSmallBinsLargeRCent10to30());}
bool cutPropagator::CheckNGenJtPtSmallBinsLargeRCent30to50(int inNGenJtPtSmallBinsLargeRCent30to50){return CheckNBins(inNGenJtPtSmallBinsLargeRCent30to50, nGenJtPtSmallBinsLargeRCent30to50, "nGenJtPtSmallBinsLargeRCent30to50");}
bool cutPropagator::CheckNGenJtPtSmallBinsLargeRCent30to50(cutPropagator inCutProp){return CheckNGenJtPtSmallBinsLargeRCent30to50(inCutProp.GetNGenJtPtSmallBinsLargeRCent30to50());}
bool cutPropagator::CheckNGenJtPtSmallBinsLargeRCent50to90(int inNGenJtPtSmallBinsLargeRCent50to90){return CheckNBins(inNGenJtPtSmallBinsLargeRCent50to90, nGenJtPtSmallBinsLargeRCent50to90, "nGenJtPtSmallBinsLargeRCent50to90");}
bool cutPropagator::CheckNGenJtPtSmallBinsLargeRCent50to90(cutPropagator inCutProp){return CheckNGenJtPtSmallBinsLargeRCent50to90(inCutProp.GetNGenJtPtSmallBinsLargeRCent50to90());}

bool cutPropagator::CheckNGenJtPtLargeBinsLargeRCent0to10(int inNGenJtPtLargeBinsLargeRCent0to10){return CheckNBins(inNGenJtPtLargeBinsLargeRCent0to10, nGenJtPtLargeBinsLargeRCent0to10, "nGenJtPtLargeBinsLargeRCent0to10");}
bool cutPropagator::CheckNGenJtPtLargeBinsLargeRCent0to10(cutPropagator inCutProp){return CheckNGenJtPtLargeBinsLargeRCent0to10(inCutProp.GetNGenJtPtLargeBinsLargeRCent0to10());}
bool cutPropagator::CheckNGenJtPtLargeBinsLargeRCent10to30(int inNGenJtPtLargeBinsLargeRCent10to30){return CheckNBins(inNGenJtPtLargeBinsLargeRCent10to30, nGenJtPtLargeBinsLargeRCent10to30, "nGenJtPtLargeBinsLargeRCent10to30");}
bool cutPropagator::CheckNGenJtPtLargeBinsLargeRCent10to30(cutPropagator inCutProp){return CheckNGenJtPtLargeBinsLargeRCent10to30(inCutProp.GetNGenJtPtLargeBinsLargeRCent10to30());}
bool cutPropagator::CheckNGenJtPtLargeBinsLargeRCent30to50(int inNGenJtPtLargeBinsLargeRCent30to50){return CheckNBins(inNGenJtPtLargeBinsLargeRCent30to50, nGenJtPtLargeBinsLargeRCent30to50, "nGenJtPtLargeBinsLargeRCent30to50");}
bool cutPropagator::CheckNGenJtPtLargeBinsLargeRCent30to50(cutPropagator inCutProp){return CheckNGenJtPtLargeBinsLargeRCent30to50(inCutProp.GetNGenJtPtLargeBinsLargeRCent30to50());}
bool cutPropagator::CheckNGenJtPtLargeBinsLargeRCent50to90(int inNGenJtPtLargeBinsLargeRCent50to90){return CheckNBins(inNGenJtPtLargeBinsLargeRCent50to90, nGenJtPtLargeBinsLargeRCent50to90, "nGenJtPtLargeBinsLargeRCent50to90");}
bool cutPropagator::CheckNGenJtPtLargeBinsLargeRCent50to90(cutPropagator inCutProp){return CheckNGenJtPtLargeBinsLargeRCent50to90(inCutProp.GetNGenJtPtLargeBinsLargeRCent50to90());}




bool cutPropagator::CheckNJtAbsEtaBins(int inNJtAbsEtaBins)
{
  bool checkVal = CheckInt(nJtAbsEtaBins, inNJtAbsEtaBins);
  if(!checkVal) std::cout << "cutPropagator check failed on nJtAbsEtaBins" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNJtAbsEtaBins(cutPropagator inCutProp){return CheckNJtAbsEtaBins(inCutProp.GetNJtAbsEtaBins());}


bool cutPropagator::CheckNID(int inNID)
{
  bool checkVal = CheckInt(nID, inNID);
  if(!checkVal) std::cout << "cutPropagator check failed on nID" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNID(cutPropagator inCutProp){return CheckNID(inCutProp.GetNID());}

bool cutPropagator::CheckRecoJtPtBinsSmallR(std::vector<double> inRecoJtPtBinsSmallR)
{
  bool checkVal = CheckVectDouble(recoJtPtBinsSmallR, inRecoJtPtBinsSmallR);
  if(!checkVal) std::cout << "cutPropagator check failed on recoJtPtBinsSmallR" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckRecoJtPtBinsSmallR(cutPropagator inCutProp){return CheckRecoJtPtBinsSmallR(inCutProp.GetRecoJtPtBinsSmallR());}

bool cutPropagator::CheckGenJtPtSmallBinsSmallR(std::vector<double> inGenJtPtSmallBinsSmallR)
{
  bool checkVal = CheckVectDouble(genJtPtSmallBinsSmallR, inGenJtPtSmallBinsSmallR);
  if(!checkVal) std::cout << "cutPropagator check failed on genJtPtSmallBinsSmallR" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckGenJtPtSmallBinsSmallR(cutPropagator inCutProp){return CheckGenJtPtSmallBinsSmallR(inCutProp.GetGenJtPtSmallBinsSmallR());}

bool cutPropagator::CheckGenJtPtLargeBinsSmallR(std::vector<double> inGenJtPtLargeBinsSmallR)
{
  bool checkVal = CheckVectDouble(genJtPtLargeBinsSmallR, inGenJtPtLargeBinsSmallR);
  if(!checkVal) std::cout << "cutPropagator check failed on genJtPtLargeBinsSmallR" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckGenJtPtLargeBinsSmallR(cutPropagator inCutProp){return CheckGenJtPtLargeBinsSmallR(inCutProp.GetGenJtPtLargeBinsSmallR());}


bool cutPropagator::CheckRecoJtPtBinsLargeR(std::vector<double> inRecoJtPtBinsLargeR)
{
  bool checkVal = CheckVectDouble(recoJtPtBinsLargeR, inRecoJtPtBinsLargeR);
  if(!checkVal) std::cout << "cutPropagator check failed on recoJtPtBinsLargeR" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckRecoJtPtBinsLargeR(cutPropagator inCutProp){return CheckRecoJtPtBinsLargeR(inCutProp.GetRecoJtPtBinsLargeR());}

bool cutPropagator::CheckGenJtPtSmallBinsLargeR(std::vector<double> inGenJtPtSmallBinsLargeR)
{
  bool checkVal = CheckVectDouble(genJtPtSmallBinsLargeR, inGenJtPtSmallBinsLargeR);
  if(!checkVal) std::cout << "cutPropagator check failed on genJtPtSmallBinsLargeR" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckGenJtPtSmallBinsLargeR(cutPropagator inCutProp){return CheckGenJtPtSmallBinsLargeR(inCutProp.GetGenJtPtSmallBinsLargeR());}

bool cutPropagator::CheckGenJtPtLargeBinsLargeR(std::vector<double> inGenJtPtLargeBinsLargeR)
{
  bool checkVal = CheckVectDouble(genJtPtLargeBinsLargeR, inGenJtPtLargeBinsLargeR);
  if(!checkVal) std::cout << "cutPropagator check failed on genJtPtLargeBinsLargeR" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckGenJtPtLargeBinsLargeR(cutPropagator inCutProp){return CheckGenJtPtLargeBinsLargeR(inCutProp.GetGenJtPtLargeBinsLargeR());}

bool cutPropagator::CheckJtAbsEtaBinsLow(std::vector<double> inJtAbsEtaBinsLow)
{
  bool checkVal = CheckVectDouble(jtAbsEtaBinsLow, inJtAbsEtaBinsLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtAbsEtaBinsLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtAbsEtaBinsLow(cutPropagator inCutProp){return CheckJtAbsEtaBinsLow(inCutProp.GetJtAbsEtaBinsLow());}

bool cutPropagator::CheckJtAbsEtaBinsHi(std::vector<double> inJtAbsEtaBinsHi)
{
  bool checkVal = CheckVectDouble(jtAbsEtaBinsHi, inJtAbsEtaBinsHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtAbsEtaBinsHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtAbsEtaBinsHi(cutPropagator inCutProp){return CheckJtAbsEtaBinsHi(inCutProp.GetJtAbsEtaBinsHi());}

bool cutPropagator::CheckIdStr(std::vector<std::string> inIDStr)
{
  bool checkVal = CheckVectString(idStr, inIDStr);
  if(!checkVal) std::cout << "cutPropagator check failed on idStr" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckIdStr(cutPropagator inCutProp){return CheckIdStr(inCutProp.GetIdStr());}

bool cutPropagator::CheckJtPfCHMFCutLow(std::vector<double> inJtPfCHMFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfCHMFCutLow, inJtPfCHMFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCHMFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCHMFCutLow(cutPropagator inCutProp){return CheckJtPfCHMFCutLow(inCutProp.GetJtPfCHMFCutLow());}

bool cutPropagator::CheckJtPfCHMFCutHi(std::vector<double> inJtPfCHMFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfCHMFCutHi, inJtPfCHMFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCHMFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCHMFCutHi(cutPropagator inCutProp){return CheckJtPfCHMFCutHi(inCutProp.GetJtPfCHMFCutHi());}

bool cutPropagator::CheckJtPfMUMFCutLow(std::vector<double> inJtPfMUMFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfMUMFCutLow, inJtPfMUMFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMUMFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfMUMFCutLow(cutPropagator inCutProp){return CheckJtPfMUMFCutLow(inCutProp.GetJtPfMUMFCutLow());}

bool cutPropagator::CheckJtPfMUMFCutHi(std::vector<double> inJtPfMUMFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfMUMFCutHi, inJtPfMUMFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMUMFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfMUMFCutHi(cutPropagator inCutProp){return CheckJtPfMUMFCutHi(inCutProp.GetJtPfMUMFCutHi());}

bool cutPropagator::CheckJtPfNHFCutLow(std::vector<double> inJtPfNHFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfNHFCutLow, inJtPfNHFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfNHFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfNHFCutLow(cutPropagator inCutProp){return CheckJtPfNHFCutLow(inCutProp.GetJtPfNHFCutLow());}

bool cutPropagator::CheckJtPfNHFCutHi(std::vector<double> inJtPfNHFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfNHFCutHi, inJtPfNHFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfNHFCutHi" << std::endl;
 return checkVal;
}
bool cutPropagator::CheckJtPfNHFCutHi(cutPropagator inCutProp){return CheckJtPfNHFCutHi(inCutProp.GetJtPfNHFCutHi());}

bool cutPropagator::CheckJtPfNEFCutLow(std::vector<double> inJtPfNEFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfNEFCutLow, inJtPfNEFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfNEFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfNEFCutLow(cutPropagator inCutProp){return CheckJtPfNEFCutLow(inCutProp.GetJtPfNEFCutLow());}

bool cutPropagator::CheckJtPfNEFCutHi(std::vector<double> inJtPfNEFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfNEFCutHi, inJtPfNEFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfNEFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfNEFCutHi(cutPropagator inCutProp){return CheckJtPfNEFCutHi(inCutProp.GetJtPfNEFCutHi());}

bool cutPropagator::CheckJtPfMUFCutLow(std::vector<double> inJtPfMUFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfMUFCutLow, inJtPfMUFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMUFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfMUFCutLow(cutPropagator inCutProp){return CheckJtPfMUFCutLow(inCutProp.GetJtPfMUFCutLow());}

bool cutPropagator::CheckJtPfMUFCutHi(std::vector<double> inJtPfMUFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfMUFCutHi, inJtPfMUFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMUFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfMUFCutHi(cutPropagator inCutProp){return CheckJtPfMUFCutHi(inCutProp.GetJtPfMUFCutHi());}

bool cutPropagator::CheckJtPfCHFCutLow(std::vector<double> inJtPfCHFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfCHFCutLow, inJtPfCHFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCHFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCHFCutLow(cutPropagator inCutProp){return CheckJtPfCHFCutLow(inCutProp.GetJtPfCHFCutLow());}

bool cutPropagator::CheckJtPfCHFCutHi(std::vector<double> inJtPfCHFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfCHFCutHi, inJtPfCHFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCHFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCHFCutHi(cutPropagator inCutProp){return CheckJtPfCHFCutHi(inCutProp.GetJtPfCHFCutHi());}

bool cutPropagator::CheckJtPfCEFCutLow(std::vector<double> inJtPfCEFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfCEFCutLow, inJtPfCEFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCEFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCEFCutLow(cutPropagator inCutProp){return CheckJtPfCEFCutLow(inCutProp.GetJtPfCEFCutLow());}

bool cutPropagator::CheckJtPfCEFCutHi(std::vector<double> inJtPfCEFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfCEFCutHi, inJtPfCEFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCEFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCEFCutHi(cutPropagator inCutProp){return CheckJtPfCEFCutHi(inCutProp.GetJtPfCEFCutHi());}

bool cutPropagator::CheckJtPfMinMult(std::vector<int> inJtPfMinMult)
{
  bool checkVal = CheckVectInt(jtPfMinMult, inJtPfMinMult);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMinMult" << std::endl;
 return checkVal;
}
bool cutPropagator::CheckJtPfMinMult(cutPropagator inCutProp){return CheckJtPfMinMult(inCutProp.GetJtPfMinMult());}

bool cutPropagator::CheckJtPfMinChgMult(std::vector<int> inJtPfMinChgMult)
{
  bool checkVal = CheckVectInt(jtPfMinChgMult, inJtPfMinChgMult);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMinChgMult" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfMinChgMult(cutPropagator inCutProp){return CheckJtPfMinChgMult(inCutProp.GetJtPfMinChgMult());}

bool cutPropagator::CheckNSyst(int inNSyst)
{
  bool checkVal = CheckInt(inNSyst, nSyst);
  if(!checkVal) std::cout << "cutPropagator check failed on nSyst" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNSyst(cutPropagator inCutProp){return CheckNSyst(inCutProp.GetNSyst());}

bool cutPropagator::CheckSystStr(std::vector<std::string> inSystStr)
{
  bool checkVal = CheckVectString(inSystStr, systStr);
  if(!checkVal) std::cout << "cutPropagator check failed on systStr" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckSystStr(cutPropagator inCutProp){return CheckSystStr(inCutProp.GetSystStr());}

bool cutPropagator::CheckNBayes(int inNBayes)
{
  bool checkVal = CheckInt(inNBayes, nBayes);
  if(!checkVal) std::cout << "cutPropagator check failed on nBayes" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNBayes(cutPropagator inCutProp){return CheckNBayes(inCutProp.GetNBayes());}

bool cutPropagator::CheckNBigBayesSymm(int inNBigBayesSymm)
{
  bool checkVal = CheckInt(inNBigBayesSymm, nBigBayesSymm);
  if(!checkVal) std::cout << "cutPropagator check failed on nBigBayesSymm" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNBigBayesSymm(cutPropagator inCutProp){return CheckNBigBayesSymm(inCutProp.GetNBigBayesSymm());}

bool cutPropagator::CheckBayesVal(std::vector<int> inBayesVal)
{
  bool checkVal = CheckVectInt(inBayesVal, bayesVal);
  if(!checkVal) std::cout << "cutPropagator check failed on bayesVal" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckBayesVal(cutPropagator inCutProp){return CheckBayesVal(inCutProp.GetBayesVal());}

bool cutPropagator::CheckNSuperBayes(int inNSuperBayes)
{
  bool checkVal = CheckInt(inNSuperBayes, nSuperBayes);
  if(!checkVal) std::cout << "cutPropagator check failed on nSuperBayes" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNSuperBayes(cutPropagator inCutProp){return CheckNSuperBayes(inCutProp.GetNSuperBayes());}

bool cutPropagator::CheckRCDiffFileName(std::string inRCDiffFileName)
{
  bool checkVal = CheckString(inRCDiffFileName, rcDiffFileName);
  if(!checkVal) std::cout << "cutPropagator check failed on rcDiffFileName" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckRCDiffFileName(cutPropagator inCutProp){return CheckRCDiffFileName(inCutProp.GetRCDiffFileName());}

bool cutPropagator::CheckPriorFlatFileName(std::string inPriorFlatFileName)
{
  bool checkVal = CheckString(inPriorFlatFileName, flatPriorFileName);
  if(!checkVal) std::cout << "cutPropagator check failed on flatPriorFileName" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckPriorFlatFileName(cutPropagator inCutProp){return CheckPriorFlatFileName(inCutProp.GetPriorFlatFileName());}

bool cutPropagator::CheckJECVarMC(double inJECVarMC)
{
  bool checkVal = CheckDouble(inJECVarMC, jecVarMC);
  if(!checkVal) std::cout << "cutPropagator check failed on jecVarMC" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJECVarMC(cutPropagator inCutProp){return CheckJECVarMC(inCutProp.GetJECVarMC());}

bool cutPropagator::CheckJERVarMC(double inJERVarMC)
{
  bool checkVal = CheckDouble(inJERVarMC, jerVarMC);
  if(!checkVal) std::cout << "cutPropagator check failed on jerVarMC" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJERVarMC(cutPropagator inCutProp){return CheckJERVarMC(inCutProp.GetJERVarMC());}

bool cutPropagator::CheckJECVarData(double inJECVarData)
{
  bool checkVal = CheckDouble(inJECVarData, jecVarData);
  if(!checkVal) std::cout << "cutPropagator check failed on jecVarData" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJECVarData(cutPropagator inCutProp){return CheckJECVarData(inCutProp.GetJECVarData());}

bool cutPropagator::CheckNResponseMod(int inNResponseMod)
{
  bool checkVal = CheckInt(inNResponseMod, nResponseMod);
  if(!checkVal) std::cout << "cutPropagator check failed on nResponseMod" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNResponseMod(cutPropagator inCutProp){return CheckNResponseMod(inCutProp.GetNResponseMod());}

bool cutPropagator::CheckResponseMod(std::vector<double> inResponseMod)
{
  bool checkVal = CheckVectDouble(inResponseMod, responseMod);
  if(!checkVal) std::cout << "cutPropagator check failed on responseMod" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckResponseMod(cutPropagator inCutProp){return CheckResponseMod(inCutProp.GetResponseMod());}

bool cutPropagator::CheckJERVarData(std::vector<double> inJERVarData)
{
  bool checkVal = CheckVectDouble(inJERVarData, jerVarData);
  if(!checkVal) std::cout << "cutPropagator check failed on jerVarData" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJERVarData(cutPropagator inCutProp){return CheckJERVarData(inCutProp.GetJERVarData());}

bool cutPropagator::CheckNPthats(int inNPthats)
{
  bool checkVal = CheckInt(inNPthats, nPthats);
  if(!checkVal) std::cout << "cutPropagator check failed on nPthats" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNPthats(cutPropagator inCutProp){return CheckNPthats(inCutProp.GetNPthats());}

bool cutPropagator::CheckPthats(std::vector<double> inPthats)
{
  bool checkVal = CheckVectDouble(inPthats, pthats);
  if(!checkVal) std::cout << "cutPropagator check failed on pthats" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckPthats(cutPropagator inCutProp){return CheckPthats(inCutProp.GetPthats());}

bool cutPropagator::CheckPthatWeights(std::vector<double> inPthatWeights)
{
  bool checkVal = CheckVectDouble(inPthatWeights, pthatWeights);
  if(!checkVal) std::cout << "cutPropagator check failed on pthatWeights" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckPthatWeights(cutPropagator inCutProp){return CheckPthatWeights(inCutProp.GetPthatWeights());}

bool cutPropagator::CheckNJtAlgos(int inNJtAlgos)
{
  bool checkVal = CheckInt(inNJtAlgos, nJtAlgos);
  if(!checkVal) std::cout << "cutPropagator check failed on nJtAlgos" << std::endl;
 return checkVal;
}
bool cutPropagator::CheckNJtAlgos(cutPropagator inCutProp){return CheckNJtAlgos(inCutProp.GetNJtAlgos());}

bool cutPropagator::CheckJtAlgos(std::vector<std::string> inJtAlgos)
{
  bool checkVal = CheckVectString(inJtAlgos, jtAlgos);
  if(!checkVal) std::cout << "cutPropagator check failed on jtAlgos" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtAlgos(cutPropagator inCutProp){return CheckJtAlgos(inCutProp.GetJtAlgos());}

bool cutPropagator::CheckMinJtPtCut(std::vector<double> inMinJtPtCut)
{
  bool checkVal = CheckVectDouble(inMinJtPtCut, minJtPtCut);
  if(!checkVal) std::cout << "cutPropagator check failed on minJtPtCut" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckMinJtPtCut(cutPropagator inCutProp){return CheckMinJtPtCut(inCutProp.GetMinJtPtCut());}

bool cutPropagator::CheckMultiJtPtCut(std::vector<double> inMultiJtPtCut)
{
  bool checkVal = CheckVectDouble(inMultiJtPtCut, multiJtPtCut);
  if(!checkVal) std::cout << "cutPropagator check failed on multiJtPtCut" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckMultiJtPtCut(cutPropagator inCutProp){return CheckMultiJtPtCut(inCutProp.GetMultiJtPtCut());}

bool cutPropagator::CheckRecoTruncPos(std::vector<int> inRecoTruncPos)
{
  bool checkVal = CheckVectInt(inRecoTruncPos, recoTruncPos);
  if(!checkVal) std::cout << "cutPropagator check failed on recoTruncPos" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckRecoTruncPos(cutPropagator inCutProp){return CheckRecoTruncPos(inCutProp.GetRecoTruncPos());}

bool cutPropagator::CheckIsPP(bool inIsPP)
{
  bool checkVal = CheckBool(inIsPP, isPP);
  if(!checkVal) std::cout << "cutPropagator check failed on isPP" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckIsPP(cutPropagator inCutProp){return CheckIsPP(inCutProp.GetIsPP());}

bool cutPropagator::CheckNCentBins(int inNCentBins)
{
  bool checkVal = CheckInt(inNCentBins, nCentBins);
  if(!checkVal) std::cout << "cutPropagator check failed on nCentBins" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNCentBins(cutPropagator inCutProp){return CheckNCentBins(inCutProp.GetNCentBins());}

bool cutPropagator::CheckCentBinsLow(std::vector<int> inCentBinsLow)
{
  bool checkVal = CheckVectInt(inCentBinsLow, centBinsLow);
  if(!checkVal) std::cout << "cutPropagator check failed on centBinsLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckCentBinsLow(cutPropagator inCutProp){return CheckCentBinsLow(inCutProp.GetCentBinsLow());}

bool cutPropagator::CheckCentBinsHi(std::vector<int> inCentBinsHi)
{
  bool checkVal = CheckVectInt(inCentBinsHi, centBinsHi);
  if(!checkVal) std::cout << "cutPropagator check failed on centBinsHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckCentBinsHi(cutPropagator inCutProp){return CheckCentBinsHi(inCutProp.GetCentBinsHi());}
//end checks


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

void cutPropagator::SetSmallRVals(int inN, const Int_t inSmallRVals[])
{
  for(int i = 0; i < inN; ++i){
    smallRVals.push_back(inSmallRVals[i]);
  }

  return;
}


void cutPropagator::SetLargeRVals(int inN, const Int_t inLargeRVals[])
{
  for(int i = 0; i < inN; ++i){
    largeRVals.push_back(inLargeRVals[i]);
  }

  return;
}


void cutPropagator::SetRecoJtPtBinsSmallR(int inN, const Double_t inRecoJtPtBinsSmallR[])
{
  for(int i = 0; i < inN; ++i){
    recoJtPtBinsSmallR.push_back(inRecoJtPtBinsSmallR[i]);
  }

  return;
}

void cutPropagator::SetGenJtPtSmallBinsSmallR(int inN, const Double_t inGenJtPtSmallBinsSmallR[])
{
  for(int i = 0; i < inN; ++i){
    genJtPtSmallBinsSmallR.push_back(inGenJtPtSmallBinsSmallR[i]);
  }

  return;
}

void cutPropagator::SetGenJtPtLargeBinsSmallR(int inN, const Double_t inGenJtPtLargeBinsSmallR[])
{
  for(int i = 0; i < inN; ++i){
    genJtPtLargeBinsSmallR.push_back(inGenJtPtLargeBinsSmallR[i]);
  }

  return;
}

void cutPropagator::SetRecoJtPtBinsLargeR(int inN, const Double_t inRecoJtPtBinsLargeR[])
{
  for(int i = 0; i < inN; ++i){
    recoJtPtBinsLargeR.push_back(inRecoJtPtBinsLargeR[i]);
  }

  return;
}

void cutPropagator::SetGenJtPtSmallBinsLargeR(int inN, const Double_t inGenJtPtSmallBinsLargeR[])
{
  for(int i = 0; i < inN; ++i){
    genJtPtSmallBinsLargeR.push_back(inGenJtPtSmallBinsLargeR[i]);
  }

  return;
}

void cutPropagator::SetGenJtPtLargeBinsLargeR(int inN, const Double_t inGenJtPtLargeBinsLargeR[])
{
  for(int i = 0; i < inN; ++i){
    genJtPtLargeBinsLargeR.push_back(inGenJtPtLargeBinsLargeR[i]);
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
