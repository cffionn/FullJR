#include <iostream>
#include <string>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TMath.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/smearingFuncs.h"

int processAndRecoPYT8(const std::string inFileName)
{
  TRandom3* randGen_p= new TRandom3(0);

  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "Given inFileName \'" << inFileName << "\' does not contain \'.root\'. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  while(outFileName.find(".root") != std::string::npos){outFileName.replace(outFileName.find(".root"), 5, "");}
  outFileName = "output/" + outFileName + "_RecoProcess_" + dateStr + ".root";

  const Int_t nRVals = 2;
  const Double_t rVals[nRVals] = {0.4, 0.8};
  const std::string rAlgos[nRVals] = {"akCs4PU3PFFlow", "akCs8PU3PFFlow"};
  fastjet::JetDefinition jetDef[nRVals];


  const fastjet::JetAlgorithm jetAlgo = fastjet::antikt_algorithm;

  const Float_t globalMinJtPt = 30.;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* outTree_p[nRVals];
  Float_t pthatOut_[nRVals];
  Float_t pthatWeightOut_[nRVals];

  const Int_t nPartMax = 2000;
  Int_t nJt_[nRVals];
  Float_t genJtPt_[nRVals][nPartMax];
  Float_t genJtPtLogLoss_[nRVals][nPartMax];
  Float_t genJtPhi_[nRVals][nPartMax];
  Float_t genJtEta_[nRVals][nPartMax];
  Float_t recoJtPt_[nRVals][nPartMax];
  Float_t recoJtPt1p1Res_[nRVals][nPartMax];
  Float_t recoJtPt0p9Res_[nRVals][nPartMax];
  
  for(Int_t rI = 0; rI < nRVals; ++rI){
    jetDef[rI] = fastjet::JetDefinition(jetAlgo, rVals[rI]);

    outTree_p[rI] = new TTree(("recoAndGenTreeR" + prettyString(rVals[rI], 1, true)).c_str(), "");
    
    outTree_p[rI]->Branch("pthat", &(pthatOut_[rI]), "pthat/F");
    outTree_p[rI]->Branch("pthatWeight", &(pthatWeightOut_[rI]), "pthatWeight/F");
    outTree_p[rI]->Branch("nJt", &(nJt_[rI]), "nJt/I");
    outTree_p[rI]->Branch("genJtPt", genJtPt_[rI], "genJtPt[nJt]/F");
    outTree_p[rI]->Branch("genJtPtLogLoss", genJtPtLogLoss_[rI], "genJtPtLogLoss[nJt]/F");
    outTree_p[rI]->Branch("genJtPhi", genJtPhi_[rI], "genJtPhi[nJt]/F");
    outTree_p[rI]->Branch("genJtEta", genJtEta_[rI], "genJtEta[nJt]/F");
    outTree_p[rI]->Branch("recoJtPt", recoJtPt_[rI], "recoJtPt[nJt]/F");
    outTree_p[rI]->Branch("recoJtPt1p1Res", recoJtPt1p1Res_[rI], "recoJtPt1p1Res[nJt]/F");
    outTree_p[rI]->Branch("recoJtPt0p9Res", recoJtPt0p9Res_[rI], "recoJtPt0p9Res[nJt]/F");
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("genTree");
  
  Float_t pthat_;
  Float_t pthatWeight_;

  Int_t nPart_;
  Float_t pt_[nPartMax];
  Float_t phi_[nPartMax];
  Float_t eta_[nPartMax];

  inTree_p->SetBranchAddress("pthat", &pthat_);
  inTree_p->SetBranchAddress("pthatWeight", &pthatWeight_);

  inTree_p->SetBranchAddress("nPart", &nPart_);
  inTree_p->SetBranchAddress("pt", &pt_);
  inTree_p->SetBranchAddress("phi", &phi_);
  inTree_p->SetBranchAddress("eta", &eta_);

  const Int_t nEntries = TMath::Min((Int_t)100000000, (Int_t)inTree_p->GetEntries());

  std::cout << "Processing Entries: " << nEntries << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;

    inTree_p->GetEntry(entry);

    std::vector<fastjet::PseudoJet> particles;

    for(Int_t i = 0; i < nPart_; ++i){
      TLorentzVector temp;
      temp.SetPtEtaPhiM(pt_[i], eta_[i], phi_[i], 0.);

      particles.push_back(fastjet::PseudoJet(temp.Px(), temp.Py(), temp.Pz(), temp.E()));
    }

    for(Int_t rI = 0; rI < nRVals; ++rI){
      fastjet::ClusterSequence* cs = new fastjet::ClusterSequence(particles, jetDef[rI]);

      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs->inclusive_jets(globalMinJtPt));
      
      pthatOut_[rI] = pthat_;
      pthatWeightOut_[rI] = pthatWeight_;
      nJt_[rI] = 0;
      for(unsigned int jI = 0; jI < jets.size(); ++jI){
	if(TMath::Abs(jets.at(jI).eta()) > 2.) continue;

	genJtPt_[rI][nJt_[rI]] = jets.at(jI).pt();
	genJtPtLogLoss_[rI][nJt_[rI]] = jets.at(jI).pt() - TMath::Log(jets.at(jI).pt());
	if(genJtPtLogLoss_[rI][nJt_[rI]] < 0.) genJtPtLogLoss_[rI][nJt_[rI]] = 0.;
	genJtPhi_[rI][nJt_[rI]] = jets.at(jI).phi();
	genJtEta_[rI][nJt_[rI]] = jets.at(jI).eta();

	recoJtPt_[rI][nJt_[rI]] = jets.at(jI).pt()*(Double_t)randGen_p->Gaus(1.0, (Double_t)getResForPtEtaCentAlgo(jets.at(jI).pt(), jets.at(jI).eta(), 0, rAlgos[rI]));
	recoJtPt1p1Res_[rI][nJt_[rI]] = jets.at(jI).pt()*(Double_t)randGen_p->Gaus(1.0, 1.1*(Double_t)getResForPtEtaCentAlgo(jets.at(jI).pt(), jets.at(jI).eta(), 0, rAlgos[rI]));
	recoJtPt0p9Res_[rI][nJt_[rI]] = jets.at(jI).pt()*(Double_t)randGen_p->Gaus(1.0, 0.9*(Double_t)getResForPtEtaCentAlgo(jets.at(jI).pt(), jets.at(jI).eta(), 0, rAlgos[rI]));
	
	++nJt_[rI];
      }

      delete cs;
    }
    
    for(Int_t rI = 0; rI < nRVals; ++rI){outTree_p[rI]->Fill();}
  }


  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  for(Int_t rI = 0; rI < nRVals; ++rI){
    outTree_p[rI]->Write("", TObject::kOverwrite);
    delete outTree_p[rI];
  }

  outFile_p->Close();
  delete outFile_p;  

  delete randGen_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/processAndRecoPYT8.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += processAndRecoPYT8(argv[1]);
  return retVal;
}
