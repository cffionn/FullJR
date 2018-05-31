#include <iostream>
#include <string>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TNamed.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/etaPhiFunc.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/smearingFuncs.h"

int processAndRecoPYT8(const std::string inFileName, const bool doRecoJets)
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

  std::string doRecoJetsStr = "DoRecoJetsTrue";
  if(!doRecoJets) doRecoJetsStr = "DoRecoJetsFalse";

  checkMakeDir("output");
  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  while(outFileName.find(".root") != std::string::npos){outFileName.replace(outFileName.find(".root"), 5, "");}
  outFileName = "output/" + outFileName + "_RecoProcess_" + doRecoJetsStr  + "_" + dateStr + ".root";

  const Int_t nRVals = 2;
  const Double_t rVals[nRVals] = {0.4, 0.8};
  const std::string rAlgos[nRVals] = {"akCs4PU3PFFlow", "akCs8PU3PFFlow"};
  fastjet::JetDefinition jetDef[nRVals];
  const Float_t cVal = .06;
  const Float_t sVal = 1.0;
  const Float_t rcWidth[nRVals] = {20., 40.};
  const fastjet::RecombinationScheme recombScheme[nRVals] = {fastjet::E_scheme, fastjet::E_scheme};
  
  const fastjet::JetAlgorithm jetAlgo = fastjet::antikt_algorithm;

  const Float_t globalMinJtPt = 30.;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* partonTree_p = new TTree("partonTree", "");

  const Int_t nEnergyFracBins = 16;
  const Float_t energyFracLow = 0.2;
  const Float_t energyFracHi = 1.0;
  Double_t energyFracBins[nEnergyFracBins+1];
  getLinBins(energyFracLow, energyFracHi, nEnergyFracBins, energyFracBins);

  const Int_t nQGMax = 2;
  Int_t nQG_;
  Float_t qgPt_[nQGMax];
  Float_t qgPhi_[nQGMax];
  Float_t qgEta_[nQGMax];
  Int_t qgID_[nQGMax];
  Float_t energyFrac_[nQGMax][nEnergyFracBins];

  Float_t pthatOut_;
  Float_t pthatWeightOut_;

  partonTree_p->Branch("nQG", &nQG_, "nQG/I");
  partonTree_p->Branch("qgPt", qgPt_, "qgPt[nQG]/F");
  partonTree_p->Branch("qgPhi", qgPhi_, "qgPhi[nQG]/F");
  partonTree_p->Branch("qgEta", qgEta_, "qgEta[nQG]/F");
  partonTree_p->Branch("qgID", qgID_, "qgID[nQG]/I");
  partonTree_p->Branch("energyFrac", energyFrac_, ("energyFrac[nQG][" + std::to_string(nEnergyFracBins) + "]/F").c_str());
  partonTree_p->Branch("pthat", &pthatOut_, "pthat/F");
  partonTree_p->Branch("pthatWeight", &pthatWeightOut_, "pthatWeight/F");

  TTree* outTree_p[nRVals];

  const Int_t nPartMax = 2000;
  Int_t nJt_[nRVals];
  Float_t genJtPt_[nRVals][nPartMax];
  Float_t genJtPhi_[nRVals][nPartMax];
  Float_t genJtEta_[nRVals][nPartMax];
  Float_t jtEnergyFrac_[nRVals][nPartMax][nEnergyFracBins];
  Float_t recoJtPt_[nRVals][nPartMax];
  Float_t recoJtPt1p1Res_[nRVals][nPartMax];
  Float_t recoJtPt1p15Res_[nRVals][nPartMax];
  Float_t recoJtPt0p85Res_[nRVals][nPartMax];
  Float_t recoJtPt0p9Res_[nRVals][nPartMax];
  
  for(Int_t rI = 0; rI < nRVals; ++rI){
    jetDef[rI] = fastjet::JetDefinition(jetAlgo, rVals[rI], recombScheme[rI]);

    std::string schemeStr = "EScheme";
    if(recombScheme[rI] == fastjet::WTA_pt_scheme) schemeStr = "WTAPtScheme";

    outTree_p[rI] = new TTree(("recoAndGenTreeR" + prettyString(rVals[rI], 1, true) + schemeStr).c_str(), "");
    
    outTree_p[rI]->Branch("nJt", &(nJt_[rI]), "nJt/I");
    outTree_p[rI]->Branch("genJtPt", genJtPt_[rI], "genJtPt[nJt]/F");
    outTree_p[rI]->Branch("genJtPhi", genJtPhi_[rI], "genJtPhi[nJt]/F");
    outTree_p[rI]->Branch("genJtEta", genJtEta_[rI], "genJtEta[nJt]/F");
    outTree_p[rI]->Branch("jtEnergyFrac", jtEnergyFrac_[rI], ("jtEnergyFrac[nJt][" + std::to_string(nEnergyFracBins) + "]/F").c_str());

    if(doRecoJets){
      outTree_p[rI]->Branch("recoJtPt", recoJtPt_[rI], "recoJtPt[nJt]/F");
      outTree_p[rI]->Branch("recoJtPt1p1Res", recoJtPt1p1Res_[rI], "recoJtPt1p1Res[nJt]/F");
      outTree_p[rI]->Branch("recoJtPt1p15Res", recoJtPt1p15Res_[rI], "recoJtPt1p15Res[nJt]/F");
      outTree_p[rI]->Branch("recoJtPt0p9Res", recoJtPt0p9Res_[rI], "recoJtPt0p9Res[nJt]/F");
      outTree_p[rI]->Branch("recoJtPt0p85Res", recoJtPt0p85Res_[rI], "recoJtPt0p85Res[nJt]/F");
    }
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("genTree");
  
  Float_t pthat_;
  Float_t pthatWeight_;

  Int_t nPart_;
  Float_t pt_[nPartMax];
  Float_t phi_[nPartMax];
  Float_t eta_[nPartMax];
  Float_t m_[nPartMax];

  inTree_p->SetBranchAddress("pthat", &pthat_);
  inTree_p->SetBranchAddress("pthatWeight", &pthatWeight_);

  inTree_p->SetBranchAddress("nQG", &nQG_);
  inTree_p->SetBranchAddress("qgPt", qgPt_);
  inTree_p->SetBranchAddress("qgPhi", qgPhi_);
  inTree_p->SetBranchAddress("qgEta", qgEta_);
  inTree_p->SetBranchAddress("qgID", qgID_);

  inTree_p->SetBranchAddress("nPart", &nPart_);
  inTree_p->SetBranchAddress("pt", pt_);
  inTree_p->SetBranchAddress("phi", phi_);
  inTree_p->SetBranchAddress("eta", eta_);
  inTree_p->SetBranchAddress("m", m_);

  const Int_t nEntries = TMath::Min((Int_t)1000000000, (Int_t)inTree_p->GetEntries());

  std::cout << "Processing Entries: " << nEntries << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;

    inTree_p->GetEntry(entry);

    std::vector<fastjet::PseudoJet> particles;

    for(Int_t qI = 0; qI < nQG_; ++qI){
      for(Int_t bI = 0; bI < nEnergyFracBins; ++bI){
	energyFrac_[qI][bI] = 0.;
      }
    }

    for(Int_t i = 0; i < nPart_; ++i){
      TLorentzVector temp;
      temp.SetPtEtaPhiM(pt_[i], eta_[i], phi_[i], m_[i]);

      for(Int_t qI = 0; qI < nQG_; ++qI){
	for(Int_t bI = 0; bI < nEnergyFracBins; ++bI){
	  
	  if(getDR(qgEta_[qI], qgPhi_[qI], eta_[i], phi_[i]) < energyFracBins[bI]){
	    energyFrac_[qI][bI] += pt_[i]/qgPt_[qI];
	  }

	}
      }
      

      particles.push_back(fastjet::PseudoJet(temp.Px(), temp.Py(), temp.Pz(), temp.E()));
    }
  
    pthatOut_ = pthat_;
    pthatWeightOut_ = pthatWeight_;

    for(Int_t rI = 0; rI < nRVals; ++rI){
      fastjet::ClusterSequence* cs = new fastjet::ClusterSequence(particles, jetDef[rI]);

      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs->inclusive_jets(globalMinJtPt));
      
      nJt_[rI] = 0;
      for(unsigned int jI = 0; jI < jets.size(); ++jI){
	if(TMath::Abs(jets.at(jI).eta()) > 2.) continue;

	genJtPt_[rI][nJt_[rI]] = jets.at(jI).pt();
	genJtPhi_[rI][nJt_[rI]] = jets.at(jI).phi();
	genJtEta_[rI][nJt_[rI]] = jets.at(jI).eta();

	if(doRecoJets){
	  //use toy values
	  Double_t sigma = TMath::Sqrt(cVal*cVal + sVal*sVal/jets.at(jI).pt() + rcWidth[rI]*rcWidth[rI]/(jets.at(jI).pt()*jets.at(jI).pt()));
	  
	  recoJtPt_[rI][nJt_[rI]] = jets.at(jI).pt()*(Double_t)randGen_p->Gaus(1.0, sigma);
	  recoJtPt1p1Res_[rI][nJt_[rI]] = jets.at(jI).pt()*(Double_t)randGen_p->Gaus(1.0, 1.1*sigma);
	  recoJtPt1p15Res_[rI][nJt_[rI]] = jets.at(jI).pt()*(Double_t)randGen_p->Gaus(1.0, 1.15*sigma);
	  recoJtPt0p9Res_[rI][nJt_[rI]] = jets.at(jI).pt()*(Double_t)randGen_p->Gaus(1.0, 0.9*sigma);
	  recoJtPt0p85Res_[rI][nJt_[rI]] = jets.at(jI).pt()*(Double_t)randGen_p->Gaus(1.0, 0.85*sigma);
	  
	  if(recoJtPt_[rI][nJt_[rI]] < 0.0) recoJtPt_[rI][nJt_[rI]] = 0.0;
	  if(recoJtPt1p1Res_[rI][nJt_[rI]] < 0.0) recoJtPt1p1Res_[rI][nJt_[rI]] = 0.0;
	  if(recoJtPt1p15Res_[rI][nJt_[rI]] < 0.0) recoJtPt1p15Res_[rI][nJt_[rI]] = 0.0;
	  if(recoJtPt0p9Res_[rI][nJt_[rI]] < 0.0) recoJtPt0p9Res_[rI][nJt_[rI]] = 0.0;
	  if(recoJtPt0p85Res_[rI][nJt_[rI]] < 0.0) recoJtPt0p85Res_[rI][nJt_[rI]] = 0.0;
	}

	++nJt_[rI];
      }
    

      for(Int_t jI = 0; jI < nJt_[rI]; ++jI){
	for(Int_t pI = 0; pI < nEnergyFracBins; ++pI){
	  jtEnergyFrac_[rI][jI][pI] = 0.;
	}
	
	for(Int_t i = 0; i < nPart_; ++i){
	  for(Int_t bI = 0; bI < nEnergyFracBins; ++bI){
	    
	    if(getDR(genJtEta_[rI][jI], genJtPhi_[rI][jI], eta_[i], phi_[i]) < energyFracBins[bI]){
	      jtEnergyFrac_[rI][jI][bI] += pt_[i]/genJtPt_[rI][jI];
	    }
	    
	  }
	}
      }

      delete cs;
    }
    
    partonTree_p->Fill();
    for(Int_t rI = 0; rI < nRVals; ++rI){outTree_p[rI]->Fill();}
  }


  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  partonTree_p->Write("", TObject::kOverwrite);
  delete partonTree_p;

  for(Int_t rI = 0; rI < nRVals; ++rI){
    outTree_p[rI]->Write("", TObject::kOverwrite);
    delete outTree_p[rI];
  }

  TDirectory* dir_p = (TDirectory*)outFile_p->mkdir("cutDir");
  dir_p->cd();

  TNamed doRecoJetsName("doRecoJets", std::to_string(doRecoJets));
  TNamed nEnergyFracBinsName("nEnergyFracBins", std::to_string(nEnergyFracBins).c_str());
  TNamed energyFracLowName("energyFracLow", std::to_string(energyFracLow).c_str());
  TNamed energyFracHiName("energyFracHi", std::to_string(energyFracHi).c_str());

  std::string energyFracBinsStr = std::to_string(energyFracBins[0]);
  for(Int_t eI = 1; eI < nEnergyFracBins+1; ++eI){
    energyFracBinsStr = energyFracBinsStr + "," + std::to_string(energyFracBins[eI]);
  }

  TNamed energyFracBinsName("energyFracBins", energyFracBinsStr.c_str());

  doRecoJetsName.Write("", TObject::kOverwrite);
  nEnergyFracBinsName.Write("", TObject::kOverwrite);
  energyFracLowName.Write("", TObject::kOverwrite);
  energyFracHiName.Write("", TObject::kOverwrite);
  energyFracBinsName.Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;  

  delete randGen_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/processAndRecoPYT8.exe <inFileName> <doRecoJets>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += processAndRecoPYT8(argv[1], std::stoi(argv[2]));
  return retVal;
}
