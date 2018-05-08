
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TNamed.h"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/pseudoTowerGeometry.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/getLogBins.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/findBinPos.h"
#include "Utility/include/etaPhiFunc.h"
#include "Utility/include/goodGlobalSelection.h"

int makeSubtractedEvent(const std::string inFileName)
{
  TDatime* date = new TDatime();

  checkMakeDir("output");
  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){
    outFileName.replace(0, outFileName.find("/")+1, "");
  }
  outFileName.replace(outFileName.find(".root"), 5, "");
  outFileName = "output/" + outFileName + "_Subtracted_" + std::to_string(date->GetDate()) + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  const Int_t nPhiBins = 72;
  Double_t phiBins[nPhiBins+1];
  getLinBins(-TMath::Pi(), TMath::Pi(), nPhiBins, phiBins);
  
  const Int_t nPhiBinsRed = 36;
  Double_t phiBinsRed[nPhiBinsRed+1];
  getLinBins(-TMath::Pi(), TMath::Pi(), nPhiBinsRed, phiBinsRed);
  
  //  const double_t hfCut = .1;
  const Double_t etaMax = 2.0;

  pseudoTowGeo towGeo;
  std::vector<double> towEtaBounds = towGeo.getEtaTowBounds();
  Int_t etaCount = 0;
  for(unsigned int eI = 0; eI < towEtaBounds.size(); ++eI){
    if(TMath::Abs(towEtaBounds.at(eI)) < etaMax) ++etaCount;
  }
  
  const Int_t nEtaBins = etaCount;
  const Int_t nEtaBinsRed = (nEtaBins-1)/2 + 1;
  Double_t etaBins[nEtaBins];
  Double_t etaBinsRed[(nEtaBins-1)/2 + 1];
  unsigned int etaPos = 0;

  int fillCount = 0;
  int currBin = 0;
  for(unsigned int eI = 0; eI < towEtaBounds.size(); ++eI){
    if(TMath::Abs(towEtaBounds.at(eI)) < etaMax){
      etaBins[etaPos] = towEtaBounds.at(eI);
      if(fillCount%2 == 0){
	etaBinsRed[currBin] = towEtaBounds.at(eI);
	currBin += 1;
      }

      ++fillCount;
      ++etaPos;
    }
  }

  for(Int_t eI = 0; eI < nEtaBins; ++eI){
    std::cout << etaBins[eI] << ",";
  }
  std::cout << std::endl;

  
  for(Int_t eI = 0; eI < nEtaBinsRed; ++eI){
    std::cout << etaBinsRed[eI] << ",";
  }
  std::cout << std::endl;
  
  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);

  const Int_t perEventNum = 80;
  const Int_t jetPerEventNum = 80;
  const Int_t histNum = 150000;
  const Int_t nCentBins = 8;
  const Int_t centBinsLow[nCentBins] = {0,  5, 10, 20, 30, 40, 50, 70};
  const Int_t centBinsHi[nCentBins] = { 5, 10, 20, 30, 40, 50, 70, 90};

  Int_t eventNum = 0;
  Int_t jetEventNum = 0;
  Int_t currNum[nCentBins];

  const Int_t nPFID = 8;

  TH1F* etaEnergy_ByID_h[nCentBins][nPFID];
  TH1F* etaEnergyCS_ByID_h[nCentBins][nPFID];
  TH1F* etaEnergyCSPU_ByID_h[nCentBins][nPFID];
  TH1F* etaEnergyCSPUFlow_ByID_h[nCentBins][nPFID];

  TH2F* etaPhiEnergy_Average_p[nCentBins];
  TH2F* etaPhiEnergyCS_Average_p[nCentBins];
  TH2F* etaPhiEnergyCSPU_Average_p[nCentBins];
  TH2F* etaPhiEnergyCSPUFlow_Average_p[nCentBins];

  TH1F* etaPhiEnergy_Average_PhiProj_p[nCentBins];
  TH1F* etaPhiEnergyCS_Average_PhiProj_p[nCentBins];
  TH1F* etaPhiEnergyCSPU_Average_PhiProj_p[nCentBins];
  TH1F* etaPhiEnergyCSPUFlow_Average_PhiProj_p[nCentBins];
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

    for(Int_t pfI = 0; pfI < nPFID; ++pfI){
      std::string pfStr = "PFIDAll";
      if(pfI != 0) pfStr = "PFID" + std::to_string(pfI);

      etaEnergy_ByID_h[cI][pfI] = new TH1F(("etaEnergy_" + centStr + "_" + pfStr +"_p").c_str(), ";#eta;Ave", nEtaBins-1, etaBins);
      etaEnergyCS_ByID_h[cI][pfI] = new TH1F(("etaEnergyCS_" + centStr + "_" + pfStr +"_p").c_str(), ";#eta;Ave", nEtaBins-1, etaBins);
      etaEnergyCSPU_ByID_h[cI][pfI] = new TH1F(("etaEnergyCSPU_" + centStr + "_" + pfStr +"_p").c_str(), ";#eta;Ave", nEtaBins-1, etaBins);
      etaEnergyCSPUFlow_ByID_h[cI][pfI] = new TH1F(("etaEnergyCSPUFlow_" + centStr + "_" + pfStr +"_p").c_str(), ";#eta;Ave", nEtaBins-1, etaBins);
    }

    etaPhiEnergy_Average_p[cI] = new TH2F(("etaPhiEnergy_Average_" + centStr + "_p").c_str(), ";#eta;#phi - #Psi_{HF,2}", nEtaBinsRed-1, etaBinsRed, nPhiBinsRed, phiBinsRed);
    etaPhiEnergyCS_Average_p[cI] = new TH2F(("etaPhiEnergyCS_Average_" + centStr + "_p").c_str(), ";#eta;#phi - #Psi_{HF,2}", nEtaBinsRed-1, etaBinsRed, nPhiBinsRed, phiBinsRed);
    etaPhiEnergyCSPU_Average_p[cI] = new TH2F(("etaPhiEnergyCSPU_Average_" + centStr + "_p").c_str(), ";#eta;#phi - #Psi_{HF,2}", nEtaBinsRed-1, etaBinsRed, nPhiBinsRed, phiBinsRed);
    etaPhiEnergyCSPUFlow_Average_p[cI] = new TH2F(("etaPhiEnergyCSPUFlow_Average_" + centStr + "_p").c_str(), ";#eta;#phi - #Psi_{HF,2}", nEtaBinsRed-1, etaBinsRed, nPhiBinsRed, phiBinsRed);
    centerTitles({etaPhiEnergy_Average_p[cI], etaPhiEnergyCS_Average_p[cI], etaPhiEnergyCSPU_Average_p[cI], etaPhiEnergyCSPUFlow_Average_p[cI]});
    setSumW2({etaPhiEnergy_Average_p[cI], etaPhiEnergyCS_Average_p[cI], etaPhiEnergyCSPU_Average_p[cI], etaPhiEnergyCSPUFlow_Average_p[cI]});
    currNum[cI] = 0;

     etaPhiEnergy_Average_p[cI]->GetXaxis()->SetTitleFont(43);
     etaPhiEnergyCS_Average_p[cI]->GetXaxis()->SetTitleFont(43);
     etaPhiEnergyCSPU_Average_p[cI]->GetXaxis()->SetTitleFont(43);
     etaPhiEnergyCSPUFlow_Average_p[cI]->GetXaxis()->SetTitleFont(43);

     etaPhiEnergy_Average_p[cI]->GetXaxis()->SetTitleSize(18);
     etaPhiEnergyCS_Average_p[cI]->GetXaxis()->SetTitleSize(18);
     etaPhiEnergyCSPU_Average_p[cI]->GetXaxis()->SetTitleSize(18);
     etaPhiEnergyCSPUFlow_Average_p[cI]->GetXaxis()->SetTitleSize(18);

     etaPhiEnergy_Average_p[cI]->GetYaxis()->SetTitleFont(43);
     etaPhiEnergyCS_Average_p[cI]->GetYaxis()->SetTitleFont(43);
     etaPhiEnergyCSPU_Average_p[cI]->GetYaxis()->SetTitleFont(43);
     etaPhiEnergyCSPUFlow_Average_p[cI]->GetYaxis()->SetTitleFont(43);

     etaPhiEnergy_Average_p[cI]->GetYaxis()->SetTitleSize(18);
     etaPhiEnergyCS_Average_p[cI]->GetYaxis()->SetTitleSize(18);
     etaPhiEnergyCSPU_Average_p[cI]->GetYaxis()->SetTitleSize(18);
     etaPhiEnergyCSPUFlow_Average_p[cI]->GetYaxis()->SetTitleSize(18);
  }


  TFile* inFile_p = TFile::Open(mntToXRootdFileString(inFileName).c_str(), "READ");
  TTree* hiEvtTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
  TTree* skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");
  TTree* pfTree_p = (TTree*)inFile_p->Get("pfcandAnalyzer/pfTree");
  TTree* pfTreeCS_p = (TTree*)inFile_p->Get("pfcandAnalyzerCS/pfTree");
  TTree* pfTreeCSPU_p = (TTree*)inFile_p->Get("pfcandAnalyzerCSR4/pfTree");
  TTree* pfTreeCSPUFlow_p = (TTree*)inFile_p->Get("pfcandAnalyzerCSFlowR4/pfTree");

  const Int_t nEntries = hiEvtTree_p->GetEntries();

  TTree* jetTree_p = (TTree*)inFile_p->Get("akCs4PU3PFFlowJetAnalyzer/t");

  unsigned int run_, lumi_;
  unsigned long long evt_;
  Int_t hiBin_;
  Float_t hiHF_;
  Float_t vz_;
  Int_t hiNevtPlane_;
  Float_t hiEvtPlanes_[29];
  
  hiEvtTree_p->SetBranchStatus("*", 0);
  hiEvtTree_p->SetBranchStatus("hiBin", 1);
  hiEvtTree_p->SetBranchStatus("hiHF", 1);
  hiEvtTree_p->SetBranchStatus("vz", 1);
  hiEvtTree_p->SetBranchStatus("hiNevtPlane", 1);
  hiEvtTree_p->SetBranchStatus("hiEvtPlanes", 1);
  hiEvtTree_p->SetBranchStatus("run", 1);
  hiEvtTree_p->SetBranchStatus("lumi", 1);
  hiEvtTree_p->SetBranchStatus("evt", 1);

  hiEvtTree_p->SetBranchAddress("hiBin", &hiBin_);
  hiEvtTree_p->SetBranchAddress("hiHF", &hiHF_);
  hiEvtTree_p->SetBranchAddress("vz", &vz_);
  hiEvtTree_p->SetBranchAddress("hiNevtPlane", &hiNevtPlane_);
  hiEvtTree_p->SetBranchAddress("hiEvtPlanes", hiEvtPlanes_);
  hiEvtTree_p->SetBranchAddress("run", &run_);
  hiEvtTree_p->SetBranchAddress("lumi", &lumi_);
  hiEvtTree_p->SetBranchAddress("evt", &evt_);

  Int_t pprimaryVertexFilter_ = -1;
  Int_t HBHENoiseFilterResultRun2Loose_ = -1;
  Int_t pclusterCompatibilityFilter_ = -1;
  Int_t phfCoincFilter3_ = -1;
  Int_t pBeamScrapingFilter_ = -1;


  skimTree_p->SetBranchStatus("*", 0);
  skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
  if(true){
    skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
    skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
    skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);
  }
  else{
    skimTree_p->SetBranchStatus("pBeamScrapingFilter", 1);
    skimTree_p->SetBranchStatus("pPAprimaryVertexFilter", 1);
  }

  skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
  if(true){
    skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
    skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
    skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);
  }
  else{
    skimTree_p->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter_);
    skimTree_p->SetBranchAddress("pPAprimaryVertexFilter", &pprimaryVertexFilter_);
  }
  
  std::vector<float>* pfPt_p=NULL;
  std::vector<float>* pfPhi_p=NULL;
  std::vector<float>* pfEta_p=NULL;
  std::vector<int>* pfId_p=NULL;

  std::vector<float>* pfPtCS_p=NULL;
  std::vector<float>* pfPhiCS_p=NULL;
  std::vector<float>* pfEtaCS_p=NULL;
  std::vector<int>* pfIdCS_p=NULL;

  std::vector<float>* pfPtCSPU_p=NULL;
  std::vector<float>* pfPhiCSPU_p=NULL;
  std::vector<float>* pfEtaCSPU_p=NULL;
  std::vector<int>* pfIdCSPU_p=NULL;

  std::vector<float>* pfPtCSPUFlow_p=NULL;
  std::vector<float>* pfPhiCSPUFlow_p=NULL;
  std::vector<float>* pfEtaCSPUFlow_p=NULL;
  std::vector<int>* pfIdCSPUFlow_p=NULL;

  pfTree_p->SetBranchStatus("*", 0);
  pfTree_p->SetBranchStatus("pfPt", 1);
  pfTree_p->SetBranchStatus("pfPhi", 1);
  pfTree_p->SetBranchStatus("pfEta", 1);
  pfTree_p->SetBranchStatus("pfId", 1);

  pfTree_p->SetBranchAddress("pfPt", &pfPt_p);
  pfTree_p->SetBranchAddress("pfPhi", &pfPhi_p);
  pfTree_p->SetBranchAddress("pfEta", &pfEta_p);
  pfTree_p->SetBranchAddress("pfId", &pfId_p);


  pfTreeCS_p->SetBranchStatus("*", 0);
  pfTreeCS_p->SetBranchStatus("pfPt", 1);
  pfTreeCS_p->SetBranchStatus("pfPhi", 1);
  pfTreeCS_p->SetBranchStatus("pfEta", 1);
  pfTreeCS_p->SetBranchStatus("pfId", 1);

  pfTreeCS_p->SetBranchAddress("pfPt", &pfPtCS_p);
  pfTreeCS_p->SetBranchAddress("pfPhi", &pfPhiCS_p);
  pfTreeCS_p->SetBranchAddress("pfEta", &pfEtaCS_p);
  pfTreeCS_p->SetBranchAddress("pfId", &pfIdCS_p);


  pfTreeCSPU_p->SetBranchStatus("*", 0);
  pfTreeCSPU_p->SetBranchStatus("pfPt", 1);
  pfTreeCSPU_p->SetBranchStatus("pfPhi", 1);
  pfTreeCSPU_p->SetBranchStatus("pfEta", 1);
  pfTreeCSPU_p->SetBranchStatus("pfId", 1);

  pfTreeCSPU_p->SetBranchAddress("pfPt", &pfPtCSPU_p);
  pfTreeCSPU_p->SetBranchAddress("pfPhi", &pfPhiCSPU_p);
  pfTreeCSPU_p->SetBranchAddress("pfEta", &pfEtaCSPU_p);
  pfTreeCSPU_p->SetBranchAddress("pfId", &pfIdCSPU_p);


  pfTreeCSPUFlow_p->SetBranchStatus("*", 0);
  pfTreeCSPUFlow_p->SetBranchStatus("pfPt", 1);
  pfTreeCSPUFlow_p->SetBranchStatus("pfPhi", 1);
  pfTreeCSPUFlow_p->SetBranchStatus("pfEta", 1);
  pfTreeCSPUFlow_p->SetBranchStatus("pfId", 1);

  pfTreeCSPUFlow_p->SetBranchAddress("pfPt", &pfPtCSPUFlow_p);
  pfTreeCSPUFlow_p->SetBranchAddress("pfPhi", &pfPhiCSPUFlow_p);
  pfTreeCSPUFlow_p->SetBranchAddress("pfEta", &pfEtaCSPUFlow_p);
  pfTreeCSPUFlow_p->SetBranchAddress("pfId", &pfIdCSPUFlow_p);

  const Int_t nMaxJets = 500;
  Int_t nref_;
  Float_t jtpt_[nMaxJets];
  Float_t jteta_[nMaxJets];
  Float_t jtphi_[nMaxJets];
  
  jetTree_p->SetBranchStatus("*", 0);
  jetTree_p->SetBranchStatus("nref", 1);
  jetTree_p->SetBranchStatus("jtpt", 1);
  jetTree_p->SetBranchStatus("jteta", 1);
  jetTree_p->SetBranchStatus("jtphi", 1);

  jetTree_p->SetBranchAddress("nref", &nref_);
  jetTree_p->SetBranchAddress("jtpt", jtpt_);
  jetTree_p->SetBranchAddress("jteta", jteta_);
  jetTree_p->SetBranchAddress("jtphi", jtphi_);

  goodGlobalSelection sel;
  sel.setIsPbPb(true);

  for(Int_t entry = 0; entry < nEntries; ++entry){
    hiEvtTree_p->GetEntry(entry);
    pfTree_p->GetEntry(entry);
    pfTreeCS_p->GetEntry(entry);
    pfTreeCSPU_p->GetEntry(entry);
    pfTreeCSPUFlow_p->GetEntry(entry);

    jetTree_p->GetEntry(entry);
    skimTree_p->GetEntry(entry);

    sel.setVz(vz_);
    sel.setHiHF(hiHF_);
    sel.setPprimaryVertexFilter(pprimaryVertexFilter_);
    sel.setPBeamScrapingFilter(pBeamScrapingFilter_);
    sel.setPhfCoincFilter3(phfCoincFilter3_);
    sel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
    sel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);

    if(!sel.isGood()) continue;

    if(hiBin_ >= 180) continue;
    if(hiNevtPlane_ < 2) continue;
    //    if(TMath::Abs(hiEvtPlanes_[8]) > hfCut) continue;


    std::vector<double> goodJetPhi_;
    std::vector<double> goodJetEta_;

    for(Int_t jI = 0; jI < nref_; ++jI){
      if(jtpt_[jI] < 100) continue;
      if(TMath::Abs(jteta_[jI]) > 2.) continue;

      goodJetPhi_.push_back(jtphi_[jI]);
      goodJetEta_.push_back(jteta_[jI]);
    }

    Int_t centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(hiBin_/2 >= centBinsLow[cI] && hiBin_/2 < centBinsHi[cI]){
	centPos = cI;
	break;
      }
    }

    if(currNum[centPos] >= histNum) continue;

    TH2F* etaPhiEnergy_p = new TH2F(("etaPhiEnergy_Run" + std::to_string(run_) + "_Lumi" + std::to_string(lumi_) + "_Evt" + std::to_string(evt_) + "_p").c_str(), ";#eta;#phi", nEtaBinsRed-1, etaBinsRed, nPhiBinsRed, phiBinsRed);
    TH2F* etaPhiEnergyCS_p = new TH2F(("etaPhiEnergyCS_Run" + std::to_string(run_) + "_Lumi" + std::to_string(lumi_) + "_Evt" + std::to_string(evt_) + "_p").c_str(), ";#eta;#phi", nEtaBinsRed-1, etaBinsRed, nPhiBinsRed, phiBinsRed);
    TH2F* etaPhiEnergyCSPU_p = new TH2F(("etaPhiEnergyCSPU_Run" + std::to_string(run_) + "_Lumi" + std::to_string(lumi_) + "_Evt" + std::to_string(evt_) + "_p").c_str(), ";#eta;#phi", nEtaBinsRed-1, etaBinsRed, nPhiBinsRed, phiBinsRed);
    TH2F* etaPhiEnergyCSPUFlow_p = new TH2F(("etaPhiEnergyCSPUFlow_Run" + std::to_string(run_) + "_Lumi" + std::to_string(lumi_) + "_Evt" + std::to_string(evt_) + "_p").c_str(), ";#eta;#phi", nEtaBinsRed-1, etaBinsRed, nPhiBinsRed, phiBinsRed);

    etaPhiEnergy_p->GetZaxis()->SetTitle("#rho = #Sigma p_{T,PF}/Area [GeV/Area]");
    etaPhiEnergyCS_p->GetZaxis()->SetTitle("#rho = #Sigma p_{T,PF}/Area [GeV/Area]");
    etaPhiEnergyCSPU_p->GetZaxis()->SetTitle("#rho = #Sigma p_{T,PF}/Area [GeV/Area]");
    etaPhiEnergyCSPUFlow_p->GetZaxis()->SetTitle("#rho = #Sigma p_{T,PF}/Area [GeV/Area]");

    

    etaPhiEnergy_p->GetZaxis()->SetLabelFont(43);
    etaPhiEnergy_p->GetZaxis()->SetLabelSize(12);
    etaPhiEnergyCS_p->GetZaxis()->SetLabelFont(43);
    etaPhiEnergyCS_p->GetZaxis()->SetLabelSize(12);
    etaPhiEnergyCSPU_p->GetZaxis()->SetLabelFont(43);
    etaPhiEnergyCSPU_p->GetZaxis()->SetLabelSize(12);
    etaPhiEnergyCSPUFlow_p->GetZaxis()->SetLabelFont(43);
    etaPhiEnergyCSPUFlow_p->GetZaxis()->SetLabelSize(12);

    etaPhiEnergy_p->GetXaxis()->SetTitleFont(43);
    etaPhiEnergy_p->GetXaxis()->SetTitleSize(12);
    etaPhiEnergyCS_p->GetXaxis()->SetTitleFont(43);
    etaPhiEnergyCS_p->GetXaxis()->SetTitleSize(12);
    etaPhiEnergyCSPU_p->GetXaxis()->SetTitleFont(43);
    etaPhiEnergyCSPU_p->GetXaxis()->SetTitleSize(12);
    etaPhiEnergyCSPUFlow_p->GetXaxis()->SetTitleFont(43);
    etaPhiEnergyCSPUFlow_p->GetXaxis()->SetTitleSize(12);

    etaPhiEnergy_p->GetZaxis()->SetTitleOffset(etaPhiEnergy_p->GetZaxis()->GetTitleOffset()*1.12);
    etaPhiEnergyCS_p->GetZaxis()->SetTitleOffset(etaPhiEnergyCS_p->GetZaxis()->GetTitleOffset()*1.12);
    etaPhiEnergyCSPU_p->GetZaxis()->SetTitleOffset(etaPhiEnergyCSPU_p->GetZaxis()->GetTitleOffset()*1.12);
    etaPhiEnergyCSPUFlow_p->GetZaxis()->SetTitleOffset(etaPhiEnergyCSPUFlow_p->GetZaxis()->GetTitleOffset()*1.12);

    centerTitles({etaPhiEnergy_p, etaPhiEnergyCS_p, etaPhiEnergyCSPU_p, etaPhiEnergyCSPUFlow_p});
    
    for(unsigned int pI = 0; pI < pfPt_p->size(); ++pI){    
      if(TMath::Abs(pfEta_p->at(pI)) > etaMax) continue;

      for(unsigned int jI = 0; jI < goodJetPhi_.size(); ++jI){
	if(getDR(pfEta_p->at(pI), pfPhi_p->at(pI), goodJetEta_.at(jI), goodJetPhi_.at(jI)) > .4) continue;

	etaEnergy_ByID_h[centPos][0]->Fill(pfEta_p->at(pI), pfPt_p->at(pI)/(Double_t(histNum)));
	etaEnergy_ByID_h[centPos][pfId_p->at(pI)]->Fill(pfEta_p->at(pI), pfPt_p->at(pI)/(Double_t(histNum)));
      }


      Double_t rotatePhi = pfPhi_p->at(pI) - hiEvtPlanes_[8];
      if(rotatePhi > TMath::Pi()) rotatePhi -= 2*TMath::Pi();
      if(rotatePhi < -TMath::Pi()) rotatePhi += 2*TMath::Pi();

      Int_t etaBinPos = findBinPos(pfEta_p->at(pI), nEtaBinsRed, etaBinsRed);
      Double_t area = TMath::Abs(etaBinsRed[etaBinPos] - etaBinsRed[etaBinPos+1])*2.*TMath::Pi()/((Double_t)nPhiBinsRed);
      //      area = 1.;

      etaPhiEnergy_p->Fill(pfEta_p->at(pI), pfPhi_p->at(pI), pfPt_p->at(pI)/area);     
      etaPhiEnergy_Average_p[centPos]->Fill(pfEta_p->at(pI), rotatePhi, pfPt_p->at(pI)/(Double_t(histNum*area)));
    }

    for(unsigned int pI = 0; pI < pfPtCS_p->size(); ++pI){
      if(TMath::Abs(pfEtaCS_p->at(pI)) > etaMax) continue;

      for(unsigned int jI = 0; jI < goodJetPhi_.size(); ++jI){
	if(getDR(pfEtaCS_p->at(pI), pfPhiCS_p->at(pI), goodJetEta_.at(jI), goodJetPhi_.at(jI)) > .4) continue;

	etaEnergyCS_ByID_h[centPos][0]->Fill(pfEtaCS_p->at(pI), pfPtCS_p->at(pI)/(Double_t(histNum)));
	etaEnergyCS_ByID_h[centPos][pfIdCS_p->at(pI)]->Fill(pfEtaCS_p->at(pI), pfPtCS_p->at(pI)/(Double_t(histNum)));
      }

      Double_t rotatePhi = pfPhiCS_p->at(pI) - hiEvtPlanes_[8];
      if(rotatePhi > TMath::Pi()) rotatePhi -= 2*TMath::Pi();
      if(rotatePhi < -TMath::Pi()) rotatePhi += 2*TMath::Pi();

      Int_t etaBinPos = findBinPos(pfEta_p->at(pI), nEtaBinsRed, etaBinsRed);
      Double_t area = TMath::Abs(etaBinsRed[etaBinPos] - etaBinsRed[etaBinPos+1])*2.*TMath::Pi()/((Double_t)nPhiBinsRed);
      //      area = 1.;

      etaPhiEnergyCS_p->Fill(pfEtaCS_p->at(pI), pfPhiCS_p->at(pI), pfPtCS_p->at(pI)/area);
      etaPhiEnergyCS_Average_p[centPos]->Fill(pfEtaCS_p->at(pI), rotatePhi, pfPtCS_p->at(pI)/(Double_t(histNum*area)));
    }


    for(unsigned int pI = 0; pI < pfPtCSPU_p->size(); ++pI){
      if(TMath::Abs(pfEtaCSPU_p->at(pI)) > etaMax) continue;

      for(unsigned int jI = 0; jI < goodJetPhi_.size(); ++jI){
	if(getDR(pfEtaCSPU_p->at(pI), pfPhiCSPU_p->at(pI), goodJetEta_.at(jI), goodJetPhi_.at(jI)) > .4) continue;

	etaEnergyCSPU_ByID_h[centPos][0]->Fill(pfEtaCSPU_p->at(pI), pfPtCSPU_p->at(pI)/(Double_t(histNum)));
	etaEnergyCSPU_ByID_h[centPos][pfIdCSPU_p->at(pI)]->Fill(pfEtaCSPU_p->at(pI), pfPtCSPU_p->at(pI)/(Double_t(histNum)));
      }

      Double_t rotatePhi = pfPhiCSPU_p->at(pI) - hiEvtPlanes_[8];
      if(rotatePhi > TMath::Pi()) rotatePhi -= 2*TMath::Pi();
      if(rotatePhi < -TMath::Pi()) rotatePhi += 2*TMath::Pi();

      Int_t etaBinPos = findBinPos(pfEta_p->at(pI), nEtaBinsRed, etaBinsRed);
      Double_t area = TMath::Abs(etaBinsRed[etaBinPos] - etaBinsRed[etaBinPos+1])*2.*TMath::Pi()/((Double_t)nPhiBinsRed);
      //      area = 1.;

      etaPhiEnergyCSPU_p->Fill(pfEtaCSPU_p->at(pI), pfPhiCSPU_p->at(pI), pfPtCSPU_p->at(pI)/area);
      etaPhiEnergyCSPU_Average_p[centPos]->Fill(pfEtaCSPU_p->at(pI), rotatePhi, pfPtCSPU_p->at(pI)/(Double_t(histNum*area)));
    }


    for(unsigned int pI = 0; pI < pfPtCSPUFlow_p->size(); ++pI){
      if(TMath::Abs(pfEtaCSPUFlow_p->at(pI)) > etaMax) continue;

      for(unsigned int jI = 0; jI < goodJetPhi_.size(); ++jI){
	if(getDR(pfEtaCSPUFlow_p->at(pI), pfPhiCSPUFlow_p->at(pI), goodJetEta_.at(jI), goodJetPhi_.at(jI)) > .4) continue;

	etaEnergyCSPUFlow_ByID_h[centPos][0]->Fill(pfEtaCSPUFlow_p->at(pI), pfPtCSPUFlow_p->at(pI)/(Double_t(histNum)));
	etaEnergyCSPUFlow_ByID_h[centPos][pfIdCSPUFlow_p->at(pI)]->Fill(pfEtaCSPUFlow_p->at(pI), pfPtCSPUFlow_p->at(pI)/(Double_t(histNum)));
      }

      Double_t rotatePhi = pfPhiCSPUFlow_p->at(pI) - hiEvtPlanes_[8];
      if(rotatePhi > TMath::Pi()) rotatePhi -= 2*TMath::Pi();
      if(rotatePhi < -TMath::Pi()) rotatePhi += 2*TMath::Pi();

      Int_t etaBinPos = findBinPos(pfEta_p->at(pI), nEtaBinsRed, etaBinsRed);
      Double_t area = TMath::Abs(etaBinsRed[etaBinPos] - etaBinsRed[etaBinPos+1])*2.*TMath::Pi()/((Double_t)nPhiBinsRed);
      //      area = 1.;

      etaPhiEnergyCSPUFlow_p->Fill(pfEtaCSPUFlow_p->at(pI), pfPhiCSPUFlow_p->at(pI), pfPtCSPUFlow_p->at(pI)/area);
      etaPhiEnergyCSPUFlow_Average_p[centPos]->Fill(pfEtaCSPUFlow_p->at(pI), rotatePhi, pfPtCSPUFlow_p->at(pI)/(Double_t(histNum*area)));
    }
  
    goodJetEta_.clear();
    goodJetPhi_.clear();
    
    std::vector<TH2*> hist_p = {etaPhiEnergy_p, etaPhiEnergyCS_p, etaPhiEnergyCSPU_p, etaPhiEnergyCSPUFlow_p};
    
    Double_t max = TMath::Max(TMath::Max(etaPhiEnergy_p->GetMaximum(), etaPhiEnergyCS_p->GetMaximum()), TMath::Max(etaPhiEnergyCSPU_p->GetMaximum(), etaPhiEnergyCSPUFlow_p->GetMaximum()));

    etaPhiEnergy_p->SetMaximum(max);
    etaPhiEnergyCS_p->SetMaximum(max);
    etaPhiEnergyCSPU_p->SetMaximum(max);
    etaPhiEnergyCSPUFlow_p->SetMaximum(max);

    etaPhiEnergy_p->SetMinimum(0.0);
    etaPhiEnergyCS_p->SetMinimum(0.0);
    etaPhiEnergyCSPU_p->SetMinimum(0.0);
    etaPhiEnergyCSPUFlow_p->SetMinimum(0.0);
    
    Bool_t isGoodJet = false;
    Int_t jtpos = -1;
    for(Int_t jI = 0; jI < nref_; ++jI){
      if(TMath::Abs(jteta_[jI]) > 1.) continue;
      if(TMath::Abs(jtphi_[jI]) > TMath::Pi()*2./3.) continue;
      if(jtpt_[jI] > 100.){
	isGoodJet = true;
	jtpos = jI;
	break;
      }
    }

    if(hiBin_ < 10 && (eventNum < perEventNum || (jetEventNum < jetPerEventNum && isGoodJet))){
      for(unsigned int hI = 0; hI < hist_p.size(); ++hI){
	TCanvas* temp_p = new TCanvas("temp_p", "", 450, 400);
	temp_p->SetTopMargin(0.01);
	temp_p->SetLeftMargin(temp_p->GetLeftMargin()*1.2);
	temp_p->SetBottomMargin(temp_p->GetLeftMargin());
	temp_p->SetBottomMargin(temp_p->GetRightMargin());
	
	hist_p.at(hI)->GetXaxis()->SetTitleOffset(hist_p.at(hI)->GetXaxis()->GetTitleOffset()*2.5);
	hist_p.at(hI)->GetYaxis()->SetTitleOffset(hist_p.at(hI)->GetYaxis()->GetTitleOffset()*2.5);
	hist_p.at(hI)->GetZaxis()->SetTitleOffset(hist_p.at(hI)->GetZaxis()->GetTitleOffset()*1.5);
	hist_p.at(hI)->DrawCopy("LEGO2 BB FB");
	gPad->SetGrid(0,0);
	gStyle->SetOptStat(0);
	
	std::string saveName = hist_p.at(hI)->GetName();
	saveName = "pdfDir/" + saveName + "_";
	if(isGoodJet) saveName = saveName + "GOODJET_";
	saveName = saveName + std::to_string(date->GetDate()) + ".pdf";

	label_p->DrawLatex(.14, .94, "#bf{CMS Preliminary} 2015 PbPb #sqrt{s_{NN}}=5.02 TeV");
	//	label_p->DrawLatex(.14, .92, "2015 PbPb #sqrt{s_{NN}}=5.02 TeV");

	if(hI == 0) label_p->DrawLatex(.14, .86, "Unsubtracted");
	else if(hI == 1) label_p->DrawLatex(.14, .86, "CS Nominal");
	else if(hI == 2) label_p->DrawLatex(.14, .86, "CS Updated");
	else if(hI == 3) label_p->DrawLatex(.14, .86, "CS Updated + Flow");

	label_p->DrawLatex(.14, .90, ("Single " + prettyString(hiBin_/2., 1, false) + "% Event").c_str());
	
	temp_p->SaveAs(saveName.c_str());
	
	if(isGoodJet) std::cout << "JET: " << jtpt_[jtpos] << ", " << jteta_[jtpos] << ", " << jtphi_[jtpos] << std::endl;

	if(isGoodJet) ++jetEventNum;
	else ++eventNum;

	delete temp_p;
      }
    }

    delete etaPhiEnergy_p;
    delete etaPhiEnergyCS_p;
    delete etaPhiEnergyCSPU_p;
    delete etaPhiEnergyCSPUFlow_p;   

    currNum[centPos]++;

    if(currNum[centPos]%1000 == 0){
      std::cout << "Centrality " << centBinsLow[centPos] << "-" << centBinsHi[centPos] << " filled: " << currNum[centPos] << "/" << histNum << std::endl;
    }

    bool doBreak = true;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(currNum[cI] < histNum){
	doBreak = false;
	break;
      }
    }

    if(doBreak) break;
  }

  inFile_p->Close();
  delete inFile_p;

  std::vector<TH2*> hist_p;
  std::vector<std::string> centLabel;
  std::vector<std::string> subLabel;
  std::vector<Double_t> maxes1;
  std::vector<Double_t> maxes2;

  std::vector<TH1*> histProjNoSub_p;
  std::vector<TH1*> histProjCS_p;
  std::vector<TH1*> histProjCSPU_p;
  std::vector<TH1*> histProjCSPUFlow_p;
  std::vector<std::string> centLabelProj;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    Double_t max = TMath::Max(TMath::Max(etaPhiEnergy_Average_p[cI]->GetMaximum(), etaPhiEnergyCS_Average_p[cI]->GetMaximum()), TMath::Max(etaPhiEnergyCSPU_Average_p[cI]->GetMaximum(), etaPhiEnergyCSPUFlow_Average_p[cI]->GetMaximum()));

    Double_t max2 = TMath::Max(etaPhiEnergyCS_Average_p[cI]->GetMaximum(), TMath::Max(etaPhiEnergyCSPU_Average_p[cI]->GetMaximum(), etaPhiEnergyCSPUFlow_Average_p[cI]->GetMaximum()));

    maxes1.push_back(max);
    maxes1.push_back(max);
    maxes1.push_back(max);
    maxes1.push_back(max);
    maxes2.push_back(max2);
    maxes2.push_back(max2);
    maxes2.push_back(max2);
    maxes2.push_back(max2);

    etaPhiEnergy_Average_p[cI]->SetMaximum(max);
    etaPhiEnergyCS_Average_p[cI]->SetMaximum(max);
    etaPhiEnergyCSPU_Average_p[cI]->SetMaximum(max);
    etaPhiEnergyCSPUFlow_Average_p[cI]->SetMaximum(max);
    etaPhiEnergy_Average_p[cI]->SetMinimum(0.0);
    etaPhiEnergyCS_Average_p[cI]->SetMinimum(0.0);
    etaPhiEnergyCSPU_Average_p[cI]->SetMinimum(0.0);
    etaPhiEnergyCSPUFlow_Average_p[cI]->SetMinimum(0.0);

    hist_p.push_back(etaPhiEnergy_Average_p[cI]);
    hist_p.push_back(etaPhiEnergyCS_Average_p[cI]);   
    hist_p.push_back(etaPhiEnergyCSPU_Average_p[cI]);   
    hist_p.push_back(etaPhiEnergyCSPUFlow_Average_p[cI]);   

    std::string nameNoSub = etaPhiEnergy_Average_p[cI]->GetName();
    std::string nameCS = etaPhiEnergyCS_Average_p[cI]->GetName();
    std::string nameCSPU = etaPhiEnergyCSPU_Average_p[cI]->GetName();
    std::string nameCSPUFlow = etaPhiEnergyCSPUFlow_Average_p[cI]->GetName();

    nameNoSub = nameNoSub + "_ProjPhi";
    nameCS = nameCS + "_ProjPhi";
    nameCSPU = nameCSPU + "_ProjPhi";
    nameCSPUFlow = nameCSPUFlow + "_ProjPhi";

    etaPhiEnergy_Average_PhiProj_p[cI] = (TH1F*)etaPhiEnergy_Average_p[cI]->ProjectionY(nameNoSub.c_str());
    etaPhiEnergyCS_Average_PhiProj_p[cI] = (TH1F*)etaPhiEnergyCS_Average_p[cI]->ProjectionY(nameCS.c_str());
    etaPhiEnergyCSPU_Average_PhiProj_p[cI] = (TH1F*)etaPhiEnergyCSPU_Average_p[cI]->ProjectionY(nameCSPU.c_str());
    etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI] = (TH1F*)etaPhiEnergyCSPUFlow_Average_p[cI]->ProjectionY(nameCSPUFlow.c_str());

    etaPhiEnergy_Average_PhiProj_p[cI]->GetXaxis()->SetTitleOffset(.2);
    etaPhiEnergy_Average_PhiProj_p[cI]->GetYaxis()->SetTitleOffset(etaPhiEnergy_Average_p[cI]->GetYaxis()->GetTitleOffset()*3./4.);
    etaPhiEnergyCS_Average_PhiProj_p[cI]->GetXaxis()->SetTitleOffset(.2);
    etaPhiEnergyCS_Average_PhiProj_p[cI]->GetYaxis()->SetTitleOffset(etaPhiEnergyCS_Average_p[cI]->GetYaxis()->GetTitleOffset()*3./4.);
    etaPhiEnergyCSPU_Average_PhiProj_p[cI]->GetXaxis()->SetTitleOffset(.2);
    etaPhiEnergyCSPU_Average_PhiProj_p[cI]->GetYaxis()->SetTitleOffset(etaPhiEnergyCSPU_Average_p[cI]->GetYaxis()->GetTitleOffset()*3./4.);
    etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->GetXaxis()->SetTitleOffset(.2);
    etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->GetYaxis()->SetTitleOffset(etaPhiEnergyCSPUFlow_Average_p[cI]->GetYaxis()->GetTitleOffset()*3./4.);

    etaPhiEnergy_Average_PhiProj_p[cI]->Scale(1./(etaMax*2.));
    etaPhiEnergyCS_Average_PhiProj_p[cI]->Scale(1./(etaMax*2.));
    etaPhiEnergyCSPU_Average_PhiProj_p[cI]->Scale(1./(etaMax*2.));
    etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->Scale(1./(etaMax*2.));

    Double_t middle = (etaPhiEnergy_Average_PhiProj_p[cI]->GetMaximum() + etaPhiEnergy_Average_PhiProj_p[cI]->GetMinimum())/2.;
    for(Int_t bI = 0; bI < etaPhiEnergy_Average_PhiProj_p[cI]->GetNbinsX(); ++bI){
      etaPhiEnergy_Average_PhiProj_p[cI]->SetBinContent(bI+1, etaPhiEnergy_Average_PhiProj_p[cI]->GetBinContent(bI+1) - middle);
    }

    middle = (etaPhiEnergyCS_Average_PhiProj_p[cI]->GetMaximum() + etaPhiEnergyCS_Average_PhiProj_p[cI]->GetMinimum())/2.;
    for(Int_t bI = 0; bI < etaPhiEnergyCS_Average_PhiProj_p[cI]->GetNbinsX(); ++bI){
      etaPhiEnergyCS_Average_PhiProj_p[cI]->SetBinContent(bI+1, etaPhiEnergyCS_Average_PhiProj_p[cI]->GetBinContent(bI+1) - middle);
    }

    middle = (etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->GetMaximum() + etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->GetMinimum())/2.;
    for(Int_t bI = 0; bI < etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->GetNbinsX(); ++bI){
      etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->SetBinContent(bI+1, etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->GetBinContent(bI+1) - middle);
    }

    middle = (etaPhiEnergyCSPU_Average_PhiProj_p[cI]->GetMaximum() + etaPhiEnergyCSPU_Average_PhiProj_p[cI]->GetMinimum())/2.;
    for(Int_t bI = 0; bI < etaPhiEnergyCSPU_Average_PhiProj_p[cI]->GetNbinsX(); ++bI){
      etaPhiEnergyCSPU_Average_PhiProj_p[cI]->SetBinContent(bI+1, etaPhiEnergyCSPU_Average_PhiProj_p[cI]->GetBinContent(bI+1) - middle);
    }

    middle = (etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->GetMaximum() + etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->GetMinimum())/2.;
    for(Int_t bI = 0; bI < etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->GetNbinsX(); ++bI){
      etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->SetBinContent(bI+1, etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->GetBinContent(bI+1) - middle);
    }

    histProjNoSub_p.push_back(etaPhiEnergy_Average_PhiProj_p[cI]);
    histProjCS_p.push_back(etaPhiEnergyCS_Average_PhiProj_p[cI]);
    histProjCSPU_p.push_back(etaPhiEnergyCSPU_Average_PhiProj_p[cI]);
    histProjCSPUFlow_p.push_back(etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]);

    const std::string centStr = "2015 PbPb " + std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%, #sqrt{s_{NN}}=5.02 TeV";
    centLabel.push_back(centStr);
    centLabel.push_back(centStr);
    centLabel.push_back(centStr);
    centLabel.push_back(centStr);
    subLabel.push_back("#bf{No subtraction}");
    subLabel.push_back("#bf{CS Nominal}");
    subLabel.push_back("#bf{CS Updated}");
    subLabel.push_back("#bf{CS Updated+Flow}");

    centLabelProj.push_back(centStr);
  }

  for(unsigned int hI = 0; hI < hist_p.size(); ++hI){
    TCanvas* temp_p = new TCanvas("temp_p", "", 450, 400);
    temp_p->SetLeftMargin(temp_p->GetLeftMargin()*1.2);
    temp_p->SetBottomMargin(temp_p->GetLeftMargin());
    temp_p->SetRightMargin(temp_p->GetLeftMargin());
    temp_p->SetTopMargin(temp_p->GetBottomMargin()*3./4.);
    
    hist_p.at(hI)->SetMaximum(maxes1.at(hI));

    hist_p.at(hI)->DrawCopy("COLZ");
    gStyle->SetOptStat(0);
   
    label_p->DrawLatex(.11, .94, (centLabel.at(hI)).c_str());
    label_p->SetTextSize(18);
    label_p->DrawLatex(.4, .3, (subLabel.at(hI)).c_str());
    label_p->DrawLatex(.4, .24, ("#LT" + std::to_string(histNum/1000) + "k MB Events#GT").c_str());
    label_p->SetTextSize(14);
    label_p->DrawLatex(.895, .96, "#LT#rho#GT");
    label_p->DrawLatex(.84, .92, "[GeV/Area]");

    std::string saveName = hist_p.at(hI)->GetName();
    saveName = "pdfDir/" + saveName + "_NoSubMax_" + std::to_string(date->GetDate()) + ".pdf";
    temp_p->SaveAs(saveName.c_str());   
    
    delete temp_p;
  }

  for(unsigned int hI = 0; hI < hist_p.size(); ++hI){
    TCanvas* temp_p = new TCanvas("temp_p", "", 450, 400);
    temp_p->SetLeftMargin(temp_p->GetLeftMargin()*1.2);
    temp_p->SetBottomMargin(temp_p->GetLeftMargin());
    temp_p->SetRightMargin(temp_p->GetLeftMargin());
    temp_p->SetTopMargin(temp_p->GetBottomMargin()*3./4.);
    
    hist_p.at(hI)->SetMaximum(maxes2.at(hI));

    hist_p.at(hI)->DrawCopy("COLZ");
    gStyle->SetOptStat(0);
   
    label_p->DrawLatex(.12, .96, "#bf{CMS Preliminary}");
    label_p->DrawLatex(.12, .92, (centLabel.at(hI)).c_str()); 
    label_p->SetTextSize(18);
    label_p->DrawLatex(.4, .3, (subLabel.at(hI)).c_str());
    label_p->DrawLatex(.4, .24, ("#LT" + std::to_string(histNum/1000) + "k MB Events#GT").c_str());
    label_p->SetTextSize(14);
    label_p->DrawLatex(.895, .96, "#LT#rho#GT");
    label_p->DrawLatex(.84, .92, "[GeV/Area]");

    std::string saveName = hist_p.at(hI)->GetName();
    saveName = "pdfDir/" + saveName + "_SubMax_" + std::to_string(date->GetDate()) + ".pdf";
    temp_p->SaveAs(saveName.c_str());   
    
    delete temp_p;
  }

  kirchnerPalette col;

  TLegend* leg_p = new TLegend(0.2, 0.70, 0.4, 0.9);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);
  leg_p->SetBorderSize(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(12);

  for(unsigned int hI = 0; hI < histProjNoSub_p.size(); ++hI){
    TCanvas* temp_p = new TCanvas("temp_p", "", 450, 400);
    temp_p->SetLeftMargin(temp_p->GetLeftMargin()*1.2);
    temp_p->SetBottomMargin(temp_p->GetLeftMargin());
    temp_p->SetRightMargin(0.01);
    temp_p->SetTopMargin(temp_p->GetBottomMargin()*3./4.);
    
    histProjNoSub_p.at(hI)->SetMarkerSize(0.8);
    histProjNoSub_p.at(hI)->SetMarkerStyle(20);
    histProjNoSub_p.at(hI)->SetMarkerColor(col.getColor(0));
    histProjNoSub_p.at(hI)->SetLineColor(col.getColor(0));

    histProjCS_p.at(hI)->SetMarkerSize(0.8);
    histProjCS_p.at(hI)->SetMarkerStyle(21);
    histProjCS_p.at(hI)->SetMarkerColor(col.getColor(1));
    histProjCS_p.at(hI)->SetLineColor(col.getColor(1));

    histProjCSPU_p.at(hI)->SetMarkerSize(0.8);
    histProjCSPU_p.at(hI)->SetMarkerStyle(33);
    histProjCSPU_p.at(hI)->SetMarkerColor(col.getColor(3));
    histProjCSPU_p.at(hI)->SetLineColor(col.getColor(3));

    histProjCSPUFlow_p.at(hI)->SetMarkerSize(0.8);
    histProjCSPUFlow_p.at(hI)->SetMarkerStyle(34);
    histProjCSPUFlow_p.at(hI)->SetMarkerColor(col.getColor(2));
    histProjCSPUFlow_p.at(hI)->SetLineColor(col.getColor(2));

    if(hI == 0){
      //      leg_p->AddEntry(histProjNoSub_p.at(hI), "Unsubtracted", "P L");
      leg_p->AddEntry(histProjCS_p.at(hI), "CS Nominal", "P L");
      leg_p->AddEntry(histProjCSPU_p.at(hI), "CS Updated", "P L");
      leg_p->AddEntry(histProjCSPUFlow_p.at(hI), "CS Updated+Flow", "P L");
    }

    histProjCS_p.at(hI)->GetXaxis()->SetTitleOffset(hist_p.at(0)->GetXaxis()->GetTitleOffset()*1.);
    histProjCS_p.at(hI)->GetYaxis()->SetTitle("#LT#rho(#phi - #Psi_{HF,2}) - #rho_{0}#GT [GeV/Area]");
    histProjCS_p.at(hI)->GetYaxis()->SetTitleOffset(hist_p.at(0)->GetYaxis()->GetTitleOffset()*3./4.);
    
    centerTitles(histProjCS_p.at(hI));

    //    histProjNoSub_p.at(hI)->DrawCopy("E1 P");
    histProjCS_p.at(hI)->SetMaximum((histProjNoSub_p.at(hI)->GetMaximum() + histProjCS_p.at(hI)->GetMaximum())/2.);
    histProjCS_p.at(hI)->SetMinimum((histProjNoSub_p.at(hI)->GetMinimum() + histProjCS_p.at(hI)->GetMinimum())/2.);

    histProjCS_p.at(hI)->DrawCopy("E1 P");
    histProjCSPU_p.at(hI)->DrawCopy("E1 P SAME");
    histProjCSPUFlow_p.at(hI)->DrawCopy("E1 P SAME");
    leg_p->Draw("SAME");
    gStyle->SetOptStat(0);
    
    label_p->DrawLatex(.12, .94, (centLabelProj.at(hI)).c_str());

    std::string saveName = histProjNoSub_p.at(hI)->GetName();
    saveName = "pdfDir/" + saveName + "_" + std::to_string(date->GetDate()) + ".pdf";
    temp_p->SaveAs(saveName.c_str());   
    delete temp_p;
  }

  delete leg_p;
  delete label_p;

  outFile_p->cd();
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    for(Int_t pfI = 0; pfI < nPFID; ++pfI){
      etaEnergy_ByID_h[cI][pfI]->Write("", TObject::kOverwrite);
      etaEnergyCS_ByID_h[cI][pfI]->Write("", TObject::kOverwrite);
      etaEnergyCSPU_ByID_h[cI][pfI]->Write("", TObject::kOverwrite);
      etaEnergyCSPUFlow_ByID_h[cI][pfI]->Write("", TObject::kOverwrite);

      delete etaEnergy_ByID_h[cI][pfI];
      delete etaEnergyCS_ByID_h[cI][pfI];
      delete etaEnergyCSPU_ByID_h[cI][pfI];
      delete etaEnergyCSPUFlow_ByID_h[cI][pfI];
    }

    etaPhiEnergy_Average_p[cI]->Write("", TObject::kOverwrite);
    etaPhiEnergyCS_Average_p[cI]->Write("", TObject::kOverwrite);
    etaPhiEnergyCSPU_Average_p[cI]->Write("", TObject::kOverwrite);
    etaPhiEnergyCSPUFlow_Average_p[cI]->Write("", TObject::kOverwrite);

    etaPhiEnergy_Average_PhiProj_p[cI]->Write("", TObject::kOverwrite);
    etaPhiEnergyCS_Average_PhiProj_p[cI]->Write("", TObject::kOverwrite);
    etaPhiEnergyCSPU_Average_PhiProj_p[cI]->Write("", TObject::kOverwrite);
    etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI]->Write("", TObject::kOverwrite);

    delete etaPhiEnergy_Average_p[cI];
    delete etaPhiEnergyCS_Average_p[cI];
    delete etaPhiEnergyCSPU_Average_p[cI];
    delete etaPhiEnergyCSPUFlow_Average_p[cI];

    delete etaPhiEnergy_Average_PhiProj_p[cI];
    delete etaPhiEnergyCS_Average_PhiProj_p[cI];
    delete etaPhiEnergyCSPU_Average_PhiProj_p[cI];
    delete etaPhiEnergyCSPUFlow_Average_PhiProj_p[cI];
  }

  outFile_p->cd();

  TNamed histNumName("histNum", std::to_string(histNum).c_str());
  histNumName.Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  delete date;


  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./makeSubtractedEvent.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0; 
  retVal += makeSubtractedEvent(argv[1]);
  return retVal;
}
