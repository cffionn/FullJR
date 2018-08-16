//cpp dependencies
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TNamed.h"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/doGlobalDebug.h"
#include "Utility/include/etaPhiFunc.h"
#include "Utility/include/findBinPos.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/getLogBins.h"
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/pseudoTowerGeometry.h"
#include "Utility/include/readJSON.h"
#include "Utility/include/returnFileList.h"
#include "Utility/include/vectorStringUtility.h"

int makeSubtractedPPEvent(const std::string inFileName)
{
  const std::string fullPath = std::getenv("FULLJRDIR");
  readJSON ppJson((fullPath + "/RandomCones/tables/Cert_262081-262273_5TeV_PromptReco_Collisions15_25ns_JSON_v2.txt").c_str());
  ppJson.Print();

  std::vector<std::string> listOfFiles;

  if(checkFile(inFileName) && inFileName.find(".root") != std::string::npos) listOfFiles.push_back(inFileName);
  else if(checkDir(inFileName)){
    listOfFiles = returnFileList(inFileName, ".root");
    removeVectorDuplicates(&listOfFiles);
    
    unsigned int pos = 0;
    while(pos < listOfFiles.size()){
      if(listOfFiles.at(pos).find("/failed/") != std::string::npos) listOfFiles.erase(listOfFiles.begin()+pos);
      else ++pos;
    }
  }

  if(listOfFiles.size() == 0){
    std::cout << "Given inFileName \'" << inFileName << "\' is not a valid root file nor a directory leading to root files. return 1" << std::endl;;
    return 1;
  }

  TDatime* date = new TDatime();

  checkMakeDir("output");
  std::string outFileName = listOfFiles.at(0);
  while(outFileName.find("/") != std::string::npos){
    outFileName.replace(0, outFileName.find("/")+1, "");
  }
  outFileName.replace(outFileName.find(".root"), 5, "");
  outFileName = "output/" + outFileName + "_PPSubtracted_" + std::to_string(date->GetDate()) + ".root";

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

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  
  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);

  const Int_t perEventNum = 5;
  const Int_t jetPerEventNum = 5;
  const Int_t histNum = 20000000;

  Int_t eventNum = 0;
  Int_t jetEventNum = 0;
  Int_t currNum = 0;

  const Int_t nPFID = 8;

  TH1F* etaEnergy_ByID_h[nPFID];

  for(Int_t pfI = 0; pfI < nPFID; ++pfI){
    std::string pfStr = "PFIDAll";
    if(pfI != 0) pfStr = "PFID" + std::to_string(pfI);
    
    etaEnergy_ByID_h[pfI] = new TH1F(("etaEnergy_PP_" + pfStr +"_p").c_str(), ";#eta;Ave", nEtaBins-1, etaBins);
  }


  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  TH2F* etaPhiEnergy_Average_p = new TH2F("etaPhiEnergy_Average_PP_p", ";#eta;#phi", nEtaBinsRed-1, etaBinsRed, nPhiBinsRed, phiBinsRed);
  centerTitles({etaPhiEnergy_Average_p});
  setSumW2({etaPhiEnergy_Average_p});

  etaPhiEnergy_Average_p->GetXaxis()->SetTitleFont(43);
  etaPhiEnergy_Average_p->GetXaxis()->SetTitleSize(18);  
  etaPhiEnergy_Average_p->GetYaxis()->SetTitleFont(43);
  etaPhiEnergy_Average_p->GetYaxis()->SetTitleSize(18);

  TTree* skimTree_p = NULL;
  TTree* pfTree_p = NULL;
  TTree* hltTree_p = NULL;
  TTree* jetTree_p = NULL;
  TTree* trackTree_p = NULL;

  int run_, lumi_;
  unsigned long long evt_;
  
  Int_t pprimaryVertexFilter_ = -1;
  Int_t HBHENoiseFilterResultRun2Loose_ = -1;
  Int_t pBeamScrapingFilter_ = -1;
  
  const Int_t nMaxPFPart = 100000;
  
  Int_t nPFpart_;
  Float_t pfPt_[nMaxPFPart];
  Float_t pfPhi_[nMaxPFPart];
  Float_t pfEta_[nMaxPFPart];
  Int_t pfId_[nMaxPFPart];

  const Int_t nMaxJets = 500;
  Int_t nref_;
  Float_t jtpt_[nMaxJets];
  Float_t jteta_[nMaxJets];
  Float_t jtphi_[nMaxJets];
  
  const Int_t nVtxMax = 100;
  Int_t nVtx_;
  Float_t zVtx_[nVtxMax];

  goodGlobalSelection sel;
  sel.setIsPbPb(true);
  
  bool doBreak = false;

  for(unsigned int fI = 0; fI < listOfFiles.size(); ++fI){
    std::cout << "Processing file " << fI << "/" << listOfFiles.size() <<  ": \'" << listOfFiles.at(fI) << "\'" << std::endl;

    TFile* inFile_p = TFile::Open(mntToXRootdFileString(listOfFiles.at(fI)).c_str(), "READ");
    skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");
    pfTree_p = (TTree*)inFile_p->Get("pfcandAnalyzer/pfTree");
    hltTree_p = (TTree*)inFile_p->Get("hltanalysis/HltTree");
    jetTree_p = (TTree*)inFile_p->Get("ak4PFJetAnalyzer/t");       
    trackTree_p = (TTree*)inFile_p->Get("ppTrack/trackTree");       

    const Int_t nEntries = skimTree_p->GetEntries();    

    if(nEntries == 0){
      inFile_p->Close();
      delete inFile_p;
      continue;
    }

    hltTree_p->SetBranchStatus("*", 0);
    hltTree_p->SetBranchStatus("Run", 1);
    hltTree_p->SetBranchStatus("LumiBlock", 1);
    hltTree_p->SetBranchStatus("Event", 1);
    
    hltTree_p->SetBranchAddress("Run", &run_);
    hltTree_p->SetBranchAddress("LumiBlock", &lumi_);
    hltTree_p->SetBranchAddress("Event", &evt_);
    
    skimTree_p->SetBranchStatus("*", 0);
    skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    skimTree_p->SetBranchStatus("pBeamScrapingFilter", 1);
    skimTree_p->SetBranchStatus("pPAprimaryVertexFilter", 1);
    
    skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
    skimTree_p->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter_);
    skimTree_p->SetBranchAddress("pPAprimaryVertexFilter", &pprimaryVertexFilter_);
         
    pfTree_p->SetBranchStatus("*", 0);
    pfTree_p->SetBranchStatus("nPFpart", 1);
    pfTree_p->SetBranchStatus("pfPt", 1);
    pfTree_p->SetBranchStatus("pfPhi", 1);
    pfTree_p->SetBranchStatus("pfEta", 1);
    pfTree_p->SetBranchStatus("pfId", 1);
    
    pfTree_p->SetBranchAddress("nPFpart", &nPFpart_);
    pfTree_p->SetBranchAddress("pfPt", pfPt_);
    pfTree_p->SetBranchAddress("pfPhi", pfPhi_);
    pfTree_p->SetBranchAddress("pfEta", pfEta_);
    pfTree_p->SetBranchAddress("pfId", pfId_);
    
    jetTree_p->SetBranchStatus("*", 0);
    jetTree_p->SetBranchStatus("nref", 1);
    jetTree_p->SetBranchStatus("jtpt", 1);
    jetTree_p->SetBranchStatus("jteta", 1);
    jetTree_p->SetBranchStatus("jtphi", 1);
    
    jetTree_p->SetBranchAddress("nref", &nref_);
    jetTree_p->SetBranchAddress("jtpt", jtpt_);
    jetTree_p->SetBranchAddress("jteta", jteta_);
    jetTree_p->SetBranchAddress("jtphi", jtphi_);
    
    trackTree_p->SetBranchStatus("*", 0);
    trackTree_p->SetBranchStatus("nVtx", 1);
    trackTree_p->SetBranchStatus("zVtx", 1);

    trackTree_p->SetBranchAddress("nVtx", &nVtx_);
    trackTree_p->SetBranchAddress("zVtx", zVtx_);

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


    for(Int_t entry = 0; entry < nEntries; ++entry){
      pfTree_p->GetEntry(entry);
      jetTree_p->GetEntry(entry);
      skimTree_p->GetEntry(entry);
      hltTree_p->GetEntry(entry);
      trackTree_p->GetEntry(entry);

      if(nVtx_ == 0){
	std::cout << "WARNING nVTX == 0. continue" << std::endl;
	continue;
      }

      if(!ppJson.isGoodEvent(run_, lumi_)) continue;

      sel.setVz(zVtx_[0]);
      sel.setHiHF(0);
      sel.setPprimaryVertexFilter(pprimaryVertexFilter_);
      sel.setPBeamScrapingFilter(pBeamScrapingFilter_);
      sel.setPhfCoincFilter3(1);
      sel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
      sel.setPclusterCompatibilityFilter(1);

      if(!sel.isGood()) continue;

      std::vector<double> goodJetPhi_;
      std::vector<double> goodJetEta_;

      for(Int_t jI = 0; jI < nref_; ++jI){
	if(jtpt_[jI] < 100) continue;
	if(TMath::Abs(jteta_[jI]) > 2.) continue;
	
	goodJetPhi_.push_back(jtphi_[jI]);
	goodJetEta_.push_back(jteta_[jI]);
      }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
      
      TH2F* etaPhiEnergy_p = new TH2F(("etaPhiEnergy_Run" + std::to_string(run_) + "_Lumi" + std::to_string(lumi_) + "_Evt" + std::to_string(evt_) + "_p").c_str(), ";#eta;#phi", nEtaBinsRed-1, etaBinsRed, nPhiBinsRed, phiBinsRed);
      
      etaPhiEnergy_p->GetZaxis()->SetTitle("#rho = #Sigma p_{T,PF}/Area [GeV/Area]");
      etaPhiEnergy_p->GetZaxis()->SetLabelFont(43);
      etaPhiEnergy_p->GetZaxis()->SetLabelSize(12);
      etaPhiEnergy_p->GetXaxis()->SetTitleFont(43);
      etaPhiEnergy_p->GetXaxis()->SetTitleSize(12);
      etaPhiEnergy_p->GetZaxis()->SetTitleOffset(etaPhiEnergy_p->GetZaxis()->GetTitleOffset()*1.12);
      centerTitles({etaPhiEnergy_p});
      
      for(int pI = 0; pI < nPFpart_; ++pI){    
	if(TMath::Abs(pfEta_[pI]) > etaMax) continue;
	
	for(unsigned int jI = 0; jI < goodJetPhi_.size(); ++jI){
	  if(getDR(pfEta_[pI], pfPhi_[pI], goodJetEta_.at(jI), goodJetPhi_.at(jI)) > .4) continue;
	  
	  etaEnergy_ByID_h[0]->Fill(pfEta_[pI], pfPt_[pI]/(Double_t(histNum)));
	  etaEnergy_ByID_h[pfId_[pI]]->Fill(pfEta_[pI], pfPt_[pI]/(Double_t(histNum)));
	}
	
	
	Double_t rotatePhi = pfPhi_[pI];
	
	Int_t etaBinPos = findBinPos(pfEta_[pI], nEtaBinsRed, etaBinsRed);
	Double_t area = TMath::Abs(etaBinsRed[etaBinPos] - etaBinsRed[etaBinPos+1])*2.*TMath::Pi()/((Double_t)nPhiBinsRed);
	//      area = 1.;
	
	etaPhiEnergy_p->Fill(pfEta_[pI], pfPhi_[pI], pfPt_[pI]/area);     
	etaPhiEnergy_Average_p->Fill(pfEta_[pI], rotatePhi, pfPt_[pI]/(Double_t(histNum*area)));
      }
      
      goodJetEta_.clear();
      goodJetPhi_.clear();
    
      std::vector<TH2*> hist_p = {etaPhiEnergy_p};
      
      etaPhiEnergy_p->SetMinimum(0.0);
      
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
      
      if(eventNum < perEventNum || (jetEventNum < jetPerEventNum && isGoodJet)){
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
	  
	  label_p->DrawLatex(.14, .90, "Single PP Event");
	  
	  quietSaveAs(temp_p, saveName);
	  
	  if(isGoodJet) std::cout << "JET: " << jtpt_[jtpos] << ", " << jteta_[jtpos] << ", " << jtphi_[jtpos] << std::endl;
	  
	  if(isGoodJet) ++jetEventNum;
	  else ++eventNum;
	  
	  delete temp_p;
	}
      }
      
      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
      
      delete etaPhiEnergy_p;
      
      currNum++;
      
      if(currNum%100000 == 0){
	std::cout << "Events filled: " << currNum << "/" << histNum << std::endl;
      }
      
      if(currNum > histNum) doBreak = true;      
      if(doBreak) break;
    }
       
    inFile_p->Close();
    delete inFile_p;

    if(doBreak) break;
  }


  delete label_p;
  
  outFile_p->cd();  
  
  for(Int_t pfI = 0; pfI < nPFID; ++pfI){
    etaEnergy_ByID_h[pfI]->Write("", TObject::kOverwrite);
    delete etaEnergy_ByID_h[pfI];
  }

  
  etaPhiEnergy_Average_p->Write("", TObject::kOverwrite);
  
  delete etaPhiEnergy_Average_p;


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
    std::cout << "Usage: ./makeSubtractedPPEvent.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0; 
  retVal += makeSubtractedPPEvent(argv[1]);
  return retVal;
}
