#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "fastjet/ClusterSequence.hh"
#include "TDatime.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TLorentzVector.h"
#include "Pythia8/Pythia.h"
#include "AMPT.h"
using namespace std;
using namespace fastjet;
using namespace Pythia8;
//--------------------------------
// Declare functions                                
//--------------------------------

float GetSphericity(vector<TLorentzVector> tracks);
int GetCharge(int id);
float Pythiadecaymass(int aPid);
bool IsFinal(int aPid);

struct ParticleInfo{
  int Id;
  TLorentzVector vec;
  int MoId;
};
auto wrap = [](double phi) {
    while (phi >  M_PI) phi -= 2 * M_PI;
    while (phi < -M_PI) phi += 2 * M_PI;
    return phi;
};
//=============================================================================
int main(int argc, char* argv[])
{
  if (argc!=3) {
	cout<<"Parameter number wrong!"<<endl;
	return 1;
  }

  TString inputFile;
  TString outputFile;
  if (argc==3) {
    inputFile = argv[1];
    outputFile = argv[2];
	cout<<"Input file:"<<inputFile<<endl;
	cout<<"Output file:"<<outputFile<<endl;
  }  

  //----------------------------------
  // Open files and add to chain
  //----------------------------------

  int fileNumber = 0;
  char fileList[512];
  TChain* chain = new TChain("AMPT");
  if (inputFile.Contains(".list"))  {
    ifstream* inputStream = new ifstream;
    inputStream->open(inputFile);
    if (!(inputStream)) {
      cout<<"can not open file list"<<endl;
      return 1;
    }//if
    while (inputStream->good()) {
      inputStream->getline(fileList, 512);
      if (inputStream->eof()) break;
      TFile *fTmp = new TFile(fileList);
      
      if (!fTmp || !(fTmp->IsOpen()) || !(fTmp->GetNkeys())) {
    cout<<"open file list error"<<endl;
    return 1;
    } else {
    cout<<"reading file "<<fileList<<endl;
    chain->Add(fileList);
    fileNumber++;
    }//if else
      delete fTmp;
  }//while
    cout<<fileNumber<<" files read in"<<endl;
  } else if (inputFile.Contains(".root")) {
    chain->Add(inputFile.Data());
  }//if else inputFile

  // Main Pythia object for managing the cascade evolution and particle decays.
  Pythia pythiaMain;
  Event& eventMain = pythiaMain.event;
  Rndm& rndm = pythiaMain.rndm;
  //double mp = pythiaMain.particleData.m0(2212);
  // Prepare to do decays but no hard processes.
  pythiaMain.readString("ProcessLevel:all = off");
  pythiaMain.readString("HadronLevel:Decay = on");
  pythiaMain.readString("ParticleDecays:limitTau0 = on");
  pythiaMain.readString("ParticleDecays:tau0Max = 10");
  pythiaMain.readString("111:mayDecay = off");
  pythiaMain.readString("211:mayDecay = off");
  pythiaMain.readString("2212:mayDecay  = off");
  pythiaMain.readString("2112:mayDecay  = off");
  pythiaMain.readString("321:mayDecay = off");
  pythiaMain.readString("130:mayDecay = off");
  pythiaMain.readString("310:mayDecay = off");
  pythiaMain.readString("3122:mayDecay = off");
  pythiaMain.readString("3312:mayDecay = off");
  pythiaMain.readString("3334:mayDecay = off");
  pythiaMain.readString("411:mayDecay = off");
  pythiaMain.readString("421:mayDecay = off");
  pythiaMain.readString("431:mayDecay = off");
  pythiaMain.readString("4122:mayDecay = off");
  pythiaMain.readString("4132:mayDecay = off");
  pythiaMain.readString("4232:mayDecay = off");
  pythiaMain.readString("4332:mayDecay = off");
  pythiaMain.readString("511:mayDecay = off");
  pythiaMain.readString("521:mayDecay = off");
  pythiaMain.readString("531:mayDecay = off");
//  pythiaMain.readString("541:mayDecay = off");
  pythiaMain.readString("5122:mayDecay = off");
  pythiaMain.readString("5132:mayDecay = off");
  pythiaMain.readString("5232:mayDecay = off");
  pythiaMain.readString("5332:mayDecay = off");

  // Initialize.
  pythiaMain.init();    



  //--------------------------------
  // define histograms
  //--------------------------------

	TList *list = new TList();

  TH1D     *hTrials = new     TH1D("hTrials", "", 1, 0., 1.); list->Add(hTrials);
  TProfile  *hXsect = new TProfile("hXsect",  "", 1, 0., 1.); list->Add(hXsect);

	//event-wise variables
	TH1D *hMultInel = new TH1D("hMultInel", "charged multiplicity |eta|<1 for INEL>0;N_{ch};P(N_{ch})", 100, 0, 200); hMultInel->Sumw2(); list->Add(hMultInel);

	TH2D *hFwdVsMid = new TH2D("hFwdVsMid", ";N_{trk}^{Mid} |#eta|<0.5; N_{trk}^{Fwd} V0 region", 500, -0.5, 499.5, 500, -0.5, 499.5); hFwdVsMid->Sumw2(); list->Add(hFwdVsMid);

	TH1D *hNtrkV0 = new TH1D("hNtrkV0", ";N_{trk}; N_{evt}", 500, -0.5, 499.5); hNtrkV0->Sumw2(); list->Add(hNtrkV0);

	TH1D *histEventCount = new TH1D("histEventCount", ";N_{trk}; N_{evt}", 100, 0, 100); list->Add(histEventCount);		
	TH1D *histPt = new TH1D("histPt", "pT Spectrum; pT (GeV/c); Counts", 100, 0, 100); histPt->Sumw2(); list->Add(histPt);	
	

	TH2D *hnMPIVsNtrkV0 = new TH2D("hnMPIVsNtrkV0", ";N_{trk}; n_{MPI}", 500, -0.5, 499.5, 100, -0.5, 99.5); hnMPIVsNtrkV0->Sumw2(); list->Add(hnMPIVsNtrkV0);

  // const int nptbins=30;
  // double ptbin[nptbins+1]={0};
  // double ptmin=0;
  // for(int ibin=0; ibin<=nptbins; ibin++){
  //   double step=0.25;
  //   if(ibin<10) 
  //     ptbin[ibin]=ptmin+ibin*step;
  //   if(ibin>=10&&ibin<20) {
  //     step=0.5;
  //     ptbin[ibin]=ptmin+0.25*10+(ibin-10)*step;
  //   }
  //   if(ibin>=20) {
  //     step=1;
  //     ptbin[ibin]=ptmin+0.25*10+0.5*10+(ibin-20)*step;
  //   }
  // }

const int nptbins = 19; // 总的 bin 数量
double ptbin[nptbins + 1] = {0};
double ptmin = 0.0;

for (int ibin = 0; ibin <= nptbins; ibin++) {
    if (ibin <= 10) { // 0.1, 0~1
        ptbin[ibin] = ptmin + ibin * 0.1;
    } else if (ibin <= 15) { // 0.2, 1~2
        ptbin[ibin] = 1. + (ibin - 10) * 0.2;
    } else if (ibin <= 17) { // 0.5, 2~3
        ptbin[ibin] = 2.0 + (ibin - 15) * 0.5;
    } else if (ibin <= 19) { // 1.0, 3~5
        ptbin[ibin] = 3.0 + (ibin - 17) * 1.0;
    } 
}

  const int nNchbins=500;
  double Nchbin[nNchbins+1]={0};
  double Nchmin=0;
  for(int ibin=0; ibin<=nNchbins; ibin++){
    double step=1;
    Nchbin[ibin]=Nchmin+ibin*step;
  }

  const int nNTbins=50;
  double NTbin[nNTbins+1]={0};
  double NTmin=-0.5;
  for(int ibin=0; ibin<=nNTbins; ibin++){
    double step=1;
    NTbin[ibin]=NTmin+ibin*step;
  }

	TH1F *nTHist = new TH1F("nTHist", "Charged Particles in Transverse Region", nNTbins, NTbin); list->Add(nTHist);
	TH1F *RT_spectrum = new TH1F("RT_spectrum", "RT Spectrum", 75, 0, 15);list->Add(RT_spectrum);

	//inclusive distribution
	TH1D *hChEta = new TH1D("hChEta", ";#eta; N_{ch}", 50, -10, 10); list->Add(hChEta);
	//TH1D *hChEta_decay = new TH1D("hChEta_decay", ";#eta; N_{ch}", 50, -10, 10); list->Add(hChEta_decay);
	TH1D *h_pion_eta = new TH1D("h_pion_eta", ";y; dN/dy", 50, -5, 5); list->Add(h_pion_eta);
	TH1D *h_pion_pt = new TH1D("h_pion_pt", ";p_{T}; dN/dp_{T}", nptbins, ptbin); list->Add(h_pion_pt);
	TH1D *h_kaon_eta = new TH1D("h_kaon_eta", ";y; dN/dy", 50, -5, 5); list->Add(h_kaon_eta);
	TH1D *h_kaon_pt = new TH1D("h_kaon_pt", ";p_{T}; dN/dp_{T}", nptbins, ptbin); list->Add(h_kaon_pt);
	TH1D *h_proton_eta = new TH1D("h_proton_eta", ";y; dN/dy", 50, -5, 5); list->Add(h_proton_eta);
	TH1D *h_proton_pt = new TH1D("h_proton_pt", ";p_{T}; dN/dp_{T}", nptbins, ptbin); list->Add(h_proton_pt);

	TH3D *hKCh_dPhi0 = new TH3D("hKCh_dPhi0", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hKCh_dPhi0->Sumw2(); list->Add(hKCh_dPhi0);
	TH3D *hKCh_dPhi1 = new TH3D("hKCh_dPhi1", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hKCh_dPhi1->Sumw2(); list->Add(hKCh_dPhi1);
	TH3D *hKCh_dPhi2 = new TH3D("hKCh_dPhi2", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hKCh_dPhi2->Sumw2(); list->Add(hKCh_dPhi2);

// ====== 新增 TR jet tree ======
TTree *tTRjet = new TTree("tTRjet", "TR jets per event");
int    tNchTR;                      // 本事件 TR 带电粒子数
int    tNjetTR;                     // 本事件 TR jet 数
vector<double> tJetPt, tJetEta, tJetPhi;
vector<int>    tJetLeadPID;         // leading constituent PID
vector<vector<int>> tJetConstPID;   // 每个 jet 内部全部 PID
vector<int> tJetConstPID_flattened;  // 扁平化存储所有 Constituent PID
vector<int> tJetConstN;              // 每个 jet 的 Constituent 数量
tJetConstPID.push_back(std::vector<int>{0});  // 占位
tTRjet->Branch("NchTR", &tNchTR);
tTRjet->Branch("NjetTR", &tNjetTR);
tTRjet->Branch("jetPt", &tJetPt);
tTRjet->Branch("jetEta", &tJetEta);
tTRjet->Branch("jetPhi", &tJetPhi);
tTRjet->Branch("jetLeadPID", &tJetLeadPID);
tTRjet->Branch("tJetConstPID_flattened", &tJetConstPID_flattened);
tTRjet->Branch("tJetConstN", &tJetConstN);

// ====== 新增 Toward jet tree ======
int    tNchTow; 
int    tNjetTow;
vector<double> tJetTowPt, tJetTowEta, tJetTowPhi;
vector<int>    tJetTowLeadPID;
vector<int> tJetTowConstPID_flattened;
vector<int> tJetTowConstN;
tTRjet->Branch("NchTow", &tNchTow);
tTRjet->Branch("NjetTow", &tNjetTow);
tTRjet->Branch("jetTowPt", &tJetTowPt);
tTRjet->Branch("jetTowEta", &tJetTowEta);
tTRjet->Branch("jetTowPhi", &tJetTowPhi);
tTRjet->Branch("jetTowLeadPID", &tJetTowLeadPID);
tTRjet->Branch("tJetTowConstPID_flattened", &tJetTowConstPID_flattened);
tTRjet->Branch("tJetTowConstN", &tJetTowConstN);


// ---- TR 喷注 ----
TH1D *hTRnJet = new TH1D("hTRnJet", ";N_{jet}^{TR};Events", 11, -0.5, 10.5); list->Add(hTRnJet);
TH1D *hTRjetPt = new TH1D("hTRjetPt", ";p_{T}^{jet} [GeV/c];dN/dp_{T}", 60, 0, 60); list->Add(hTRjetPt);
TH1D *hTRjetEta = new TH1D("hTRjetEta", ";#eta^{jet};dN/d#eta", 50, -1, 1); list->Add(hTRjetEta);
TH1D *hTRjetPhi = new TH1D("hTRjetPhi", ";#phi^{jet} [rad];dN/d#phi", 70, -3.5, 3.5); list->Add(hTRjetPhi);
TH1D *hJetPid = new TH1D("hJetPid", ";particle ID;Counts", 8001, -4000.5, 4000.5); list->Add(hJetPid);
  
	TProfile *pKAvgPt_NT0 = new TProfile("pKAvgPt_NT0", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKAvgPt_NT1 = new TProfile("pKAvgPt_NT1", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKAvgPt_NT2 = new TProfile("pKAvgPt_NT2", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKAvgPt_NT0_ptcut = new TProfile("pKAvgPt_NT0_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKAvgPt_NT1_ptcut = new TProfile("pKAvgPt_NT1_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKAvgPt_NT2_ptcut = new TProfile("pKAvgPt_NT2_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pKAvgPt_NT0->Sumw2(); list->Add(pKAvgPt_NT0);
	pKAvgPt_NT1->Sumw2(); list->Add(pKAvgPt_NT1);
	pKAvgPt_NT2->Sumw2(); list->Add(pKAvgPt_NT2);
	pKAvgPt_NT0_ptcut->Sumw2(); list->Add(pKAvgPt_NT0_ptcut);
	pKAvgPt_NT1_ptcut->Sumw2(); list->Add(pKAvgPt_NT1_ptcut);
	pKAvgPt_NT2_ptcut->Sumw2(); list->Add(pKAvgPt_NT2_ptcut);

	TH3D *hKCh_dPhi1min = new TH3D("hKCh_dPhi1min", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hKCh_dPhi1min->Sumw2(); list->Add(hKCh_dPhi1min);
	TProfile *pKAvgPt_NT1min = new TProfile("pKAvgPt_NT1min", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKAvgPt_NT1min_ptcut = new TProfile("pKAvgPt_NT1min_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pKAvgPt_NT1min->Sumw2(); list->Add(pKAvgPt_NT1min);
	pKAvgPt_NT1min_ptcut->Sumw2(); list->Add(pKAvgPt_NT1min_ptcut);
	TH3D *hKCh_dPhi1max = new TH3D("hKCh_dPhi1max", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hKCh_dPhi1max->Sumw2(); list->Add(hKCh_dPhi1max);
	TProfile *pKAvgPt_NT1max = new TProfile("pKAvgPt_NT1max", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKAvgPt_NT1max_ptcut = new TProfile("pKAvgPt_NT1max_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pKAvgPt_NT1max->Sumw2(); list->Add(pKAvgPt_NT1max);
	pKAvgPt_NT1max_ptcut->Sumw2(); list->Add(pKAvgPt_NT1max_ptcut);


	TH3D *hPiCh_dPhi0 = new TH3D("hPiCh_dPhi0", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hPiCh_dPhi0->Sumw2(); list->Add(hPiCh_dPhi0);
	TH3D *hPiCh_dPhi1 = new TH3D("hPiCh_dPhi1", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hPiCh_dPhi1->Sumw2(); list->Add(hPiCh_dPhi1);
	TH3D *hPiCh_dPhi2 = new TH3D("hPiCh_dPhi2", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hPiCh_dPhi2->Sumw2(); list->Add(hPiCh_dPhi2);

	TProfile *pPiAvgPt_NT0 = new TProfile("pPiAvgPt_NT0", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPiAvgPt_NT1 = new TProfile("pPiAvgPt_NT1", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPiAvgPt_NT2 = new TProfile("pPiAvgPt_NT2", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPiAvgPt_NT0_ptcut = new TProfile("pPiAvgPt_NT0_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPiAvgPt_NT1_ptcut = new TProfile("pPiAvgPt_NT1_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPiAvgPt_NT2_ptcut = new TProfile("pPiAvgPt_NT2_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pPiAvgPt_NT0->Sumw2(); list->Add(pPiAvgPt_NT0);
	pPiAvgPt_NT1->Sumw2(); list->Add(pPiAvgPt_NT1);
	pPiAvgPt_NT2->Sumw2(); list->Add(pPiAvgPt_NT2);
	pPiAvgPt_NT0_ptcut->Sumw2(); list->Add(pPiAvgPt_NT0_ptcut);
	pPiAvgPt_NT1_ptcut->Sumw2(); list->Add(pPiAvgPt_NT1_ptcut);
	pPiAvgPt_NT2_ptcut->Sumw2(); list->Add(pPiAvgPt_NT2_ptcut);

	TH3D *hPiCh_dPhi1min = new TH3D("hPiCh_dPhi1min", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hPiCh_dPhi1min->Sumw2(); list->Add(hPiCh_dPhi1min);
	TProfile *pPiAvgPt_NT1min = new TProfile("pPiAvgPt_NT1min", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPiAvgPt_NT1min_ptcut = new TProfile("pPiAvgPt_NT1min_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pPiAvgPt_NT1min->Sumw2(); list->Add(pPiAvgPt_NT1min);
	pPiAvgPt_NT1min_ptcut->Sumw2(); list->Add(pPiAvgPt_NT1min_ptcut);
	TH3D *hPiCh_dPhi1max = new TH3D("hPiCh_dPhi1max", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hPiCh_dPhi1max->Sumw2(); list->Add(hPiCh_dPhi1max);
	TProfile *pPiAvgPt_NT1max = new TProfile("pPiAvgPt_NT1max", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPiAvgPt_NT1max_ptcut = new TProfile("pPiAvgPt_NT1max_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pPiAvgPt_NT1max->Sumw2(); list->Add(pPiAvgPt_NT1max);
	pPiAvgPt_NT1max_ptcut->Sumw2(); list->Add(pPiAvgPt_NT1max_ptcut);


	TH3D *hKshort_dPhi0 = new TH3D("hKshort_dPhi0", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hKshort_dPhi0->Sumw2(); list->Add(hKshort_dPhi0);
	TH3D *hKshort_dPhi1 = new TH3D("hKshort_dPhi1", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hKshort_dPhi1->Sumw2(); list->Add(hKshort_dPhi1);
	TH3D *hKshort_dPhi2 = new TH3D("hKshort_dPhi2", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hKshort_dPhi2->Sumw2(); list->Add(hKshort_dPhi2);

	TProfile *pKSAvgPt_NT0 = new TProfile("pKSAvgPt_NT0", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKSAvgPt_NT1 = new TProfile("pKSAvgPt_NT1", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKSAvgPt_NT2 = new TProfile("pKSAvgPt_NT2", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKSAvgPt_NT0_ptcut = new TProfile("pKSAvgPt_NT0_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKSAvgPt_NT1_ptcut = new TProfile("pKSAvgPt_NT1_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pKSAvgPt_NT2_ptcut = new TProfile("pKSAvgPt_NT2_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pKSAvgPt_NT0->Sumw2(); list->Add(pKSAvgPt_NT0);
	pKSAvgPt_NT1->Sumw2(); list->Add(pKSAvgPt_NT1);
	pKSAvgPt_NT2->Sumw2(); list->Add(pKSAvgPt_NT2);
	pKSAvgPt_NT0_ptcut->Sumw2(); list->Add(pKSAvgPt_NT0_ptcut);
	pKSAvgPt_NT1_ptcut->Sumw2(); list->Add(pKSAvgPt_NT1_ptcut);
	pKSAvgPt_NT2_ptcut->Sumw2(); list->Add(pKSAvgPt_NT2_ptcut);

	TH3D *hLambda_dPhi0 = new TH3D("hLambda_dPhi0", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hLambda_dPhi0->Sumw2(); list->Add(hLambda_dPhi0);
	TH3D *hLambda_dPhi1 = new TH3D("hLambda_dPhi1", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hLambda_dPhi1->Sumw2(); list->Add(hLambda_dPhi1);
	TH3D *hLambda_dPhi2 = new TH3D("hLambda_dPhi2", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hLambda_dPhi2->Sumw2(); list->Add(hLambda_dPhi2);

	TProfile *pLAvgPt_NT0 = new TProfile("pLAvgPt_NT0", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pLAvgPt_NT1 = new TProfile("pLAvgPt_NT1", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pLAvgPt_NT2 = new TProfile("pLAvgPt_NT2", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pLAvgPt_NT0_ptcut = new TProfile("pLAvgPt_NT0_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pLAvgPt_NT1_ptcut = new TProfile("pLAvgPt_NT1_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pLAvgPt_NT2_ptcut = new TProfile("pLAvgPt_NT2_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pLAvgPt_NT0->Sumw2(); list->Add(pLAvgPt_NT0);
	pLAvgPt_NT1->Sumw2(); list->Add(pLAvgPt_NT1);
	pLAvgPt_NT2->Sumw2(); list->Add(pLAvgPt_NT2);
	pLAvgPt_NT0_ptcut->Sumw2(); list->Add(pLAvgPt_NT0_ptcut);
	pLAvgPt_NT1_ptcut->Sumw2(); list->Add(pLAvgPt_NT1_ptcut);
	pLAvgPt_NT2_ptcut->Sumw2(); list->Add(pLAvgPt_NT2_ptcut);

	TH3D *hProton_dPhi0 = new TH3D("hProton_dPhi0", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hProton_dPhi0->Sumw2(); list->Add(hProton_dPhi0);
	TH3D *hProton_dPhi1 = new TH3D("hProton_dPhi1", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hProton_dPhi1->Sumw2(); list->Add(hProton_dPhi1);
	TH3D *hProton_dPhi2 = new TH3D("hProton_dPhi2", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hProton_dPhi2->Sumw2(); list->Add(hProton_dPhi2);

	TProfile *pPAvgPt_NT0 = new TProfile("pPAvgPt_NT0", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPAvgPt_NT1 = new TProfile("pPAvgPt_NT1", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPAvgPt_NT2 = new TProfile("pPAvgPt_NT2", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPAvgPt_NT0_ptcut = new TProfile("pPAvgPt_NT0_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPAvgPt_NT1_ptcut = new TProfile("pPAvgPt_NT1_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPAvgPt_NT2_ptcut = new TProfile("pPAvgPt_NT2_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pPAvgPt_NT0->Sumw2(); list->Add(pPAvgPt_NT0);
	pPAvgPt_NT1->Sumw2(); list->Add(pPAvgPt_NT1);
	pPAvgPt_NT2->Sumw2(); list->Add(pPAvgPt_NT2);
	pPAvgPt_NT0_ptcut->Sumw2(); list->Add(pPAvgPt_NT0_ptcut);
	pPAvgPt_NT1_ptcut->Sumw2(); list->Add(pPAvgPt_NT1_ptcut);
	pPAvgPt_NT2_ptcut->Sumw2(); list->Add(pPAvgPt_NT2_ptcut);

	TH3D *hProton_dPhi1min = new TH3D("hProton_dPhi1min", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hProton_dPhi1min->Sumw2(); list->Add(hProton_dPhi1min);
	TProfile *pPAvgPt_NT1min = new TProfile("pPAvgPt_NT1min", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPAvgPt_NT1min_ptcut = new TProfile("pPAvgPt_NT1min_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pPAvgPt_NT1min->Sumw2(); list->Add(pPAvgPt_NT1min);
	pPAvgPt_NT1min_ptcut->Sumw2(); list->Add(pPAvgPt_NT1min_ptcut);
	TH3D *hProton_dPhi1max = new TH3D("hProton_dPhi1max", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hProton_dPhi1max->Sumw2(); list->Add(hProton_dPhi1max);
	TProfile *pPAvgPt_NT1max = new TProfile("pPAvgPt_NT1max", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pPAvgPt_NT1max_ptcut = new TProfile("pPAvgPt_NT1max_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pPAvgPt_NT1max->Sumw2(); list->Add(pPAvgPt_NT1max);
	pPAvgPt_NT1max_ptcut->Sumw2(); list->Add(pPAvgPt_NT1max_ptcut);


	TH3D *hD0_dPhi0 = new TH3D("hD0_dPhi0", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hD0_dPhi0->Sumw2(); list->Add(hD0_dPhi0);
	TH3D *hD0_dPhi1 = new TH3D("hD0_dPhi1", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hD0_dPhi1->Sumw2(); list->Add(hD0_dPhi1);
	TH3D *hD0_dPhi2 = new TH3D("hD0_dPhi2", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hD0_dPhi2->Sumw2(); list->Add(hD0_dPhi2);

	TProfile *pD0AvgPt_NT0 = new TProfile("pD0AvgPt_NT0", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pD0AvgPt_NT1 = new TProfile("pD0AvgPt_NT1", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pD0AvgPt_NT2 = new TProfile("pD0AvgPt_NT2", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pD0AvgPt_NT0_ptcut = new TProfile("pD0AvgPt_NT0_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pD0AvgPt_NT1_ptcut = new TProfile("pD0AvgPt_NT1_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pD0AvgPt_NT2_ptcut = new TProfile("pD0AvgPt_NT2_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pD0AvgPt_NT0->Sumw2(); list->Add(pD0AvgPt_NT0);
	pD0AvgPt_NT1->Sumw2(); list->Add(pD0AvgPt_NT1);
	pD0AvgPt_NT2->Sumw2(); list->Add(pD0AvgPt_NT2);
	pD0AvgPt_NT0_ptcut->Sumw2(); list->Add(pD0AvgPt_NT0_ptcut);
	pD0AvgPt_NT1_ptcut->Sumw2(); list->Add(pD0AvgPt_NT1_ptcut);
	pD0AvgPt_NT2_ptcut->Sumw2(); list->Add(pD0AvgPt_NT2_ptcut);

	TH3D *hDs_dPhi0 = new TH3D("hDs_dPhi0", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hDs_dPhi0->Sumw2(); list->Add(hDs_dPhi0);
	TH3D *hDs_dPhi1 = new TH3D("hDs_dPhi1", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hDs_dPhi1->Sumw2(); list->Add(hDs_dPhi1);
	TH3D *hDs_dPhi2 = new TH3D("hDs_dPhi2", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hDs_dPhi2->Sumw2(); list->Add(hDs_dPhi2);

	TProfile *pDsAvgPt_NT0 = new TProfile("pDsAvgPt_NT0", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pDsAvgPt_NT1 = new TProfile("pDsAvgPt_NT1", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pDsAvgPt_NT2 = new TProfile("pDsAvgPt_NT2", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pDsAvgPt_NT0_ptcut = new TProfile("pDsAvgPt_NT0_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pDsAvgPt_NT1_ptcut = new TProfile("pDsAvgPt_NT1_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pDsAvgPt_NT2_ptcut = new TProfile("pDsAvgPt_NT2_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pDsAvgPt_NT0->Sumw2(); list->Add(pDsAvgPt_NT0);
	pDsAvgPt_NT1->Sumw2(); list->Add(pDsAvgPt_NT1);
	pDsAvgPt_NT2->Sumw2(); list->Add(pDsAvgPt_NT2);
	pDsAvgPt_NT0_ptcut->Sumw2(); list->Add(pDsAvgPt_NT0_ptcut);
	pDsAvgPt_NT1_ptcut->Sumw2(); list->Add(pDsAvgPt_NT1_ptcut);
	pDsAvgPt_NT2_ptcut->Sumw2(); list->Add(pDsAvgPt_NT2_ptcut);

	TH3D *hLambc_dPhi0 = new TH3D("hLambc_dPhi0", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hLambc_dPhi0->Sumw2(); list->Add(hLambc_dPhi0);
	TH3D *hLambc_dPhi1 = new TH3D("hLambc_dPhi1", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hLambc_dPhi1->Sumw2(); list->Add(hLambc_dPhi1);
	TH3D *hLambc_dPhi2 = new TH3D("hLambc_dPhi2", ";p_{T} [GeV];N_{ch}^{fwd}; S", nptbins, ptbin, nNchbins, Nchbin, nNTbins, NTbin);
	hLambc_dPhi2->Sumw2(); list->Add(hLambc_dPhi2);

	TProfile *pLcAvgPt_NT0 = new TProfile("pLcAvgPt_NT0", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pLcAvgPt_NT1 = new TProfile("pLcAvgPt_NT1", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pLcAvgPt_NT2 = new TProfile("pLcAvgPt_NT2", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pLcAvgPt_NT0_ptcut = new TProfile("pLcAvgPt_NT0_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pLcAvgPt_NT1_ptcut = new TProfile("pLcAvgPt_NT1_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	TProfile *pLcAvgPt_NT2_ptcut = new TProfile("pLcAvgPt_NT2_ptcut", ";N_{T} ;<p_{T}> [GeV]", nNTbins, NTbin);
	pLcAvgPt_NT0->Sumw2(); list->Add(pLcAvgPt_NT0);
	pLcAvgPt_NT1->Sumw2(); list->Add(pLcAvgPt_NT1);
	pLcAvgPt_NT2->Sumw2(); list->Add(pLcAvgPt_NT2);
	pLcAvgPt_NT0_ptcut->Sumw2(); list->Add(pLcAvgPt_NT0_ptcut);
	pLcAvgPt_NT1_ptcut->Sumw2(); list->Add(pLcAvgPt_NT1_ptcut);
	pLcAvgPt_NT2_ptcut->Sumw2(); list->Add(pLcAvgPt_NT2_ptcut);

//=============================================================================

	//cuts
	//alice cut |eta|<0.5, |y|<0.5
  const double dcEtaCut  = 1.0;
  const double dfEtaMin  = 2.;
  const double dfEtaMax  = 5.;
  const double dPtCut    = 0.25;

	const double dPtTrMax = 3;
	const double dPtTrMin = 1;
	const double dPtAsMax = 3;
	const double dPtAsMin = 1;

	const double dHighMult = 80;
	const double dLowMult = 20;
	//temperary containers
	TLorentzVector vec_tmp;
	TLorentzVector vec_LeadJet;
	vector<TLorentzVector> sph_container;
	vector<ParticleInfo> decayList;


  ParticleInfo pInfo;

	//jet related
	vector<PseudoJet> jetParticles;
	double R=0.4;
	JetDefinition jet_def(antikt_algorithm, R);
	cout << "Clustering with " << jet_def.description() << endl;
   
  //--------------------------------
  // loop events
  //--------------------------------

long nEvents = (int)chain->GetEntries();
cout << "Total events : " << nEvents << endl;
AMPT* mctree = new AMPT(chain);
  double nEvents_Inel=0;
const long maxEvents = 100e6; // 设置最大事件数为一百万
for (int iEvent = 0; iEvent < min(nEvents, maxEvents); ++iEvent) {
	  if (iEvent%10000==0) cout<<"Proccessing event # "<<iEvent<<endl;
	  mctree->GetEntry(iEvent);

	  int nTracks = mctree->Event_multiplicity;
    double dnMPI = mctree->Event_ncoll;

		//event wise initializations
    double dFwdCh = 0.;
    double dMidCh = 0.;
		double dNch05 = 0.;
		double dNch08 = 0.;
		double dNtrkV0 = 0.;
		double dEtrkV0 = 0.;
		double dNtrkV0A = 0.;
		double dNtrkV0C = 0.;
    double dNchCMS = 0;
    int nEvents_1 = 0; 
		sph_container.clear();
    jetParticles.clear();
    decayList.clear();
    
    double maxPt = 0.0;
		//bool keepEvent = false;
		//int maxPtIndex = -1; 
		double maxPhi = 0.0;
	  //track loop
	  for( int iTrack = 0; iTrack < nTracks; ++iTrack ){
		  int id = mctree->id[iTrack];
		  float px = mctree->px[iTrack];
		  float py = mctree->py[iTrack];
		  float pz = mctree->pz[iTrack];
		  float mass = mctree->m[iTrack];
			float E = sqrt(mass*mass + px*px + py*py + pz*pz);
		  vec_tmp.SetPxPyPzE(px, py, pz, E);
			double Pt = sqrt(px * px + py * py);

      if( !IsFinal(abs(id)) ){
        //cout<<id<<" "<<vec_tmp.M()<<endl;
        eventMain.clear();
        //eventMain.append(id,11,0,0,vec_tmp.Px(), vec_tmp.Py(),vec_tmp.Pz(),vec_tmp.E(), vec_tmp.M());
        eventMain.append(id,11,0,0,vec_tmp.Px(), vec_tmp.Py(),vec_tmp.Pz(),vec_tmp.E(), mass);
        pythiaMain.moreDecays();
        //eventMain.list();
        //loop from the decay product, initial index is 1 instead of 0
        for(int iTrk = 1; iTrk<eventMain.size(); ++iTrk)
        {
            if(eventMain[iTrk].isFinal())
            {
                int dauId = eventMain[iTrk].id();
                //cout<<iTrk<<" "<<dauId<<endl;
                //assign the 4-mom of this decay daughter to vDaughter vector
                TLorentzVector vDaughter(eventMain[iTrk].px(),eventMain[iTrk].py(),eventMain[iTrk].pz(),eventMain[iTrk].e());
                pInfo.Id=dauId;
                pInfo.vec=vDaughter;
                pInfo.MoId=id;
                decayList.push_back(pInfo);
            }
        }
        continue;
      }

			double dPt = sqrt(px*px+py*py);
      double dRap = vec_tmp.Rapidity();
		  //if ( dPt<1E-3 || vec_tmp.P()<1e-3 ) continue;

			double dEta = vec_tmp.Eta();
			
			//V0 charge amplitude
			if( ((dEta>2.8&&dEta<5.1) || (dEta>-3.7&&dEta<-1.7)) && GetCharge(id)!=0 ){
				dNtrkV0++;
        dEtrkV0+=E;
				if(dEta>2.8&&dEta<5.1) dNtrkV0A++;
				if(dEta>-3.7&&dEta<-1.7) dNtrkV0C++;
			}

      double dEtaAbs = TMath::Abs(dEta);

			if( GetCharge(id)!=0 ){
				if ((dEtaAbs>dfEtaMin) && (dEtaAbs<=dfEtaMax)) {
					dFwdCh += 1.;
				}

				if( dEtaAbs<dcEtaCut ) {
					if(dPt>0) dMidCh++; //used for INEL>0 cut
					if( dPt>dPtCut ) {
            sph_container.push_back(vec_tmp);
            jetParticles.push_back(PseudoJet(px, py, pz, E));
          }
					if( dEtaAbs<0.5 ) dNch05++;
					if( dEtaAbs<0.8 && dPt>0.2&& dPt<3 ) dNch08++; //used for ALICE flow ana event classifier
				}
			}//charge fi

			if(fabs(dEta)<0.8 && GetCharge(id) != 0 && Pt > maxPt) {
        maxPt = Pt;
        //maxPtIndex = iTrack;
        maxPhi = atan2(py, px);
      }
		}//track loop

    histPt->Fill(maxPt);
		
    /*
		if (maxPt >= 5.0 && maxPt < 40.0) {
    nEvents_1++;
		} else {
		    continue;
		}
    */

    for(int idecay=0; idecay<decayList.size(); idecay++){
      int id = decayList.at(idecay).Id;
      int MoId = decayList.at(idecay).MoId;
      vec_tmp = decayList.at(idecay).vec;
			double dEta = vec_tmp.Eta();
      double dPt = vec_tmp.Pt();
      double E = vec_tmp.E();
			//V0 charge amplitude
      if( GetCharge(id)!=0 ){
        if( ((dEta>2.8&&dEta<5.1) || (dEta>-3.7&&dEta<-1.7)) ){
          dNtrkV0++;
          dEtrkV0+=E;
          if(dEta>2.8&&dEta<5.1) dNtrkV0A++;
          if(dEta>-3.7&&dEta<-1.7) dNtrkV0C++;
        }
				if ((fabs(dEta)>dfEtaMin) && (fabs(dEta)<=dfEtaMax)) {
					dFwdCh += 1.;
				}
        if( fabs(dEta)<dcEtaCut ){
          dMidCh++;
					if( dPt>dPtCut ) {
            sph_container.push_back(vec_tmp);
            jetParticles.push_back(PseudoJet(vec_tmp.Px(), vec_tmp.Py(), vec_tmp.Pz(), vec_tmp.E()));
          }
          if( fabs(dEta)<0.5 ) dNch05++;
          if( fabs(dEta)<0.8 && dPt>0.2 && dPt<3 ) dNch08++; 
        }
      }

      if(fabs(dEta)<0.8 && GetCharge(id) != 0 && dPt > maxPt) {
        maxPt = dPt;
        //maxPtIndex = iTrack;
        maxPhi = atan2(vec_tmp.Py(), vec_tmp.Py());
      }

    }//decay loop

		if (maxPt < 5.0 || maxPt >= 40.0) continue;
    nEvents_1++;


		if(dMidCh>0||(dNtrkV0>0)){ //INEL cut

      nEvents_Inel++;

      //INEL events analysis loop
      for (int iTrack=0; iTrack<nTracks; iTrack++) {
        int id = mctree->id[iTrack];

        if( !IsFinal(abs(id)) ) continue;

        float px = mctree->px[iTrack];
        float py = mctree->py[iTrack];
        float pz = mctree->pz[iTrack];
        float mass = mctree->m[iTrack];
        float Energy = sqrt(mass*mass+px*px+py*py+pz*pz);
        vec_tmp.SetPxPyPzE(px,py,pz,Energy);
        double dPt = sqrt(px*px + py*py);
        double dEta = vec_tmp.Eta();
        double dRap = vec_tmp.Rapidity();
        double dPhi = vec_tmp.Phi();
        if( GetCharge(id)!=0 ){
  //			if(abs(id)==211||abs(id)==321||abs(id)==2212){
          hChEta->Fill(dEta);
          if(abs(id)==211) h_pion_eta->Fill(dRap);
          if(abs(id)==321) h_kaon_eta->Fill(dRap);
          if(abs(id)==2212) h_proton_eta->Fill(dRap);

        }
      }
    }
		//if(dMidCh<1&&dNtrkV0A<1&&dNtrkV0C<1) continue; //INEL>0 cut
		if(dNtrkV0A<1||dNtrkV0C<1) continue; //MB cut

		nEvents_Inel++;

		hNtrkV0->Fill(dNtrkV0);

		histEventCount->Fill(nEvents_1);

    hnMPIVsNtrkV0->Fill(dNtrkV0, dnMPI);

		hMultInel->Fill(dMidCh); 

		//hFwdVsMid->Fill(dMidCh, dNtrkV0);

		hFwdVsMid->Fill(dNch05, dNtrkV0);

		double dSph = GetSphericity(sph_container);

//=============================================================================
		double dNchTransI = 0.; // transverse I区域的带电粒子数
    double dNchTransII = 0.; // transverse II区域的带电粒子数
		double dNchTrans = 0.; // transverse I区域的带电粒子数

		for (int iTrack=0; iTrack<nTracks; iTrack++){
      int id = mctree->id[iTrack];
      float px = mctree->px[iTrack];
      float py = mctree->py[iTrack];
      float pz = mctree->pz[iTrack];
			float m = mctree->m[iTrack];
			float E = sqrt(px*px+py*py+pz*pz+m*m);
			vec_tmp.SetPxPyPzE(px,py,pz,E);
      double dPt = sqrt(px*px + py*py);
			double dEta = vec_tmp.Eta();
      double dRap = vec_tmp.Rapidity();

      //take both particle and anti-particle
      int aId=abs(id);

      //double dEtaAbs = TMath::Abs(dEta); 
			if (fabs(dEta)>=0.8) continue;
      
      if (GetCharge(id) != 0 && dPt > 0.15 && dPt < 5) {
        double Phi = vec_tmp.Phi();
        double dPhi = Phi - maxPhi; // maxPhi为触发粒子的phi角
        dPhi = TVector2::Phi_mpi_pi(dPhi); // 将dPhi限制在-pi到pi之间

        // transverse I区域：pi/3 < deltaPhi < 2pi/3
        if (dPhi > TMath::Pi() / 3 && dPhi < 2 * TMath::Pi() / 3) {
            dNchTransI++;
        }
        // transverse II区域：-2pi/3 < deltaPhi < -pi/3
        else if (dPhi < -TMath::Pi() / 3 && dPhi > -2 * TMath::Pi() / 3) {
            dNchTransII++;
        }
      }
    }//track loop end

    for(int idecay=0; idecay<decayList.size(); idecay++){
      int id = decayList.at(idecay).Id;
      int MoId = decayList.at(idecay).MoId;
      vec_tmp = decayList.at(idecay).vec;
      double dEta = vec_tmp.Eta();
      double dRap = vec_tmp.Rapidity();
      double dPt = vec_tmp.Pt();
      double E = vec_tmp.E();

      if (fabs(dEta)>=0.8) continue;

      if (GetCharge(id) != 0 && dPt > 0.15 && dPt < 5) {
          double Phi = vec_tmp.Phi();
          double dPhi = Phi - maxPhi; // maxPhi为触发粒子的phi角
          dPhi = TVector2::Phi_mpi_pi(dPhi); // 将dPhi限制在-pi到pi之间

          // transverse I区域：pi/3 < deltaPhi < 2pi/3
          if (dPhi > TMath::Pi() / 3 && dPhi < 2 * TMath::Pi() / 3) {
              dNchTransI++;
          }
          // transverse II区域：-2pi/3 < deltaPhi < -pi/3
          else if (dPhi < -TMath::Pi() / 3 && dPhi > -2 * TMath::Pi() / 3) {
              dNchTransII++;
          }
      }
    }

    // 确定trans-max和trans-min区域
    double dNchTransMax = (dNchTransI >= dNchTransII) ? dNchTransI : dNchTransII;
    double dNchTransMin = (dNchTransI <= dNchTransII) ? dNchTransI : dNchTransII;

    //cout<<dNchTransMax<<" "<<dNchTransMin<<endl;

    dNchTrans=dNchTransI+dNchTransII;
    //dNchTrans=dNchTransMin;
    //nTHist->Fill(dNchTransMin);
    nTHist->Fill(dNchTrans);
    
    /**************************************************************************
     * 
     *  ====================== NEW JET FINDING LOGIC START ======================
     *  
     *  步骤 1: 收集事件中所有 |eta|<0.8, pT>0.25 的粒子.
     *  步骤 2: 对所有收集到的粒子运行一次 anti-kT 算法, 找到事件中的所有喷注.
     *  步骤 3: 遍历找到的喷注, 根据它们相对于领头粒子(maxPhi)的方位角,
     *          将它们分类到 Toward 和 Transverse 区域.
     *  步骤 4: 使用正确分类的喷注集合来填充直方图和 TTree.
     * 
     **************************************************************************/

    // 步骤 1: 收集整个事件中符合条件的粒子
    vector<PseudoJet> allParticlesForJets;
    for (int iTrk = 0; iTrk < nTracks; ++iTrk) {
        if (!IsFinal(abs(mctree->id[iTrk]))) continue; // 只使用末态粒子

        TLorentzVector v;
        float px = mctree->px[iTrk];
        float py = mctree->py[iTrk];
        float pz = mctree->pz[iTrk];
        float m  = mctree->m[iTrk];
        v.SetPxPyPzE(px, py, pz, sqrt(m*m + px*px + py*py + pz*pz));

        // 应用统一的粒子运动学接收度
        if (fabs(v.Eta()) < 0.8 && v.Pt() > 0.25) {
            allParticlesForJets.emplace_back(v.Px(), v.Py(), v.Pz(), v.E());
            // 保存原始粒子在 mctree 中的索引, 以便稍后获取PID等信息
            allParticlesForJets.back().set_user_index(iTrk);
        }
    }

    // 步骤 2: 对所有粒子运行一次喷注寻找
    double R_jet = 0.4;
    JetDefinition jet_def_global(antikt_algorithm, R_jet);
    ClusterSequence cs_global(allParticlesForJets, jet_def_global);
    // 寻找 pT > 0 GeV/c 的所有喷注 (您也可以设置一个阈值, 如 2.0 GeV/c)
    vector<PseudoJet> all_jets_found = sorted_by_pt(cs_global.inclusive_jets(0.0)); 

    // 步骤 3: 准备容器并对喷注进行分类
    vector<PseudoJet> toward_jets;
    vector<PseudoJet> transverse_jets;

    for (const PseudoJet& jet : all_jets_found) {
        double dPhi_jet = TVector2::Phi_mpi_pi(jet.phi_std() - maxPhi);
        double dPhi_jet_abs = fabs(dPhi_jet);

        if (dPhi_jet_abs < TMath::Pi() / 3.0) {
            toward_jets.push_back(jet);
        } else if (dPhi_jet_abs >= TMath::Pi() / 3.0 && dPhi_jet_abs < 2.0 * TMath::Pi() / 3.0) {
            transverse_jets.push_back(jet);
        }
        // 其他区域 (Away) 的喷注在此被忽略
    }

    // 步骤 4: 使用正确分类的喷注填充直方图和 TTree

    // --- 填充 Transverse 区域的喷注信息 ---
    hTRnJet->Fill(transverse_jets.size());
    for (const PseudoJet &j : transverse_jets) {
        hTRjetPt->Fill(j.pt());
        hTRjetEta->Fill(j.eta());
        hTRjetPhi->Fill(j.phi_std());
        for (const PseudoJet &cons : j.constituents()) {
            int idx = cons.user_index();
            int pid = mctree->id[idx];
            hJetPid->Fill(pid);
        }
    }
    
    // --- 计算 Toward 和 Transverse 区域的带电粒子数 (用于 TTree) ---
    int nchTow = 0;
    int nchTR = 0;
    // (这段逻辑和您原来的 dNchTrans 计算类似, 但现在可以精确地定义区域)
    double tow_phi_min_wrapped = wrap(maxPhi - M_PI / 3.0);
    double tow_phi_max_wrapped = wrap(maxPhi + M_PI / 3.0);
    Selector tow_phi_selector = (tow_phi_min_wrapped < tow_phi_max_wrapped) ?
                                SelectorPhiRange(tow_phi_min_wrapped, tow_phi_max_wrapped) :
                                (SelectorPhiRange(tow_phi_min_wrapped, M_PI) || SelectorPhiRange(-M_PI, tow_phi_max_wrapped));
    Selector tow_selector_final = SelectorAbsEtaMax(0.8) && tow_phi_selector;
    
    // Transverse 区域的选择器定义 (来自您原始代码)
    double phi1_min_wrapped = wrap(maxPhi + M_PI / 3.0);
    double phi1_max_wrapped = wrap(maxPhi + 2.0 * M_PI / 3.0);
    Selector tr_sel1 = (phi1_min_wrapped < phi1_max_wrapped) ? 
                       SelectorPhiRange(phi1_min_wrapped, phi1_max_wrapped) : 
                       (SelectorPhiRange(phi1_min_wrapped, M_PI) || SelectorPhiRange(-M_PI, phi1_max_wrapped));
    double phi2_min_wrapped = wrap(maxPhi - 2.0 * M_PI / 3.0);
    double phi2_max_wrapped = wrap(maxPhi - M_PI / 3.0);
    Selector tr_sel2 = (phi2_min_wrapped < phi2_max_wrapped) ? 
                       SelectorPhiRange(phi2_min_wrapped, phi2_max_wrapped) :
                       (SelectorPhiRange(phi2_min_wrapped, M_PI) || SelectorPhiRange(-M_PI, phi2_max_wrapped));
    Selector tr_selector_final = SelectorAbsEtaMax(0.8) && (tr_sel1 || tr_sel2);

    for (int iTrk = 0; iTrk < nTracks; ++iTrk) {
        if (!IsFinal(abs(mctree->id[iTrk]))) continue;
        if (GetCharge(mctree->id[iTrk]) == 0) continue; // 只统计带电粒子
        
        TLorentzVector v;
        v.SetPxPyPzE(mctree->px[iTrk], mctree->py[iTrk], mctree->pz[iTrk], sqrt(pow(mctree->m[iTrk],2) + pow(mctree->px[iTrk],2) + pow(mctree->py[iTrk],2) + pow(mctree->pz[iTrk],2)));
        
        if (fabs(v.Eta()) > 0.8 || v.Pt() < 0.25) continue;
        
        PseudoJet pj(v.Px(), v.Py(), v.Pz(), v.E());
        if (tow_selector_final(pj)) nchTow++;
        if (tr_selector_final(pj)) nchTR++;
    }

    tNchTow = nchTow;
    tNchTR = nchTR;

    // --- 填充 TTree ---
    
    // Transverse jets
    tNjetTR = transverse_jets.size();
    tJetPt.clear(); tJetEta.clear(); tJetPhi.clear();
    tJetLeadPID.clear(); tJetConstPID_flattened.clear(); tJetConstN.clear();
    for (const PseudoJet &j : transverse_jets) {
        tJetPt.push_back(j.pt());
        tJetEta.push_back(j.eta());
        tJetPhi.push_back(j.phi_std());

        double maxPtConst = 0;
        int leadPID = 0;
        vector<int> oneJetPIDs;
        for (const PseudoJet &cons : j.constituents()) {
            int idx = cons.user_index();
            int pid = mctree->id[idx];
            oneJetPIDs.push_back(pid);
            if (cons.pt() > maxPtConst) {
                maxPtConst = cons.pt();
                leadPID = pid;
            }
        }
        tJetLeadPID.push_back(leadPID);
        tJetConstN.push_back(oneJetPIDs.size());
        for (int pid : oneJetPIDs) {
            tJetConstPID_flattened.push_back(pid);
        }
    }

    // Toward jets
    tNjetTow = toward_jets.size();
    tJetTowPt.clear(); tJetTowEta.clear(); tJetTowPhi.clear();
    tJetTowLeadPID.clear(); tJetTowConstPID_flattened.clear(); tJetTowConstN.clear();
    for (const PseudoJet &j : toward_jets) {
        tJetTowPt.push_back(j.pt());
        tJetTowEta.push_back(j.eta());
        tJetTowPhi.push_back(j.phi_std());
        
        double maxPtConst = 0;
        int leadPID = 0;
        vector<int> oneJetPIDs;
        for (const PseudoJet &cons : j.constituents()) {
            int idx = cons.user_index();
            int pid = mctree->id[idx];
            oneJetPIDs.push_back(pid);
            if (cons.pt() > maxPtConst) {
                maxPtConst = cons.pt();
                leadPID = pid;
            }
        }
        tJetTowLeadPID.push_back(leadPID);
        tJetTowConstN.push_back(oneJetPIDs.size());
        for (int pid : oneJetPIDs) {
            tJetTowConstPID_flattened.push_back(pid);
        }
    }

    tTRjet->Fill(); // 每个事件填充一次 TTree
    
    // ======================= NEW JET FINDING LOGIC END =======================

		for (int iTrack=0; iTrack<nTracks; iTrack++){
      int id = mctree->id[iTrack];
      float px = mctree->px[iTrack];
      float py = mctree->py[iTrack];
      float pz = mctree->pz[iTrack];
			float m = mctree->m[iTrack];
			float E = sqrt(px*px+py*py+pz*pz+m*m);
			vec_tmp.SetPxPyPzE(px,py,pz,E);
      double dPt = sqrt(px*px + py*py);
			double dEta = vec_tmp.Eta();
      double dRap = vec_tmp.Rapidity();

      //take both particle and anti-particle
      int aId=abs(id);

      //double dEtaAbs = TMath::Abs(dEta); 
			if (fabs(dEta)>=0.8) continue;
	    double Phi = vec_tmp.Phi();
	    double dPhi = Phi - maxPhi;
	    double dPhiAbs = fabs(TVector2::Phi_mpi_pi(dPhi));

      bool isTransMin = false;

	    int dPhi_segment = 0;
		  if (dPhiAbs >= 0 && dPhiAbs < TMath::Pi() / 3) {
	      dPhi_segment = 0; //to
	      //cout << "Particle " << iTrack << " in 'to' region - phi: " << Phi* 180.0 / M_PI << ", dPhi: " << dPhi* 180.0 / M_PI << ", maxPhi: " << maxPhi* 180.0 / M_PI << endl;
      }else if (dPhiAbs >= 2 * TMath::Pi() / 3 && dPhiAbs < TMath::Pi()) {
	      dPhi_segment = 2; //away
	      //cout << "Particle " << iTrack << " in 'away' region - phi: " << Phi* 180.0 / M_PI << ", dPhi: " << dPhi* 180.0 / M_PI << ", maxPhi: " << maxPhi* 180.0 / M_PI << endl;
	    } else {
	      dPhi_segment = 1; //tr
	      //cout << "Particle " << iTrack << " in 'tr' region - phi: " << Phi* 180.0 / M_PI << ", dPhi: " << dPhi* 180.0 / M_PI << ", maxPhi: " << maxPhi* 180.0 / M_PI << endl;
        
        //if TransI==TransII, define dPhi>0 as TransMin
        if (dNchTransMin==dNchTransMax){
          if (dPhi>0) isTransMin=true;
          //isTransMin=false;
        }
        //if TransI!=TransII, simply choose the min region
        else{
          if (dNchTransII == dNchTransMin && dPhi < 0) isTransMin=true;
          if (dNchTransI == dNchTransMin && dPhi > 0) isTransMin=true;
        }
	    }

			if (GetCharge(id) != 0) {
				//if (maxPtIndex != -1) {
					if(fabs(dRap)<0.5){
	        if(abs(id)==321) {
           	if (dPhi_segment == 0) {
              if(vec_tmp.Pt()<5) pKAvgPt_NT0_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pKAvgPt_NT0->Fill(dNchTrans, vec_tmp.Pt());
            	hKCh_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
              if(vec_tmp.Pt()<5) pKAvgPt_NT1_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pKAvgPt_NT1->Fill(dNchTrans, vec_tmp.Pt());
              hKCh_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              if(isTransMin){
                if(vec_tmp.Pt()<5) pKAvgPt_NT1min_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pKAvgPt_NT1min->Fill(dNchTrans, vec_tmp.Pt());
                hKCh_dPhi1min->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
              else{
                if(vec_tmp.Pt()<5) pKAvgPt_NT1max_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pKAvgPt_NT1max->Fill(dNchTrans, vec_tmp.Pt());
                hKCh_dPhi1max->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
            } else if (dPhi_segment == 2) {
              if(vec_tmp.Pt()<5) pKAvgPt_NT2_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pKAvgPt_NT2->Fill(dNchTrans, vec_tmp.Pt());
              hKCh_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
	        }
	        if(abs(id)==211) {
            if (dPhi_segment == 0) {
              if(vec_tmp.Pt()<5) pPiAvgPt_NT0_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPiAvgPt_NT0->Fill(dNchTrans, vec_tmp.Pt());
            	hPiCh_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
              if(vec_tmp.Pt()<5) pPiAvgPt_NT1_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPiAvgPt_NT1->Fill(dNchTrans, vec_tmp.Pt());
              hPiCh_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              if(isTransMin){
                if(vec_tmp.Pt()<5) pPiAvgPt_NT1min_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pPiAvgPt_NT1min->Fill(dNchTrans, vec_tmp.Pt());
                hPiCh_dPhi1min->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
              else{
                if(vec_tmp.Pt()<5) pPiAvgPt_NT1max_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pPiAvgPt_NT1max->Fill(dNchTrans, vec_tmp.Pt());
                hPiCh_dPhi1max->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
            } else if (dPhi_segment == 2) {
              if(vec_tmp.Pt()<5) pPiAvgPt_NT2_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPiAvgPt_NT2->Fill(dNchTrans, vec_tmp.Pt());
              hPiCh_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
					}
					
					if(abs(id)==2212) {
            if (dPhi_segment == 0) {
              if(vec_tmp.Pt()<5) pPAvgPt_NT0_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPAvgPt_NT0->Fill(dNchTrans, vec_tmp.Pt());
            	hProton_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
              if(vec_tmp.Pt()<5) pPAvgPt_NT1_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPAvgPt_NT1->Fill(dNchTrans, vec_tmp.Pt());
              hProton_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              if(isTransMin){
                if(vec_tmp.Pt()<5) pPAvgPt_NT1min_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pPAvgPt_NT1min->Fill(dNchTrans, vec_tmp.Pt());
                hProton_dPhi1min->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
              else{
                if(vec_tmp.Pt()<5) pPAvgPt_NT1max_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pPAvgPt_NT1max->Fill(dNchTrans, vec_tmp.Pt());
                hProton_dPhi1max->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
            } else if (dPhi_segment == 2) {
              if(vec_tmp.Pt()<5) pPAvgPt_NT2_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPAvgPt_NT2->Fill(dNchTrans, vec_tmp.Pt());
              hProton_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
					}

					if(abs(id)==431) {
            if (dPhi_segment == 0) {
              if(vec_tmp.Pt()<5) pDsAvgPt_NT0_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pDsAvgPt_NT0->Fill(dNchTrans, vec_tmp.Pt());
            	hDs_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
              if(vec_tmp.Pt()<5) pDsAvgPt_NT1_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pDsAvgPt_NT1->Fill(dNchTrans, vec_tmp.Pt());
              hDs_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 2) {
              if(vec_tmp.Pt()<5) pDsAvgPt_NT2_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pDsAvgPt_NT2->Fill(dNchTrans, vec_tmp.Pt());
              hDs_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
					}
					if(abs(id)==4122) {
            if (dPhi_segment == 0) {
              if(vec_tmp.Pt()<5) pLcAvgPt_NT0_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pLcAvgPt_NT0->Fill(dNchTrans, vec_tmp.Pt());
            	hLambc_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
              if(vec_tmp.Pt()<5) pLcAvgPt_NT1_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pLcAvgPt_NT1->Fill(dNchTrans, vec_tmp.Pt());
              hLambc_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 2) {
              if(vec_tmp.Pt()<5) pLcAvgPt_NT2_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pLcAvgPt_NT2->Fill(dNchTrans, vec_tmp.Pt());
              hLambc_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
    			}
    			}
    		//}
    	}
    	if(fabs(dRap)<0.5){
    		  if(abs(id)==3122) {
            if (dPhi_segment == 0) {
              if(vec_tmp.Pt()<5) pLAvgPt_NT0_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pLAvgPt_NT0->Fill(dNchTrans, vec_tmp.Pt());
            	hLambda_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
              if(vec_tmp.Pt()<5) pLAvgPt_NT1_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pLAvgPt_NT1->Fill(dNchTrans, vec_tmp.Pt());
              hLambda_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 2) {
              if(vec_tmp.Pt()<5) pLAvgPt_NT2_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pLAvgPt_NT2->Fill(dNchTrans, vec_tmp.Pt());
              hLambda_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
					}
					if(abs(id)==411) {
            if (dPhi_segment == 0) {
              if(vec_tmp.Pt()<5) pD0AvgPt_NT0_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pD0AvgPt_NT0->Fill(dNchTrans, vec_tmp.Pt());
            	hD0_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
              if(vec_tmp.Pt()<5) pD0AvgPt_NT1_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pD0AvgPt_NT1->Fill(dNchTrans, vec_tmp.Pt());
              hD0_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 2) {
              if(vec_tmp.Pt()<5) pD0AvgPt_NT2_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pD0AvgPt_NT2->Fill(dNchTrans, vec_tmp.Pt());
              hD0_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
					}
					if(abs(id)==310) {
            if (dPhi_segment == 0) {
              if(vec_tmp.Pt()<5) pKSAvgPt_NT0_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pKSAvgPt_NT0->Fill(dNchTrans, vec_tmp.Pt());
            	hKshort_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
              if(vec_tmp.Pt()<5) pKSAvgPt_NT1_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pKSAvgPt_NT1->Fill(dNchTrans, vec_tmp.Pt());
              hKshort_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 2) {
              if(vec_tmp.Pt()<5) pKSAvgPt_NT2_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pKSAvgPt_NT2->Fill(dNchTrans, vec_tmp.Pt());
              hKshort_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
					}
    	}
    }//track loop end

    for(int idecay=0; idecay<decayList.size(); idecay++){
      int id = decayList.at(idecay).Id;
      int MoId = decayList.at(idecay).MoId;
      vec_tmp = decayList.at(idecay).vec;
      double dEta = vec_tmp.Eta();
      double dRap = vec_tmp.Rapidity();
      double dPt = vec_tmp.Pt();
      double E = vec_tmp.E();

			if (fabs(dEta)>=0.8) continue;
    	
	    double Phi = vec_tmp.Phi();
	    double dPhi = Phi - maxPhi;
	    double dPhiAbs = fabs(TVector2::Phi_mpi_pi(dPhi));
      bool isTransMin = false;

	    int dPhi_segment = 0;
		  if (dPhiAbs >= 0 && dPhiAbs < TMath::Pi() / 3) {
	      dPhi_segment = 0; //to
	      //cout << "Particle " << iTrack << " in 'to' region - phi: " << Phi* 180.0 / M_PI << ", dPhi: " << dPhi* 180.0 / M_PI << ", maxPhi: " << maxPhi* 180.0 / M_PI << endl;
      }else if (dPhiAbs >= 2 * TMath::Pi() / 3 && dPhiAbs < TMath::Pi()) {
	      dPhi_segment = 2; //away
	      //cout << "Particle " << iTrack << " in 'away' region - phi: " << Phi* 180.0 / M_PI << ", dPhi: " << dPhi* 180.0 / M_PI << ", maxPhi: " << maxPhi* 180.0 / M_PI << endl;
	    } else {
	      dPhi_segment = 1; //tr
	      //cout << "Particle " << iTrack << " in 'tr' region - phi: " << Phi* 180.0 / M_PI << ", dPhi: " << dPhi* 180.0 / M_PI << ", maxPhi: " << maxPhi* 180.0 / M_PI << endl;
        
        //if TransI==TransII, define dPhi>0 as TransMin
        if (dNchTransMin==dNchTransMax){
          if (dPhi>0) isTransMin=true;
        }
        //if TransI!=TransII, simply choose the min region
        else{
          if (dNchTransII == dNchTransMin && dPhi < 0) isTransMin=true;
          if (dNchTransI == dNchTransMin && dPhi > 0) isTransMin=true;
        }
	    }

			if (GetCharge(id) != 0) {
        if(fabs(dRap)<0.5){
	        if(abs(id)==321) {
           	if (dPhi_segment == 0) {
              if(vec_tmp.Pt()<5) pKAvgPt_NT0_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pKAvgPt_NT0->Fill(dNchTrans, vec_tmp.Pt());
            	hKCh_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
              if(vec_tmp.Pt()<5) pKAvgPt_NT1_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pKAvgPt_NT1->Fill(dNchTrans, vec_tmp.Pt());
              hKCh_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              if(isTransMin){
                if(vec_tmp.Pt()<5) pKAvgPt_NT1min_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pKAvgPt_NT1min->Fill(dNchTrans, vec_tmp.Pt());
                hKCh_dPhi1min->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
              else{
                if(vec_tmp.Pt()<5) pKAvgPt_NT1max_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pKAvgPt_NT1max->Fill(dNchTrans, vec_tmp.Pt());
                hKCh_dPhi1max->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
            } else if (dPhi_segment == 2) {
              if(vec_tmp.Pt()<5) pKAvgPt_NT2_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pKAvgPt_NT2->Fill(dNchTrans, vec_tmp.Pt());
              hKCh_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
	        }
	        if(abs(id)==211) {
            if (dPhi_segment == 0) {
              if(vec_tmp.Pt()<5) pPiAvgPt_NT0_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPiAvgPt_NT0->Fill(dNchTrans, vec_tmp.Pt());
            	hPiCh_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
              if(vec_tmp.Pt()<5) pPiAvgPt_NT1_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPiAvgPt_NT1->Fill(dNchTrans, vec_tmp.Pt());
              hPiCh_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              if(isTransMin){
                if(vec_tmp.Pt()<5) pPiAvgPt_NT1min_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pPiAvgPt_NT1min->Fill(dNchTrans, vec_tmp.Pt());
                hPiCh_dPhi1min->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
              else{
                if(vec_tmp.Pt()<5) pPiAvgPt_NT1max_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pPiAvgPt_NT1max->Fill(dNchTrans, vec_tmp.Pt());
                hPiCh_dPhi1max->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
            } else if (dPhi_segment == 2) {
              if(vec_tmp.Pt()<5) pPiAvgPt_NT2_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPiAvgPt_NT2->Fill(dNchTrans, vec_tmp.Pt());
              hPiCh_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
					}
					
					if(abs(id)==2212) {
            if (dPhi_segment == 0) {
              if(vec_tmp.Pt()<5) pPAvgPt_NT0_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPAvgPt_NT0->Fill(dNchTrans, vec_tmp.Pt());
            	hProton_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
              if(vec_tmp.Pt()<5) pPAvgPt_NT1_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPAvgPt_NT1->Fill(dNchTrans, vec_tmp.Pt());
              hProton_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              if(isTransMin){
                if(vec_tmp.Pt()<5) pPAvgPt_NT1min_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pPAvgPt_NT1min->Fill(dNchTrans, vec_tmp.Pt());
                hProton_dPhi1min->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
              else{
                if(vec_tmp.Pt()<5) pPAvgPt_NT1max_ptcut->Fill(dNchTrans, vec_tmp.Pt());
                pPAvgPt_NT1max->Fill(dNchTrans, vec_tmp.Pt());
                hProton_dPhi1max->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
              }
            } else if (dPhi_segment == 2) {
              if(vec_tmp.Pt()<5) pPAvgPt_NT2_ptcut->Fill(dNchTrans, vec_tmp.Pt());
              pPAvgPt_NT2->Fill(dNchTrans, vec_tmp.Pt());
              hProton_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
					}

          if(abs(id)==431) {
            if (dPhi_segment == 0) {
                hDs_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
                hDs_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 2) {
                hDs_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
          }
          if(abs(id)==4122) {
            if (dPhi_segment == 0) {
                hLambc_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 1) {
                hLambc_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            } else if (dPhi_segment == 2) {
                hLambc_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
            }
          }
        }
    	}
    	if(fabs(dRap)<0.5){
    		  if(abs(id)==3122) {
                if (dPhi_segment == 0) {
                    hLambda_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
                } else if (dPhi_segment == 1) {
                    hLambda_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
                } else if (dPhi_segment == 2) {
                    hLambda_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
                }
					}
					if(abs(id)==411) {
                if (dPhi_segment == 0) {
                    hD0_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
                } else if (dPhi_segment == 1) {
                    hD0_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
                } else if (dPhi_segment == 2) {
                    hD0_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
                }
					}
					if(abs(id)==310) {
                if (dPhi_segment == 0) {
                    hKshort_dPhi0->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
                } else if (dPhi_segment == 1) {
                    hKshort_dPhi1->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
                } else if (dPhi_segment == 2) {
                    hKshort_dPhi2->Fill(vec_tmp.Pt(), dNtrkV0, dNchTrans);
                }
					}
    	}
    }//decay loop
  }
//=============================================================================

	
  hTrials->Fill(0.5, nEvents);
	hChEta->Scale(1./nEvents_Inel, "width");

  //--------------------------------
  // Write Hists and exit                               
  //--------------------------------

  TFile *file = new TFile(outputFile,"recreate"); 
  file->cd(); 
	list->Write(); 
tTRjet->Write();
	file->Close();

  return 0;
}


//--------------------------------
// Define functions here                               
//--------------------------------


float GetSphericity(vector<TLorentzVector> tracks)
{
  // a critical check
  if(tracks.size()<=1) return 0;
  
	double sumPt = 0;
  // first fill the momentum tensor
//  TMatrixDSym MomentumTensor(3);
  TMatrixDSym MomentumTensor(2);
  for(vector<TLorentzVector>::iterator itTrack = tracks.begin(); itTrack!=tracks.end(); ++itTrack) {
//  std::vector<double> momentum(3);
  std::vector<double> momentum(2);
  momentum[0] = itTrack->X();
  momentum[1] = itTrack->Y();
	double trackPt = itTrack->Pt();
	sumPt+=trackPt;
//  momentum[2] = itTrack->Z();
//    for(unsigned int i=0;i<3;i++)
    for(unsigned int i=0;i<2;i++)
      for(unsigned int j=0;j<=i;j++) {
        MomentumTensor[i][j] += momentum[i]*momentum[j]/trackPt;
      }
  }
//  MomentumTensor*=1/(MomentumTensor[0][0]+MomentumTensor[1][1]+MomentumTensor[2][2]);
  MomentumTensor*=1/(sumPt);
  // find the eigen values
  TMatrixDSymEigen eigen(MomentumTensor);
  TVectorD eigenvals = eigen.GetEigenValues();
//  vector<float> eigenvaluess(3);
  vector<float> eigenvaluess(2);
  eigenvaluess[0] = eigenvals[0];
  eigenvaluess[1] = eigenvals[1];
//  eigenvaluess[2] = eigenvals[2];
	//sort the eigen value from low to high
  sort(eigenvaluess.begin(),eigenvaluess.end());
  // compute spericity
//  float sph = ( 1.5*(1-eigenvaluess[2]));
  float sph = 2*eigenvaluess[0]/(eigenvaluess[0]+eigenvaluess[1]);
  return sph;
}

int GetCharge(int id)
{
  int charge = 0;

  // QUARKS
  if      (id==1)    charge=-1/3.; // d
  else if (id==-1)   charge=1/3.;
  else if (id==2)    charge=2/3.;  // u
  else if (id==-2)   charge=-2/3.;
  else if (id==3)    charge=-1/3.; // s
  else if (id==-3)   charge=1/3.;
  else if (id==4)    charge=2/3.;  // c
  else if (id==-4)   charge=-2/3.;
  else if (id==5)    charge=-1/3.; // b
  else if (id==-5)   charge=1/3.;
  else if (id==6)    charge=2/3.;  // t
  else if (id==-6)   charge=-2/3.;
  // LEPTONS
  if      (id==11)   charge=-1;    // e-
  else if (id==-11)  charge=1;
  // LIGHT I = 1 MESONS
  else if (id==211)  charge=1;     // π+
  else if (id==213)  charge=1;     // ρ(770)+
  else if (id==-211) charge=-1;
  else if (id==-213) charge=-1;
  // STRANGE MESONS
  else if (id==321)  charge=1;     // Κ+
  else if (id==323)  charge=1;     // Κ*(892)+
  else if (id==-321) charge=-1;
  else if (id==-323) charge=-1;
  // CHARMED MESONS
  else if (id==411)  charge=1;     // D+
  else if (id==413)  charge=1;     // D*(2010)+
  else if (id==431)  charge=1;     // Ds+
  else if (id==433)  charge=1;     // Ds*+
  else if (id==-411) charge=-1;
  else if (id==-413) charge=-1;
  else if (id==-431) charge=-1;
  else if (id==-433) charge=-1;
  // BOTTOM MESONS
  else if (id==521)  charge=1;     // B+
  else if (id==523)  charge=1;     // B*+
  else if (id==541)  charge=1;     // Bc+
  else if (id==543)  charge=1;     // Bc*+
  else if (id==-521) charge=-1;
  else if (id==-523) charge=-1;
  else if (id==-541)  charge=-1;     
  else if (id==-543)  charge=-1;     

  // LIGHT BARYONS
  else if (id==2212)   charge=1;    // p+
  else if (id==2224)   charge=2;    // Δ++
  else if (id==2214)   charge=1;    // Δ+
  else if (id==1114)   charge=-1;   // Δ-
  else if (id==-2212)  charge=-1;
  else if (id==-2224)  charge=-2;
  else if (id==-2214)  charge=-1;
  else if (id==-1114)  charge=1;
  // STRANGE BARYONS
  else if (id==3222) charge=1;    // Σ+
  else if (id==3112) charge=-1;   // Σ-
  else if (id==3224) charge=1;    // Σ*+
  else if (id==3114) charge=-1;   // Σ*-
  else if (id==3312) charge=-1;   // Ξ-
  else if (id==3314) charge=-1;   // Ξ*-
  else if (id==3334) charge=-1;   // Ω-
  else if (id==-3222) charge=-1;
  else if (id==-3112) charge=1;
  else if (id==-3224) charge=-1;
  else if (id==-3114) charge=1;
  else if (id==-3312) charge=1;
  else if (id==-3314) charge=1;
  else if (id==-3334) charge=1;
  // CHARMED BARYONS
  else if (id==4122) charge=1;    // Λc+
  else if (id==4222) charge=2;    // Σc++
  else if (id==4212) charge=1;    // Σc+
  else if (id==4224) charge=2;    // Σc*++
  else if (id==4214) charge=1;    // Σc*+
  else if (id==4232) charge=1;    // Ξc+
  else if (id==4322) charge=1;    // Ξ'c+
  else if (id==4324) charge=1;    // Ξc*+
  else if (id==4412) charge=1;    // Ξcc+
  else if (id==4422) charge=2;    // Ξcc++
  else if (id==4414) charge=1;    // Ξcc*+
  else if (id==4424) charge=2;    // Ξcc*++
  else if (id==4432) charge=1;    //Ωcc+
  else if (id==4434) charge=1;    //Ωcc*+
  else if (id==4444) charge=2;    //Ωccc++
  else if (id==-4122) charge=-1;
  else if (id==-4222) charge=-2;
  else if (id==-4212) charge=-1;
  else if (id==-4224) charge=-2;
  else if (id==-4214) charge=-1;
  else if (id==-4232) charge=-1;
  else if (id==-4322) charge=-1;
  else if (id==-4324) charge=-1;
  else if (id==4412) charge=-1;
  else if (id==4422) charge=-2;
  else if (id==4414) charge=-1;
  else if (id==4424) charge=-2;
  else if (id==4432) charge=-1;
  else if (id==4434) charge=-1;
  else if (id==4444) charge=-2;
  // BOTTOMED BARYONS
  else if (id==5112) charge=-1;    // Σb-
  else if (id==5222) charge=1;    // Σb+
  else if (id==5114) charge=-1;    // Σb*-
  else if (id==5224) charge=1;    // Σb*+
  else if (id==5132) charge=-1;    // Ξb-
  else if (id==5312) charge=-1;    // Ξ'b-
  else if (id==5314) charge=-1;    // Ξb*-
  else if (id==5332) charge=-1;    //Ωb-
  else if (id==5334) charge=-1;    //Ωb*-
  else if (id==5242) charge=1;    // Ξbc+
  else if (id==5422) charge=1;    // Ξ'bc+
  else if (id==5424) charge=1;    // Ξ'bc*+
  else if (id==5442) charge=1;    //Ωbcc+
  else if (id==5444) charge=1;    //Ωbcc*+
  else if (id==5512) charge=-1;    //Ξbb-
  else if (id==5514) charge=-1;    //Ξbb*-
  else if (id==5532) charge=-1;    //Ωbb-
  else if (id==5534) charge=-1;    //Ωbb*-
  else if (id==5554) charge=-1;    //Ωbbb-
  else if (id==-5112) charge=1; 
  else if (id==-5222) charge=-1;  
  else if (id==-5114) charge=1; 
  else if (id==-5224) charge=-1;  
  else if (id==-5132) charge=1; 
  else if (id==-5312) charge=1; 
  else if (id==-5314) charge=1; 
  else if (id==-5332) charge=1; 
  else if (id==-5334) charge=1; 
  else if (id==-5242) charge=-1;  
  else if (id==-5422) charge=-1;  
  else if (id==-5424) charge=-1;  
  else if (id==-5442) charge=-1;  
  else if (id==-5444) charge=-1;  
  else if (id==-5512) charge=1; 
  else if (id==-5514) charge=1; 
  else if (id==-5532) charge=1; 
  else if (id==-5534) charge=1; 
  else if (id==-5554) charge=1; 

  return charge;
}

float Pythiadecaymass(int aPid) 
{

    double mass=0;

    if(aPid==421) mass=1.8645;
    if(aPid==411) mass=1.8693;
    if(aPid==431) mass=1.9685;
    if(aPid==423) mass=2.0067;
    if(aPid==433) mass=2.1124;
    if(aPid==3112) mass=1.19744;
    if(aPid==3222) mass=1.18937;
    if(aPid==4112) mass=2.4521;
    if(aPid==4122) mass=2.2849;
    if(aPid==4212) mass=2.4535;
    if(aPid==4222) mass=2.4529;
    if(aPid==4132) mass=2.4703;
    if(aPid==4332) mass=2.704;
    if(aPid==4232) mass=2.4656;
    if(aPid==4432) mass=3.78663;

    return mass;
}


bool IsFinal(int aPid){

  bool finalFlag=false;

  const int N=24;

  int finalList[N]={ 111, 211, 2212, 2112, 
    321, 130, 310, 3122, 3312, 3334,
    411, 421, 431, 4122, 4132, 4232, 4332,
    511, 521, 531, 5122, 5132, 5232, 5332
  };

  for(int i=0; i<N; i++){
    if(aPid==finalList[i]) finalFlag=true;
  }

  return finalFlag;
}
