//--------------------------------
// Load header                              
//--------------------------------

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PlotFile;
#endif

#ifndef __CINT__
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "math.h"
#include "string.h"

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TComplex.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"


using namespace std;

#endif

#include "AMPT.h"

//--------------------------------
// Declare functions                                
//--------------------------------

int GetCentrality(AMPT* ampt);
int GetCharge(int id);

float GetSphericity(vector<TVector3> tracks);

void FillDPhiDEta(TVector3 vTr, TVector3 vAs, TH2D *hout);
void FillDPhiDEta3D(TVector3 vTr, TVector3 vAs, TH3D *hout);
void PrintArray(double array[], int size);

//--------------------------------
// Main function starts here                               
//--------------------------------

int main(int argc, char** argv)
{
  if (argc!=3) {
	cout<<"Parameter number wrong!"<<endl;
	return 1;
  }

  //   char* inputFile;
  //   char* outputFile;
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


  //--------------------------------
  // define histograms
  //--------------------------------
  const double PI = 3.1415926;
  TH1::SetDefaultSumw2();
	TH1D *hNChMid = new TH1D("hNChMid", ";N_{ch}; counts", 2000, 0, 4000);

	TH2D *hDEtaDPhiSameEventTPCFMD12 = new TH2D("hDEtaDPhiSameEventTPCFMD12", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
  TH2D *hDEtaDPhiSameEventTPCFMD3 = new TH2D("hDEtaDPhiSameEventTPCFMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
  TH2D *hDEtaDPhiSameEventFMD12FMD3 = new TH2D("hDEtaDPhiSameEventFMD12FMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiMixEventTPCFMD12 = new TH2D("hDEtaDPhiMixEventTPCFMD12", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
  TH2D *hDEtaDPhiMixEventTPCFMD3 = new TH2D("hDEtaDPhiMixEventTPCFMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
  TH2D *hDEtaDPhiMixEventFMD12FMD3 = new TH2D("hDEtaDPhiMixEventFMD12FMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

	TH2D *hDEtaDPhiSameEventLowMidTPCFMD12 = new TH2D("hDEtaDPhiSameEventLowMidTPCFMD12", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiSameEventLowMidTPCFMD3 = new TH2D("hDEtaDPhiSameEventLowMidTPCFMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiSameEventLowMidFMD12FMD3 = new TH2D("hDEtaDPhiSameEventLowMidFMD12FMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiMixEventLowMidTPCFMD12 = new TH2D("hDEtaDPhiMixEventLowMidTPCFMD12", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiMixEventLowMidTPCFMD3 = new TH2D("hDEtaDPhiMixEventLowMidTPCFMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiMixEventLowMidFMD12FMD3 = new TH2D("hDEtaDPhiMixEventLowMidFMD12FMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

	TH2D *hDEtaDPhiSameEventHighMidTPCFMD12 = new TH2D("hDEtaDPhiSameEventHighMidTPCFMD12", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiSameEventHighMidTPCFMD3 = new TH2D("hDEtaDPhiSameEventHighMidTPCFMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiSameEventHighMidFMD12FMD3 = new TH2D("hDEtaDPhiSameEventHighMidFMD12FMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiMixEventHighMidTPCFMD12 = new TH2D("hDEtaDPhiMixEventHighMidTPCFMD12", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiMixEventHighMidTPCFMD3 = new TH2D("hDEtaDPhiMixEventHighMidTPCFMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiMixEventHighMidFMD12FMD3 = new TH2D("hDEtaDPhiMixEventHighMidFMD12FMD3", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

	TH3D *hDEtaDPhiTrigEtaSameEventLowMidTPCFMD12 = new TH3D("hDEtaDPhiTrigEtaSameEventLowMidTPCFMD12", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaSameEventLowMidTPCFMD3 = new TH3D("hDEtaDPhiTrigEtaSameEventLowMidTPCFMD3", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaSameEventLowMidFMD12FMD3 = new TH3D("hDEtaDPhiTrigEtaSameEventLowMidFMD12FMD3", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaMixEventLowMidTPCFMD12 = new TH3D("hDEtaDPhiTrigEtaMixEventLowMidTPCFMD12", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaMixEventLowMidTPCFMD3 = new TH3D("hDEtaDPhiTrigEtaMixEventLowMidTPCFMD3", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaMixEventLowMidFMD12FMD3 = new TH3D("hDEtaDPhiTrigEtaMixEventLowMidFMD12FMD3", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);

	TH3D *hDEtaDPhiTrigEtaSameEventHighMidTPCFMD12 = new TH3D("hDEtaDPhiTrigEtaSameEventHighMidTPCFMD12", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaSameEventHighMidTPCFMD3 = new TH3D("hDEtaDPhiTrigEtaSameEventHighMidTPCFMD3", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaSameEventHighMidFMD12FMD3 = new TH3D("hDEtaDPhiTrigEtaSameEventHighMidFMD12FMD3", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaMixEventHighMidTPCFMD12 = new TH3D("hDEtaDPhiTrigEtaMixEventHighMidTPCFMD12", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaMixEventHighMidTPCFMD3 = new TH3D("hDEtaDPhiTrigEtaMixEventHighMidTPCFMD3", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaMixEventHighMidFMD12FMD3 = new TH3D("hDEtaDPhiTrigEtaMixEventHighMidFMD12FMD3", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);

  TH1D *hTrigPtTPCFMD12 = new TH1D("hTrigPtTPCFMD12", ";p_{T}", 50, 0, 10);
  TH1D *hTrigPtTPCFMD3 = new TH1D("hTrigPtTPCFMD3", ";p_{T}", 50, 0, 10);
  TH1D *hTrigPtFMD12FMD3 = new TH1D("hTrigPtFMD12FMD3", ";p_{T}", 50, 0, 10);
  TH1D *hTrigPtLowTPCFMD12 = new TH1D("hTrigPtLowTPCFMD12", ";p_{T}", 50, 0, 10);
  TH1D *hTrigPtLowTPCFMD3 = new TH1D("hTrigPtLowTPCFMD3", ";p_{T}", 50, 0, 10);
  TH1D *hTrigPtLowFMD12FMD3 = new TH1D("hTrigPtLowFMD12FMD3", ";p_{T}", 50, 0, 10);
  TH1D *hTrigPtHighTPCFMD12 = new TH1D("hTrigPtHighTPCFMD12", ";p_{T}", 50, 0, 10);
  TH1D *hTrigPtHighTPCFMD3 = new TH1D("hTrigPtHighTPCFMD3", ";p_{T}", 50, 0, 10);
  TH1D *hTrigPtHighFMD12FMD3 = new TH1D("hTrigPtHighFMD12FMD3", ";p_{T}", 50, 0, 10);
  TH2D *hTrigPtEtaLowTPCFMD12 = new TH2D("hTrigPtEtaLowTPCFMD12", ";p_{T};#eta_{tr}", 50, 0, 10, 60, -3, 3);
  TH2D *hTrigPtEtaLowTPCFMD3 = new TH2D("hTrigPtEtaLowTPCFMD3", ";p_{T};#eta_{tr}", 50, 0, 10, 60, -3, 3);
  TH2D *hTrigPtEtaLowFMD12FMD3 = new TH2D("hTrigPtEtaLowFMD12FMD3", ";p_{T};#eta_{tr}", 50, 0, 10, 60, -3, 3);
  TH2D *hTrigPtEtaHighTPCFMD12 = new TH2D("hTrigPtEtaHighTPCFMD12", ";p_{T};#eta_{tr}", 50, 0, 10, 60, -3, 3);
  TH2D *hTrigPtEtaHighTPCFMD3 = new TH2D("hTrigPtEtaHighTPCFMD3", ";p_{T};#eta_{tr}", 50, 0, 10, 60, -3, 3);
  TH2D *hTrigPtEtaHighFMD12FMD3 = new TH2D("hTrigPtEtaHighFMD12FMD3", ";p_{T};#eta_{tr}", 50, 0, 10, 60, -3, 3);

  //set 9 event classes, 10 Nsel bin edge
  //include 0-10 bin for consistency check
  const int N_cent=9;
  double Nsel_min=0;
  double Nsel_bin=10;
  double Nsel_cut[N_cent+1]={0};
  TH2D *hDEtaDPhiSameEvent_CentTPCFMD12[N_cent];
  TH2D *hDEtaDPhiSameEvent_CentTPCFMD3[N_cent];
  TH2D *hDEtaDPhiSameEvent_CentFMD12FMD3[N_cent];
  TH2D *hDEtaDPhiMixEvent_CentTPCFMD12[N_cent];
  TH2D *hDEtaDPhiMixEvent_CentTPCFMD3[N_cent];
  TH2D *hDEtaDPhiMixEvent_CentFMD12FMD3[N_cent];
  TH1D *hTrigPt_CentTPCFMD12[N_cent];
  TH1D *hTrigPt_CentTPCFMD3[N_cent];
  TH1D *hTrigPt_CentFMD12FMD3[N_cent];
  
  TH3D *hDEtaDPhiTrigEtaSameEvent_CentTPCFMD12[N_cent];
  TH3D *hDEtaDPhiTrigEtaSameEvent_CentTPCFMD3[N_cent];
  TH3D *hDEtaDPhiTrigEtaSameEvent_CentFMD12FMD3[N_cent];
  TH3D *hDEtaDPhiTrigEtaMixEvent_CentTPCFMD12[N_cent];
  TH3D *hDEtaDPhiTrigEtaMixEvent_CentTPCFMD3[N_cent];
  TH3D *hDEtaDPhiTrigEtaMixEvent_CentFMD12FMD3[N_cent];
  TH2D *hTrigPtEta_CentTPCFMD12[N_cent];
  TH2D *hTrigPtEta_CentTPCFMD3[N_cent];
  TH2D *hTrigPtEta_CentFMD12FMD3[N_cent];

  for(int i=0; i<N_cent; i++) {
    Nsel_cut[i]=Nsel_min+i*Nsel_bin;
    if(i==N_cent-1) Nsel_cut[N_cent]=Nsel_cut[i]+Nsel_bin;
    hDEtaDPhiSameEvent_CentTPCFMD12[i] = new TH2D(Form("hDEtaDPhiSameEvent_CentTPCFMD12%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiSameEvent_CentTPCFMD3[i] = new TH2D(Form("hDEtaDPhiSameEvent_CentTPCFMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiSameEvent_CentFMD12FMD3[i] = new TH2D(Form("hDEtaDPhiSameEvent_CentFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEvent_CentTPCFMD12[i] = new TH2D(Form("hDEtaDPhiMixEvent_CentTPCFMD12%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEvent_CentTPCFMD3[i] = new TH2D(Form("hDEtaDPhiMixEvent_CentTPCFMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEvent_CentFMD12FMD3[i] = new TH2D(Form("hDEtaDPhiMixEvent_CentFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

    hDEtaDPhiTrigEtaSameEvent_CentTPCFMD12[i] = new TH3D(Form("hDEtaDPhiTrigEtaSameEvent_CentTPCFMD12%d",i),"; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
    hDEtaDPhiTrigEtaSameEvent_CentTPCFMD3[i] = new TH3D(Form("hDEtaDPhiTrigEtaSameEvent_CentTPCFMD3%d",i),"; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
    hDEtaDPhiTrigEtaSameEvent_CentFMD12FMD3[i] = new TH3D(Form("hDEtaDPhiTrigEtaSameEvent_CentFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
    hDEtaDPhiTrigEtaMixEvent_CentTPCFMD12[i] = new TH3D(Form("hDEtaDPhiTrigEtaMixEvent_CentTPCFMD12%d",i),"; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
    hDEtaDPhiTrigEtaMixEvent_CentTPCFMD3[i] = new TH3D(Form("hDEtaDPhiTrigEtaMixEvent_CentTPCFMD3%d",i),"; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
    hDEtaDPhiTrigEtaMixEvent_CentFMD12FMD3[i] = new TH3D(Form("hDEtaDPhiTrigEtaMixEvent_CentFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);

    hTrigPt_CentTPCFMD12[i]= new TH1D(Form("hTrigPt_CentTPCFMD12%d",i), ";p_{T}", 50, 0, 10);
    hTrigPt_CentTPCFMD3[i]= new TH1D(Form("hTrigPt_CentTPCFMD3%d",i), ";p_{T}", 50, 0, 10);
    hTrigPt_CentFMD12FMD3[i]= new TH1D(Form("hTrigPt_CentFMD12FMD3%d",i), ";p_{T}", 50, 0, 10);
    hTrigPtEta_CentTPCFMD12[i]= new TH2D(Form("hTrigPtEta_CentTPCFMD12%d",i), ";p_{T}", 50, 0, 10, 60, -3, 3);
    hTrigPtEta_CentTPCFMD3[i]= new TH2D(Form("hTrigPtEta_CentTPCFMD3%d",i), ";p_{T}", 50, 0, 10, 60, -3, 3);
    hTrigPtEta_CentFMD12FMD3[i]= new TH2D(Form("hTrigPtEta_CentFMD12FMD3%d",i), ";p_{T}", 50, 0, 10, 60, -3, 3);
  }
  cout<<"cent bin edge:"<<endl;
  PrintArray(Nsel_cut, N_cent+1);


  //8 pt bin interval, 9 pt bin edge
  //const int N_ptbin=8;
  //10 pt bin interval, 11 pt bin edge
  const int N_ptbin=10;
  double pt_min=0;
  double pt_bin=0.5;
  //double ptbin_cut[N_ptbin+1]={0};
  double ptbin_cut[N_ptbin+1]={0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2.0, 2.5, 3.5, 4.5};
  //keep N_ptbin hist in N_ptbin ptbin intervals
  //associate in the fixed range e.g.(0,3), trigger in different bins
  TH2D *hDEtaDPhiSameEventHigh_ptbinTPCFMD12[N_ptbin];
  TH2D *hDEtaDPhiMixEventHigh_ptbinTPCFMD12[N_ptbin];
  TH2D *hDEtaDPhiSameEventLow_ptbinTPCFMD12[N_ptbin];
  TH2D *hDEtaDPhiMixEventLow_ptbinTPCFMD12[N_ptbin];
  TH1D *hTrigPtHigh_ptbinTPCFMD12[N_ptbin];
  TH1D *hTrigPtLow_ptbinTPCFMD12[N_ptbin];

  TH2D *hDEtaDPhiSameEventHigh_ptbinTPCFMD3[N_ptbin];
  TH2D *hDEtaDPhiMixEventHigh_ptbinTPCFMD3[N_ptbin];
  TH2D *hDEtaDPhiSameEventLow_ptbinTPCFMD3[N_ptbin];
  TH2D *hDEtaDPhiMixEventLow_ptbinTPCFMD3[N_ptbin];
  TH1D *hTrigPtHigh_ptbinTPCFMD3[N_ptbin];
  TH1D *hTrigPtLow_ptbinTPCFMD3[N_ptbin];

  TH2D *hDEtaDPhiSameEventHigh_ptbinFMD12FMD3[N_ptbin];
  TH2D *hDEtaDPhiMixEventHigh_ptbinFMD12FMD3[N_ptbin];
  TH2D *hDEtaDPhiSameEventLow_ptbinFMD12FMD3[N_ptbin];
  TH2D *hDEtaDPhiMixEventLow_ptbinFMD12FMD3[N_ptbin];
  TH1D *hTrigPtHigh_ptbinFMD12FMD3[N_ptbin];
  TH1D *hTrigPtLow_ptbinFMD12FMD3[N_ptbin];
  //keep trig/asso in the same pt range
  TH2D *hDEtaDPhiSameEventHigh_ptbinTASameTPCFMD12[N_ptbin];
  TH2D *hDEtaDPhiMixEventHigh_ptbinTASameTPCFMD12[N_ptbin];
  TH2D *hDEtaDPhiSameEventLow_ptbinTASameTPCFMD12[N_ptbin];
  TH2D *hDEtaDPhiMixEventLow_ptbinTASameTPCFMD12[N_ptbin];
  TH2D *hDEtaDPhiSameEventHigh_ptbinTASameTPCFMD3[N_ptbin];
  TH2D *hDEtaDPhiMixEventHigh_ptbinTASameTPCFMD3[N_ptbin];
  TH2D *hDEtaDPhiSameEventLow_ptbinTASameTPCFMD3[N_ptbin];
  TH2D *hDEtaDPhiMixEventLow_ptbinTASameTPCFMD3[N_ptbin];
  TH2D *hDEtaDPhiSameEventHigh_ptbinTASameFMD12FMD3[N_ptbin];
  TH2D *hDEtaDPhiMixEventHigh_ptbinTASameFMD12FMD3[N_ptbin];
  TH2D *hDEtaDPhiSameEventLow_ptbinTASameFMD12FMD3[N_ptbin];
  TH2D *hDEtaDPhiMixEventLow_ptbinTASameFMD12FMD3[N_ptbin];

  for(int i=0; i<N_ptbin; i++){
//    ptbin_cut[i]=pt_min+i*pt_bin;
//    if(i>4) ptbin_cut[i]=ptbin_cut[i-1]+(2*pt_bin);
//    if(i==N_ptbin-1) ptbin_cut[N_ptbin]=ptbin_cut[i]+2*pt_bin;
    hDEtaDPhiSameEventHigh_ptbinTPCFMD12[i] = new TH2D(Form("hDEtaDPhiSameEventHigh_ptbinTPCFMD12%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventHigh_ptbinTPCFMD12[i] = new TH2D(Form("hDEtaDPhiMixEventHigh_ptbinTPCFMD12%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hTrigPtHigh_ptbinTPCFMD12[i]= new TH1D(Form("hTrigPtHigh_ptbinTPCFMD12%d",i), ";p_{T}", 500, 0, 10);
    hDEtaDPhiSameEventLow_ptbinTPCFMD12[i] = new TH2D(Form("hDEtaDPhiSameEventLow_ptbinTPCFMD12%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventLow_ptbinTPCFMD12[i] = new TH2D(Form("hDEtaDPhiMixEventLow_ptbinTPCFMD12%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hTrigPtLow_ptbinTPCFMD12[i]= new TH1D(Form("hTrigPtLow_ptbinTPCFMD12%d",i), ";p_{T}", 500, 0, 10);

    hDEtaDPhiSameEventHigh_ptbinTASameTPCFMD12[i] = new TH2D(Form("hDEtaDPhiSameEventHigh_ptbinTASameTPCFMD12%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventHigh_ptbinTASameTPCFMD12[i] = new TH2D(Form("hDEtaDPhiMixEventHigh_ptbinTASameTPCFMD12%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiSameEventLow_ptbinTASameTPCFMD12[i] = new TH2D(Form("hDEtaDPhiSameEventLow_ptbinTASameTPCFMD12%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventLow_ptbinTASameTPCFMD12[i] = new TH2D(Form("hDEtaDPhiMixEventLow_ptbinTASameTPCFMD12%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

    hDEtaDPhiSameEventHigh_ptbinTPCFMD3[i] = new TH2D(Form("hDEtaDPhiSameEventHigh_ptbinTPCFMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventHigh_ptbinTPCFMD3[i] = new TH2D(Form("hDEtaDPhiMixEventHigh_ptbinTPCFMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hTrigPtHigh_ptbinTPCFMD3[i]= new TH1D(Form("hTrigPtHigh_ptbinTPCFMD3%d",i), ";p_{T}", 500, 0, 10);
    hDEtaDPhiSameEventLow_ptbinTPCFMD3[i] = new TH2D(Form("hDEtaDPhiSameEventLow_ptbinTPCFMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventLow_ptbinTPCFMD3[i] = new TH2D(Form("hDEtaDPhiMixEventLow_ptbinTPCFMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hTrigPtLow_ptbinTPCFMD3[i]= new TH1D(Form("hTrigPtLow_ptbinTPCFMD3%d",i), ";p_{T}", 500, 0, 10);

    hDEtaDPhiSameEventHigh_ptbinTASameTPCFMD3[i] = new TH2D(Form("hDEtaDPhiSameEventHigh_ptbinTASameTPCFMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventHigh_ptbinTASameTPCFMD3[i] = new TH2D(Form("hDEtaDPhiMixEventHigh_ptbinTASameTPCFMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiSameEventLow_ptbinTASameTPCFMD3[i] = new TH2D(Form("hDEtaDPhiSameEventLow_ptbinTASameTPCFMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventLow_ptbinTASameTPCFMD3[i] = new TH2D(Form("hDEtaDPhiMixEventLow_ptbinTASameTPCFMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

    hDEtaDPhiSameEventHigh_ptbinFMD12FMD3[i] = new TH2D(Form("hDEtaDPhiSameEventHigh_ptbinFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventHigh_ptbinFMD12FMD3[i] = new TH2D(Form("hDEtaDPhiMixEventHigh_ptbinFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hTrigPtHigh_ptbinFMD12FMD3[i]= new TH1D(Form("hTrigPtHigh_ptbinFMD12FMD3%d",i), ";p_{T}", 500, 0, 10);
    hDEtaDPhiSameEventLow_ptbinFMD12FMD3[i] = new TH2D(Form("hDEtaDPhiSameEventLow_ptbinFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventLow_ptbinFMD12FMD3[i] = new TH2D(Form("hDEtaDPhiMixEventLow_ptbinFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hTrigPtLow_ptbinFMD12FMD3[i]= new TH1D(Form("hTrigPtLow_ptbinFMD12FMD3%d",i), ";p_{T}", 500, 0, 10);

    hDEtaDPhiSameEventHigh_ptbinTASameFMD12FMD3[i] = new TH2D(Form("hDEtaDPhiSameEventHigh_ptbinTASameFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventHigh_ptbinTASameFMD12FMD3[i] = new TH2D(Form("hDEtaDPhiMixEventHigh_ptbinTASameFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiSameEventLow_ptbinTASameFMD12FMD3[i] = new TH2D(Form("hDEtaDPhiSameEventLow_ptbinTASameFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventLow_ptbinTASameFMD12FMD3[i] = new TH2D(Form("hDEtaDPhiMixEventLow_ptbinTASameFMD12FMD3%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

  }
  cout<<"pt bin edge:"<<endl;
  PrintArray(ptbin_cut, N_ptbin+1);


	TH3D *hNchDEtaDPhiSameEventTPCFMD12 = new TH3D("hNchDEtaDPhiSameEventTPCFMD12", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 300, -0.5, 299.5);
	TH3D *hNchDEtaDPhiMixEventTPCFMD12 = new TH3D("hNchDEtaDPhiMixEventTPCFMD12", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 300, -0.5, 299.5);
	TH3D *hNchDEtaDPhiSameEventTPCFMD3 = new TH3D("hNchDEtaDPhiSameEventTPCFMD3", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 300, -0.5, 299.5);
	TH3D *hNchDEtaDPhiMixEventTPCFMD3 = new TH3D("hNchDEtaDPhiMixEventTPCFMD3", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 300, -0.5, 299.5);
	TH3D *hNchDEtaDPhiSameEventFMD12FMD3 = new TH3D("hNchDEtaDPhiSameEventFMD12FMD3", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 300, -0.5, 299.5);
	TH3D *hNchDEtaDPhiMixEventFMD12FMD3 = new TH3D("hNchDEtaDPhiMixEventFMD12FMD3", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 300, -0.5, 299.5);

	TH2D *hSphericity = new TH2D("hSphericity", ";N_{ch}^{mid}; S", 2000, 0, 4000, 50, 0, 1);

   
  //--------------------------------
  // loop events
  //--------------------------------

  int nEvents = (int)chain->GetEntries();
  cout<<"Total events : "<<nEvents<<endl;

//  if(nEvents>8e6) nEvents=8e6;
  AMPT* ampt = new AMPT(chain);

	//cuts
  //double accept_cut = 2.4;
  double accept_cut = 6;
	double highMultCut = 80;
	double lowMultCut = 10;

	
 	//temperary containers
	TVector3 vec_tmp;
	TVector3 vec_trig;
	vector<int> iTrIndexTPCFMD12;
	vector<int> iTrPtIdxTPCFMD12;
	vector<int> iAsIndexTPCFMD12;
  vector<int> iTrIndexTPCFMD3;
	vector<int> iTrPtIdxTPCFMD3;
	vector<int> iAsIndexTPCFMD3;
  vector<int> iTrIndexFMD12FMD3;
	vector<int> iTrPtIdxFMD12FMD3;
	vector<int> iAsIndexFMD12FMD3;
	vector<TVector3> vRefArrayTPCFMD12;
  vector<TVector3> vRefArrayTPCFMD3;
  vector<TVector3> vRefArrayFMD12FMD3;


	const int iMixSize=5;
	vector< vector<TVector3> > vRefArrayPoolTPCFMD12;
	vector< vector<TVector3> > vRefArrayPoolTPCFMD3;
	vector< vector<TVector3> > vRefArrayPoolFMD12FMD3;

	vector< vector<TVector3> > vRefArrayPool_centTPCFMD12[10];
  vector< vector<TVector3> > vRefArrayPool_centTPCFMD3[10];
	vector< vector<TVector3> > vRefArrayPool_centFMD12FMD3[10];

	vector< vector<TVector3> > vRefArrayPoolLowTPCFMD12;
	vector< vector<TVector3> > vRefArrayPoolLowTPCFMD3;
	vector< vector<TVector3> > vRefArrayPoolLowFMD12FMD3;
  
	vector< vector<TVector3> > vRefArrayPoolHighTPCFMD12;
	vector< vector<TVector3> > vRefArrayPoolHighTPCFMD3;
	vector< vector<TVector3> > vRefArrayPoolHighFMD12FMD3;

	vector<TVector3> sph_container;

  double NtrigHighTPCFMD12=0;
  double NtrigLowTPCFMD12=0;
  double NtrigHighTPCFMD3=0;
  double NtrigLowTPCFMD3=0;
  double NtrigHighFMD12FMD3=0;
  double NtrigLowFMD12FMD3=0;
  
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (iEvent%1000==0) cout << "Processing event # " << iEvent << endl;
    ampt->GetEntry(iEvent);

    int nTracks = ampt->Event_multiplicity;

    //--------------------------------
    // loop tracks                               
    //--------------------------------

		//random angle to be implemented for every track
		//to eliminate the phi dependence in event plane
		double rand_phi = gRandom->Uniform(-1,1)*PI;

		iTrIndexTPCFMD12.clear();
		iTrPtIdxTPCFMD12.clear();
		iAsIndexTPCFMD12.clear();
    iTrIndexTPCFMD3.clear();
		iTrPtIdxTPCFMD3.clear();
		iAsIndexTPCFMD3.clear();
    iTrIndexFMD12FMD3.clear();
		iTrPtIdxFMD12FMD3.clear();
		iAsIndexFMD12FMD3.clear();
		vRefArrayTPCFMD12.clear();
    vRefArrayTPCFMD3.clear();
    vRefArrayFMD12FMD3.clear();
		sph_container.clear();

    int dMidCh=0;

    bool foundEast=false;
    bool foundWest=false;

    //first loop to find the trigger/assc index
    for (int iTrack = 0; iTrack < nTracks; ++iTrack) {
      int id = ampt->id[iTrack];
      float px = ampt->px[iTrack];
      float py = ampt->py[iTrack];
      float pz = ampt->pz[iTrack];
			double mass = ampt->m[iTrack];
			float energy = sqrt(px*px+py*py+pz*pz+mass*mass);
      vec_tmp.SetXYZ(px,py,pz);
      double pt = sqrt(px*px + py*py);

      if(vec_tmp.Eta()>2.8&&vec_tmp.Eta()<5.1&&energy>3) foundEast=true;
      if(vec_tmp.Eta()>-3.7&&vec_tmp.Eta()<-1.7&&energy>3) foundWest=true;
			
      //acceptance cut
			if (fabs(vec_tmp.Eta()) > accept_cut) continue;

			//randomize the event plane (b direction)
			//which eliminates the delta-phi dependence
			//in the mixed event sample
			vec_tmp.RotateZ(rand_phi);

      bool pass_id=false;
      bool pass_trig_id=false;
      if(abs(id)==211||abs(id)==321||abs(id)==2212) pass_id=true;
      //if(GetCharge(id)!=0) pass_id=true;
      //if(abs(id)==310) pass_trig_id=true; //KS
      //if(abs(id)==3122) pass_trig_id=true; //Lambda0
      //if(abs(id)==2212) pass_trig_id=true; //p 
      //if(abs(id)==321) pass_trig_id=true; //K+
      //if(abs(id)==211) pass_trig_id=true; //pi
      if(abs(id)==211||abs(id)==321||abs(id)==2212) pass_trig_id=true; //primary charge

      //if(pass_id&&pt>0.4) dMidCh++;//Nsel
      if(pass_id&&pt>0.2&&fabs(vec_tmp.Eta())<0.8) dMidCh++;//Nch ALICE

      //reference particle, fixed for charge between 0.3~3
			//if(GetCharge(id)!=0&&vec_tmp.Eta()<3.2&&vec_tmp.Eta()>2.6)
      //TPCFMD12
			if(pass_id&&vec_tmp.Eta()<3.2&&vec_tmp.Eta()>2.6) {
        sph_container.push_back(vec_tmp);
        iAsIndexTPCFMD12.push_back(iTrack);
        //randomize ref particle container for mix event
        vRefArrayTPCFMD12.push_back(vec_tmp);
      }
      //TPCFMD3
      if(pass_id&&vec_tmp.Eta()>-3.0&&vec_tmp.Eta()<-2.0) {
        sph_container.push_back(vec_tmp);
        iAsIndexTPCFMD3.push_back(iTrack);
        //randomize ref particle container for mix event
        vRefArrayTPCFMD3.push_back(vec_tmp);
      }
      //FMD12FMD
      if(pass_id&&vec_tmp.Eta()>-3.0&&vec_tmp.Eta()<-2.0) {
        sph_container.push_back(vec_tmp);
        iAsIndexFMD12FMD3.push_back(iTrack);
        //randomize ref particle container for mix event
        vRefArrayFMD12FMD3.push_back(vec_tmp);
      }

      //get trig particle index (0.3~6) and do selection in second loop
      //if need to include other Pid trig hadron, remove charge cut here
      //and do Pid cut later in the second loop
			//if(GetCharge(id)!=0&&pt>0.3&&pt<6) {
			//if(pass_id&&pt>0.001&&pt<6) {
			if(pass_trig_id&&pt>0.001&&pt<5) {
				iTrIndexTPCFMD12.push_back(iTrack);
      }
      if(pass_trig_id&&pt>0.001&&pt<5) {
				iTrIndexTPCFMD3.push_back(iTrack);
      }
      if(pass_trig_id&&pt>0.001&&pt<5) {
				iTrIndexFMD12FMD3.push_back(iTrack);
      }

    }//for first loop

    //must found at least 1 track with E>3 at each side of detector
    if(!foundEast||!foundWest) continue;
			
		//sphericity analysis
//		double dSph = GetSphericity(sph_container);
//		if(dMidCh>0) 
//			hSphericity->Fill(dMidCh, dSph);

    //if(dMidCh<10) continue;
    if(dMidCh>highMultCut) NtrigHighTPCFMD12+=iTrIndexTPCFMD12.size();
    if(dMidCh<lowMultCut) NtrigLowTPCFMD12+=iTrIndexTPCFMD12.size();
    if(dMidCh>highMultCut) NtrigHighTPCFMD3+=iTrIndexTPCFMD3.size();
    if(dMidCh<lowMultCut) NtrigLowTPCFMD3+=iTrIndexTPCFMD3.size();
    if(dMidCh>highMultCut) NtrigHighFMD12FMD3+=iTrIndexFMD12FMD3.size();
    if(dMidCh<lowMultCut) NtrigLowFMD12FMD3+=iTrIndexFMD12FMD3.size();

		hNChMid->Fill(dMidCh);
    
    //find current event centrality
    int icent=-1;
    for(int i=0; i<N_cent; i++){
      if(dMidCh>=Nsel_cut[i]&&dMidCh<Nsel_cut[i+1]){
        icent=i;
        break;
      }
    }

		//do correlation in the same event
		//reject trig and asso with the same index
    //TPCFMD12
		for (unsigned int iTr=0; iTr<iTrIndexTPCFMD12.size(); iTr++ ) {
      int iIdxTrig = iTrIndexTPCFMD12[iTr];
      vec_trig.SetXYZ(ampt->px[iIdxTrig], ampt->py[iIdxTrig], ampt->pz[iIdxTrig] );

      //fill trig info 0~3 and |Eta|<0.8
			if(vec_trig.Pt()>0.2&&vec_trig.Pt()<3&&fabs(vec_trig.Eta())<0.8) {
        hTrigPtTPCFMD12->Fill(vec_trig.Pt());
        //Nsel dependence
        if(icent>=0) {
          hTrigPt_CentTPCFMD12[icent]->Fill(vec_trig.Pt());
          hTrigPtEta_CentTPCFMD12[icent]->Fill(vec_trig.Pt(), vec_trig.Eta());
        }

        if(dMidCh<lowMultCut){
          hTrigPtLowTPCFMD12->Fill(vec_trig.Pt());
          hTrigPtEtaLowTPCFMD12->Fill(vec_trig.Pt(), vec_trig.Eta());
        }
        else if(dMidCh>highMultCut){
          hTrigPtHighTPCFMD12->Fill(vec_trig.Pt());
          hTrigPtEtaHighTPCFMD12->Fill(vec_trig.Pt(), vec_trig.Eta());
        }
      }

      //trig info pt dependence
      int findPtBin=-1;
      if(dMidCh<lowMultCut||dMidCh>highMultCut){
        for(int iptbin=0; iptbin<N_ptbin; iptbin++){
          if(vec_trig.Pt()>ptbin_cut[iptbin]&&vec_trig.Pt()<ptbin_cut[iptbin+1]){
            findPtBin=iptbin;
            if(dMidCh<lowMultCut)
              hTrigPtLow_ptbinTPCFMD12[iptbin]->Fill(vec_trig.Pt());
            else if(dMidCh>highMultCut)
              hTrigPtHigh_ptbinTPCFMD12[iptbin]->Fill(vec_trig.Pt());
            break;
          }
        }
        iTrPtIdxTPCFMD12.push_back(findPtBin);
      }
      if(findPtBin>=N_ptbin||(findPtBin<0&&findPtBin!=-1)) {
        cout<<"findPtBin="<<findPtBin<<" Error!"<<endl;
        exit(1);
      }

			for (unsigned int iAs=0; iAs<iAsIndexTPCFMD12.size(); iAs++ ) {
				if( iAsIndexTPCFMD12[iAs]!=iTrIndexTPCFMD12[iTr] ) {
					int iIdxAs = iAsIndexTPCFMD12[iAs];
					vec_tmp.SetXYZ(ampt->px[iIdxAs], ampt->py[iIdxAs], ampt->pz[iIdxAs] );

					double dDPhi = vec_tmp.DeltaPhi(vec_trig);
					if(dDPhi<-0.5*PI) dDPhi=dDPhi + 2*PI;
					if(dDPhi>1.5*PI) dDPhi=dDPhi - 2*PI;
					double dDEta = vec_tmp.Eta() - vec_trig.Eta();

          //minbias or centrality dependence
					if(vec_trig.Pt()>0.3&&vec_trig.Pt()<3) {
            FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventTPCFMD12);

            //Nsel dependence
            if(icent>=0) {
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEvent_CentTPCFMD12[icent]);
              FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEvent_CentTPCFMD12[icent]);
            }

						if(dMidCh<lowMultCut){
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLowMidTPCFMD12);
              FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEventLowMidTPCFMD12);
            }
						else if(dMidCh>highMultCut){
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHighMidTPCFMD12);
              FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEventHighMidTPCFMD12);
            }
					}

          //pt dependence
          if(findPtBin>=0&&findPtBin<N_ptbin){
            if(dMidCh<lowMultCut){
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLow_ptbinTPCFMD12[findPtBin]);
              if(vec_tmp.Pt()>ptbin_cut[findPtBin]&&vec_tmp.Pt()<ptbin_cut[findPtBin+1]) FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLow_ptbinTASameTPCFMD12[findPtBin]);
            }
            else if(dMidCh>highMultCut) {
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHigh_ptbinTPCFMD12[findPtBin]);
              if(vec_tmp.Pt()>ptbin_cut[findPtBin]&&vec_tmp.Pt()<ptbin_cut[findPtBin+1]) FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHigh_ptbinTASameTPCFMD12[findPtBin]);
            }
          }

				}//fi asso indx different from trig

			}//asso for
		}//trig for

    //TPCFMD3
		for (unsigned int iTr=0; iTr<iTrIndexTPCFMD3.size(); iTr++ ) {
      int iIdxTrig = iTrIndexTPCFMD3[iTr];
      vec_trig.SetXYZ(ampt->px[iIdxTrig], ampt->py[iIdxTrig], ampt->pz[iIdxTrig] );

      //fill trig info 0~3 and |Eta|<0.8
      if(vec_trig.Pt()>0.2&&vec_trig.Pt()<3&&fabs(vec_trig.Eta())<0.8) {
        hTrigPtTPCFMD3->Fill(vec_trig.Pt());
        //Nsel dependence
        if(icent>=0) {
          hTrigPt_CentTPCFMD3[icent]->Fill(vec_trig.Pt());
          hTrigPtEta_CentTPCFMD3[icent]->Fill(vec_trig.Pt(), vec_trig.Eta());
        }

        if(dMidCh<lowMultCut){
          hTrigPtLowTPCFMD3->Fill(vec_trig.Pt());
          hTrigPtEtaLowTPCFMD3->Fill(vec_trig.Pt(), vec_trig.Eta());
        }
        else if(dMidCh>highMultCut){
          hTrigPtHighTPCFMD3->Fill(vec_trig.Pt());
          hTrigPtEtaHighTPCFMD3->Fill(vec_trig.Pt(), vec_trig.Eta());
        }
      }

      //trig info pt dependence
      int findPtBin=-1;
      if(dMidCh<lowMultCut||dMidCh>highMultCut){
        for(int iptbin=0; iptbin<N_ptbin; iptbin++){
          if(vec_trig.Pt()>ptbin_cut[iptbin]&&vec_trig.Pt()<ptbin_cut[iptbin+1]){
            findPtBin=iptbin;
            if(dMidCh<lowMultCut)
              hTrigPtLow_ptbinTPCFMD3[iptbin]->Fill(vec_trig.Pt());
            else if(dMidCh>highMultCut)
              hTrigPtHigh_ptbinTPCFMD3[iptbin]->Fill(vec_trig.Pt());
            break;
          }
        }
        iTrPtIdxTPCFMD3.push_back(findPtBin);
      }
      if(findPtBin>=N_ptbin||(findPtBin<0&&findPtBin!=-1)) {
        cout<<"findPtBin="<<findPtBin<<" Error!"<<endl;
        exit(1);
      }

			for (unsigned int iAs=0; iAs<iAsIndexTPCFMD3.size(); iAs++ ) {
				if( iAsIndexTPCFMD3[iAs]!=iTrIndexTPCFMD3[iTr] ) {
					int iIdxAs = iAsIndexTPCFMD3[iAs];
					vec_tmp.SetXYZ(ampt->px[iIdxAs], ampt->py[iIdxAs], ampt->pz[iIdxAs] );

					double dDPhi = vec_tmp.DeltaPhi(vec_trig);
					if(dDPhi<-0.5*PI) dDPhi=dDPhi + 2*PI;
					if(dDPhi>1.5*PI) dDPhi=dDPhi - 2*PI;
					double dDEta = vec_tmp.Eta() - vec_trig.Eta();

      //minbias or centrality dependence
      if(vec_trig.Pt()>0.3&&vec_trig.Pt()<3) {
        FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventTPCFMD3);

        //Nsel dependence
        if(icent>=0) {
          FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEvent_CentTPCFMD3[icent]);
          FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEvent_CentTPCFMD3[icent]);
        }

        if(dMidCh<lowMultCut){
          FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLowMidTPCFMD3);
          FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEventLowMidTPCFMD3);
        }
        else if(dMidCh>highMultCut){
          FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHighMidTPCFMD3);
          FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEventHighMidTPCFMD3);
        }
      }

          //pt dependence
          if(findPtBin>=0&&findPtBin<N_ptbin){
            if(dMidCh<lowMultCut){
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLow_ptbinTPCFMD3[findPtBin]);
              if(vec_tmp.Pt()>ptbin_cut[findPtBin]&&vec_tmp.Pt()<ptbin_cut[findPtBin+1]) FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLow_ptbinTASameTPCFMD3[findPtBin]);
            }
            else if(dMidCh>highMultCut) {
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHigh_ptbinTPCFMD3[findPtBin]);
              if(vec_tmp.Pt()>ptbin_cut[findPtBin]&&vec_tmp.Pt()<ptbin_cut[findPtBin+1]) FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHigh_ptbinTASameTPCFMD3[findPtBin]);
            }
          }

				}//fi asso indx different from trig

			}//asso for
		}//trig for

    //FMD12FMD3
    for (unsigned int iTr=0; iTr<iTrIndexFMD12FMD3.size(); iTr++ ) {
      int iIdxTrig = iTrIndexFMD12FMD3[iTr];
      vec_trig.SetXYZ(ampt->px[iIdxTrig], ampt->py[iIdxTrig], ampt->pz[iIdxTrig] );


      //fill trig info 0~3
      if(vec_trig.Eta()>2.6&&vec_trig.Eta()<3.2) {
        hTrigPtFMD12FMD3->Fill(vec_trig.Pt());
        //Nsel dependence
        if(icent>=0) {
          hTrigPt_CentFMD12FMD3[icent]->Fill(vec_trig.Pt());
          hTrigPtEta_CentFMD12FMD3[icent]->Fill(vec_trig.Pt(), vec_trig.Eta());
        }

        if(dMidCh<lowMultCut){
          hTrigPtLowFMD12FMD3->Fill(vec_trig.Pt());
          hTrigPtEtaLowFMD12FMD3->Fill(vec_trig.Pt(), vec_trig.Eta());
        }
        else if(dMidCh>highMultCut){
          hTrigPtHighFMD12FMD3->Fill(vec_trig.Pt());
          hTrigPtEtaHighFMD12FMD3->Fill(vec_trig.Pt(), vec_trig.Eta());
        }
      }

      //trig info pt dependence
      int findPtBin=-1;
      if(dMidCh<lowMultCut||dMidCh>highMultCut){
        for(int iptbin=0; iptbin<N_ptbin; iptbin++){
          if(vec_trig.Pt()>ptbin_cut[iptbin]&&vec_trig.Pt()<ptbin_cut[iptbin+1]){
            findPtBin=iptbin;
            if(dMidCh<lowMultCut)
              hTrigPtLow_ptbinFMD12FMD3[iptbin]->Fill(vec_trig.Pt());
            else if(dMidCh>highMultCut)
              hTrigPtHigh_ptbinFMD12FMD3[iptbin]->Fill(vec_trig.Pt());
            break;
          }
        }
        iTrPtIdxFMD12FMD3.push_back(findPtBin);
      }
      if(findPtBin>=N_ptbin||(findPtBin<0&&findPtBin!=-1)) {
        cout<<"findPtBin="<<findPtBin<<" Error!"<<endl;
        exit(1);
      }

			for (unsigned int iAs=0; iAs<iAsIndexFMD12FMD3.size(); iAs++ ) {
				if( iAsIndexFMD12FMD3[iAs]!=iTrIndexFMD12FMD3[iTr] ) {
					int iIdxAs = iAsIndexFMD12FMD3[iAs];
					vec_tmp.SetXYZ(ampt->px[iIdxAs], ampt->py[iIdxAs], ampt->pz[iIdxAs] );

					double dDPhi = vec_tmp.DeltaPhi(vec_trig);
					if(dDPhi<-0.5*PI) dDPhi=dDPhi + 2*PI;
					if(dDPhi>1.5*PI) dDPhi=dDPhi - 2*PI;
					double dDEta = vec_tmp.Eta() - vec_trig.Eta();

          //minbias or centrality dependence
					if(vec_trig.Pt()>0.3&&vec_trig.Pt()<3) {
            FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventFMD12FMD3);

            //Nsel dependence
            if(icent>=0) {
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEvent_CentFMD12FMD3[icent]);
              FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEvent_CentFMD12FMD3[icent]);
            }

						if(dMidCh<lowMultCut){
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLowMidFMD12FMD3);
              FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEventLowMidFMD12FMD3);
            }
						else if(dMidCh>highMultCut){
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHighMidFMD12FMD3);
              FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEventHighMidFMD12FMD3);
            }
					}

          //pt dependence
          if(findPtBin>=0&&findPtBin<N_ptbin){
            if(dMidCh<lowMultCut){
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLow_ptbinFMD12FMD3[findPtBin]);
              if(vec_tmp.Pt()>ptbin_cut[findPtBin]&&vec_tmp.Pt()<ptbin_cut[findPtBin+1]) FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLow_ptbinTASameFMD12FMD3[findPtBin]);
            }
            else if(dMidCh>highMultCut) {
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHigh_ptbinFMD12FMD3[findPtBin]);
              if(vec_tmp.Pt()>ptbin_cut[findPtBin]&&vec_tmp.Pt()<ptbin_cut[findPtBin+1]) FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHigh_ptbinTASameFMD12FMD3[findPtBin]);
            }
          }

				}//fi asso indx different from trig

			}//asso for
		}//trig for


		//do correlation in the mix event
    //loop current found trigger particles
    //TPCFMD12
		for(unsigned int iTr=0; iTr<iTrIndexTPCFMD12.size(); iTr++ ) {
      int iIdxTrig = iTrIndexTPCFMD12[iTr];
      vec_trig.SetXYZ(ampt->px[iIdxTrig], ampt->py[iIdxTrig], ampt->pz[iIdxTrig] );
      //minbias or centrality event mixing
      if(vec_trig.Pt()>0.2&&vec_trig.Pt()<3&&fabs(vec_trig.Eta())<0.8) {
        for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolTPCFMD12.size(); iRefEvt++) if(vRefArrayPoolTPCFMD12.size()==iMixSize) {
          for(unsigned int iAs=0; iAs<vRefArrayPoolTPCFMD12.at(iRefEvt).size(); iAs++){
            FillDPhiDEta(vec_trig, vRefArrayPoolTPCFMD12.at(iRefEvt).at(iAs), hDEtaDPhiMixEventTPCFMD12);
          }//for asso track in each ref pool evt
        }//for ref pool loop minbias

        if(icent>=0){
          for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPool_centTPCFMD12[icent].size(); iRefEvt++) if(vRefArrayPool_centTPCFMD12[icent].size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPool_centTPCFMD12[icent].at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPool_centTPCFMD12[icent].at(iRefEvt).at(iAs), hDEtaDPhiMixEvent_CentTPCFMD12[icent]);
              FillDPhiDEta3D(vec_trig, vRefArrayPool_centTPCFMD12[icent].at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEvent_CentTPCFMD12[icent]);
            }
          }
        }

        if(dMidCh<lowMultCut){
          for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolLowTPCFMD12.size(); iRefEvt++) if(vRefArrayPoolLowTPCFMD12.size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPoolLowTPCFMD12.at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPoolLowTPCFMD12.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLowMidTPCFMD12);
              FillDPhiDEta3D(vec_trig, vRefArrayPoolLowTPCFMD12.at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEventLowMidTPCFMD12);
            }//for asso track in each ref pool evt
          }//for ref pool loop lowMult
        }//fi low Mult
        else if(dMidCh>highMultCut){
          for(unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolHighTPCFMD12.size(); iRefEvt++) if(vRefArrayPoolHighTPCFMD12.size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPoolHighTPCFMD12.at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPoolHighTPCFMD12.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHighMidTPCFMD12);
              FillDPhiDEta3D(vec_trig, vRefArrayPoolHighTPCFMD12.at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEventHighMidTPCFMD12);
            }//for asso track in each ref pool evt
          }//for ref pool loop highMult
        }//else fi high Mult
      }//fi minbias or centrality trig cut 0~3

      //pt dependence
      if(dMidCh>highMultCut||dMidCh<lowMultCut){
        int iTrPt=iTrPtIdxTPCFMD12.at(iTr);
        if(iTrPt>=0&&iTrPt<N_ptbin){
          if(dMidCh>highMultCut){
            if(iTrPtIdxTPCFMD12.size()!=iTrIndexTPCFMD12.size()) {
              cout<<"Trig pt range idx size not match trig array size!"<<endl;
              exit(1);
            }
            for(unsigned  int iRefEvt=0; iRefEvt<vRefArrayPoolHighTPCFMD12.size(); iRefEvt++) if(vRefArrayPoolHighTPCFMD12.size()==iMixSize) {
              for(unsigned int iAs=0; iAs<vRefArrayPoolHighTPCFMD12.at(iRefEvt).size(); iAs++){
                FillDPhiDEta(vec_trig, vRefArrayPoolHighTPCFMD12.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHigh_ptbinTPCFMD12[iTrPt]);
                if(vRefArrayPoolHighTPCFMD12.at(iRefEvt).at(iAs).Pt()>ptbin_cut[iTrPt]&&vRefArrayPoolHighTPCFMD12.at(iRefEvt).at(iAs).Pt()<ptbin_cut[iTrPt+1]) FillDPhiDEta(vec_trig, vRefArrayPoolHighTPCFMD12.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHigh_ptbinTASameTPCFMD12[iTrPt]);

              }//for asso track in each ref pool evt
            }//for ref pool loop highMult
          }//fi high Mult
          else if(dMidCh<lowMultCut){
            if(iTrPtIdxTPCFMD12.size()!=iTrIndexTPCFMD12.size()) {
              cout<<"Trig pt range idx size not match trig array size!"<<endl;
              exit(1);
            }
            for(unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolLowTPCFMD12.size(); iRefEvt++) if(vRefArrayPoolLowTPCFMD12.size()==iMixSize) {
              for(unsigned int iAs=0; iAs<vRefArrayPoolLowTPCFMD12.at(iRefEvt).size(); iAs++){
                FillDPhiDEta(vec_trig, vRefArrayPoolLowTPCFMD12.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLow_ptbinTPCFMD12[iTrPt]);
                if(vRefArrayPoolLowTPCFMD12.at(iRefEvt).at(iAs).Pt()>ptbin_cut[iTrPt]&&vRefArrayPoolLowTPCFMD12.at(iRefEvt).at(iAs).Pt()<ptbin_cut[iTrPt+1]) FillDPhiDEta(vec_trig, vRefArrayPoolLowTPCFMD12.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLow_ptbinTASameTPCFMD12[iTrPt]);
              }//for asso track in each ref pool evt
            }//for ref pool loop lowMult
          }//fi low Mult
        }
      }

    }//trig for

    //TPCFMD3
    for(unsigned int iTr=0; iTr<iTrIndexTPCFMD3.size(); iTr++ ) {
      int iIdxTrig = iTrIndexTPCFMD3[iTr];
      vec_trig.SetXYZ(ampt->px[iIdxTrig], ampt->py[iIdxTrig], ampt->pz[iIdxTrig] );
      //minbias or centrality event mixing
      if(vec_trig.Pt()>0.2&&vec_trig.Pt()<3&&fabs(vec_trig.Eta())<0.8) {
        for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolTPCFMD3.size(); iRefEvt++) if(vRefArrayPoolTPCFMD3.size()==iMixSize) {
          for(unsigned int iAs=0; iAs<vRefArrayPoolTPCFMD3.at(iRefEvt).size(); iAs++){
            FillDPhiDEta(vec_trig, vRefArrayPoolTPCFMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventTPCFMD3);
          }//for asso track in each ref pool evt
        }//for ref pool loop minbias
  
        if(icent>=0){
          for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPool_centTPCFMD3[icent].size(); iRefEvt++) if(vRefArrayPool_centTPCFMD3[icent].size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPool_centTPCFMD3[icent].at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPool_centTPCFMD3[icent].at(iRefEvt).at(iAs), hDEtaDPhiMixEvent_CentTPCFMD3[icent]);
              FillDPhiDEta3D(vec_trig, vRefArrayPool_centTPCFMD3[icent].at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEvent_CentTPCFMD3[icent]);
            }
          }
        }
  
        if(dMidCh<lowMultCut){
          for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolLowTPCFMD3.size(); iRefEvt++) if(vRefArrayPoolLowTPCFMD3.size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPoolLowTPCFMD3.at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPoolLowTPCFMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLowMidTPCFMD3);
              FillDPhiDEta3D(vec_trig, vRefArrayPoolLowTPCFMD3.at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEventLowMidTPCFMD3);
            }//for asso track in each ref pool evt
          }//for ref pool loop lowMult
        }//fi low Mult
        else if(dMidCh>highMultCut){
          for(unsigned  int iRefEvt=0; iRefEvt<vRefArrayPoolHighTPCFMD3.size(); iRefEvt++) if(vRefArrayPoolHighTPCFMD3.size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPoolHighTPCFMD3.at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPoolHighTPCFMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHighMidTPCFMD3);
              FillDPhiDEta3D(vec_trig, vRefArrayPoolHighTPCFMD3.at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEventHighMidTPCFMD3);
            }//for asso track in each ref pool evt
          }//for ref pool loop highMult
        }//else fi high Mult
      }//fi minbias or centrality trig cut 0~3

      //pt dependence
      if(dMidCh>highMultCut||dMidCh<lowMultCut){
        int iTrPt=iTrPtIdxTPCFMD3.at(iTr);
        if(iTrPt>=0&&iTrPt<N_ptbin){
          if(dMidCh>highMultCut){
            if(iTrPtIdxTPCFMD3.size()!=iTrIndexTPCFMD3.size()) {
              cout<<"Trig pt range idx size not match trig array size!"<<endl;
              exit(1);
            }
            for(unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolHighTPCFMD3.size(); iRefEvt++) if(vRefArrayPoolHighTPCFMD3.size()==iMixSize) {
              for(unsigned int iAs=0; iAs<vRefArrayPoolHighTPCFMD3.at(iRefEvt).size(); iAs++){
                FillDPhiDEta(vec_trig, vRefArrayPoolHighTPCFMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHigh_ptbinTPCFMD3[iTrPt]);
                if(vRefArrayPoolHighTPCFMD3.at(iRefEvt).at(iAs).Pt()>ptbin_cut[iTrPt]&&vRefArrayPoolHighTPCFMD3.at(iRefEvt).at(iAs).Pt()<ptbin_cut[iTrPt+1]) FillDPhiDEta(vec_trig, vRefArrayPoolHighTPCFMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHigh_ptbinTASameTPCFMD3[iTrPt]);

              }//for asso track in each ref pool evt
            }//for ref pool loop highMult
          }//fi high Mult
          else if(dMidCh<lowMultCut){
            if(iTrPtIdxTPCFMD3.size()!=iTrIndexTPCFMD3.size()) {
              cout<<"Trig pt range idx size not match trig array size!"<<endl;
              exit(1);
            }
            for(unsigned  int iRefEvt=0; iRefEvt<vRefArrayPoolLowTPCFMD3.size(); iRefEvt++) if(vRefArrayPoolLowTPCFMD3.size()==iMixSize) {
              for(unsigned int iAs=0; iAs<vRefArrayPoolLowTPCFMD3.at(iRefEvt).size(); iAs++){
                FillDPhiDEta(vec_trig, vRefArrayPoolLowTPCFMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLow_ptbinTPCFMD3[iTrPt]);
                if(vRefArrayPoolLowTPCFMD3.at(iRefEvt).at(iAs).Pt()>ptbin_cut[iTrPt]&&vRefArrayPoolLowTPCFMD3.at(iRefEvt).at(iAs).Pt()<ptbin_cut[iTrPt+1]) FillDPhiDEta(vec_trig, vRefArrayPoolLowTPCFMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLow_ptbinTASameTPCFMD3[iTrPt]);
              }//for asso track in each ref pool evt
            }//for ref pool loop lowMult
          }//fi low Mult
        }
      }

    }//trig for

    //FMD12FMD3
    for(unsigned int iTr=0; iTr<iTrIndexFMD12FMD3.size(); iTr++ ) {
      int iIdxTrig = iTrIndexFMD12FMD3[iTr];
      vec_trig.SetXYZ(ampt->px[iIdxTrig], ampt->py[iIdxTrig], ampt->pz[iIdxTrig] );
      //minbias or centrality event mixing
      if(vec_trig.Eta()>2.6&&vec_trig.Eta()<3.2) {
        for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolFMD12FMD3.size(); iRefEvt++) if(vRefArrayPoolFMD12FMD3.size()==iMixSize) {
          for(unsigned int iAs=0; iAs<vRefArrayPoolFMD12FMD3.at(iRefEvt).size(); iAs++){
            FillDPhiDEta(vec_trig, vRefArrayPoolFMD12FMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventFMD12FMD3);
          }//for asso track in each ref pool evt
        }//for ref pool loop minbias
  
        if(icent>=0){
          for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPool_centFMD12FMD3[icent].size(); iRefEvt++) if(vRefArrayPool_centFMD12FMD3[icent].size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPool_centFMD12FMD3[icent].at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPool_centFMD12FMD3[icent].at(iRefEvt).at(iAs), hDEtaDPhiMixEvent_CentFMD12FMD3[icent]);
              FillDPhiDEta3D(vec_trig, vRefArrayPool_centFMD12FMD3[icent].at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEvent_CentFMD12FMD3[icent]);
            }
          }
        }
  
        if(dMidCh<lowMultCut){
          for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolLowFMD12FMD3.size(); iRefEvt++) if(vRefArrayPoolLowFMD12FMD3.size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPoolLowFMD12FMD3.at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPoolLowFMD12FMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLowMidFMD12FMD3);
              FillDPhiDEta3D(vec_trig, vRefArrayPoolLowFMD12FMD3.at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEventLowMidFMD12FMD3);
            }//for asso track in each ref pool evt
          }//for ref pool loop lowMult
        }//fi low Mult
        else if(dMidCh>highMultCut){
          for(unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolHighFMD12FMD3.size(); iRefEvt++) if(vRefArrayPoolHighFMD12FMD3.size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPoolHighFMD12FMD3.at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPoolHighFMD12FMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHighMidFMD12FMD3);
              FillDPhiDEta3D(vec_trig, vRefArrayPoolHighFMD12FMD3.at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEventHighMidFMD12FMD3);
            }//for asso track in each ref pool evt
          }//for ref pool loop highMult
        }//else fi high Mult
      }//fi minbias or centrality trig cut 0~3

      //pt dependence
      if(dMidCh>highMultCut||dMidCh<lowMultCut){
        int iTrPt=iTrPtIdxFMD12FMD3.at(iTr);
        if(iTrPt>=0&&iTrPt<N_ptbin){
          if(dMidCh>highMultCut){
            if(iTrPtIdxFMD12FMD3.size()!=iTrIndexFMD12FMD3.size()) {
              cout<<"Trig pt range idx size not match trig array size!"<<endl;
              exit(1);
            }
            for(unsigned  int iRefEvt=0; iRefEvt<vRefArrayPoolHighFMD12FMD3.size(); iRefEvt++) if(vRefArrayPoolHighFMD12FMD3.size()==iMixSize) {
              for(unsigned int iAs=0; iAs<vRefArrayPoolHighFMD12FMD3.at(iRefEvt).size(); iAs++){
                FillDPhiDEta(vec_trig, vRefArrayPoolHighFMD12FMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHigh_ptbinFMD12FMD3[iTrPt]);
                if(vRefArrayPoolHighFMD12FMD3.at(iRefEvt).at(iAs).Pt()>ptbin_cut[iTrPt]&&vRefArrayPoolHighFMD12FMD3.at(iRefEvt).at(iAs).Pt()<ptbin_cut[iTrPt+1]) FillDPhiDEta(vec_trig, vRefArrayPoolHighFMD12FMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHigh_ptbinTASameFMD12FMD3[iTrPt]);

              }//for asso track in each ref pool evt
            }//for ref pool loop highMult
          }//fi high Mult
          else if(dMidCh<lowMultCut){
            if(iTrPtIdxFMD12FMD3.size()!=iTrIndexFMD12FMD3.size()) {
              cout<<"Trig pt range idx size not match trig array size!"<<endl;
              exit(1);
            }
            for(unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolLowFMD12FMD3.size(); iRefEvt++) if(vRefArrayPoolLowFMD12FMD3.size()==iMixSize) {
              for(unsigned int iAs=0; iAs<vRefArrayPoolLowFMD12FMD3.at(iRefEvt).size(); iAs++){
                FillDPhiDEta(vec_trig, vRefArrayPoolLowFMD12FMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLow_ptbinFMD12FMD3[iTrPt]);
                if(vRefArrayPoolLowFMD12FMD3.at(iRefEvt).at(iAs).Pt()>ptbin_cut[iTrPt]&&vRefArrayPoolLowFMD12FMD3.at(iRefEvt).at(iAs).Pt()<ptbin_cut[iTrPt+1]) FillDPhiDEta(vec_trig, vRefArrayPoolLowFMD12FMD3.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLow_ptbinTASameFMD12FMD3[iTrPt]);
              }//for asso track in each ref pool evt
            }//for ref pool loop lowMult
          }//fi low Mult
        }
      }

    }//trig for


		//found 1 accept ref list, save to the
		//event mixing pool
		//if mix pool exceeds limit erases the 0th element
    //TPCFMD12
		if(vRefArrayTPCFMD12.size()>=1) {
			vRefArrayPoolTPCFMD12.push_back(vRefArrayTPCFMD12);
		  if(vRefArrayPoolTPCFMD12.size()>iMixSize) vRefArrayPoolTPCFMD12.erase(vRefArrayPoolTPCFMD12.begin());

      //Nsel dependence
      if(icent>=0){
        vRefArrayPool_centTPCFMD12[icent].push_back(vRefArrayTPCFMD12);
        if(vRefArrayPool_centTPCFMD12[icent].size()>iMixSize) vRefArrayPool_centTPCFMD12[icent].erase(vRefArrayPool_centTPCFMD12[icent].begin());
      }

			if(dMidCh<lowMultCut) {
				vRefArrayPoolLowTPCFMD12.push_back(vRefArrayTPCFMD12);
        if(vRefArrayPoolLowTPCFMD12.size()>iMixSize) vRefArrayPoolLowTPCFMD12.erase(vRefArrayPoolLowTPCFMD12.begin());
      }
			else if(dMidCh>highMultCut){
				vRefArrayPoolHighTPCFMD12.push_back(vRefArrayTPCFMD12);
        if(vRefArrayPoolHighTPCFMD12.size()>iMixSize) vRefArrayPoolHighTPCFMD12.erase(vRefArrayPoolHighTPCFMD12.begin());
      }
		}

    //TPCFMD3
    if(vRefArrayTPCFMD3.size()>=1) {
			vRefArrayPoolTPCFMD3.push_back(vRefArrayTPCFMD3);
		  if(vRefArrayPoolTPCFMD3.size()>iMixSize) vRefArrayPoolTPCFMD3.erase(vRefArrayPoolTPCFMD3.begin());

      //Nsel dependence
      if(icent>=0){
        vRefArrayPool_centTPCFMD3[icent].push_back(vRefArrayTPCFMD3);
        if(vRefArrayPool_centTPCFMD3[icent].size()>iMixSize) vRefArrayPool_centTPCFMD3[icent].erase(vRefArrayPool_centTPCFMD3[icent].begin());
      }

			if(dMidCh<lowMultCut) {
				vRefArrayPoolLowTPCFMD3.push_back(vRefArrayTPCFMD3);
        if(vRefArrayPoolLowTPCFMD3.size()>iMixSize) vRefArrayPoolLowTPCFMD3.erase(vRefArrayPoolLowTPCFMD3.begin());
      }
			else if(dMidCh>highMultCut){
				vRefArrayPoolHighTPCFMD3.push_back(vRefArrayTPCFMD3);
        if(vRefArrayPoolHighTPCFMD3.size()>iMixSize) vRefArrayPoolHighTPCFMD3.erase(vRefArrayPoolHighTPCFMD3.begin());
      }
		}

    //FMD12FMD3
    if(vRefArrayFMD12FMD3.size()>=1) {
			vRefArrayPoolFMD12FMD3.push_back(vRefArrayFMD12FMD3);
		  if(vRefArrayPoolFMD12FMD3.size()>iMixSize) vRefArrayPoolFMD12FMD3.erase(vRefArrayPoolFMD12FMD3.begin());

      //Nsel dependence
      if(icent>=0){
        vRefArrayPool_centFMD12FMD3[icent].push_back(vRefArrayFMD12FMD3);
        if(vRefArrayPool_centFMD12FMD3[icent].size()>iMixSize) vRefArrayPool_centFMD12FMD3[icent].erase(vRefArrayPool_centFMD12FMD3[icent].begin());
      }

			if(dMidCh<lowMultCut) {
				vRefArrayPoolLowFMD12FMD3.push_back(vRefArrayFMD12FMD3);
        if(vRefArrayPoolLowFMD12FMD3.size()>iMixSize) vRefArrayPoolLowFMD12FMD3.erase(vRefArrayPoolLowFMD12FMD3.begin());
      }
			else if(dMidCh>highMultCut){
				vRefArrayPoolHighFMD12FMD3.push_back(vRefArrayFMD12FMD3);
        if(vRefArrayPoolHighFMD12FMD3.size()>iMixSize) vRefArrayPoolHighFMD12FMD3.erase(vRefArrayPoolHighFMD12FMD3.begin());
      }
		}
 
  }//for event loop


  //--------------------------------
  // Write Hists and exit                               
  //--------------------------------

  TFile* f = new TFile(outputFile, "RECREATE");
  f->cd();

	hNChMid->Write();

	hDEtaDPhiSameEventTPCFMD12->Write();
  hDEtaDPhiSameEventTPCFMD3->Write();
  hDEtaDPhiSameEventFMD12FMD3->Write();
	hDEtaDPhiMixEventTPCFMD12->Write();
  hDEtaDPhiMixEventTPCFMD3->Write();
  hDEtaDPhiMixEventFMD12FMD3->Write();
	hDEtaDPhiSameEventLowMidTPCFMD12->Write();
  hDEtaDPhiSameEventLowMidTPCFMD3->Write();
  hDEtaDPhiSameEventLowMidFMD12FMD3->Write();
	hDEtaDPhiMixEventLowMidTPCFMD12->Write();
  hDEtaDPhiMixEventLowMidTPCFMD3->Write();
  hDEtaDPhiMixEventLowMidFMD12FMD3->Write();
	hDEtaDPhiSameEventHighMidTPCFMD12->Write();
  hDEtaDPhiSameEventHighMidTPCFMD3->Write();
  hDEtaDPhiSameEventHighMidFMD12FMD3->Write();
	hDEtaDPhiMixEventHighMidTPCFMD12->Write();
  hDEtaDPhiMixEventHighMidTPCFMD3->Write();
  hDEtaDPhiMixEventHighMidFMD12FMD3->Write();
	hDEtaDPhiTrigEtaSameEventLowMidTPCFMD12->Write();
  hDEtaDPhiTrigEtaSameEventLowMidTPCFMD3->Write();
  hDEtaDPhiTrigEtaSameEventLowMidFMD12FMD3->Write();
	hDEtaDPhiTrigEtaMixEventLowMidTPCFMD12->Write();
  hDEtaDPhiTrigEtaMixEventLowMidTPCFMD3->Write();
	hDEtaDPhiTrigEtaMixEventLowMidFMD12FMD3->Write();
	hDEtaDPhiTrigEtaSameEventHighMidTPCFMD12->Write();
  hDEtaDPhiTrigEtaSameEventHighMidTPCFMD3->Write();
  hDEtaDPhiTrigEtaSameEventHighMidFMD12FMD3->Write();
	hDEtaDPhiTrigEtaMixEventHighMidTPCFMD12->Write();
  hDEtaDPhiTrigEtaMixEventHighMidTPCFMD3->Write();
  hDEtaDPhiTrigEtaMixEventHighMidFMD12FMD3->Write();

  hTrigPtTPCFMD12->Write();
  hTrigPtLowTPCFMD12->Write();
  hTrigPtHighTPCFMD12->Write();
  hTrigPtEtaLowTPCFMD12->Write();
  hTrigPtEtaHighTPCFMD12->Write();
  hTrigPtTPCFMD3->Write();
  hTrigPtLowTPCFMD3->Write();
  hTrigPtHighTPCFMD3->Write();
  hTrigPtEtaLowTPCFMD3->Write();
  hTrigPtEtaHighTPCFMD3->Write();
  hTrigPtFMD12FMD3->Write();
  hTrigPtLowFMD12FMD3->Write();
  hTrigPtHighFMD12FMD3->Write();
  hTrigPtEtaLowFMD12FMD3->Write();
  hTrigPtEtaHighFMD12FMD3->Write();
//	hSphericity->Write();

  for(int i=0; i<N_cent; i++) {
    hDEtaDPhiSameEvent_CentTPCFMD12[i]->Write();
    hDEtaDPhiSameEvent_CentTPCFMD3[i]->Write();
    hDEtaDPhiSameEvent_CentFMD12FMD3[i]->Write();
    hDEtaDPhiMixEvent_CentTPCFMD12[i]->Write();
    hDEtaDPhiMixEvent_CentTPCFMD3[i]->Write();
    hDEtaDPhiMixEvent_CentFMD12FMD3[i]->Write();
    hTrigPt_CentTPCFMD12[i]->Write();
    hTrigPt_CentTPCFMD3[i]->Write();
    hTrigPt_CentFMD12FMD3[i]->Write();
    hDEtaDPhiTrigEtaSameEvent_CentTPCFMD12[i]->Write();
    hDEtaDPhiTrigEtaSameEvent_CentTPCFMD3[i]->Write();
    hDEtaDPhiTrigEtaSameEvent_CentFMD12FMD3[i]->Write();
    hDEtaDPhiTrigEtaMixEvent_CentTPCFMD12[i]->Write();
    hDEtaDPhiTrigEtaMixEvent_CentTPCFMD3[i]->Write();
    hDEtaDPhiTrigEtaMixEvent_CentFMD12FMD3[i]->Write();
    hTrigPtEta_CentTPCFMD12[i]->Write();
    hTrigPtEta_CentTPCFMD3[i]->Write();
    hTrigPtEta_CentFMD12FMD3[i]->Write();
  }

  for(int i=0; i<N_ptbin; i++) {
    hDEtaDPhiSameEventHigh_ptbinTPCFMD12[i]->Write();
    hDEtaDPhiMixEventHigh_ptbinTPCFMD12[i]->Write();
    hTrigPtHigh_ptbinTPCFMD12[i]->Write();
    hDEtaDPhiSameEventLow_ptbinTPCFMD12[i]->Write();
    hDEtaDPhiMixEventLow_ptbinTPCFMD12[i]->Write();
    hTrigPtLow_ptbinTPCFMD12[i]->Write();
    hDEtaDPhiSameEventHigh_ptbinTASameTPCFMD12[i]->Write();
    hDEtaDPhiMixEventHigh_ptbinTASameTPCFMD12[i]->Write();
    hDEtaDPhiSameEventLow_ptbinTASameTPCFMD12[i]->Write();
    hDEtaDPhiMixEventLow_ptbinTASameTPCFMD12[i]->Write();

    hDEtaDPhiSameEventHigh_ptbinTPCFMD3[i]->Write();
    hDEtaDPhiMixEventHigh_ptbinTPCFMD3[i]->Write();
    hTrigPtHigh_ptbinTPCFMD3[i]->Write();
    hDEtaDPhiSameEventLow_ptbinTPCFMD3[i]->Write();
    hDEtaDPhiMixEventLow_ptbinTPCFMD3[i]->Write();
    hTrigPtLow_ptbinTPCFMD3[i]->Write();
    hDEtaDPhiSameEventHigh_ptbinTASameTPCFMD3[i]->Write();
    hDEtaDPhiMixEventHigh_ptbinTASameTPCFMD3[i]->Write();
    hDEtaDPhiSameEventLow_ptbinTASameTPCFMD3[i]->Write();
    hDEtaDPhiMixEventLow_ptbinTASameTPCFMD3[i]->Write();

    hDEtaDPhiSameEventHigh_ptbinFMD12FMD3[i]->Write();
    hDEtaDPhiMixEventHigh_ptbinFMD12FMD3[i]->Write();
    hTrigPtHigh_ptbinFMD12FMD3[i]->Write();
    hDEtaDPhiSameEventLow_ptbinFMD12FMD3[i]->Write();
    hDEtaDPhiMixEventLow_ptbinFMD12FMD3[i]->Write();
    hTrigPtLow_ptbinFMD12FMD3[i]->Write();
    hDEtaDPhiSameEventHigh_ptbinTASameFMD12FMD3[i]->Write();
    hDEtaDPhiMixEventHigh_ptbinTASameFMD12FMD3[i]->Write();
    hDEtaDPhiSameEventLow_ptbinTASameFMD12FMD3[i]->Write();
    hDEtaDPhiMixEventLow_ptbinTASameFMD12FMD3[i]->Write();
  }


  f->Close();

  delete chain;
  return 0;
}

//--------------------------------
// Define functions here                               
//--------------------------------

int GetCentrality(AMPT* ampt) // different with StRefMultCorr
{
  float b = ampt->Event_b;

  int cent=-1;
  for (int i=1 ; i<=10 ; ++i) {
    if (b<=sqrt(10*i/100.0)*2.0*pow(197, 1.0/3.0)*1.2) { // ...*1.124
      cent = i-1;
      return cent;
    }
  }

  return cent;
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

float GetSphericity(vector<TVector3> tracks)
{
  // a critical check
  if(tracks.size()==0) return 0;
  
  // first fill the momentum tensor
  TMatrixDSym MomentumTensor(3);
  for(vector<TVector3>::iterator itTrack = tracks.begin(); itTrack!=tracks.end(); ++itTrack) {
  std::vector<double> momentum(3);
  momentum[0] = itTrack->X();
  momentum[1] = itTrack->Y();
  momentum[2] = itTrack->Z();
    for(unsigned int i=0;i<3;i++)
      for(unsigned int j=0;j<=i;j++) {
        MomentumTensor[i][j] += momentum[i]*momentum[j];
      }
  }
  MomentumTensor*=1/(MomentumTensor[0][0]+MomentumTensor[1][1]+MomentumTensor[2][2]);
  // find the eigen values
  TMatrixDSymEigen eigen(MomentumTensor);
  TVectorD eigenvals = eigen.GetEigenValues();
  vector<float> eigenvaluess(3);
  eigenvaluess[0] = eigenvals[0];
  eigenvaluess[1] = eigenvals[1];
  eigenvaluess[2] = eigenvals[2];
  sort(eigenvaluess.begin(),eigenvaluess.end());
  // compute spericity
  float sph = ( 1.5*(1-eigenvaluess[2]));
  return sph;
}

void FillDPhiDEta(TVector3 vTr, TVector3 vAs, TH2D *hout){
  double dDPhi = vTr.DeltaPhi(vAs);
  const double PI = 3.1415926;
  if(dDPhi<-0.5*PI) dDPhi=dDPhi + 2*PI;
  if(dDPhi>1.5*PI) dDPhi=dDPhi - 2*PI;
  double dDEta = vTr.Eta() - vAs.Eta();
  hout->Fill(dDPhi, dDEta);
}

void FillDPhiDEta3D(TVector3 vTr, TVector3 vAs, TH3D *hout){
  double dDPhi = vTr.DeltaPhi(vAs);
  const double PI = 3.1415926;
  if(dDPhi<-0.5*PI) dDPhi=dDPhi + 2*PI;
  if(dDPhi>1.5*PI) dDPhi=dDPhi - 2*PI;
  double dDEta = vTr.Eta() - vAs.Eta();
  hout->Fill(dDPhi, dDEta, vTr.Eta());
}


void PrintArray(double array[], int size){
  for(int i=0; i<size; i++){
    cout<<array[i]<<" ";
  }
  cout<<endl;
}

