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

	TH2D *hDEtaDPhiSameEvent = new TH2D("hDEtaDPhiSameEvent", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiMixEvent = new TH2D("hDEtaDPhiMixEvent", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

	TH2D *hDEtaDPhiSameEventLowMid = new TH2D("hDEtaDPhiSameEventLowMid", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiMixEventLowMid = new TH2D("hDEtaDPhiMixEventLowMid", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

	TH2D *hDEtaDPhiSameEventHighMid = new TH2D("hDEtaDPhiSameEventHighMid", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
	TH2D *hDEtaDPhiMixEventHighMid = new TH2D("hDEtaDPhiMixEventHighMid", "; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

	TH3D *hDEtaDPhiTrigEtaSameEventLowMid = new TH3D("hDEtaDPhiTrigEtaSameEventLowMid", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaMixEventLowMid = new TH3D("hDEtaDPhiTrigEtaMixEventLowMid", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);

	TH3D *hDEtaDPhiTrigEtaSameEventHighMid = new TH3D("hDEtaDPhiTrigEtaSameEventHighMid", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
	TH3D *hDEtaDPhiTrigEtaMixEventHighMid = new TH3D("hDEtaDPhiTrigEtaMixEventHighMid", "; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);


  TH1D *hTrigPt = new TH1D("hTrigPt", ";p_{T}", 50, 0, 10);
  TH1D *hTrigPtLow = new TH1D("hTrigPtLow", ";p_{T}", 50, 0, 10);
  TH1D *hTrigPtHigh = new TH1D("hTrigPtHigh", ";p_{T}", 50, 0, 10);
  TH2D *hTrigPtEtaLow = new TH2D("hTrigPtEtaLow", ";p_{T};#eta_{tr}", 50, 0, 10, 60, -3, 3);
  TH2D *hTrigPtEtaHigh = new TH2D("hTrigPtEtaHigh", ";p_{T};#eta_{tr}", 50, 0, 10, 60, -3, 3);

  //set 9 event classes, 10 Nsel bin edge
  //include 0-10 bin for consistency check
  const int N_cent=9;
  double Nsel_min=0;
  double Nsel_bin=10;
  //double Nsel_cut[N_cent+1]={0};
  double Nsel_cut[N_cent+1]={10, 30, 40, 50, 60, 80, 100, 120, 150, 200};
  TH2D *hDEtaDPhiSameEvent_Cent[N_cent];
  TH2D *hDEtaDPhiMixEvent_Cent[N_cent];
  TH1D *hTrigPt_Cent[N_cent];

  TH3D *hDEtaDPhiTrigEtaSameEvent_Cent[N_cent];
  TH3D *hDEtaDPhiTrigEtaMixEvent_Cent[N_cent];
  TH2D *hTrigPtEta_Cent[N_cent];

  for(int i=0; i<N_cent; i++) {
    //Nsel_cut[i]=Nsel_min+i*Nsel_bin;
    //if(i==N_cent-1) Nsel_cut[N_cent]=Nsel_cut[i]+Nsel_bin;
    hDEtaDPhiSameEvent_Cent[i] = new TH2D(Form("hDEtaDPhiSameEvent_Cent%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEvent_Cent[i] = new TH2D(Form("hDEtaDPhiMixEvent_Cent%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

    hDEtaDPhiTrigEtaSameEvent_Cent[i] = new TH3D(Form("hDEtaDPhiTrigEtaSameEvent_Cent%d",i),"; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);
    hDEtaDPhiTrigEtaMixEvent_Cent[i] = new TH3D(Form("hDEtaDPhiTrigEtaMixEvent_Cent%d",i),"; #Delta#phi; #Delta#eta; #eta_{tr}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 60, -3, 3);

    hTrigPt_Cent[i]= new TH1D(Form("hTrigPt_Cent%d",i), ";p_{T}", 50, 0, 10);
    hTrigPtEta_Cent[i]= new TH2D(Form("hTrigPtEta_Cent%d",i), ";p_{T}", 50, 0, 10, 60, -3, 3);
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
  TH2D *hDEtaDPhiSameEventHigh_ptbin[N_ptbin];
  TH2D *hDEtaDPhiMixEventHigh_ptbin[N_ptbin];
  TH2D *hDEtaDPhiSameEventLow_ptbin[N_ptbin];
  TH2D *hDEtaDPhiMixEventLow_ptbin[N_ptbin];
  TH1D *hTrigPtHigh_ptbin[N_ptbin];
  TH1D *hTrigPtLow_ptbin[N_ptbin];

  //keep trig/asso in the same pt range
  TH2D *hDEtaDPhiSameEventHigh_ptbinTASame[N_ptbin];
  TH2D *hDEtaDPhiMixEventHigh_ptbinTASame[N_ptbin];
  TH2D *hDEtaDPhiSameEventLow_ptbinTASame[N_ptbin];
  TH2D *hDEtaDPhiMixEventLow_ptbinTASame[N_ptbin];


  for(int i=0; i<N_ptbin; i++){
//    ptbin_cut[i]=pt_min+i*pt_bin;
//    if(i>4) ptbin_cut[i]=ptbin_cut[i-1]+(2*pt_bin);
//    if(i==N_ptbin-1) ptbin_cut[N_ptbin]=ptbin_cut[i]+2*pt_bin;
    hDEtaDPhiSameEventHigh_ptbin[i] = new TH2D(Form("hDEtaDPhiSameEventHigh_ptbin%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventHigh_ptbin[i] = new TH2D(Form("hDEtaDPhiMixEventHigh_ptbin%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hTrigPtHigh_ptbin[i]= new TH1D(Form("hTrigPtHigh_ptbin%d",i), ";p_{T}", 500, 0, 10);
    hDEtaDPhiSameEventLow_ptbin[i] = new TH2D(Form("hDEtaDPhiSameEventLow_ptbin%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventLow_ptbin[i] = new TH2D(Form("hDEtaDPhiMixEventLow_ptbin%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hTrigPtLow_ptbin[i]= new TH1D(Form("hTrigPtLow_ptbin%d",i), ";p_{T}", 500, 0, 10);

    hDEtaDPhiSameEventHigh_ptbinTASame[i] = new TH2D(Form("hDEtaDPhiSameEventHigh_ptbinTASame%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventHigh_ptbinTASame[i] = new TH2D(Form("hDEtaDPhiMixEventHigh_ptbinTASame%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiSameEventLow_ptbinTASame[i] = new TH2D(Form("hDEtaDPhiSameEventLow_ptbinTASame%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);
    hDEtaDPhiMixEventLow_ptbinTASame[i] = new TH2D(Form("hDEtaDPhiMixEventLow_ptbinTASame%d",i),"; #Delta#phi; #Delta#eta", 40, -0.5*PI, 1.5*PI, 200, -12, 12);

  }
  cout<<"pt bin edge:"<<endl;
  PrintArray(ptbin_cut, N_ptbin+1);


	TH3D *hNchDEtaDPhiSameEvent = new TH3D("hNchDEtaDPhiSameEvent", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 300, -0.5, 299.5);
	TH3D *hNchDEtaDPhiMixEvent = new TH3D("hNchDEtaDPhiMixEvent", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*PI, 1.5*PI, 200, -12, 12, 300, -0.5, 299.5);

	TH2D *hSphericity = new TH2D("hSphericity", ";N_{ch}^{mid}; S", 2000, 0, 4000, 50, 0, 1);

   
  //--------------------------------
  // loop events
  //--------------------------------

  int nEvents = (int)chain->GetEntries();
  cout<<"Total events : "<<nEvents<<endl;

//  if(nEvents>8e6) nEvents=8e6;
  AMPT* ampt = new AMPT(chain);

	//cuts
  double accept_cut = 2.5;
  //double accept_cut = 6;
	double highMultCut = 60;
	double lowMultCut = 30;

	
 	//temperary containers
	TVector3 vec_tmp;
	TVector3 vec_trig;
	vector<int> iTrIndex;
	vector<int> iTrPtIdx;
	vector<int> iAsIndex;
	vector<TVector3> vRefArray;
	const int iMixSize=5;
	vector< vector<TVector3> > vRefArrayPool;
	vector< vector<TVector3> > vRefArrayPool_cent[10];
	vector< vector<TVector3> > vRefArrayPoolLow;
	vector< vector<TVector3> > vRefArrayPoolHigh;
	vector<TVector3> sph_container;

  double NtrigHigh=0;
  double NtrigLow=0;
 
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

		iTrIndex.clear();
		iTrPtIdx.clear();
		iAsIndex.clear();
		vRefArray.clear();
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

      if(vec_tmp.Eta()>3&&vec_tmp.Eta()<5&&energy>3) foundEast=true;
      if(vec_tmp.Eta()>-5&&vec_tmp.Eta()<-3&&energy>3) foundWest=true;
			

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

      //reference particle, using forward FCal in ATLAS
			if(pt>0.5&&pt<5&&vec_tmp.Eta()>4&&vec_tmp.Eta()<4.9) {
        iAsIndex.push_back(iTrack);
        //randomize ref particle container for mix event
        vRefArray.push_back(vec_tmp);
      }


      //mid-rapidity acceptance cut
			if (fabs(vec_tmp.Eta()) > accept_cut) continue;

      if(pass_id&&pt>0.4) dMidCh++;//Nsel
      //if(pass_id&&pt>0.2&&fabs(vec_tmp.Eta())<0.8) dMidCh++;//Nch ALICE

      /*
      //reference particle, fixed for charge between 0.3~3
			//if(GetCharge(id)!=0&&pt>0.3&&pt<3) {
			if(pass_id&&pt>0.3&&pt<5) {
        sph_container.push_back(vec_tmp);
        iAsIndex.push_back(iTrack);
        //randomize ref particle container for mix event
        vRefArray.push_back(vec_tmp);
      }
      */

      //get trig particle index (0.3~6) and do selection in second loop
      //if need to include other Pid trig hadron, remove charge cut here
      //and do Pid cut later in the second loop
			//if(GetCharge(id)!=0&&pt>0.3&&pt<6) {
			//if(pass_id&&pt>0.001&&pt<6) {
			if(pass_trig_id&&pt>0.3&&pt<5) {
				iTrIndex.push_back(iTrack);
      }

    }//for first loop

    //must found at least 1 track with E>3 at each side of detector
    if(!foundEast||!foundWest) continue;
			
		//sphericity analysis
//		double dSph = GetSphericity(sph_container);
//		if(dMidCh>0) 
//			hSphericity->Fill(dMidCh, dSph);

    //ATLAS nsd events
    if(dMidCh<10) continue;

    if(dMidCh>highMultCut) NtrigHigh+=iTrIndex.size();
    if(dMidCh<lowMultCut) NtrigLow+=iTrIndex.size();

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
		for (unsigned int iTr=0; iTr<iTrIndex.size(); iTr++ ) {
      int iIdxTrig = iTrIndex[iTr];
      vec_trig.SetXYZ(ampt->px[iIdxTrig], ampt->py[iIdxTrig], ampt->pz[iIdxTrig] );

      //fill trig info
			if(vec_trig.Pt()>0.3&&vec_trig.Pt()<5) {
        hTrigPt->Fill(vec_trig.Pt());
        //Nsel dependence
        if(icent>=0) {
          hTrigPt_Cent[icent]->Fill(vec_trig.Pt());
          hTrigPtEta_Cent[icent]->Fill(vec_trig.Pt(), vec_trig.Eta());
        }

        if(dMidCh<lowMultCut){
          hTrigPtLow->Fill(vec_trig.Pt());
          hTrigPtEtaLow->Fill(vec_trig.Pt(), vec_trig.Eta());
        }
        else if(dMidCh>highMultCut){
          hTrigPtHigh->Fill(vec_trig.Pt());
          hTrigPtEtaHigh->Fill(vec_trig.Pt(), vec_trig.Eta());
        }
      }

      //trig info pt dependence
      int findPtBin=-1;
      if(dMidCh<lowMultCut||dMidCh>highMultCut){
        for(int iptbin=0; iptbin<N_ptbin; iptbin++){
          if(vec_trig.Pt()>ptbin_cut[iptbin]&&vec_trig.Pt()<ptbin_cut[iptbin+1]){
            findPtBin=iptbin;
            if(dMidCh<lowMultCut)
              hTrigPtLow_ptbin[iptbin]->Fill(vec_trig.Pt());
            else if(dMidCh>highMultCut)
              hTrigPtHigh_ptbin[iptbin]->Fill(vec_trig.Pt());
            break;
          }
        }
        iTrPtIdx.push_back(findPtBin);
      }
      if(findPtBin>=N_ptbin||(findPtBin<0&&findPtBin!=-1)) {
        cout<<"findPtBin="<<findPtBin<<" Error!"<<endl;
        exit(1);
      }

			for (unsigned int iAs=0; iAs<iAsIndex.size(); iAs++ ) {
				if( iAsIndex[iAs]!=iTrIndex[iTr] ) {
					int iIdxAs = iAsIndex[iAs];
					vec_tmp.SetXYZ(ampt->px[iIdxAs], ampt->py[iIdxAs], ampt->pz[iIdxAs] );

					double dDPhi = vec_tmp.DeltaPhi(vec_trig);
					if(dDPhi<-0.5*PI) dDPhi=dDPhi + 2*PI;
					if(dDPhi>1.5*PI) dDPhi=dDPhi - 2*PI;
					double dDEta = vec_tmp.Eta() - vec_trig.Eta();

          //minbias or centrality dependence
					if(vec_trig.Pt()>0.3&&vec_trig.Pt()<5) {
            FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEvent);

            //Nsel dependence
            if(icent>=0) {
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEvent_Cent[icent]);
              FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEvent_Cent[icent]);
            }

						if(dMidCh<lowMultCut){
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLowMid);
              FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEventLowMid);
            }
						else if(dMidCh>highMultCut){
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHighMid);
              FillDPhiDEta3D(vec_trig, vec_tmp, hDEtaDPhiTrigEtaSameEventHighMid);
            }
					}

          //pt dependence
          if(findPtBin>=0&&findPtBin<N_ptbin){
            if(dMidCh<lowMultCut){
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLow_ptbin[findPtBin]);
              if(vec_tmp.Pt()>ptbin_cut[findPtBin]&&vec_tmp.Pt()<ptbin_cut[findPtBin+1]) FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventLow_ptbinTASame[findPtBin]);
            }
            else if(dMidCh>highMultCut) {
              FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHigh_ptbin[findPtBin]);
              if(vec_tmp.Pt()>ptbin_cut[findPtBin]&&vec_tmp.Pt()<ptbin_cut[findPtBin+1]) FillDPhiDEta(vec_trig, vec_tmp, hDEtaDPhiSameEventHigh_ptbinTASame[findPtBin]);
            }
          }

				}//fi asso indx different from trig

			}//asso for
		}//trig for

		//do correlation in the mix event
    //loop current found trigger particles
		for(unsigned int iTr=0; iTr<iTrIndex.size(); iTr++ ) {
      int iIdxTrig = iTrIndex[iTr];
      vec_trig.SetXYZ(ampt->px[iIdxTrig], ampt->py[iIdxTrig], ampt->pz[iIdxTrig] );
      //minbias or centrality event mixing
      if(vec_trig.Pt()>0.3&&vec_trig.Pt()<5) {
        for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPool.size(); iRefEvt++) if(vRefArrayPool.size()==iMixSize) {
          for(unsigned int iAs=0; iAs<vRefArrayPool.at(iRefEvt).size(); iAs++){
            FillDPhiDEta(vec_trig, vRefArrayPool.at(iRefEvt).at(iAs), hDEtaDPhiMixEvent);
          }//for asso track in each ref pool evt
        }//for ref pool loop minbias

        if(icent>=0){
          for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPool_cent[icent].size(); iRefEvt++) if(vRefArrayPool_cent[icent].size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPool_cent[icent].at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPool_cent[icent].at(iRefEvt).at(iAs), hDEtaDPhiMixEvent_Cent[icent]);
              FillDPhiDEta3D(vec_trig, vRefArrayPool_cent[icent].at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEvent_Cent[icent]);
            }
          }
        }

        if(dMidCh<lowMultCut){
          for( unsigned int iRefEvt=0; iRefEvt<vRefArrayPoolLow.size(); iRefEvt++) if(vRefArrayPoolLow.size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPoolLow.at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPoolLow.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLowMid);
              FillDPhiDEta3D(vec_trig, vRefArrayPoolLow.at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEventLowMid);
            }//for asso track in each ref pool evt
          }//for ref pool loop lowMult
        }//fi low Mult
        else if(dMidCh>highMultCut){
          for(unsigned  int iRefEvt=0; iRefEvt<vRefArrayPoolHigh.size(); iRefEvt++) if(vRefArrayPoolHigh.size()==iMixSize) {
            for(unsigned int iAs=0; iAs<vRefArrayPoolHigh.at(iRefEvt).size(); iAs++){
              FillDPhiDEta(vec_trig, vRefArrayPoolHigh.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHighMid);
              FillDPhiDEta3D(vec_trig, vRefArrayPoolHigh.at(iRefEvt).at(iAs), hDEtaDPhiTrigEtaMixEventHighMid);
            }//for asso track in each ref pool evt
          }//for ref pool loop highMult
        }//else fi high Mult
      }//fi minbias or centrality trig cut 0~3

      //pt dependence
      if(dMidCh>highMultCut||dMidCh<lowMultCut){
        int iTrPt=iTrPtIdx.at(iTr);
        if(iTrPt>=0&&iTrPt<N_ptbin){
          if(dMidCh>highMultCut){
            if(iTrPtIdx.size()!=iTrIndex.size()) {
              cout<<"Trig pt range idx size not match trig array size!"<<endl;
              exit(1);
            }
            for(unsigned  int iRefEvt=0; iRefEvt<vRefArrayPoolHigh.size(); iRefEvt++) if(vRefArrayPoolHigh.size()==iMixSize) {
              for(unsigned int iAs=0; iAs<vRefArrayPoolHigh.at(iRefEvt).size(); iAs++){
                FillDPhiDEta(vec_trig, vRefArrayPoolHigh.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHigh_ptbin[iTrPt]);
                if(vRefArrayPoolHigh.at(iRefEvt).at(iAs).Pt()>ptbin_cut[iTrPt]&&vRefArrayPoolHigh.at(iRefEvt).at(iAs).Pt()<ptbin_cut[iTrPt+1]) FillDPhiDEta(vec_trig, vRefArrayPoolHigh.at(iRefEvt).at(iAs), hDEtaDPhiMixEventHigh_ptbinTASame[iTrPt]);

              }//for asso track in each ref pool evt
            }//for ref pool loop highMult
          }//fi high Mult
          else if(dMidCh<lowMultCut){
            if(iTrPtIdx.size()!=iTrIndex.size()) {
              cout<<"Trig pt range idx size not match trig array size!"<<endl;
              exit(1);
            }
            for(unsigned  int iRefEvt=0; iRefEvt<vRefArrayPoolLow.size(); iRefEvt++) if(vRefArrayPoolLow.size()==iMixSize) {
              for(unsigned int iAs=0; iAs<vRefArrayPoolLow.at(iRefEvt).size(); iAs++){
                FillDPhiDEta(vec_trig, vRefArrayPoolLow.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLow_ptbin[iTrPt]);
                if(vRefArrayPoolLow.at(iRefEvt).at(iAs).Pt()>ptbin_cut[iTrPt]&&vRefArrayPoolLow.at(iRefEvt).at(iAs).Pt()<ptbin_cut[iTrPt+1]) FillDPhiDEta(vec_trig, vRefArrayPoolLow.at(iRefEvt).at(iAs), hDEtaDPhiMixEventLow_ptbinTASame[iTrPt]);
              }//for asso track in each ref pool evt
            }//for ref pool loop lowMult
          }//fi low Mult
        }
      }

    }//trig for

		//found 1 accept ref list, save to the
		//event mixing pool
		//if mix pool exceeds limit erases the 0th element
		if(vRefArray.size()>=1) {
			vRefArrayPool.push_back(vRefArray);
		  if(vRefArrayPool.size()>iMixSize) vRefArrayPool.erase(vRefArrayPool.begin());

      //Nsel dependence
      if(icent>=0){
        vRefArrayPool_cent[icent].push_back(vRefArray);
        if(vRefArrayPool_cent[icent].size()>iMixSize) vRefArrayPool_cent[icent].erase(vRefArrayPool_cent[icent].begin());
      }

			if(dMidCh<lowMultCut) {
				vRefArrayPoolLow.push_back(vRefArray);
        if(vRefArrayPoolLow.size()>iMixSize) vRefArrayPoolLow.erase(vRefArrayPoolLow.begin());
      }
			else if(dMidCh>highMultCut){
				vRefArrayPoolHigh.push_back(vRefArray);
        if(vRefArrayPoolHigh.size()>iMixSize) vRefArrayPoolHigh.erase(vRefArrayPoolHigh.begin());
      }
		}
 
  }//for event loop


  //--------------------------------
  // Write Hists and exit                               
  //--------------------------------

  TFile* f = new TFile(outputFile, "RECREATE");
  f->cd();

	hNChMid->Write();

	hDEtaDPhiSameEvent->Write();
	hDEtaDPhiMixEvent->Write();
	hDEtaDPhiSameEventLowMid->Write();
	hDEtaDPhiMixEventLowMid->Write();
	hDEtaDPhiSameEventHighMid->Write();
	hDEtaDPhiMixEventHighMid->Write();
	hDEtaDPhiTrigEtaSameEventLowMid->Write();
	hDEtaDPhiTrigEtaMixEventLowMid->Write();
	hDEtaDPhiTrigEtaSameEventHighMid->Write();
	hDEtaDPhiTrigEtaMixEventHighMid->Write();

  hTrigPt->Write();
  hTrigPtLow->Write();
  hTrigPtHigh->Write();
  hTrigPtEtaLow->Write();
  hTrigPtEtaHigh->Write();


//	hSphericity->Write();

  for(int i=0; i<N_cent; i++) {
    hDEtaDPhiSameEvent_Cent[i]->Write();
    hDEtaDPhiMixEvent_Cent[i]->Write();
    hTrigPt_Cent[i]->Write();

    hDEtaDPhiTrigEtaSameEvent_Cent[i]->Write();
    hDEtaDPhiTrigEtaMixEvent_Cent[i]->Write();
    hTrigPtEta_Cent[i]->Write();
  }

  for(int i=0; i<N_ptbin; i++) {
    hDEtaDPhiSameEventHigh_ptbin[i]->Write();
    hDEtaDPhiMixEventHigh_ptbin[i]->Write();
    hTrigPtHigh_ptbin[i]->Write();
    hDEtaDPhiSameEventLow_ptbin[i]->Write();
    hDEtaDPhiMixEventLow_ptbin[i]->Write();
    hTrigPtLow_ptbin[i]->Write();

    hDEtaDPhiSameEventHigh_ptbinTASame[i]->Write();
    hDEtaDPhiMixEventHigh_ptbinTASame[i]->Write();
    hDEtaDPhiSameEventLow_ptbinTASame[i]->Write();
    hDEtaDPhiMixEventLow_ptbinTASame[i]->Write();
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

