/*************************************************************************
    > File Name: main.C
    > Author: Liang Zheng
    > Mail: liangzhphy@gmail.com 
    > Created Time: Sun Sep 27 10:51:49 2020
 ************************************************************************/

#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"

#include "Pythia8/Pythia.h"

// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;

//event level tree info
typedef struct {
  int   multiplicity;		// multiplicity
  float b; 			// impact parameter
  int   npart1;		// proj npart
  int   npart2;		// targ npart
  int   ncoll;		// targ npart
} Event_tree;

//Minimum size Track structure.
//(There is a limit for the struct size
//if too large MAXTRK or too many variables
//stored in Track Struct, the program crashes
//when running.)
const int MAXTRK = 200000;

//track level tree info
typedef struct { 
  //PID
  int   id[MAXTRK];
  float m[MAXTRK];
  //momentum 
  float px[MAXTRK];
  float py[MAXTRK];
  float pz[MAXTRK];
  //rapidity space-time
  float x[MAXTRK];
  float y[MAXTRK];
  float z[MAXTRK];
  float t[MAXTRK];
  float eta[MAXTRK];
} Track;

int GetRandomSeed();


int main(int argc, char* argv[]) {

  TString outFileName="output.root";

  if(argc>2) {
    cout<<"wrong input parameter number:"<<argc<<endl;
    exit(1);
  }
  else if(argc==2){
    outFileName=argv[1];
  }

	Pythia pythia;

  pythia.settings.addMode("Main:triggerNcut", 0, true, true, 0, 1e9);

	pythia.readFile("input.dat");

  int kSeed = GetRandomSeed();
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d",kSeed));


//	pythia.readString("SoftQCD:nonDiffractive = off");
//	pythia.readString("SoftQCD:singleDiffractive = on");
//	pythia.readString("SoftQCD:doubleDiffractive = on");

//	pythia.readString("HadronLevel:all = on");
	//switch off decay to let resonance go for further evolution
//	pythia.readString("HadronLevel:Decay = off");
//	pythia.readString("Fragmentation:setVertices = off");
//	pythia.readString("HadronVertex:smearOn = off");
//	pythia.readString("PartonVertex:setVertex = on");
//	pythia.readString("PartonVertex:modeVertex = 4"); //D=2 Gaussion, 1 for uniform proton disc
//	pythia.readString("PartonVertex:randomPlane = off");
//	pythia.readString("PartonVertex:ProtonRadius = 0.7"); //D=0.7 at 14 TeV
//	pythia.readString("PartonVertex:EmissionWidth = 0.1"); //D=0.1 
	
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10.");
  pythia.readString("310:mayDecay = off"); //KS
  pythia.readString("111:mayDecay = off"); //pi0
  pythia.readString("3122:mayDecay = off"); //Lambda
  pythia.readString("3222:mayDecay = off"); //Sigma+
  pythia.readString("3212:mayDecay = off"); //Sigma0
  pythia.readString("3112:mayDecay = off"); //Sigma-
  pythia.readString("3322:mayDecay = off"); //Xi0
  pythia.readString("3312:mayDecay = off"); //Xi-
  pythia.readString("3334:mayDecay = off"); //Omega

  // Setup the beams.
	// AA collisions
//  pythia.readString("Beams:idA = 1000822080");
//  pythia.readString("Beams:idB = 1000822080"); // The lead ion.
//  pythia.readString("Beams:idA = 1000060120");
//  pythia.readString("Beams:idB = 1000060120"); // The carbon ion.
//  pythia.readString("Beams:idA = 1000030060");
//  pythia.readString("Beams:idB = 1000030060"); // The carbon ion.

	//pp collisions
//  pythia.readString("Beams:idA = 2212");
//  pythia.readString("Beams:idB = 2212"); 
//
//	// This forces the HeavyIons model to be used even for pp collisons.
//	pythia.readString("HeavyIon:mode = 2"); //force Angantyr model

  // Initialize the Angantyr model to fit the total and semi-includive
  // cross sections in Pythia within some tolerance.
  pythia.readString("HeavyIon:SigFitErr = "
                    "0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
  // These parameters are typicall suitable for sqrt(S_NN)=5TeV
  pythia.readString("HeavyIon:SigFitDefPar = "
                    "17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
  // A simple genetic algorithm is run for 20 generations to fit the
  // parameters.
  pythia.readString("HeavyIon:SigFitNGen = 5");

	//choose Monash 2013 tune = 14
	//pythia.readString("Tune:pp = 5");
	//pythia.readString("BeamRemnants:remnantMode = 1");
	//pythia.readString("ColourReconnection:mode = 1");
  // Enabling flavour ropes, setting model parameters.
  // The model is still untuned. These parameter values
  // are choosen for illustrative purposes.
//  pythia.readString("Ropewalk:RopeHadronization = on");
//  pythia.readString("Ropewalk:doShoving = on");
//  pythia.readString("Ropewalk:doFlavour = on");
//  pythia.readString("Ropewalk:r0 = 0.5");
//  pythia.readString("Ropewalk:m0 = 0.2");
//  pythia.readString("Ropewalk:beta = 0.2");
//	pythia.readString("Ropewalk:setFixedKappa = on");
//	pythia.readString("Ropewalk:presetKappa = 2");
	
	//pythia.readString("StringFlav:probQQtoQ = 0.04");
//	pythia.readString("StringZ:aLund = 1.5"); //D=0.68
//	pythia.readString("StringZ:bLund = 0.9"); //D=0.98
	//pythia.readString("StringFlav:probSQtoQQ = 0.19");

//	pythia.readString("MultipartonInteractions:pT0Ref = 1.90");

	//tree construction
//  TFile* treefile = new TFile("tree_AA.root","RECREATE");
  TTree* tree = new TTree("AMPT","AMPT Tree");

  //assign branch address
  Event_tree evt;
  //PID
  int   *id = new int[MAXTRK];
  float *m = new float[MAXTRK];
  //momentum 
  float *px = new float[MAXTRK];
  float *py = new float[MAXTRK];
  float *pz = new float[MAXTRK];
  //rapidity space-time
  float *x = new float[MAXTRK];
  float *y = new float[MAXTRK];
  float *z = new float[MAXTRK];
  float *t = new float[MAXTRK];
  float *eta = new float[MAXTRK];

  //event wise branch
  tree->Branch("Event", &evt, "multiplicity/I:b/F:npart1/I:npart2/I:ncoll/I");
  //PID branch
  tree->Branch("id",id,"id[multiplicity]/I");
  tree->Branch("m",m,"m[multiplicity]/F");
  //moment branch
  tree->Branch("px",px,"px[multiplicity]/F");
  tree->Branch("py",py,"py[multiplicity]/F");
  tree->Branch("pz",pz,"pz[multiplicity]/F");  
  //rapidity space-time branch
  tree->Branch("eta",eta,"eta[multiplicity]/F");  
  tree->Branch("x",x,"x[multiplicity]/F");
  tree->Branch("y",y,"y[multiplicity]/F");
  tree->Branch("z",z,"z[multiplicity]/F");  
  tree->Branch("time",t,"t[multiplicity]/F");

  //set up temp variables
//  int	mult   	= 0;
//  float	b      	= 0;
//  int	npart_1   	= 0;
//  int	npart_2   	= 0;
  int	ncoll_tmp   	= 0;


  // Book histogram for the centrality measure.
	TProfile *hpTVsNch = new TProfile("hpTVsNch", ";N_{ch}; <p_{T}>", 1000, 0, 5000, 0, 50);
  TH1D *hChargeEta = new TH1D("hChargeEta", ";#eta;N_{ch}", 40, -8, 8);

  // Sum up the weights of all generated events.
  double sumw = 0.0;

  // Initialise Pythia.
  pythia.init();

  bool doHeavyIons = HeavyIons::isHeavyIon(pythia.settings) || pythia.settings.mode("HeavyIon:mode") == 2;

  if(doHeavyIons) 
    cout<<"do heavyion set!"<<endl;
  else
    cout<<"No heavyion set!"<<endl;

  // Loop over events.
//  int nEvents = 1000;
  int nEvents = pythia.mode("Main:numberOfEvents");

  //cut to veto events with too few FS particles
  int Ncut = pythia.mode("Main:triggerNcut");

  if(nEvents<0) {
    nEvents=100;
    cout<<"Main:numberOfEvents undefined in input.dat"<<endl;
  }

//  shared_ptr<Angantyr> hiPtr = make_shared<Angantyr>(pythia);
  //Angantyr(pythia.getHeavyIonsPtr()).projBegin();

//  for ( int iEvent = 0; iEvent < nEvents; ++iEvent ) {
  int iEvent=0;
  while ( iEvent < nEvents ) {
    if ( !pythia.next() ) continue;

    //cut to veto events if not enough final state charge particles
    if ( pythia.event.nFinal(true)<Ncut ) continue;

//		if(pythia.info.code()==102) continue; //remove elastic for pp
		if(iEvent%100==0) cout<<iEvent<<endl;

		if(iEvent<3){
//			pythia.event.list(true, true);
//			pythia.partonSystems.list();

		}

    //cout<<"projectile number:"<<hiPtr->projectile.size()<<endl;
//    cout<<"projectile number:"<<hiPtr->HADRON<<endl;
//      for(vector<Nucleon>::iterator iter = hiPtr->projBegin();iter != hiPtr->projEnd(); ++iter){
//        cout << iter->id() << endl;
//      }


    if(doHeavyIons) {
		//used in version before 8243
		
//		evt.b = pythia.info.hiinfo->b();
//		evt.npart1 = pythia.info.hiinfo->nPartProj();
//		evt.npart2 = pythia.info.hiinfo->nPartTarg();
//		evt.ncoll = pythia.info.hiinfo->nCollTot();

		//used in version since 8301
		  evt.b = pythia.info.hiInfo->b();
		  evt.npart1 = pythia.info.hiInfo->nPartProj();
		  evt.npart2 = pythia.info.hiInfo->nPartTarg();
//	  	evt.ncoll = pythia.info.hiInfo->nCollTot();

//      cout<<evt.b<<endl;
      //cout<<pythia.EventInfo.code<<endl;
//      cout<<typeid(*pythia.info.hiInfo->subCollisionsPtr()).name()<<endl;
      ncoll_tmp=0;
      multiset<SubCollision>* collPtr = pythia.info.hiInfo->subCollisionsPtr();
      for ( multiset<SubCollision>::iterator cit = collPtr->begin(); cit!=collPtr->end(); ++cit){
//        cout<<"col nucleon type:"<<cit->nucleons()<<" b:"<<cit->b<<" proj idx:"<<cit->proj->index()<<" x:"<<cit->proj->bPos()<<" tgt idx:"<<cit->targ->index()<<" x:"<<cit->targ->bPos()<<endl;
//        cout<<cit->proj->bPos().pT()<<" "<<cit->targ->bPos().pT()<<endl;
        ncoll_tmp++;
      }
//      cout<<"b:"<<evt.b<<" tot col num:"<<ncoll_tmp<<" ncoll="<<pythia.info.hiInfo->nCollTot()<<endl;

      evt.ncoll=ncoll_tmp;
    }
    else {
    //when not in HeavyIon mode, set dumb values to these variables
      evt.b = 0;
      evt.npart1 = 0;
      evt.npart2 = 0;
      evt.ncoll = 0;
    }

    if(abs(pythia.info.idA())<3000&&abs(pythia.info.idB())<3000) evt.ncoll=pythia.info.nMPI();

//    cout<<pythia.info.code()<<" "<<pythia.info.name()<<" "<<pythia.info.nMPI()<<endl;

//    if(pythia.info.code()==101||pythia.info.nMPI()==0) continue;

    int nc=0; //charge particle number in an event
		
    //momentum conservation check
		double sum_x = 0;
		double sum_y = 0;
		double sum_z = 0;

		int itrk=0; //track counter

    int n_large_z = 0;

    //track-wise loop
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() ) {

				//assign the tree variables
				id[itrk] = p.id();
				m[itrk] = p.m();
				px[itrk] = p.px();
				py[itrk] = p.py();
				pz[itrk] = p.pz();
				eta[itrk] = p.eta();
				sum_x+=px[itrk];
				sum_y+=py[itrk];
				sum_z+=pz[itrk];

//				if(iEvent<3) cout<<i<<"id="<<p.id()<<" x="<<p.xProd()<<" px="<<p.px()<<endl;
				//convert vertex info from mm to fm
				x[itrk] = p.xProd()*1e12;
				y[itrk] = p.yProd()*1e12;
				z[itrk] = p.zProd()*1e12;
				t[itrk] = p.tProd()*1e12;

        if(z[itrk]>1) n_large_z++;

        sum_x+=p.px();
        sum_y+=p.py();
        sum_z+=p.pz();

        if(p.isCharged()) nc++;

				itrk++;
      }
    }
		if(iEvent<3) cout<<sum_x<<" "<<sum_y<<" "<<sum_z<<endl;

//    if(n_large_z>0) {
//      cout<<"Event:"<<iEvent<<" with large z, n="<<n_large_z<<endl;
//			pythia.event.list(true, true);
//    }

		evt.multiplicity = itrk;

		//fill the tree branches
		tree->Fill();

    // Keep track of the sum of waights
    double weight = pythia.info.weight();
    sumw += weight;

    // Go through the event again and fill the eta distributions.
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];

      if ( p.isFinal() && p.isCharged() && p.pT() > 0.1 ) {
        hChargeEta->Fill(p.eta());

				if(fabs(p.eta())<0.5)
					hpTVsNch->Fill(nc, p.pT());
      }
    }

    iEvent++;
  }

  // The run is over, so we write out some statistics.
  pythia.stat();
	
	TFile *outfile = new TFile(outFileName,"recreate");
	outfile->cd();
	tree->Write();
	hpTVsNch->Write();
  hChargeEta->Write();
	outfile->Close();

  return 0;
}

int GetRandomSeed()
{
//=============================================================================

  TDatime adt;
  new TSystem();
  UInt_t kTime = (UInt_t)adt.Get();
  UInt_t kProc = (UInt_t)gSystem->GetPid();//Get process id.
  UInt_t kInit = kTime - kProc*10;
  int kSeed = (int)kInit - (int)(((int)(kInit/1e8))*1e8);

  printf("Proc ID     = %d \n", kProc);
  printf("System time = %d \n", kTime);
  printf("Random number seed = %d \n", kSeed);

  return kSeed;
}


