/*************************************************************************
    > File Name: main.C
    > Author: Liang Zheng
    > Mail: liangzhphy@gmail.com 
    > Created Time: Wed Nov 11 19:59:58 2020
    > main steering file to coordinate the running for pythia8 based 
    > mutliphase transport code
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

#include "amptWrapper.h"

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

class myImpactGenerator : public ImpactParameterGenerator {
public:
  Vec4 generate(double &weight) const{
    //double b = 1.0;
    double b = sqrt(-2.0*log(rndPtr->flat()))*width();;
    double phi = 0.5*M_PI;
    if(bmax<998||bmin>0.001) {
      while(b>bmax||b<bmin)
        b = sqrt(-2.0*log(rndPtr->flat()))*width();;
    }
    weight = 2.0*M_PI*width()*width()*exp(0.5*b*b/(width()*width()));
//    cout<<"using my generator"<<endl;
    return Vec4(b*sin(phi), b*cos(phi), 0.0, 0.0);
  }

  void SetMaxImp(double bin) { bmax=bin; }
  void SetMinImp(double bin) { bmin=bin; }

private:
  double bmin;
  double bmax;
};

class MyHIUserHooks : public HIUserHooks {
public:

  MyHIUserHooks() { 
    impactGen = new myImpactGenerator(); 
    impactGen->SetMaxImp(999);
    impactGen->SetMinImp(0);
  }
  ~MyHIUserHooks() { delete impactGen; }

  virtual bool hasImpactParameterGenerator() const { return true; }
  virtual ImpactParameterGenerator * impactParameterGenerator() const {
    return impactGen;
  }

  void setImpactMax(double bin) {
    if(bin<999) cout<<"New bmax in Angantyr:"<<bin<<endl;
    impactGen->SetMaxImp(bin);
  }

  void setImpactMin(double bin) {
    if(bin>0) cout<<"New bmin in Angantyr:"<<bin<<endl;
    impactGen->SetMinImp(bin);
  }

private:
  myImpactGenerator *impactGen;

};


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

  pythia.settings.addParm("AMPT:xmu",3.203, true, true, 0, 1e9);
  pythia.settings.addParm("AMPT:alpha",0.33, true, true, 0, 1);
  pythia.settings.addParm("AMPT:NTMAX",150, true, true, 1, 1000);
  pythia.settings.addMode("AMPT:isoft",4, true, true, 1, 4);
  pythia.settings.addMode("AMPT:ICmode",1, true, true, 1, 2); //initial spatial condition for excited string
  pythia.settings.addParm("Angantyr:bmax",999, true, true, 0, 999);
  pythia.settings.addParm("Angantyr:bmin",0, true, true, 0, 999);

  pythia.settings.addMode("Main:triggerNcut", 0, true, true, 0, 10000);

  int kSeed = GetRandomSeed();
  //int kSeed = 8993;
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d",kSeed));


	pythia.readFile("input.dat");

  para2.xmu=pythia.settings.parm("AMPT:xmu");
  para2.alpha=pythia.settings.parm("AMPT:alpha");
  input2.ntmax=int(pythia.settings.parm("AMPT:ntmax"));
  anim.isoft=pythia.mode("AMPT:isoft");

  int geoICmode = pythia.mode("AMPT:ICmode");

  cout<<para2.xmu<<" "<<para2.alpha<<" "<<input2.ntmax<<" "<<anim.isoft<<endl;


//	pythia.readString("SoftQCD:all = off");
//	pythia.readString("HardQCD:all = off");
	pythia.readString("SoftQCD:nonDiffractive = on");
//  pythia.readString("HardQCD:hardccbar = on");
////  pythia.readString("HardQCD:hardbbbar = on");
//
//	pythia.readString("PartonLevel:MPI = off");
//	pythia.readString("PartonLevel:ISR = off");
//	pythia.readString("PartonLevel:FSR = off");
//
//	pythia.readString("HadronLevel:all = on");
	pythia.readString("HadronLevel:Hadronize = on");
	//switch off decay to let resonance go for further evolution
	pythia.readString("HadronLevel:Decay = off");
	pythia.readString("Fragmentation:setVertices = off");
//	pythia.readString("HadronVertex:mode = -1");
//	pythia.readString("HadronVertex:smearOn = off");
	pythia.readString("PartonVertex:setVertex = off");
	pythia.readString("PartonVertex:modeVertex = 1"); //D=2 Gaussion, 1 for uniform proton disc with b constraint
	pythia.readString("PartonVertex:randomPlane = off");
	pythia.readString("PartonVertex:ProtonRadius = 0.7"); //D=0.7 at 14 TeV
	pythia.readString("PartonVertex:EmissionWidth = 0.1"); //D=0.1 


	/*
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10.");
  pythia.readString("310:mayDecay = off");
  pythia.readString("111:mayDecay = off");
  pythia.readString("411:mayDecay = off");
  pythia.readString("421:mayDecay = off");
	*/

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
	pythia.readString("HeavyIon:mode = 2"); //force Angantyr model
//
//  pythia.readString("Beams:eCM = 5020.0"); //c.m.s frame
//  pythia.readString("Beams:frameType = 1"); //chose c.m.s frame
//

//  pythia.readString("Angantyr:GlauberOnly = on");
//
	//pA collisions asymmetric lab frame
//  pythia.readString("Beams:idA = 2212");
//  pythia.readString("Beams:idB = 1000822080"); // The lead ion.
//  pythia.readString("Beams:eA = 4000");
//  pythia.readString("Beams:eB = 1570");
//  pythia.readString("Beams:frameType = 1");
//  pythia.readString("Beams:eCM = 7000.0"); //c.m.s frame

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
	//pythia.readString("MultipartonInteractions:pT0Ref = 1.90");

  auto myHIUserHooks = make_shared<MyHIUserHooks>();
  pythia.settings.parm("AMPT:xmu");
  //myHIUserHooks->setImpactMax(1);
  myHIUserHooks->setImpactMax(pythia.settings.parm("Angantyr:bmax"));
  myHIUserHooks->setImpactMin(pythia.settings.parm("Angantyr:bmin"));
  pythia.setHIHooks(myHIUserHooks);

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

  cout<<"pzA="<<pythia.info.pzA()<<" "<<pythia.info.pzB()<<endl;
  cout<<"ezA="<<pythia.info.eA()<<" "<<pythia.info.eB()<<endl;
  cout<<"pzA="<<pythia.settings.parm("Beams:eA")<<" "<<pythia.settings.parm("Beams:eB")<<endl;
  float pzA=pythia.settings.parm("Beams:eA");
  float pzB=pythia.settings.parm("Beams:eB");
  //initialize ampt data
  call_ini_ampt(pzA, pzB, kSeed);

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
  //		evt.ncoll = pythia.info.hiInfo->nCollTot();
  //
  //    cout<<evt.b<<endl;
      //cout<<pythia.EventInfo.code<<endl;
  //    cout<<typeid(*pythia.info.hiInfo->subCollisionsPtr()).name()<<endl;
      ncoll_tmp=0;
      multiset<SubCollision>* collPtr = pythia.info.hiInfo->subCollisionsPtr();
      for ( multiset<SubCollision>::iterator cit = collPtr->begin(); cit!=collPtr->end(); ++cit){
  //      cout<<"col nucleon type:"<<cit->nucleons()<<" b:"<<cit->b<<" proj idx:"<<cit->proj->index()<<" x:"<<cit->proj->bPos()<<" tgt idx:"<<cit->targ->index()<<" x:"<<cit->targ->bPos()<<endl;
  //      cout<<cit->proj->bPos().pT()<<" "<<cit->targ->bPos().pT()<<endl;
        ncoll_tmp++;
      }
  //    cout<<"b:"<<evt.b<<" tot col num:"<<ncoll<<" ncoll="<<pythia.info.hiInfo->nCollTot()<<endl;
      evt.ncoll=ncoll_tmp;
    }
    else {
    //when not in HeavyIon mode, set dumb values to these variables
      evt.b = 0;
      evt.npart1 = 0;
      evt.npart2 = 0;
      evt.ncoll = 0;
    }

    //rewrite pp ncoll as the MPI coll number instead of n-n level coll number
    if(abs(pythia.info.idA())<3000&&abs(pythia.info.idB())<3000) evt.ncoll=pythia.info.nMPI();

		if(iEvent<3){
      cout<<"b:"<<evt.b<<" ncoll:"<<evt.ncoll<<endl;
//			pythia.event.list(true, true);
//			pythia.partonSystems.list();
		}


    int nc=0;
		

		int itrk=0;

    //momentum conservation check
		double sum_x = 0;
		double sum_y = 0;
		double sum_z = 0;
		double sum_E = 0;

    //initialize ampt data array ready to take pythia particle list
    hmain1.natt=0;
    long natt_tmp=0;

    //track-wise loop
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() ) {

        if(p.isCharged()) nc++;

				//assign the tree variables
//				id[itrk] = p.id();
//				m[itrk] = p.m();
//				px[itrk] = p.px();
//				py[itrk] = p.py();
//				pz[itrk] = p.pz();
//				eta[itrk] = p.eta();

//        if(p.status()>=81&&p.status()<=89){
//          Vec4 mo1=pythia.event[p.mother1()].vProd();
//          Vec4 mo2=pythia.event[p.mother2()].vProd();
//          Vec4 vtmp=mo1+(mo2-mo1)*(p.y()-pythia.event[p.mother1()].y())/(pythia.event[p.mother2()].y()-pythia.event[p.mother1()].y());
//          p.vProd(vtmp.px(), vtmp.py(), vtmp.pz(), vtmp.e());
//        }

        if(p.status()>=81&&p.status()<=89){
          double bHalf=evt.b*0.5;
          
          //assign hadron vertex in the overlap region with a minimum 
          //overlap region at large impact parameter
          if(geoICmode==1){
            double rProton=0.8;
            if(bHalf>0.95*rProton) bHalf=0.95*rProton;
            double xmax=rProton-bHalf;
            double ymax=sqrt(rProton*rProton-bHalf*bHalf);
            bool accept=false;
            while(!accept){
              double xrand=(2*pythia.rndm.flat()-1)*xmax;
              double yrand=(2*pythia.rndm.flat()-1)*ymax;
              double rA2=pow2(xrand-bHalf)+yrand*yrand;
              double rB2=pow2(xrand+bHalf)+yrand*yrand;
              if(rA2<rProton*rProton&&rB2<rProton*rProton){
                p.vProd(xrand*1e-12, yrand*1e-12, 0, 0);
                accept=true;
              }
            }
          }
          //assign hadron vertex based on ordered rapidity with
          //beam rapidity and pos sit at two ends of impact parameter
          else if(geoICmode==2){
            double rap_max=pythia.event[1].y();
            double rap_min=pythia.event[2].y();
            double rap_p=p.y();
            double x_tmp = (bHalf-(-bHalf))*(rap_p-rap_min)/(rap_max-rap_min) +(-bHalf);
            //cout<<"rap_max="<<rap_max<<" rap_min="<<rap_min<<" bHalf="<<bHalf<<" x_tmp:"<<x_tmp<<" rap_p"<<rap_p<<endl;
            pair<double, double> xy=pythia.rndm.gauss2();
            double r_smear=0.01;
            //p.vProd(x_tmp*1e-12, 0, 0, 0);
            p.vProd((xy.first+x_tmp)*r_smear*1e-12, xy.second*r_smear*1e-12, 0, 0);
          }
          else
            p.vProd(0, 0, 0, 0);
        }
        //non hard string fragmentation parts assigned to beam vertex
        else{
          if(p.y()>0) 
            p.vProd(evt.b*0.5, 0, 0, 0);
          else
            p.vProd(-evt.b*0.5, 0, 0, 0);
          cout<<"other than 81-89 string hadronizaton KS:"<<p.status()<<endl;
        }


//        if(p.status()>=81&&p.status()<=84){
////          cout<<pythia.event[p.mother1()].xProd()*1e12<<" "<<pythia.event[p.mother1()].yProd()*1e12<<endl;
////          cout<<pythia.event[p.mother2()].xProd()*1e12<<" "<<pythia.event[p.mother2()].yProd()*1e12<<endl;
//        }
//        else if(p.status()>=85&&p.status()<=89){
//          cout<<"junction hadron:"<<i<<" id:"<<p.id()<<" mother1:"<<pythia.event[p.mother1()].id()<<" mother2:"<<pythia.event[p.mother2()].id()<<endl;
//        }
//        else {
//          cout<<"other hadron:"<<i<<" status:"<<p.status()<<" id:"<<p.id()<<" mother1:"<<pythia.event[p.mother1()].id()<<" mother2:"<<pythia.event[p.mother2()].id()<<endl;
//        }

//				if(iEvent<3) cout<<i<<"id="<<p.id()<<" x="<<p.xProd()<<" px="<<p.px()<<endl;
				//convert vertex info from mm to fm
//				x[itrk] = p.xProd()*1e12;
//				y[itrk] = p.yProd()*1e12;
//				z[itrk] = p.zProd()*1e12;
//				t[itrk] = p.tProd()*1e12;
//        cout<<itrk<<" "<<id[itrk]<<" "<<x[itrk]<<" "<<y[itrk]<<" "<<z[itrk]<<" "<<t[itrk]<<endl;
//        cout<<itrk<<" "<<id[itrk]<<" "<<p.xProd()*1e12<<" "<<p.yProd()*1e12<<" "<<p.zProd()*1e12<<" "<<p.tProd()*1e12<<endl;

        //copy pythia particle list to ampt ready for string melting
        hmain2.katt[0][natt_tmp]=p.id();
        hmain2.katt[1][natt_tmp]=20;     

        hmain2.katt[2][natt_tmp]=1; //status number
        hmain2.katt[3][natt_tmp]=0;
//       ****** identify the mother particle
        hmain2.patt[0][natt_tmp]=p.px();
        hmain2.patt[1][natt_tmp]=p.py();
        hmain2.patt[2][natt_tmp]=p.pz();
        hmain2.patt[3][natt_tmp]=p.e();
        hmain1.eatt+=p.e();
        arprc.gxar[natt_tmp] = p.xProd()*1e12;
        arprc.gyar[natt_tmp] = p.yProd()*1e12;
        //assume all hadrons are produced at 0 point longitudinally
        arprc.gzar[natt_tmp] = 0;
        arprc.ftar[natt_tmp] = 0;
        //arprc.gzar[natt_tmp] = p.zProd()*1e12;
        //arprc.ftar[natt_tmp] = p.tProd()*1e12;
        arprc.itypar[natt_tmp] = p.id();
        arprc.pxar[natt_tmp] = p.px();
        arprc.pyar[natt_tmp] = p.py();
        arprc.pzar[natt_tmp] = p.pz();
        arprc.pear[natt_tmp] = p.e();
        arprc.xmar[natt_tmp] = p.m();
        precpa.xstrg0[natt_tmp]=p.xProd()*1e12;
        precpa.ystrg0[natt_tmp]=p.yProd()*1e12;
        precpa.istrg0[natt_tmp]=0;


//        cout<<p.px()<<" "<<hmain2.patt[1][natt_tmp]<<endl;

        //for photon or leptonic particles, assign status number=40
        if((p.id()>=22&&p.id()<=37)||(abs(p.id())>=11&&abs(p.id())<=18)) hmain2.katt[2][natt_tmp]=40; 
        natt_tmp++;


//				cout << setiosflags(ios::scientific);
//				if(t[itrk]<0)  cout<<iEvent<<" "<<i<<" "<<itrk<<" "<<id[itrk]<<" "<<z[itrk]<<" "<<t[itrk]<<endl;
//				if(itrk<100) cout<<i<<" "<<itrk<<" "<<id[itrk]<<" "<<x[itrk]<<" "<<y[itrk]<<" "<<z[itrk]<<" "<<t[itrk]<<endl;
				itrk++;
      }
    }

    //save the number of tracks to be used in zpc
    hmain1.natt=natt_tmp;
//    cout<<"this event natt="<<hmain1.natt<<endl;
    call_run_ampt(iEvent,evt.ncoll, evt.b);

		evt.multiplicity = rootdt.nline;

    // Keep track of the sum of waights
    double weight = pythia.info.weight();
    sumw += weight;

    // Go through the event again and fill the tree
    for (int i = 0; i < evt.multiplicity; ++i) {

      id[i]=rootdt.id[i];
      px[i]=rootdt.proot[0][i];
      py[i]=rootdt.proot[1][i];
      pz[i]=rootdt.proot[2][i];
      m[i]=rootdt.proot[3][i];
      x[i]=rootdt.xroot[0][i];
      y[i]=rootdt.xroot[1][i];
      z[i]=rootdt.xroot[2][i];
      t[i]=rootdt.xroot[3][i];
      double p=sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
      double pt=sqrt(px[i]*px[i]+py[i]*py[i]);
      eta[i]=0.5*log((p+pz[i])/(p-pz[i]));

      sum_x+=px[i];
      sum_y+=py[i];
      sum_z+=pz[i];
      sum_E+=sqrt(p*p+m[i]*m[i]);

      if (  pt > 0.1 ) hChargeEta->Fill(eta[i]);
      if(fabs(eta[i])<0.5) hpTVsNch->Fill(evt.multiplicity, pt);
    }

		if(iEvent<3) cout<<"mom conservation:"<<sum_x<<" "<<sum_y<<" "<<sum_z<<" "<<sum_E<<endl;

		//fill the tree branches
		tree->Fill();


    iEvent++;
  }

  // The run is over, so we write out some statistics.

  // Befor we end we print out some statistics. Also, we want to check
  // that our generated centrality classes were the same as we
  // guessed.
  pythia.stat();
  // And we're done!

	
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
  int kSeed = (int)kInit - (int)(((int)(kInit/1e6))*1e6);

  printf("Proc ID     = %d \n", kProc);
  printf("System time = %d \n", kTime);
  printf("Random number seed = %d \n", kSeed);

  return kSeed;
}


