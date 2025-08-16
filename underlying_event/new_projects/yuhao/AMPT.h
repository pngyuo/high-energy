//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 10 13:42:29 2016 by ROOT version 6.01/02
// from TTree AMPT/AMPT Tree
// found on file: ../create_tree/pp_lhc_out.root
//////////////////////////////////////////////////////////

#ifndef AMPT_h
#define AMPT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class AMPT {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Event_multiplicity;
   Float_t         Event_b;
   Int_t           Event_npart1; //new variables added after 2018/10/19
   Int_t           Event_npart2; //new variables added after 2018/10/19
   Int_t           Event_ncoll; //new variables added after 2018/10/19
   Int_t           id[100000];   //[multiplicity]
   Float_t         m[100000];   //[multiplicity]
   Float_t         px[100000];   //[multiplicity]
   Float_t         py[100000];   //[multiplicity]
   Float_t         pz[100000];   //[multiplicity]
   Float_t         eta[100000];   //[multiplicity]
   Float_t         x[100000];   //[multiplicity]
   Float_t         y[100000];   //[multiplicity]
   Float_t         z[100000];   //[multiplicity]
   Float_t         time[100000];   //[multiplicity]

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_id;   //!
   TBranch        *b_m;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_time;   //!

   AMPT(TTree *tree=0);
   virtual ~AMPT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AMPT_cxx
AMPT::AMPT(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../create_tree/pp_lhc_out.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../create_tree/pp_lhc_out.root");
      }
      f->GetObject("AMPT",tree);

   }
   Init(tree);
}

AMPT::~AMPT()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AMPT::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AMPT::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AMPT::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_multiplicity, &b_Event);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("m", m, &b_m);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("time", time, &b_time);
   Notify();
}

Bool_t AMPT::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AMPT::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AMPT::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AMPT_cxx
