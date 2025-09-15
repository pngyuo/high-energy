#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveText.h>

void draw(){
    // ---------- 打开文件 ----------
    TFile *f = TFile::Open("output.root");
    if (!f){ cout << "open output.root failed!\n"; return; }

    // ---------- 平均 <N^TR> ----------
    TH1D *hNchTR = (TH1D*)f->Get("chargedParticlesHist");
    if (!hNchTR){ cout << "chargedParticlesHist not found!\n"; return; }
    double sumw = 0, sumn = 0;
    for (int i = 1; i <= hNchTR->GetNbinsX(); ++i){
        double ntr = hNchTR->GetBinCenter(i);
        double w   = hNchTR->GetBinContent(i);
        sumw += w * ntr;
        sumn += w;
    }
    double avgNTR = (sumn ? sumw / sumn : 0);
    cout << "Average N^TR = " << avgNTR << endl;

    // ---------- R_T 区间 ----------
    double rtLo[2] = {0.0, 0.5};
    double rtHi[2] = {0.5, 2.5};

    // ---------- 新建输出 hist ----------
    TH1D *hNj[2];
    hNj[0] = new TH1D("hNj05", "0#leqR_{T}<0.5;N_{jet}^{TR};Events",11,-0.5,10.5);
    hNj[1] = new TH1D("hNj525","0.5#leqR_{T}<2.5;N_{jet}^{TR};Events",11,-0.5,10.5);

    // ---------- 取 TTree ----------
    TTree *t = (TTree*)f->Get("tTRjet");
    if (!t){ cout << "tTRjet not found! 请先跑修改后的 analyse_decayer.cxx\n"; return; }

    int  NchTR, NjetTR;
    vector<double> *jetPt=0, *jetEta=0, *jetPhi=0;
    vector<int>    *jetLeadPID=0;
    vector<vector<int>> *jetConstPID=0;

    t->SetBranchAddress("NchTR",&NchTR);
    t->SetBranchAddress("NjetTR",&NjetTR);
    t->SetBranchAddress("jetPt",&jetPt);
    t->SetBranchAddress("jetEta",&jetEta);
    t->SetBranchAddress("jetPhi",&jetPhi);
    t->SetBranchAddress("jetLeadPID",&jetLeadPID);
    t->SetBranchAddress("jetConstPID",&jetConstPID);

    // ---------- 事件循环 ----------
    for (Long64_t ev=0; ev<t->GetEntries(); ++ev){
        t->GetEntry(ev);
        double rt = (avgNTR>0 ? double(NchTR)/avgNTR : -1);
        int idx = -1;
        if (rt>=rtLo[0] && rt<rtHi[0]) idx = 0;
        if (rt>=rtLo[1] && rt<rtHi[1]) idx = 1;
        if (idx<0) continue;

        hNj[idx]->Fill(NjetTR);

        cout << "Event " << ev << "  R_T=" << rt << "  NjetTR=" << NjetTR << endl;
        for (int j=0; j<NjetTR; ++j){
            cout << "  jet#" << j
                 << "  pt=" << jetPt->at(j)
                 << "  eta=" << jetEta->at(j)
                 << "  phi=" << jetPhi->at(j)
                 << "  leadPID=" << jetLeadPID->at(j)
                 << "  const={ ";
            for (int pid : jetConstPID->at(j)) cout << pid << " ";
            cout << "}\n";
        }
    }

    // ---------- 简单画图 ----------
    TCanvas c("c","TR jet multiplicity",800,400);
    c.Divide(2,1);
    for (int i=0;i<2;++i){ c.cd(i+1); hNj[i]->Draw(); }
    c.SaveAs("TRjetMult_vs_RT.png");
}