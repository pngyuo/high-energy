#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>

void drawConnectedPoints(TH1D* hist, Color_t color) {
  int n = hist->GetNbinsX();
  double* x = new double[n];
  double* y = new double[n];
  for (int i = 1; i <= n; ++i) {
    x[i-1] = hist->GetBinCenter(i);
    y[i-1] = hist->GetBinContent(i);
  }
  TGraph* graph = new TGraph(n, x, y);
  graph->SetLineColor(color);
  graph->SetLineWidth(2);
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(color);
  graph->Draw("l,same");
}
void addText(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.035);  // 设置字体大小
    latex->Draw();
  }
void doSliceHist(TH2D *hin, double val_mid[], double cut_fwd[], int Nbins){

  TH1D *hfwd = (TH1D*)hin->ProjectionY("hfwd");
  double Ntotal = hin->GetEntries();
  int start_bin = hfwd->FindFirstBinAbove(0);
  int end_bin = hfwd->FindLastBinAbove(0);
  cut_fwd[0] = hfwd->GetBinCenter(start_bin);
  //cut_fwd[Nbins] = hfwd->GetBinCenter(start_bin);
  cut_fwd[Nbins] = hfwd->GetBinCenter(end_bin);

  double fraction[9] = {0.92/100., 4.6/100., 9.2/100., 13.8/100., 18.4/100., 27.6/100., 36.8/100., 46./100., 64.5/100.};
  //double fraction[9] = {0.92/100., 10/100., 20./100., 30./100., 40./100., 50./100., 60./100., 70./100., 80./100.};
  double counts=0;
  int i_mark = 0; //tagger for current cen bin
  for(int i=end_bin; i>=start_bin; i--){
    if(hfwd->GetBinContent(i)>0){
      counts+=hfwd->GetBinContent(i);
      if(counts/Ntotal>fraction[i_mark]){
        i_mark++;
        cut_fwd[Nbins-i_mark]=hfwd->GetBinCenter(i);
        cout<<counts/Ntotal*100<<"%: "<<cut_fwd[Nbins-i_mark]<<endl;
      }
      if(i_mark==Nbins-1) break;
    }
  }

  for(int i=0; i<Nbins; i++){
    int low_bin = hfwd->FindBin(cut_fwd[i]);
    int high_bin = hfwd->FindBin(cut_fwd[i+1]);
    if(i!=Nbins-1) high_bin=high_bin-1;
    val_mid[i] = hin->ProjectionX(Form("hx_%d",i), low_bin, high_bin)->GetMean();
  }
}


//dN/dptdy
TH1D *getCentPt(TString filename, TString histname){
  const int Nbins_centrality=10;
  TFile *infile = new TFile(filename);
  TH2D *hFwdVsMid = (TH2D*)infile->Get("hFwdVsMid");
  TH1D *hNtrkV0 = (TH1D*)infile->Get("hNtrkV0");
  double val_mid[Nbins_centrality];
  double cut_fwd[Nbins_centrality+1];
  doSliceHist(hFwdVsMid, val_mid, cut_fwd, Nbins_centrality);
  TH3D *hist3D = (TH3D*)infile->Get(histname);
  double yield[Nbins_centrality];
  TH1D *hist;
  for(int i=0; i<Nbins_centrality; i++){
    if(i==0||i==5||i==9) {
      hist = hist3D->ProjectionX(Form("%s_%i",histname.Data(),i),hist3D->GetYaxis()->FindBin(cut_fwd[i]), hist3D->GetYaxis()->FindBin(cut_fwd[i+1])); 
      hist->Scale(1./hNtrkV0->Integral(hNtrkV0->FindBin(cut_fwd[i]),hNtrkV0->FindBin(cut_fwd[i+1])), "width");
      hist->SetLineWidth(0);
      hist->Draw("same");
      if(i==9){
        hist->SetMarkerStyle(20);
        hist->SetMarkerColor(kBlack);
        drawConnectedPoints(hist, kBlack);
      }
      if(i==5){
        hist->SetMarkerStyle(21);
        hist->SetMarkerColor(kRed);
        drawConnectedPoints(hist, kRed);
      }
      if(i==0){
        hist->SetMarkerStyle(22);
        hist->SetMarkerColor(kBlue);
        drawConnectedPoints(hist, kBlue);
      }
    }
  }
  return hist;
}


void draw(){
  TCanvas *c1 = new TCanvas();
  TH2D *hframe = new TH2D("hframe", ";p_{T}(GeV/c); dN/(dp_{T}dy)(GeV/c)^{-1}", 100, 0, 6, 1000, 5e-5, 1e1);
  hframe->Draw();
  //————————————————————
  //这里就是提取下载的数据，
  TFile *infile1 = new TFile("K_class1.root");
  TGraph *h1 = (TGraph*)infile1->Get("Pt");
  TFile *infile2 = new TFile("K_class5.root");
  TGraph *h2 = (TGraph*)infile2->Get("Pt");
  TFile *infile3 = new TFile("K_class10.root");
  TGraph *h3 = (TGraph*)infile3->Get("Pt");
  //——————————————————
  TH1D *hist = getCentPt("hist_output.root", "hKCh");

  h1->SetMarkerStyle(24);
  h1->SetMarkerSize(1);
  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(h1->GetLineColor());
  h1->Draw("eP");

  h2->SetMarkerStyle(25);
  h2->SetMarkerSize(1);
  h2->SetLineColor(kRed);
  h2->SetMarkerColor(h2->GetLineColor());
  h2->Draw("eP");

  h3->SetMarkerStyle(26);
  h3->SetMarkerSize(1);
  h3->SetLineColor(kBlue);
  h3->SetMarkerColor(h3->GetLineColor());
  h3->Draw("eP");
  c1->SetLogy(1);

TLegend *legend = new TLegend(0.59, 0.65, 0.61, 0.8);
legend->SetNColumns(1);
legend->SetFillStyle(0);

TGraphErrors *marker1 = new TGraphErrors();
marker1->SetMarkerStyle(24); // 黑色空心圆标记
marker1->SetMarkerColor(kBlack);
marker1->SetMarkerSize(1);
marker1->SetLineWidth(2); 
marker1->SetLineColor(kBlack);
legend->AddEntry(marker1, " ", "ep");

TGraphErrors *marker3 = new TGraphErrors();
marker3->SetMarkerStyle(25); // 红色空心方块
marker3->SetMarkerColor(kRed);
marker3->SetMarkerSize(1);
marker3->SetLineWidth(2); 
marker3->SetLineColor(kRed);
legend->AddEntry(marker3, " ", "ep");

TGraphErrors *marker5 = new TGraphErrors();
marker5->SetMarkerStyle(26); // 蓝色空心三角形
marker5->SetMarkerColor(kBlue);
marker5->SetMarkerSize(1);
marker5->SetLineWidth(2); 
marker5->SetLineColor(kBlue);
legend->AddEntry(marker5, " ", "ep");

TLegend *legend2 = new TLegend(0.70, 0.65, 1, 0.8);
legend2->SetNColumns(1);
legend2->SetFillStyle(0);

TGraphErrors *marker2 = new TGraphErrors();
marker2->SetMarkerStyle(20); // 黑色实心圆标记
marker2->SetMarkerColor(kBlack);
marker2->SetMarkerSize(1);
marker2->SetLineWidth(2); 
marker2->SetLineColor(kBlack);
legend2->AddEntry(marker2, " ", "lp");

TGraphErrors *marker4 = new TGraphErrors();
marker4->SetMarkerStyle(21); // 红色实心方块
marker4->SetMarkerColor(kRed);
marker4->SetMarkerSize(1);
marker4->SetLineWidth(2); 
marker4->SetLineColor(kRed);
legend2->AddEntry(marker4, " ", "lp");

TGraphErrors *marker6 = new TGraphErrors();
marker6->SetMarkerStyle(22); // 蓝色实心三角形
marker6->SetMarkerColor(kBlue);
marker6->SetMarkerSize(1);
marker6->SetLineWidth(2); 
marker6->SetLineColor(kBlue);
legend2->AddEntry(marker6, " ", "lp");

gStyle->SetOptStat(0); // 隐藏统计框
legend->Draw("same");
legend2->Draw("same");
addText(0.56, 0.81, "ALICE");
addText(0.65, 0.81, "PYTHIA8 Monash");
addText(0.46, 0.76, "0-0.92%");
addText(0.46, 0.71, "13.8-18.4%");
addText(0.46, 0.66, "64.5-100%");
TLatex *latex = new TLatex(2, 2, "\\Kappa^{+}+\\Kappa^{-}");
//TLatex *latex = new TLatex(2, 20, "\\pi^{+}+\\pi^{-}");
//TLatex *latex = new TLatex(2, 2, "P+\\bar{P}");
TLatex *latex1 = new TLatex(0.5, 0.0005, "pp,\\sqrt{s}=13TeV,|y|<0.5");
latex->SetTextSize(0.05);
latex->Draw();
latex1->SetTextSize(0.04);
latex1->Draw();
c1->SaveAs("K_Pt.png");
}


