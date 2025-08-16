#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>

void addText(Double_t x, Double_t y, const char* text)
  {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.035);  // 设置字体大小
    latex->Draw();
  }

void drawtest(){
  TCanvas *c1 = new TCanvas();
  TH2D *hframe = new TH2D("hframe", ";p_{T}(GeV/c); dN/(dp_{T}dy)(GeV/c)^{-1}", 100, 0, 6, 1000, 5e-5, 1e1);
  hframe->Draw();
  //————————————————————
  //这里就是提取下载的数据，这里是用了三个数据图
  TFile *infile1 = new TFile("K_class1.root");//修改名字
  TGraph *h1 = (TGraph*)infile1->Get("Pt");//这就是提取数据时修改五填的
  TFile *infile2 = new TFile("K_class5.root");
  TGraph *h2 = (TGraph*)infile2->Get("Pt");
  TFile *infile3 = new TFile("K_class10.root");
  TGraph *h3 = (TGraph*)infile3->Get("Pt");
  //——————————————————

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


gStyle->SetOptStat(0); // 隐藏统计框
legend->Draw("same");

addText(0.56, 0.81, "ALICE");
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
//c1->SaveAs("K_Pt.png");//保存为图片
}


