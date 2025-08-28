#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TLine.h>
#include <TLegend.h>
#include <vector>
#include <map>

const double PI = 3.1415926;
const int N_cent = 8;  // 将N_cent定义为常量
double gYIntegral = 0;

// 全局变量用于存储所有文件的rn数据
std::map<int, std::vector<TGraphErrors*>> an_f2_graphs; // 按中心性存储所有文件的an_f2图形
std::map<int, std::vector<TGraphErrors*>> cn_f2_graphs; // 按中心性存储所有文件的cn_f2图形
std::vector<TString> file_names; // 存储文件名

void divideHist(TH2D *hin) {
    for (int i = 1; i <= hin->GetNbinsX(); ++i) {
        for (int j = 1; j <= hin->GetNbinsY(); ++j) {
            double dphi = hin->GetXaxis()->GetBinCenter(i);
            double deta = hin->GetYaxis()->GetBinCenter(j);
            double ratio = (1 - 1e-4) * (4.75 - fabs(deta)) / 4.75 + 1e-4;
            double val = hin->GetBinContent(i, j);
            val = val / ratio;
            if (fabs(deta) > 4.75) {
                val = 0;
            }
            hin->SetBinContent(i, j, val);
        }
    }
}

void removeBkg(TH1D *hin) {
    double bkg = hin->GetMinimum();
    for (int i = 1; i <= hin->GetNbinsX(); ++i) {
        hin->SetBinContent(i, hin->GetBinContent(i) - bkg);
    }
}

TH1D *sumRidge(TH2D *hin, TString name){
    TH1D *hout = new TH1D(name, "", hin->GetNbinsX(), hin->GetXaxis()->GetXmin(), hin->GetXaxis()->GetXmax()); 
    TH1D *htmp1 = (TH1D*)hin->ProjectionX("htmp1", hin->GetYaxis()->FindBin(-11.0), hin->GetYaxis()->FindBin(-1.0));
    TH1D *htmp2 = (TH1D*)hin->ProjectionX("htmp2", hin->GetYaxis()->FindBin(1.0), hin->GetYaxis()->FindBin(11.0)); 
    hout->Add(htmp1, htmp2);
    delete htmp1;
    delete htmp2;
    return hout;
}

TH1D *getY(TH1D *hSame, TH1D *hMix, TString name, double Ntrig) {
    TH1D *hout = new TH1D(name, "", hSame->GetNbinsX(), hSame->GetXaxis()->GetXmin(), hSame->GetXaxis()->GetXmax());
    hout->Divide(hSame, hMix);
    double Btot = hMix->Integral();
    hout->Scale(1. / Ntrig / 2 / PI * Btot);
    return hout;
}

TH1D *getYFromHist(TH2D *hSame, TH2D *hMix, TString name, double Ntrig) {
    TH1D *hSameRidge = sumRidge(hSame, "same" + name);
    TH1D *hMixRidge = sumRidge(hMix, "mix" + name);
    TH1D *hY = getY(hSameRidge, hMixRidge, "Y" + name, Ntrig);
    delete hSameRidge;
    delete hMixRidge;
    return hY;
}

double funTemplate(double *x, double *par) {
    static TF1 *funRidge = nullptr;
    if (!funRidge) {
        funRidge = new TF1("funRidge", "[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))", -0.5*PI, 1.5*PI);
        funRidge->SetParameters(par[0], par[1], par[2], par[3], par[4]);
    } else {
        funRidge->SetParameters(par[0], par[1], par[2], par[3], par[4]);
    }
    return funRidge->Eval(x[0]);
}
  
void plotRnVsEta(TGraphErrors* graph, const TString& name, const TString& title) {
    TCanvas c("c", "c", 800, 600);
    graph->SetTitle(title);
    graph->GetXaxis()->SetTitle("|#eta^{a}|");
    graph->GetYaxis()->SetTitle("r_{n}");
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue);
    graph->Draw("AP");
    
    // 添加参考线 y=1
    TLine line(graph->GetXaxis()->GetXmin(), 1, 
               graph->GetXaxis()->GetXmax(), 1);
    line.SetLineStyle(2);
    line.SetLineColor(kRed);
    line.Draw();
    
    c.SaveAs(Form("rn_plots/%s.png", name.Data()));
}

// 反对称相除拟合函数
void fitRnAndSave(TGraphErrors* graph, const TString& name, double& F_n, double& F_n_err) {
    TCanvas c;
    graph->Draw("AP");
    
    // 定义反对称相除拟合函数：r_n(|\eta^a|) = 1 - 2F_n|\eta^a| + 2F_n^2|\eta^a|^2
    TF1 *fitFunc = new TF1("fitFunc", "1 - 2*[0]*x + 2*[0]*[0]*x*x", 
                          graph->GetXaxis()->GetXmin(), 
                          graph->GetXaxis()->GetXmax());
    fitFunc->SetParameter(0, 0.03); // 初始参数
    graph->Fit("fitFunc", "Q");
    
    // 提取参数
    F_n = fitFunc->GetParameter(0);
    F_n_err = fitFunc->GetParError(0);
    
    // 绘图设置
    TLatex text;
    text.SetTextSize(0.04);
    text.DrawLatexNDC(0.25, 0.85, Form("F_{n} = %.4f #pm %.4f", F_n, F_n_err));
    
    c.SaveAs(Form("fit_results/%s.png", name.Data()));
    delete fitFunc;
}

// 计算反对称相除值 r_n(|\eta^a|) = v_n(-|\eta^a|)/v_n(|\eta^a|)
void getRnFromVn(const double* vn, const double* vn_err, double* rn, double* rn_err, int nBins) {
    // 假设 eta bins 是对称的，从负到正
    int halfBins = nBins / 2;
    
    for (int k = 0; k < halfBins; k++) {
        int negIdx = halfBins - 1 - k; // 负 eta 索引
        int posIdx = halfBins + k;     // 正 eta 索引
        
        double vn_neg = vn[negIdx];
        double vn_pos = vn[posIdx];
        double vn_neg_err = vn_err[negIdx];
        double vn_pos_err = vn_err[posIdx];
        
        if (vn_pos != 0) {
            rn[k] = vn_neg / vn_pos;
            // 误差传播：r = a/b, dr = r * sqrt( (da/a)^2 + (db/b)^2 )
            rn_err[k] = rn[k] * sqrt(pow(vn_neg_err/vn_neg, 2) + pow(vn_pos_err/vn_pos, 2));
        } else {
            rn[k] = 0;
            rn_err[k] = 0;
        }
    }
}

void SavePlot(TH1D *hY, TString name) {
    double par[5];
    for (int i = 0; i < 5; i++) {
        par[i] = hY->GetFunction("fitFun")->GetParameter(i);
    }
    double chi2 = hY->GetFunction("fitFun")->GetChisquare();
    int ndf = hY->GetFunction("fitFun")->GetNDF();
    TLatex *text = new TLatex();
    text->SetTextSize(0.04);
    text->DrawLatexNDC(0.25, 0.85, Form("G=%f", par[0]));
    text->DrawLatexNDC(0.25, 0.80, Form("a1=%f", par[1]));
    text->DrawLatexNDC(0.25, 0.75, Form("a2=%f", par[2]));
    text->DrawLatexNDC(0.25, 0.70, Form("a3=%f", par[3]));
    text->DrawLatexNDC(0.25, 0.65, Form("a4=%f", par[4]));
    text->DrawLatexNDC(0.25, 0.60, Form("chi2=%.2f, ndf=%d", chi2, ndf));
}

double getA2(TH1D *hY, double &a1, double &a2, double &a3, double &a4, double &G, double &chi2, int &ndf, double &a2err) {
    gYIntegral = hY->Integral("width");
    TF1 *fitFun = new TF1("fitFun", funTemplate, -0.5 * PI, 1.5 * PI, 5); 
    fitFun->SetParameters(0.1,0.1, 0.1, 0.1, 0.1); 
    hY->Fit("fitFun");
    a1 = fitFun->GetParameter(1);
    a2 = fitFun->GetParameter(2);
    a3 = fitFun->GetParameter(3);
    a4 = fitFun->GetParameter(4);
    a2err=fitFun->GetParError(2);
    G = gYIntegral / (2.0 * PI);
    chi2 = fitFun->GetChisquare();
    ndf = fitFun->GetNDF();
    delete fitFun;
    return a2;
}

double getA3(TH1D *hY, double &a1, double &a2, double &a3, double &a4, double &G, double &chi2, int &ndf, double &a3err) {
    gYIntegral = hY->Integral("width");
    TF1 *fitFun = new TF1("fitFun", funTemplate, -0.5 * PI, 1.5 * PI, 5); 
    fitFun->SetParameters(0.1,0.1, 0.1, 0.1, 0.1); 
    hY->Fit("fitFun");
    a1 = fitFun->GetParameter(1);
    a2 = fitFun->GetParameter(2);
    a3 = fitFun->GetParameter(3);
    a4 = fitFun->GetParameter(4);
    a3err=fitFun->GetParError(3);
    G = gYIntegral / (2.0 * PI);
    chi2 = fitFun->GetChisquare();
    ndf = fitFun->GetNDF();
    delete fitFun;
    return a3;
}

// 绘制四个文件的rn比较图
void plotCombinedRn(const std::vector<TGraphErrors*>& graphs, const TString& name, 
                   const TString& title, const std::vector<TString>& labels) {
    TCanvas c("c", "c", 800, 600);
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    
    bool first = true;
    int colors[] = {kRed, kBlue, kGreen+2, kMagenta}; // 添加第四个颜色
    int markers[] = {20, 21, 22, 23}; // 添加第四个标记
    
    for (int i = 0; i < graphs.size(); i++) {
        if (!graphs[i]) continue;
        
        graphs[i]->SetMarkerStyle(markers[i]);
        graphs[i]->SetMarkerColor(colors[i]);
        graphs[i]->SetLineColor(colors[i]);
        
        if (first) {
            graphs[i]->SetTitle(title);
            graphs[i]->GetXaxis()->SetTitle("|#eta^{a}|");
            graphs[i]->GetYaxis()->SetTitle("r_{n}");
            graphs[i]->GetYaxis()->SetRangeUser(0.3, 1.8);
            graphs[i]->Draw("APL");
            first = false;
        } else {
            graphs[i]->Draw("PL SAME");
        }
        
        leg->AddEntry(graphs[i], labels[i], "p");
    }
    
    // 添加参考线 y=1
    TLine line(graphs[0]->GetXaxis()->GetXmin(), 1, 
               graphs[0]->GetXaxis()->GetXmax(), 1);
    line.SetLineStyle(2);
    line.SetLineColor(kBlack);
    line.Draw();
    
    leg->Draw();
    
    c.SaveAs(Form("rn_plots/%s.png", name.Data()));
    delete leg;
}

// 处理单个文件的函数
void processFile(const char* filename, const double* Nsel_val, const double* Nsel_bin_width, 
                double* cent_an_f2, double* cent_cn_f2, double* err_an_f2, double* err_cn_f2, int file_index) {
    
    TFile *file_1 = new TFile(filename);
    TFile *file_2 = new TFile(filename);

    if (!file_1 || file_1->IsZombie()) {
        std::cout << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    TH1D *hYHigh = nullptr;
    TH1D *hYLow = nullptr;

    double G_high, G_low;
    double a1_high, a2_high, a3_high, a4_high;
    double a1_low, a2_low, a3_low, a4_low;
    double chi2_high, chi2_low;
    int ndf_high, ndf_low;

    TH1D *hTrigPtLow = (TH1D *)file_1->Get("hTrigPtLow");
    TH1D *hTrigPtHigh = (TH1D *)file_1->Get("hTrigPtHigh");
    TH1D *hNChMid = (TH1D *)file_1->Get("hNChMid");

    if (!hTrigPtLow || !hTrigPtHigh || !hNChMid) {
        std::cout << "Error: Cannot find required histograms in " << filename << std::endl;
        file_1->Close();
        file_2->Close();
        return;
    }

    int nEtaBins = 50;
    double etaMin = -2.5;
    double etaMax = 2.5;
    double etaValues[nEtaBins];
    double etaBinWidths[nEtaBins];
    double a2_Eta_high[50], a2_Eta_err[50];
    double a3_Eta_high[50], a3_Eta_err[50];
    double c2_Eta[50], c2_Eta_err[50];
    double c3_Eta[50], c3_Eta_err[50];

    double Nsel_cut[N_cent + 1] = {10, 30, 40, 50, 60, 80, 100, 120, 150};

    TH3D *hsame_cent[N_cent];
    TH3D *hmixed_cent[N_cent];
    TH1D *hTrigPt_cent[N_cent];
    TH3D *hsame_centlow = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaSameEvent_Cent0");
    TH3D *hmixed_centlow = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaMixEvent_Cent0");

    if (!hsame_centlow || !hmixed_centlow) {
        std::cout << "Error: Cannot find centrality 0 histograms in " << filename << std::endl;
        file_1->Close();
        file_2->Close();
        return;
    }

    for (int icent = 1; icent < N_cent; icent++) {
        double centMin = Nsel_cut[icent];
        double centMax = Nsel_cut[icent + 1];
        
        if (Nsel_val[icent] < 85) {
            if (icent == 0) {
                hTrigPt_cent[icent] = (TH1D *)file_1->Get(Form("hTrigPt_Cent%d", icent));
            } else {
                hsame_cent[icent] = (TH3D *)file_1->Get(Form("hDEtaDPhiTrigEtaSameEvent_Cent%d", icent));
                hmixed_cent[icent] = (TH3D *)file_1->Get(Form("hDEtaDPhiTrigEtaMixEvent_Cent%d", icent));
                hTrigPt_cent[icent] = (TH1D *)file_1->Get(Form("hTrigPt_Cent%d", icent));
            }
        } else {
            if (icent == 0) {
                hTrigPt_cent[icent] = (TH1D *)file_2->Get(Form("hTrigPt_Cent%d", icent));
            } else {
                hsame_cent[icent] = (TH3D *)file_2->Get(Form("hDEtaDPhiTrigEtaSameEvent_Cent%d", icent));
                hmixed_cent[icent] = (TH3D *)file_2->Get(Form("hDEtaDPhiTrigEtaMixEvent_Cent%d", icent));
                hTrigPt_cent[icent] = (TH1D *)file_2->Get(Form("hTrigPt_Cent%d", icent));
            }
        }
        
        if (!hsame_cent[icent] || !hmixed_cent[icent] || !hTrigPt_cent[icent]) {
            std::cout << "Warning: Missing histograms for centrality " << icent << " in " << filename << std::endl;
            continue;
        }
        
        double totalChargedParticleslow_cent = 0;
        double nEventslow_cent = 0;
        for (int bin = 1; bin <= hNChMid->GetNbinsX(); bin++) {
            double binCenter = hNChMid->GetXaxis()->GetBinCenter(bin);
            double binContent = hNChMid->GetBinContent(bin);
            if (binCenter >= 10 && binCenter < 30) {
                totalChargedParticleslow_cent += binCenter * binContent;
                nEventslow_cent += binContent;
            }
        }
        double averageNTlow_cent = (nEventslow_cent > 0) ? totalChargedParticleslow_cent / nEventslow_cent : 0;
        
        double totalChargedParticleshigh_cent = 0;
        double nEventshigh_cent = 0;
        for (int bin = 1; bin <= hNChMid->GetNbinsX(); bin++) {
            double binCenter = hNChMid->GetXaxis()->GetBinCenter(bin);
            double binContent = hNChMid->GetBinContent(bin);
            if (binCenter >= centMin && binCenter < centMax) {
                totalChargedParticleshigh_cent += binCenter * binContent;
                nEventshigh_cent += binContent;
            }
        }
        double averageNThigh_cent = (nEventshigh_cent > 0) ? totalChargedParticleshigh_cent / nEventshigh_cent : 0;
        
        double k = (averageNThigh_cent > 0) ? averageNTlow_cent / averageNThigh_cent : 0;

        for (int iEta = 0; iEta < nEtaBins; iEta++) {
            double eta = etaMin + iEta * (etaMax - etaMin) / nEtaBins;
            etaBinWidths[iEta] = (etaMax - etaMin) / nEtaBins;
            etaValues[iEta] = eta;
            int zBin = iEta + 6;
            
            hsame_centlow->GetZaxis()->SetRange(zBin, zBin);
            hmixed_centlow->GetZaxis()->SetRange(zBin, zBin);
            hsame_cent[icent]->GetZaxis()->SetRange(zBin, zBin);
            hmixed_cent[icent]->GetZaxis()->SetRange(zBin, zBin);

            TH2D *hsame_eta_high = (TH2D *)hsame_cent[icent]->Project3D("yx");
            TH2D *hmix_eta_high = (TH2D *)hmixed_cent[icent]->Project3D("yx");
            TH2D *hsame_eta_low = (TH2D *)hsame_centlow->Project3D("yx");
            TH2D *hmix_eta_low = (TH2D *)hmixed_centlow->Project3D("yx");

            if (hsame_eta_high->GetEntries() > 0 && hmix_eta_high->GetEntries() > 0) {
                hYHigh = getYFromHist(hsame_eta_high, hmix_eta_high, Form("YHigh_eta%d", iEta), hTrigPtHigh->GetEntries());
                hYLow = getYFromHist(hsame_eta_low, hmix_eta_low, Form("YLow_eta%d", iEta), hTrigPtLow->GetEntries());
                
                double a2higherr = 0, a2lowerr = 0;
                double a2_high = getA2(hYHigh, a1_high, a2_high, a3_high, a4_high, G_high, chi2_high, ndf_high, a2higherr);
                double a2_low = getA2(hYLow, a1_low, a2_low, a3_low, a4_low, G_low, chi2_low, ndf_low, a2lowerr);
                
                a2_Eta_high[iEta] = a2_high;
                a2_Eta_err[iEta] = a2higherr;
                
                double c2 = a2_high - a2_low * k;
                c2_Eta[iEta] = c2;
                double c2_err = sqrt(pow(a2higherr, 2) + pow(a2lowerr * k, 2));
                c2_Eta_err[iEta] = c2_err;

                delete hYHigh;
                delete hYLow;
            }
            delete hsame_eta_high;
            delete hmix_eta_high;
            delete hsame_eta_low;
            delete hmix_eta_low;
        }

        const int halfBins = nEtaBins / 2;
        double absEtaValues[halfBins];
        
        for (int k = 0; k < halfBins; k++) {
            absEtaValues[k] = 0.1 * (k + 1);
        }
        
        // 对an_f2进行反对称相除处理
        double r_an_f2[halfBins], r_an_f2_err[halfBins];
        getRnFromVn(a2_Eta_high, a2_Eta_err, r_an_f2, r_an_f2_err, nEtaBins);
        TGraphErrors *gr_an_f2 = new TGraphErrors(halfBins, absEtaValues, r_an_f2, 0, r_an_f2_err);
        
        // 保存an_f2的rn图形到rn_plots目录
        plotRnVsEta(gr_an_f2, Form("an_f2_%s_cent%d", filename, icent), 
                   Form("a_{n} f2 r_{n} vs |#eta^{a}|, cent %d", icent));
        
        // 存储an_f2图形用于后续合并
        if (an_f2_graphs.find(icent) == an_f2_graphs.end()) {
            an_f2_graphs[icent] = std::vector<TGraphErrors*>(4, nullptr);
        }
        an_f2_graphs[icent][file_index] = (TGraphErrors*)gr_an_f2->Clone();
        
        fitRnAndSave(gr_an_f2, Form("an_f2_%s_cent%d", filename, icent), cent_an_f2[icent], err_an_f2[icent]);
        delete gr_an_f2;
        
        // 对cn_f2进行反对称相除处理
        double r_cn_f2[halfBins], r_cn_f2_err[halfBins];
        getRnFromVn(c2_Eta, c2_Eta_err, r_cn_f2, r_cn_f2_err, nEtaBins);
        TGraphErrors *gr_cn_f2 = new TGraphErrors(halfBins, absEtaValues, r_cn_f2, 0, r_cn_f2_err);
        
        // 保存cn_f2的rn图形到rn_plots目录
        plotRnVsEta(gr_cn_f2, Form("cn_f2_%s_cent%d", filename, icent), 
                   Form("c_{n} f2 r_{n} vs |#eta^{a}|, cent %d", icent));
        
        // 存储cn_f2图形用于后续合并
        if (cn_f2_graphs.find(icent) == cn_f2_graphs.end()) {
            cn_f2_graphs[icent] = std::vector<TGraphErrors*>(4, nullptr);
        }
        cn_f2_graphs[icent][file_index] = (TGraphErrors*)gr_cn_f2->Clone();
        
        fitRnAndSave(gr_cn_f2, Form("cn_f2_%s_cent%d", filename, icent), cent_cn_f2[icent], err_cn_f2[icent]);
        delete gr_cn_f2;
    }

    file_1->Close();
    file_2->Close();
    delete file_1;
    delete file_2;
}

void calc_flow_fn_longRange_new_files() {
    double Nsel_cut[N_cent + 1] = {10, 30, 40, 50, 60, 80, 100, 120, 150};
    double Nsel_val[N_cent] = {0};
    double Nsel_bin_width[N_cent] = {0};
    
    // 计算 centrality 值
    for (int icent = 1; icent < N_cent; icent++) {
        Nsel_val[icent] = (Nsel_cut[icent] + Nsel_cut[icent + 1]) / 2.0;
        Nsel_bin_width[icent] = (Nsel_cut[icent + 1] - Nsel_cut[icent]) / 2.0;
    }
    
    // 四个文件的数据存储
    double cent_an_f2_file1[N_cent], cent_cn_f2_file1[N_cent];
    double err_an_f2_file1[N_cent], err_cn_f2_file1[N_cent];
    
    double cent_an_f2_file2[N_cent], cent_cn_f2_file2[N_cent];
    double err_an_f2_file2[N_cent], err_cn_f2_file2[N_cent];
    
    double cent_an_f2_file3[N_cent], cent_cn_f2_file3[N_cent];
    double err_an_f2_file3[N_cent], err_cn_f2_file3[N_cent];
    
    double cent_an_f2_file4[N_cent], cent_cn_f2_file4[N_cent];
    double err_an_f2_file4[N_cent], err_cn_f2_file4[N_cent];
    
    // 初始化数组
    for (int i = 0; i < N_cent; i++) {
        cent_an_f2_file1[i] = 0; cent_cn_f2_file1[i] = 0;
        err_an_f2_file1[i] = 0; err_cn_f2_file1[i] = 0;
        cent_an_f2_file2[i] = 0; cent_cn_f2_file2[i] = 0;
        err_an_f2_file2[i] = 0; err_cn_f2_file2[i] = 0;
        cent_an_f2_file3[i] = 0; cent_cn_f2_file3[i] = 0;
        err_an_f2_file3[i] = 0; err_cn_f2_file3[i] = 0;
        cent_an_f2_file4[i] = 0; cent_cn_f2_file4[i] = 0;
        err_an_f2_file4[i] = 0; err_cn_f2_file4[i] = 0;
    }
    
    // 存储文件名
    file_names.push_back("1.5mb");
    file_names.push_back("0.15mb");
    file_names.push_back("1.5mb a0.8 b0.4");
    file_names.push_back("0.15mb a0.8 b0.4");
    
    // 处理四个文件
    std::cout << "Processing file 1: hist_ampt_normal_1.5mb_Decorr_yuhao.root" << std::endl;
    processFile("hist_ampt_normal_1.5mb_Decorr_yuhao.root", Nsel_val, Nsel_bin_width,
               cent_an_f2_file1, cent_cn_f2_file1, err_an_f2_file1, err_cn_f2_file1, 0);
    
    std::cout << "Processing file 2: hist_ampt_normal_0.15mb_Decorr_yuhao.root" << std::endl;
    processFile("hist_ampt_normal_0.15mb_Decorr_yuhao.root", Nsel_val, Nsel_bin_width,
               cent_an_f2_file2, cent_cn_f2_file2, err_an_f2_file2, err_cn_f2_file2, 1);
    
    std::cout << "Processing file 3: hist_ampt_normal_1.5mb_a_0.8_b_0.4_Decorr_yuhao.root" << std::endl;
    processFile("hist_ampt_normal_1.5mb_a_0.8_b_0.4_Decorr_yuhao.root", Nsel_val, Nsel_bin_width,
               cent_an_f2_file3, cent_cn_f2_file3, err_an_f2_file3, err_cn_f2_file3, 2);
               
    std::cout << "Processing file 4: hist_ampt_normal_0.15mb_a_0.8_b_0.4_Decorr_yuhao.root" << std::endl;
    processFile("hist_ampt_normal_0.15mb_a_0.8_b_0.4_Decorr_yuhao.root", Nsel_val, Nsel_bin_width,
               cent_an_f2_file4, cent_cn_f2_file4, err_an_f2_file4, err_cn_f2_file4, 3);
    // 绘制合并的rn图形
    for (int icent = 1; icent < N_cent; icent++) {
        if (an_f2_graphs.find(icent) != an_f2_graphs.end()) {
            plotCombinedRn(an_f2_graphs[icent], 
                          Form("an_f2_combined_cent%d", icent),
                          Form("a_{n} f2 r_{n} vs |#eta^{a}|, cent %d", icent),
                          file_names);
        }
        
        if (cn_f2_graphs.find(icent) != cn_f2_graphs.end()) {
            plotCombinedRn(cn_f2_graphs[icent], 
                          Form("cn_f2_combined_cent%d", icent),
                          Form("c_{n} f2 r_{n} vs |#eta^{a}|, cent %d", icent),
                          file_names);
        }
    }

    // 创建并绘制综合图表
TCanvas *cCombined = new TCanvas("cCombined", "Fn vs Nch - Four Files Comparison", 1200, 800);
cCombined->SetLogy();

// 创建八个图表对象（每个文件两条线）
TGraphErrors *gr_f2_raw_file1 = new TGraphErrors(N_cent, Nsel_val, cent_an_f2_file1, Nsel_bin_width, err_an_f2_file1);
TGraphErrors *gr_f2_sub_file1 = new TGraphErrors(N_cent, Nsel_val, cent_cn_f2_file1, Nsel_bin_width, err_cn_f2_file1);

TGraphErrors *gr_f2_raw_file2 = new TGraphErrors(N_cent, Nsel_val, cent_an_f2_file2, Nsel_bin_width, err_an_f2_file2);
TGraphErrors *gr_f2_sub_file2 = new TGraphErrors(N_cent, Nsel_val, cent_cn_f2_file2, Nsel_bin_width, err_cn_f2_file2);

TGraphErrors *gr_f2_raw_file3 = new TGraphErrors(N_cent, Nsel_val, cent_an_f2_file3, Nsel_bin_width, err_an_f2_file3);
TGraphErrors *gr_f2_sub_file3 = new TGraphErrors(N_cent, Nsel_val, cent_cn_f2_file3, Nsel_bin_width, err_cn_f2_file3);

TGraphErrors *gr_f2_raw_file4 = new TGraphErrors(N_cent, Nsel_val, cent_an_f2_file4, Nsel_bin_width, err_an_f2_file4);
TGraphErrors *gr_f2_sub_file4 = new TGraphErrors(N_cent, Nsel_val, cent_cn_f2_file4, Nsel_bin_width, err_cn_f2_file4);

// 设置样式 - File 1 (1.5mb)
gr_f2_raw_file1->SetMarkerStyle(20); gr_f2_raw_file1->SetMarkerColor(kRed); gr_f2_raw_file1->SetLineColor(kRed);
gr_f2_sub_file1->SetMarkerStyle(24); gr_f2_sub_file1->SetMarkerColor(kRed); gr_f2_sub_file1->SetLineColor(kRed);

// 设置样式 - File 2 (0.15mb)
gr_f2_raw_file2->SetMarkerStyle(21); gr_f2_raw_file2->SetMarkerColor(kBlue); gr_f2_raw_file2->SetLineColor(kBlue);
gr_f2_sub_file2->SetMarkerStyle(25); gr_f2_sub_file2->SetMarkerColor(kBlue); gr_f2_sub_file2->SetLineColor(kBlue);

// 设置样式 - File 3 (1.5mb_a_0.8_b_0.4)
gr_f2_raw_file3->SetMarkerStyle(22); gr_f2_raw_file3->SetMarkerColor(kGreen+2); gr_f2_raw_file3->SetLineColor(kGreen+2);
gr_f2_sub_file3->SetMarkerStyle(26); gr_f2_sub_file3->SetMarkerColor(kGreen+2); gr_f2_sub_file3->SetLineColor(kGreen+2);

// 设置样式 - File 4 (0.15mb_a_0.8_b_0.4)
gr_f2_raw_file4->SetMarkerStyle(23); gr_f2_raw_file4->SetMarkerColor(kMagenta); gr_f2_raw_file4->SetLineColor(kMagenta);
gr_f2_sub_file4->SetMarkerStyle(27); gr_f2_sub_file4->SetMarkerColor(kMagenta); gr_f2_sub_file4->SetLineColor(kMagenta);

// 绘制所有图表
gr_f2_raw_file1->Draw("AP");
gr_f2_sub_file1->Draw("P SAME");
gr_f2_raw_file2->Draw("P SAME");
gr_f2_sub_file2->Draw("P SAME");
gr_f2_raw_file3->Draw("P SAME");
gr_f2_sub_file3->Draw("P SAME");
gr_f2_raw_file4->Draw("P SAME");
gr_f2_sub_file4->Draw("P SAME");

// 设置坐标轴标题和范围
gr_f2_raw_file1->GetXaxis()->SetTitle("N_{ch}");
gr_f2_raw_file1->GetYaxis()->SetTitle("F_{n}");
gr_f2_raw_file1->SetTitle("Long-range flow coefficients - Four Files Comparison");
gr_f2_raw_file1->GetYaxis()->SetRangeUser(-0.01, 0.8);

// 创建图例
TLegend *leg = new TLegend(0.15, 0.15, 0.4, 0.5);
leg->AddEntry(gr_f2_raw_file1, "f2^{raw} (1.5mb)", "p");
leg->AddEntry(gr_f2_sub_file1, "f2^{sub} (1.5mb)", "p");
leg->AddEntry(gr_f2_raw_file2, "f2^{raw} (0.15mb)", "p");
leg->AddEntry(gr_f2_sub_file2, "f2^{sub} (0.15mb)", "p");
leg->AddEntry(gr_f2_raw_file3, "f2^{raw} (1.5mb_a0.8_b0.4)", "p");
leg->AddEntry(gr_f2_sub_file3, "f2^{sub} (1.5mb_a0.8_b0.4)", "p");
leg->AddEntry(gr_f2_raw_file4, "f2^{raw} (0.15mb_a0.8_b0.4)", "p");
leg->AddEntry(gr_f2_sub_file4, "f2^{sub} (0.15mb_a0.8_b0.4)", "p");
leg->Draw();

// 保存图片
cCombined->SaveAs("Fn_vs_Nch_FourFiles.png");

// 保存数据到ROOT文件
TFile *outfile = new TFile("longRange_flow_FourFiles.root", "recreate");

// File 1 数据
gr_f2_raw_file1->Write("gr_f2_raw_1p5mb");
gr_f2_sub_file1->Write("gr_f2_sub_1p5mb");

// File 2 数据
gr_f2_raw_file2->Write("gr_f2_raw_0p15mb");
gr_f2_sub_file2->Write("gr_f2_sub_0p15mb");

// File 3 数据
gr_f2_raw_file3->Write("gr_f2_raw_1p5mb_a0p8_b0p4");
gr_f2_sub_file3->Write("gr_f2_sub_1p5mb_a0p8_b0p4");

// File 4 数据
gr_f2_raw_file4->Write("gr_f2_raw_0p15mb_a0p8_b0p4");
gr_f2_sub_file4->Write("gr_f2_sub_0p15mb_a0p8_b0p4");

outfile->Close();

std::cout << "Analysis completed. Results saved to:" << std::endl;
std::cout << "  - Fn_vs_Nch_FourFiles.png" << std::endl;
std::cout << "  - longRange_flow_FourFiles.root" << std::endl;
std::cout << "  - rn_plots/ directory contains individual and combined r_n plots" << std::endl;
    
// 清理内存
for (auto& pair : an_f2_graphs) {
    for (auto graph : pair.second) {
        if (graph) delete graph;
    }
}
for (auto& pair : cn_f2_graphs) {
    for (auto graph : pair.second) {
        if (graph) delete graph;
    }
}

// 删除图表对象
delete gr_f2_raw_file1; delete gr_f2_sub_file1;
delete gr_f2_raw_file2; delete gr_f2_sub_file2;
delete gr_f2_raw_file3; delete gr_f2_sub_file3;
delete gr_f2_raw_file4; delete gr_f2_sub_file4;
delete cCombined;
}