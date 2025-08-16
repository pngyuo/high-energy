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

const double PI = 3.1415926;
double gYIntegral = 0;

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

void calc_flow_fn_longRange_new() {
    TFile *file_1 = new TFile("hist_ampt_3q_1.5mb_Decorr_yuhao.root");
    TFile *file_2 = new TFile("hist_ampt_3q_1.5mb_Decorr_yuhao.root");

    TH1D *hYHigh = nullptr;
    TH1D *hYLow = nullptr;

    double G_high, G_low;
    double a1_high, a2_high, a3_high, a4_high;
    double a1_low, a2_low, a3_low, a4_low;
    double chi2_high, chi2_low;
    int ndf_high, ndf_low;
    double totalChargedParticleslow, totalChargedParticleshigh;
    double nEventslow, nEventshigh;
    double averageNThigh, averageNTlow;

    TH1D *hTrigPtLow = (TH1D *)file_1->Get("hTrigPtLow");
    TH1D *hTrigPtHigh = (TH1D *)file_1->Get("hTrigPtHigh");

    TH1D *hNChMid = (TH1D *)file_1->Get("hNChMid");

    int nEtaBins = 50;
    double etaMin = -2.5;
    double etaMax = 2.5;
    double etaValues[nEtaBins];
    double c2_Eta[50];
    double c2_Eta_err[50];
    double a2_Eta_high[50];
    double a2_Eta_err[50];
    double a3_Eta_high[50];      // 存储a3原始值
    double a3_Eta_err[50];       // 存储a3误差
    double c3_Eta[50];           // 存储修正后的c3值
    double c3_Eta_err[50];       // 存储c3误差
    double etaBinWidths[nEtaBins];

    const int N_cent = 8;
    double Nsel_cut[N_cent + 1] = {10, 30, 40, 50, 60, 80, 100, 120, 150}; // 修复数组初始化
    double Nsel_val[N_cent] = {0};
    double Nsel_bin_width[N_cent] = {0}; // 定义为数组

    TH3D *hsame_cent[N_cent];  // 修复类型声明
    TH3D *hmixed_cent[N_cent];
    TH1D *hTrigPt_cent[N_cent];
    double cent_an_f2[N_cent] = {0}, cent_an_f3[N_cent] = {0};
    double cent_cn_f2[N_cent] = {0}, cent_cn_f3[N_cent] = {0};
    double err_an_f2[N_cent] = {0}, err_an_f3[N_cent] = {0};
    double err_cn_f2[N_cent] = {0}, err_cn_f3[N_cent] = {0};
    TH3D *hsame_centlow = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaSameEvent_Cent0");
    TH3D *hmixed_centlow = (TH3D *)file_1->Get("hDEtaDPhiTrigEtaMixEvent_Cent0");

    TH1D *hRefPtCount_low = (TH1D*)file_1->Get("hRefPtCount0");
    double totalChargedParticleslow_cent = 0;
    double nEventslow_cent = 0;
    for (int bin = 1; bin <= hRefPtCount_low->GetNbinsX(); bin++) {
        double binCenter = hRefPtCount_low->GetXaxis()->GetBinCenter(bin);
        double binContent = hRefPtCount_low->GetBinContent(bin);
            totalChargedParticleslow_cent += binCenter * binContent;
            nEventslow_cent += binContent;
    }
    double averageNTlow_cent = totalChargedParticleslow_cent / nEventslow_cent;
    std::cout << "averageNTlow_cent: " << averageNTlow_cent << std::endl;

    for (int icent = 1; icent < N_cent; icent++) {
        // 动态计算每个 centrality bin 的中心值和宽度
        Nsel_val[icent] = (Nsel_cut[icent] + Nsel_cut[icent + 1]) / 2.0;
        Nsel_bin_width[icent] = (Nsel_cut[icent + 1] - Nsel_cut[icent]) / 2.0;
            double centMin = Nsel_cut[icent];
            double centMax = Nsel_cut[icent + 1];
        // 从文件中获取对应的 hist
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
        
            TH1D *hRefPtCount_cent = (TH1D*)file_1->Get(Form("hRefPtCount%d", icent));
            // 计算 averageNThigh（基于当前 cent）
            double totalChargedParticleshigh_cent = 0;
            double nEventshigh_cent = 0;
            for (int bin = 1; bin <= hRefPtCount_cent->GetNbinsX(); bin++) {
                double binCenter = hRefPtCount_cent->GetXaxis()->GetBinCenter(bin);
                double binContent = hRefPtCount_cent->GetBinContent(bin);
                    totalChargedParticleshigh_cent += binCenter * binContent;
                    nEventshigh_cent += binContent;
            }
            double averageNThigh_cent = totalChargedParticleshigh_cent / nEventshigh_cent;
            std::cout << "averageNThigh_cent: " << averageNThigh_cent << std::endl;
        
            // 计算 k，使用当前 cent 的 averageNThigh_cent 和固定的 averageNTlow_cent
            double k = averageNTlow_cent / averageNThigh_cent;
            std::cout << "k = " << k << std::endl;

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

                double a3higherr = 0, a3lowerr = 0;
                double a3_high = getA3(hYHigh, a1_high, a2_high, a3_high, a4_high, G_high, chi2_high, ndf_high, a3higherr);
                double a3_low = getA3(hYLow, a1_low, a2_low, a3_low, a4_low, G_low, chi2_low, ndf_low, a3lowerr);
                
                a3_Eta_high[iEta] = a3_high;
                a3_Eta_err[iEta] = a3higherr;
                
                double c3 = a3_high - a3_low * k;
                c3_Eta[iEta] = c3;
                double c3_err = sqrt(pow(a3higherr, 2) + pow(a3lowerr * k, 2));
                c3_Eta_err[iEta] = c3_err;
                
                delete hYHigh;
                delete hYLow;
            }
            delete hsame_eta_high;
            delete hmix_eta_high;
            delete hsame_eta_low;
            delete hmix_eta_low;
        }

        // 准备反对称相除的数据
        const int halfBins = nEtaBins / 2;
        double absEtaValues[halfBins];  // |eta| 值
        
        // 计算 |eta| 值 (0.1, 0.2, ..., 2.5)
        for (int k = 0; k < halfBins; k++) {
            absEtaValues[k] = 0.1 * (k + 1);
        }
        
        // 对an_f2进行反对称相除处理
        double r_an_f2[halfBins], r_an_f2_err[halfBins];
        getRnFromVn(a2_Eta_high, a2_Eta_err, r_an_f2, r_an_f2_err, nEtaBins);
        TGraphErrors *gr_an_f2 = new TGraphErrors(halfBins, absEtaValues, r_an_f2, 0, r_an_f2_err);
        fitRnAndSave(gr_an_f2, Form("an_f2_cent%d", icent), cent_an_f2[icent], err_an_f2[icent]);
        delete gr_an_f2;
        
        // 对an_f3进行反对称相除处理
        double r_an_f3[halfBins], r_an_f3_err[halfBins];
        getRnFromVn(a3_Eta_high, a3_Eta_err, r_an_f3, r_an_f3_err, nEtaBins);
        TGraphErrors *gr_an_f3 = new TGraphErrors(halfBins, absEtaValues, r_an_f3, 0, r_an_f3_err);
        fitRnAndSave(gr_an_f3, Form("an_f3_cent%d", icent), cent_an_f3[icent], err_an_f3[icent]);
        delete gr_an_f3;
        
        // 对cn_f2进行反对称相除处理
        double r_cn_f2[halfBins], r_cn_f2_err[halfBins];
        getRnFromVn(c2_Eta, c2_Eta_err, r_cn_f2, r_cn_f2_err, nEtaBins);
        TGraphErrors *gr_cn_f2 = new TGraphErrors(halfBins, absEtaValues, r_cn_f2, 0, r_cn_f2_err);
        fitRnAndSave(gr_cn_f2, Form("cn_f2_cent%d", icent), cent_cn_f2[icent], err_cn_f2[icent]);
        delete gr_cn_f2;
        
        // 对cn_f3进行反对称相除处理
        double r_cn_f3[halfBins], r_cn_f3_err[halfBins];
        getRnFromVn(c3_Eta, c3_Eta_err, r_cn_f3, r_cn_f3_err, nEtaBins);
        TGraphErrors *gr_cn_f3 = new TGraphErrors(halfBins, absEtaValues, r_cn_f3, 0, r_cn_f3_err);
        fitRnAndSave(gr_cn_f3, Form("cn_f3_cent%d", icent), cent_cn_f3[icent], err_cn_f3[icent]);
        delete gr_cn_f3;
    }

    // 创建并绘制四个数据集
    TCanvas *cCombined = new TCanvas("cCombined", "Fn vs Nch", 800, 600);
    cCombined->SetLogy();
    TGraphErrors *gr_f2_raw = new TGraphErrors(N_cent, Nsel_val, cent_an_f2, Nsel_bin_width, err_an_f2);
    TGraphErrors *gr_f2_sub = new TGraphErrors(N_cent, Nsel_val, cent_cn_f2, Nsel_bin_width, err_cn_f2); 
    TGraphErrors *gr_f3_raw = new TGraphErrors(N_cent-1, &Nsel_val[1], &cent_an_f3[1], 0, &err_an_f3[1]);
    TGraphErrors *gr_f3_sub = new TGraphErrors(N_cent-1, &Nsel_val[1], &cent_cn_f3[1], 0, &err_cn_f3[1]);

    gr_f2_raw->SetMarkerStyle(20); gr_f2_raw->SetMarkerColor(kRed);
    gr_f3_raw->SetMarkerStyle(21); gr_f3_raw->SetMarkerColor(kBlue);
    gr_f2_sub->SetMarkerStyle(22); gr_f2_sub->SetMarkerColor(kGreen+2);
    gr_f3_sub->SetMarkerStyle(23); gr_f3_sub->SetMarkerColor(kMagenta);

    gr_f2_raw->Draw("AP");
    gr_f3_raw->Draw("P SAME");
    gr_f2_sub->Draw("P SAME");
    gr_f3_sub->Draw("P SAME");

    gr_f2_raw->GetXaxis()->SetTitle("N_{ch}");
    gr_f2_raw->GetYaxis()->SetTitle("F_{n}");
    gr_f2_raw->SetTitle("Long-range flow coefficients");
    gr_f2_raw->GetYaxis()->SetRangeUser(-0.01, 0.8);

    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(gr_f2_raw, "f2^{raw}", "p");
    leg->AddEntry(gr_f3_raw, "f3^{raw}", "p");
    leg->AddEntry(gr_f2_sub, "f2^{sub}", "p");
    leg->AddEntry(gr_f3_sub, "f3^{sub}", "p");
    leg->Draw();

    cCombined->SaveAs("Fn_vs_Nch.png");

    TFile *outfile = new TFile("longRange_flow.root", "recreate");
    gr_f2_raw->Write("gr_f2_raw");
    gr_f3_raw->Write("gr_f3_raw");
    gr_f2_sub->Write("gr_f2_sub");
    gr_f3_sub->Write("gr_f3_sub");
    outfile->Close();
}