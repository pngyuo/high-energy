#include <iostream>
#include <vector>
#include <utility>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <cmath>
#include <TF1.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <map>
#include <TRandom3.h> 

std::vector<double> globalRtLowerEdges = {0, 0.25, 0.5, 0.7, 1.1,1.5, 2, 2.8};
TF1 *fitfun = new TF1("fitfun", "x*[0]*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt(x*x+[3]*[3])-[3])/[1]/[2], -[1])", 0, 20);

void addText(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.04);
    latex->Draw();
}

void addText2(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.08);
    latex->Draw();
}

void addText3(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.06);
    latex->Draw();
}

void addText4(Double_t x, Double_t y, const char* text) {
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(0.15);
    latex->Draw();
}

// 获取粒子类型字符串
TString getParticleTypeFromHistName(const TString& histname) {
    if (histname.Contains("PiCh")) return "Pi";
    if (histname.Contains("KCh")) return "K";
    if (histname.Contains("Proton")) return "Proton";
    return "Unknown";
}

// 获取dPhi区域字符串
TString getRegionFromHistName(const TString& histname) {
    if (histname.Contains("dPhi0")) return "Toward";
    if (histname.Contains("dPhi1")) return "Transverse";
    return "Unknown";
}

// 修改函数以使用全局RT范围
std::pair<std::vector<TH1D*>, double> getCentPtProjections(const TString& filename, const TString& histname) {
    TFile* infile = new TFile(filename, "READ");
    TH3D* hist3D = (TH3D*)infile->Get(histname);
    TH1D* chargedParticlesCount = (TH1D*)infile->Get("nTHist");

    double totalChargedParticles = 0;
    double nEvents = 0;
    for (int bin = 1; bin <= chargedParticlesCount->GetNbinsX(); bin++) {
        double binCenter = chargedParticlesCount->GetXaxis()->GetBinCenter(bin);
        double binContent = chargedParticlesCount->GetBinContent(bin);
        totalChargedParticles += binCenter * binContent;
        nEvents += binContent;
    }
    double averageNT = totalChargedParticles / nEvents;

    std::vector<TH1D*> projections;
    int nBins = globalRtLowerEdges.size();

    for (int i = 0; i < nBins; i++) {
        double rtLow = globalRtLowerEdges[i];
        double rtUp;
        if (i < nBins - 1) {
            rtUp = globalRtLowerEdges[i+1];
        } else {
            rtUp = 30.0 / averageNT;
        }

        double ntRangeLow = rtLow * averageNT;
        double ntRangeUp = rtUp * averageNT;

        int binLow = hist3D->GetZaxis()->FindBin(ntRangeLow);
        int binUp = hist3D->GetZaxis()->FindBin(ntRangeUp);
        binUp -= 1;

        TH1D* histProjX = hist3D->ProjectionX(Form("%s_histProjX_%i", histname.Data(), i + 1), 1, hist3D->GetNbinsY(), binLow, binUp);

        double totalEventsInRange = 0;
        for (int binIndex = binLow; binIndex <= binUp; ++binIndex) {
            totalEventsInRange += chargedParticlesCount->GetBinContent(binIndex);
        }

        histProjX->Scale(1.0 / totalEventsInRange, "width");
        TH1D* histProjXTemp = (TH1D*)histProjX->Clone(Form("%s_temp_%i", histname.Data(), i + 1));
        histProjXTemp->GetXaxis()->SetRangeUser(0, 5.0);
        projections.push_back(histProjXTemp);
    }
    return std::make_pair(projections, averageNT);
}

TGraphErrors* drawAveragePt(const std::vector<double>& avgPt, const std::vector<double>& avgPtErr, 
    double averageNT, Color_t color, Style_t markerStyle, double maxRT = 100.0) {
    int nBins = globalRtLowerEdges.size();
    std::vector<double> rtCenters;

    for (int i = 0; i < nBins; i++) {
        double rtLow = globalRtLowerEdges[i];
        double rtUp;
        if (i < nBins - 1) {
            rtUp = globalRtLowerEdges[i+1];
        } else {
            rtUp = 30.0 / averageNT;
        }
        rtCenters.push_back((rtLow + rtUp) / 2.0);
    }

    int count = 0;
    for (int i = 0; i < avgPt.size(); ++i) {
        if (rtCenters[i] < maxRT) {
            count++;
        }
    }

    double* x = new double[count];
    double* y = new double[count];
    double* ex = new double[count];
    double* ey = new double[count];

    int index = 0;
    for (int i = 0; i < avgPt.size(); ++i) {
        if (rtCenters[i] < maxRT) {
            x[index] = rtCenters[i];
            y[index] = avgPt[i];
            ex[index] = 0.0;
            ey[index] = avgPtErr[i];
            index++;
        }
    }

    TGraphErrors* gr = new TGraphErrors(count, x, y, ex, ey);
    gr->SetTitle("");
    gr->SetMarkerStyle(markerStyle);
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);
    gr->SetLineWidth(2);

    delete[] x;
    delete[] y;
    delete[] ex;
    delete[] ey;

    return gr;
}

// 获取粒子质量
double getMassFromHistName(const TString& histname) {
    if (histname.Contains("hPiCh")) return 0.134;
    if (histname.Contains("hKCh")) return 0.493;
    if (histname.Contains("hProton")) return 0.938;
    return 0.0;
}

TGraphErrors* analyzeAndDrawpi(const std::vector<TH1D*>& projections, TVirtualPad* pad, 
    double averageNT, int modelIndex, Color_t color, Style_t markerStyle, 
    double maxRT, double mass) {
    pad->cd();
    std::vector<double> avgPt;
    std::vector<double> avgPtErr;
    int i = 0;
    
    TString modelName;
    if (modelIndex == 0) modelName = "noFSI";
    else if (modelIndex == 1) modelName = "nohFSI";
    else if (modelIndex == 2) modelName = "allFSI";
    else modelName = "Unknown";
    
    for (const auto& h : projections) {
        double meanPt = 0.0;
        double chi2 = 0.0;
        int ndf = 0;
        int status = -1;
        double par[4] = {0};
        
        // 创建唯一文件名
        TString particleType = getParticleTypeFromHistName(h->GetName());
        TString region = getRegionFromHistName(h->GetName());
        TString canvasName = Form("fit_results/%s_%s_%s_RTbin%d.png", modelName.Data(), particleType.Data(), region.Data(), i);
        
        // 创建临时画布用于拟合
        // TCanvas c(Form("c_%s", h->GetName()), "Fit Result", 800, 600);
        // h->Draw();
        fitfun->SetParameters(1.1,300, 0.1, mass);
        TFitResultPtr fitResult = h->Fit("fitfun", "SQ0");
        status = fitResult->Status();
        meanPt = fitfun->Mean(0, 5);
        chi2 = fitfun->GetChisquare();
        ndf = fitfun->GetNDF();
        for (int iParam = 0; iParam < 4; iParam++) {
            par[iParam] = fitfun->GetParameter(iParam);
        }
        double avgPtinjetErr = 0.0;
        const TMatrixDSym& covMatrix = fitResult->GetCovarianceMatrix();
        double integralError = fitfun->IntegralError(0, 5, fitResult->GetParams(), covMatrix.GetMatrixArray());
        avgPtinjetErr = integralError / (5.0 - 0.0);  // 归一化到积分范围
        avgPt.push_back(meanPt);
        avgPtErr.push_back(avgPtinjetErr);
        
        // 添加详细的拟合信息
        // TLatex text;
        // text.SetTextSize(0.04);
        // text.SetNDC(kTRUE);
        // text.DrawLatex(0.25, 0.85, Form("Model: %s", modelName.Data()));
        // text.DrawLatex(0.25, 0.80, Form("Particle: %s", particleType.Data()));
        // text.DrawLatex(0.25, 0.75, Form("Region: %s", region.Data()));
        // text.DrawLatex(0.25, 0.70, Form("RT bin: %d", i));
        // text.DrawLatex(0.25, 0.65, Form("p0 = %.4f", par[0]));
        // text.DrawLatex(0.25, 0.60, Form("p1 = %.4f", par[1]));
        // text.DrawLatex(0.25, 0.55, Form("p2 = %.4f", par[2]));
        // text.DrawLatex(0.25, 0.50, Form("p3 = %.4f", par[3]));
        // text.DrawLatex(0.25, 0.45, Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf));
        // text.DrawLatex(0.25, 0.40, Form("Mean p_{T}: %.4f #pm %.4f", meanPt, avgPtinjetErr));
        // text.DrawLatex(0.25, 0.35, Form("Fit Status: %d", status));
        
        // if (i < globalRtLowerEdges.size() - 1) {
        //     text.DrawLatex(0.25, 0.30, Form("RT: [%.2f, %.2f]", globalRtLowerEdges[i], globalRtLowerEdges[i+1]));
        // } else {
        //     text.DrawLatex(0.25, 0.30, Form("RT: [%.2f, MAX]", globalRtLowerEdges[i]));
        // }

        // c.SaveAs(canvasName);
        // i++;
    }
    
    if (!avgPt.empty()) {
        return drawAveragePt(avgPt, avgPtErr, averageNT, color, markerStyle, maxRT);
    }
    return nullptr;
}

TGraphErrors*  drawInJetRegion(const std::vector<TH1D*>& projectionsToward, 
    const std::vector<TH1D*>& projectionsTransverse, 
    double averageNT, 
    int modelIndex,
    Color_t color, 
    Style_t markerStyle) {
    std::vector<double> avgPtDiff;
    std::vector<double> avgPtDiffErr;
    std::vector<double> rtCenters;
    TFitResultPtr fitResult; // 声明fitResult变量
    std::vector<double> avgPt;
    std::vector<double> avgPtErr;
    double meanPt;
    double avgPtinjetErr;
    
    TString modelName;
    if (modelIndex == 0) modelName = "noFSI";
    else if (modelIndex == 1) modelName = "nohFSI";
    else if (modelIndex == 2) modelName = "allFSI";
    else modelName = "Unknown";

    for (size_t j = 0; j < projectionsToward.size(); ++j) {
        TH1D* histToward = projectionsToward[j];
        TH1D* histTransverse = projectionsTransverse[j];
        TH1D* diffHist = (TH1D*)histToward->Clone(Form("diffHist_%zu", j));
        diffHist->Add(histTransverse, -1);

        double mass = getMassFromHistName(histTransverse->GetName());
        fitfun->SetParameters(1.1,500, 0.05, mass);
        
        // 为差谱创建详细拟合图
        TString particleType = getParticleTypeFromHistName(histTransverse->GetName());
        TString canvasName = Form("fit_results/diff_%s_%s_RTbin%d.png", modelName.Data(), particleType.Data(), j);
        
        // TCanvas c(Form("c_diff_%s_%d", particleType.Data(), j), "Diff Fit Result", 800, 600);
        // diffHist->Draw();
        fitResult = diffHist->Fit("fitfun", "SQ0");
        meanPt = fitfun->Mean(0,5);
        
        // 添加详细的拟合信息
        // TLatex text;
        // text.SetTextSize(0.04);
        // text.SetNDC(kTRUE);
        // text.DrawLatex(0.25, 0.85, Form("Model: %s", modelName.Data()));
        // text.DrawLatex(0.25, 0.80, Form("Particle: %s", particleType.Data()));
        // text.DrawLatex(0.25, 0.75, Form("RT bin: %d", j));
        // text.DrawLatex(0.25, 0.70, "Region: In-Jet (Toward - Transverse)");
        // text.DrawLatex(0.25, 0.65, Form("p0 = %.4f", fitfun->GetParameter(0)));
        // text.DrawLatex(0.25, 0.60, Form("p1 = %.4f", fitfun->GetParameter(1)));
        // text.DrawLatex(0.25, 0.55, Form("p2 = %.4f", fitfun->GetParameter(2)));
        // text.DrawLatex(0.25, 0.50, Form("p3 = %.4f", fitfun->GetParameter(3)));
        // text.DrawLatex(0.25, 0.45, Form("#chi^{2}/ndf = %.2f/%d", fitfun->GetChisquare(), fitfun->GetNDF()));
        // text.DrawLatex(0.25, 0.40, Form("Mean p_{T}: %.4f", meanPt));
        // text.DrawLatex(0.25, 0.35, Form("Fit Status: %d", fitResult->Status()));
        
        // if (j < globalRtLowerEdges.size() - 1) {
        //     text.DrawLatex(0.25, 0.30, Form("RT: [%.2f, %.2f]", globalRtLowerEdges[j], globalRtLowerEdges[j+1]));
        // } else {
        //     text.DrawLatex(0.25, 0.30, Form("RT: [%.2f, MAX]", globalRtLowerEdges[j]));
        // }
        
        // c.SaveAs(canvasName);
        
        const TMatrixDSym& covMatrix = fitResult->GetCovarianceMatrix();
        double integralError = fitfun->IntegralError(0, 5, fitResult->GetParams(), covMatrix.GetMatrixArray());
        avgPtinjetErr = integralError / (5.0 - 0.0);  // 归一化到积分范围
        avgPt.push_back(meanPt);
        avgPtErr.push_back(avgPtinjetErr);

        avgPtDiff.push_back(meanPt);
        avgPtDiffErr.push_back(avgPtinjetErr);

        // 计算RT中心点
        double rtLow = globalRtLowerEdges[j];
        double rtUp = (j < globalRtLowerEdges.size() - 1) ? globalRtLowerEdges[j+1] : 30.0 / averageNT;
        rtCenters.push_back((rtLow + rtUp) / 2.0);
    }

    std::vector<double> filteredRtCenters;
    std::vector<double> filteredAvgPtDiff;
    std::vector<double> filteredAvgPtDiffErr;

    for (size_t i = 0; i < rtCenters.size(); i++) {
        if (rtCenters[i] < 1.2) {  // 只保留RT<1.2的点
            filteredRtCenters.push_back(rtCenters[i]);
            filteredAvgPtDiff.push_back(avgPtDiff[i]);
            filteredAvgPtDiffErr.push_back(avgPtDiffErr[i]);
        }
    }

    // 创建TGraphErrors并绘制（使用过滤后的数据）
    if (!filteredRtCenters.empty()) {
        TGraphErrors* graph = new TGraphErrors(filteredRtCenters.size(), 
                             filteredRtCenters.data(), 
                             filteredAvgPtDiff.data(), 
                             nullptr, 
                             filteredAvgPtDiffErr.data());
        graph->SetMarkerStyle(markerStyle);
        graph->SetMarkerColor(color);
        graph->SetLineColor(color);
        graph->SetLineWidth(2);
        graph->Draw("LP");
        return graph; 
    }
   return nullptr; 
}


void drawPi() {
    // 设置全局RT范围
    globalRtLowerEdges = {0, 0.25, 0.5, 0.7, 1.1,1.5, 2, 2.5, 2.8};
    
    // 确保输出目录存在
    gSystem->mkdir("fit_results", kTRUE);
    
    const TString filenames[3] = {
        "hist_outputnohFSI_01mb.root",
        "hist_outputnohFSI_new.root",
        "hist_outputnohFSI_02mb.root"
    };

    TCanvas *c1 = new TCanvas("c1", "c1", 1500, 550);
    c1->SetTicks(1, 1);
    c1->SetMargin(0.17, 0.10, 0.2, 0.2);
    c1->Divide(3, 1, 0, 0);
    gStyle->SetOptStat(kFALSE);
    
    // 通用样式设置
    std::vector<Color_t> colors = {kRed, kBlue, kMagenta};
    std::vector<Style_t> markerStyles = {20, 21, 22};

    TString filenames1[2]= {"pi_toward.root",  "pi_transverse.root"};
    TString filenames2[2]= {"K_toward.root", "K_transverse.root"};
    TString filenames3[2]= {"proton_toward.root", "proton_transverse.root"};
            // 创建输出ROOT文件
            TFile *outFile = new TFile("pi_results.root", "RECREATE");
    
            // 创建目录结构
            outFile->mkdir("Toward");
            outFile->mkdir("Transverse");
            outFile->mkdir("InJet");
        
            // 存储图形的容器
            std::vector<TGraphErrors*> graphsToward;
            std::vector<TGraphErrors*> graphsTransverse;
            std::vector<TGraphErrors*> graphsInJet;
    // 第一子图: Toward区域(π)
    c1->cd(1);
    {
        gPad->SetTicks(1,1);
        TH2D *axisFrame1 = new TH2D("axisFrame1", "", 100, -0.1, 5.1, 100, 0.3, 1.7);
        axisFrame1->GetYaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetLabelSize(0.080);
        axisFrame1->GetYaxis()->SetLabelSize(0.080);
        axisFrame1->Draw("axis");

        TLatex* label1 = new TLatex();
        label1->SetTextSize(0.1);
        label1->SetTextAngle(90);
        label1->SetTextAlign(22);
        label1->DrawLatexNDC(0.05, 0.54, "<p_{T}>(GeV/c)");

        for (int i = 0; i < 3; i++) {
            TString histname = "hPiCh_dPhi0";
            double mass = getMassFromHistName(histname);
            auto result = getCentPtProjections(filenames[i], histname);
            std::vector<TH1D*> projections = result.first;
            double averageNT = result.second;
            TGraphErrors* graph = analyzeAndDrawpi(projections, gPad, averageNT, i, colors[i], markerStyles[i], 5.0, mass);
            graphsToward.push_back(graph);  // 存储图形
            if (i == 0) {
                graph->Draw("LP");
            } else {
                graph->Draw("LP same");
            }

            TFile* infile = new TFile(filenames1[0], "READ");
            if (infile && !infile->IsZombie()) {
                TGraphErrors *graph1 = (TGraphErrors*)infile->Get("Pt");
                if (graph1) {
                    graph1->SetMarkerStyle(24);
                    graph1->SetMarkerSize(1);
                    graph1->SetLineColor(kBlack);
                    graph1->SetMarkerColor(kBlack);
                    graph1->Draw("P same");
                    graphsToward.push_back(graph1);  // 存储图形
                }
                delete infile;
            }
        }
        addText2(0.51, 0.9, "Toward");
    }
    
    // 创建图例
    TLegend *legend = new TLegend(0.2, 0.57, 0.85, 0.87);
    legend->SetNColumns(1);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.07); 

    TGraphErrors *marker1 = new TGraphErrors();
    marker1->SetMarkerStyle(24);
    marker1->SetMarkerColor(kBlack);
    marker1->SetMarkerSize(1);
    marker1->SetLineWidth(2); 
    marker1->SetLineColor(kBlack);
    legend->AddEntry(marker1, "ALICE", "lp");

    TGraphErrors *marker2 = new TGraphErrors();
    marker2->SetMarkerStyle(20);
    marker2->SetMarkerColor(kRed);
    marker2->SetMarkerSize(1);
    marker2->SetLineWidth(2); 
    marker2->SetLineColor(kRed);
    legend->AddEntry(marker2, "noFSI", "lp");

    TGraphErrors *marker3 = new TGraphErrors();
    marker3->SetMarkerStyle(21);
    marker3->SetMarkerColor(kBlue);
    marker3->SetMarkerSize(1);
    marker3->SetLineWidth(2); 
    marker3->SetLineColor(kBlue);
    legend->AddEntry(marker3, "nohFSI", "lp");

    TGraphErrors *marker4 = new TGraphErrors();
    marker4->SetMarkerStyle(22);
    marker4->SetMarkerColor(kMagenta);
    marker4->SetMarkerSize(1);
    marker4->SetLineWidth(2); 
    marker4->SetLineColor(kMagenta);
    legend->AddEntry(marker4, "allFSI", "lp");
    legend->Draw(); 

    // 第二子图: Transverse区域(π)
    c1->cd(2);
    {
        gPad->SetTicks(1,1);
        TH2D *axisFrame1 = new TH2D("axisFrame1", "", 100, -0.1, 5.1, 100,  0.3, 1.7);
        axisFrame1->GetYaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetLabelSize(0.080);
        axisFrame1->GetYaxis()->SetLabelSize(0.080);
        axisFrame1->Draw("axis");

        TPaveText *text = new TPaveText(0.4, 0.01, 0.6, 0.1, "NDC");
        text->AddText("R_{T}");
        text->SetFillColor(0);
        text->SetTextAlign(22); 
        text->SetTextSize(0.12);
        text->Draw("same");

        for (int i = 0; i < 3; i++) {
            TString histname = "hPiCh_dPhi1";
            double mass = getMassFromHistName(histname);
            auto result = getCentPtProjections(filenames[i], histname);
            std::vector<TH1D*> projections = result.first;
            double averageNT = result.second;
            TGraphErrors* graph = analyzeAndDrawpi(projections, gPad, averageNT, i, colors[i], markerStyles[i], 5.0, mass);
            graphsTransverse.push_back(graph);  // 存储图形
            if (i == 0) {
                graph->Draw("LP");
            } else {
                graph->Draw("LP same");
            }
            
            TFile* infile = new TFile(filenames1[1], "READ");
            if (infile && !infile->IsZombie()) {
                TGraphErrors *graph1 = (TGraphErrors*)infile->Get("Pt");
                if (graph1) {
                    graph1->SetMarkerStyle(24);
                    graph1->SetMarkerSize(1);
                    graph1->SetLineColor(kBlack);
                    graph1->SetMarkerColor(kBlack);
                    graph1->Draw("P same");
                    graphsTransverse.push_back(graph1);  // 存储图形
                }
                delete infile;
            }
        }
        addText2(0.35, 0.9, "Transverse");
        addText4(0.45, 0.25, "#pi^{+}+#pi^{-}");
    }

    // 第三子图: In-Jet区域(π) - 使用拟合方法计算差谱
    c1->cd(3);
    {
        gPad->SetTicks(1,1);
        gPad->SetRightMargin(0.01);
        TH2D *axisFrame1 = new TH2D("axisFrame1", "", 50, -0.02, 1.15, 100, 0.3, 1.7);
        axisFrame1->GetYaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetTickLength(0.02);
        axisFrame1->GetXaxis()->SetLabelSize(0.080);
        axisFrame1->GetYaxis()->SetLabelSize(0.080);
        axisFrame1->Draw("axis");
        
        for (int i = 0; i < 3; i++) {
            TString histnameToward = "hPiCh_dPhi0";
            TString histnameTransverse = "hPiCh_dPhi1";
            auto resultToward = getCentPtProjections(filenames[i], histnameToward);
            auto resultTransverse = getCentPtProjections(filenames[i], histnameTransverse);
            std::vector<TH1D*> projectionsToward = resultToward.first;
            std::vector<TH1D*> projectionsTransverse = resultTransverse.first;
            double averageNT = resultToward.second;
        
            TGraphErrors *graph=drawInJetRegion(projectionsToward, projectionsTransverse, averageNT, i, colors[i], markerStyles[i]);
            graphsInJet.push_back(graph);
                }
        
        addText2(0.4, 0.9, "In-Jet");
    }
    outFile->cd("Toward");
    for (size_t i = 0; i < graphsToward.size(); ++i) {
        TString name = TString::Format("Pi_Toward_%zu", i);
        graphsToward[i]->Write(name);
    }
    
    outFile->cd("Transverse");
    for (size_t i = 0; i < graphsTransverse.size(); ++i) {
        TString name = TString::Format("Pi_Transverse_%zu", i);
        graphsTransverse[i]->Write(name);
    }
    
    outFile->cd("InJet");
    for (size_t i = 0; i < graphsInJet.size(); ++i) {
        TString name = TString::Format("Pi_InJet_%zu", i);
        graphsInJet[i]->Write(name);
    }

    // 关闭输出文件
    outFile->Close();
    delete outFile;
    c1->SaveAs("pi.png");
    c1->SaveAs("pi.pdf");
}