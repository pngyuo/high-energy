void addText(Double_t x, Double_t y, const char* text)
{
    TLatex* latex = new TLatex(x, y, text);
    latex->SetNDC(kTRUE);  // 将坐标转换为百分比
    latex->SetTextSize(0.04);  // 设置字体大小
    latex->Draw();
}

vector<TH1D*> getCentPtProjections(const TString& filename, const TString& histname, const vector<pair<double, double>>& RTcuts = {{0.0, 1.0}}) {

   	TFile* infile = new TFile(filename, "READ");
    TH3D* hist3D = (TH3D*)infile->Get(histname);

    TH1D *nTHist = (TH1D*)infile->Get("nTHist");

    double totalChargedParticles = 0; // 加权总和
    double nEvents = 0;  // 总的事件数


    double averageNT = nTHist->GetMean();
    cout<<"averageNT: "<<averageNT<<endl;
    double nEvents_tot = nTHist->GetEntries();

    vector<TH1D*> hist_list;
    
    // Loop through each RT cut pair
    for (const auto& RTcut : RTcuts) {
        double RTmin = RTcut.first;
        double RTmax = RTcut.second;
        
        double NT_min = averageNT * RTmin;
        double NT_max = averageNT * RTmax;
        int bin_NT_min = nTHist->FindBin(NT_min);
        int bin_NT_max = nTHist->FindBin(NT_max);
        
        if (RTmax < 0) {
            bin_NT_min = 1;
            bin_NT_max = nTHist->GetNbinsX();
        }

        // Project X axis histogram for this RT cut
        TH1D* histProjX = hist3D->ProjectionX(Form("histProjX_%g_%g", RTmin, RTmax), 
                                             1, hist3D->GetNbinsY(), 
                                             bin_NT_min, bin_NT_max);

        double totalEventsInRange = nTHist->Integral(bin_NT_min, bin_NT_max);
        histProjX->Scale(1.0 / totalEventsInRange, "width");
        hist_list.push_back(histProjX);
    }
    
    return hist_list;
}

void draw_pt_ratio(){
  gStyle->SetOptStat(kFALSE);
  TH1::SetDefaultSumw2();
  TString filenames[3] = {"hist_outputallFSI_liang.root", "hist_outputnoFSI_liang.root", "hist_outputnohFSI_liang.root"};

  // Histogram names for Pi, P, and K
  TString histPiTow = "hPiCh_dPhi0";
  TString histPiTrans = "hPiCh_dPhi1";
  TString histPiTransMax = "hPiCh_dPhi1max";
  TString histPiTransMin = "hPiCh_dPhi1min";

  TString histPTow = "hProton_dPhi0";
  TString histPTrans = "hProton_dPhi1";
  TString histPTransMax = "hProton_dPhi1max";
  TString histPTransMin = "hProton_dPhi1min";

  // Add K histograms
  TString histKTow = "hKCh_dPhi0";
  TString histKTrans = "hKCh_dPhi1";
  TString histKTransMax = "hKCh_dPhi1max";
  TString histKTransMin = "hKCh_dPhi1min";

  // Define RT cuts
  //vector<pair<double, double>> RTcuts = {{0.0, 1.0}, {0.0, 0.5}, {0.5, 1.0}};
  vector<pair<double, double>> RTcuts = {{0.0, 0.5}, {0.5, 1.5}};

  // Get projections with RT cuts for Pi, P, and K
  vector<TH1D*> projectionsPiTow = getCentPtProjections(filenames[0], histPiTow, RTcuts);
  vector<TH1D*> projectionsPiTrans = getCentPtProjections(filenames[0], histPiTrans, RTcuts);
  vector<TH1D*> projectionsPiTransMax = getCentPtProjections(filenames[0], histPiTransMax, RTcuts);
  vector<TH1D*> projectionsPiTransMin = getCentPtProjections(filenames[0], histPiTransMin, RTcuts);

  vector<TH1D*> projectionsPTow = getCentPtProjections(filenames[0], histPTow, RTcuts);
  vector<TH1D*> projectionsPTrans = getCentPtProjections(filenames[0], histPTrans, RTcuts);
  vector<TH1D*> projectionsPTransMax = getCentPtProjections(filenames[0], histPTransMax, RTcuts);
  vector<TH1D*> projectionsPTransMin = getCentPtProjections(filenames[0], histPTransMin, RTcuts);

  // Add K projections
  vector<TH1D*> projectionsKTow = getCentPtProjections(filenames[0], histKTow, RTcuts);
  vector<TH1D*> projectionsKTrans = getCentPtProjections(filenames[0], histKTrans, RTcuts);
  vector<TH1D*> projectionsKTransMax = getCentPtProjections(filenames[0], histKTransMax, RTcuts);
  vector<TH1D*> projectionsKTransMin = getCentPtProjections(filenames[0], histKTransMin, RTcuts);

  // Create vectors to store ratios for different RT cuts
  vector<TH1D*> hPOverPiTow_cuts, hPOverPiTrans_cuts, hPOverPiInjet_cuts;
  vector<TH1D*> hPTransMaxOverPiTransMax_cuts, hPTransMinOverPiTransMin_cuts;
  // Add K/Pi ratio vectors
  vector<TH1D*> hKOverPiTow_cuts, hKOverPiTrans_cuts, hKOverPiInjet_cuts;
  vector<TH1D*> hKTransMaxOverPiTransMax_cuts, hKTransMinOverPiTransMin_cuts;

  // Process each RT cut
  for (int i = 0; i < RTcuts.size(); i++) {
    // Clone Pi histograms
    TH1D *hPiTow = (TH1D*)projectionsPiTow[i]->Clone(Form("hPiTow_%d", i));
    TH1D *hPiTrans = (TH1D*)projectionsPiTrans[i]->Clone(Form("hPiTrans_%d", i));
    TH1D *hPiTransMax = (TH1D*)projectionsPiTransMax[i]->Clone(Form("hPiTransMax_%d", i));
    TH1D *hPiTransMin = (TH1D*)projectionsPiTransMin[i]->Clone(Form("hPiTransMin_%d", i));

    // Clone P histograms
    TH1D *hPTow = (TH1D*)projectionsPTow[i]->Clone(Form("hPTow_%d", i));
    TH1D *hPTrans = (TH1D*)projectionsPTrans[i]->Clone(Form("hPTrans_%d", i));
    TH1D *hPTransMax = (TH1D*)projectionsPTransMax[i]->Clone(Form("hPTransMax_%d", i));
    TH1D *hPTransMin = (TH1D*)projectionsPTransMin[i]->Clone(Form("hPTransMin_%d", i));

    // Clone K histograms
    TH1D *hKTow = (TH1D*)projectionsKTow[i]->Clone(Form("hKTow_%d", i));
    TH1D *hKTrans = (TH1D*)projectionsKTrans[i]->Clone(Form("hKTrans_%d", i));
    TH1D *hKTransMax = (TH1D*)projectionsKTransMax[i]->Clone(Form("hKTransMax_%d", i));
    TH1D *hKTransMin = (TH1D*)projectionsKTransMin[i]->Clone(Form("hKTransMin_%d", i));

    // Calculate P/Pi ratios
    TH1D *hPOverPiTow = (TH1D*)hPTow->Clone(Form("hPOverPiTow_%d", i));
    TH1D *hPOverPiTrans = (TH1D*)hPTrans->Clone(Form("hPOverPiTrans_%d", i));
    TH1D *hPOverPiInjet = (TH1D*)hPTow->Clone(Form("hPOverPiInjet_%d", i));

    // Calculate K/Pi ratios
    TH1D *hKOverPiTow = (TH1D*)hKTow->Clone(Form("hKOverPiTow_%d", i));
    TH1D *hKOverPiTrans = (TH1D*)hKTrans->Clone(Form("hKOverPiTrans_%d", i));
    TH1D *hKOverPiInjet = (TH1D*)hKTow->Clone(Form("hKOverPiInjet_%d", i));

    // Calculate ratios for P/Pi
    hPOverPiTow->Divide(hPiTow);
    hPOverPiTrans->Divide(hPiTrans);

    // Calculate ratios for K/Pi
    hKOverPiTow->Divide(hPiTow);
    hKOverPiTrans->Divide(hPiTrans);

    // Calculate injet contributions
    TH1D *hPInjet = (TH1D*)hPTow->Clone(Form("hPInjet_%d", i));
    TH1D *hKInjet = (TH1D*)hKTow->Clone(Form("hKInjet_%d", i));
    TH1D *hPiInjet = (TH1D*)hPiTow->Clone(Form("hPiInjet_%d", i));

    hPInjet->Add(hPTrans, -1);
    hKInjet->Add(hKTrans, -1);
    hPiInjet->Add(hPiTrans, -1);

    hPOverPiInjet->Add(hPTrans, -1);
    hPOverPiInjet->Divide(hPiInjet);

    hKOverPiInjet->Add(hKTrans, -1);
    hKOverPiInjet->Divide(hPiInjet);

    // Store all ratios
    hPOverPiTow_cuts.push_back(hPOverPiTow);
    hPOverPiTrans_cuts.push_back(hPOverPiTrans);
    hPOverPiInjet_cuts.push_back(hPOverPiInjet);

    hKOverPiTow_cuts.push_back(hKOverPiTow);
    hKOverPiTrans_cuts.push_back(hKOverPiTrans);
    hKOverPiInjet_cuts.push_back(hKOverPiInjet);

    // Process TransMax and TransMin for both P/Pi and K/Pi
    TH1D *hPTransMaxRatio = (TH1D*)hPTransMax->Clone(Form("hPTransMaxRatio_%d", i));
    TH1D *hPTransMinRatio = (TH1D*)hPTransMin->Clone(Form("hPTransMinRatio_%d", i));
    TH1D *hKTransMaxRatio = (TH1D*)hKTransMax->Clone(Form("hKTransMaxRatio_%d", i));
    TH1D *hKTransMinRatio = (TH1D*)hKTransMin->Clone(Form("hKTransMinRatio_%d", i));

    hPTransMaxRatio->Divide(hPiTransMax);
    hPTransMinRatio->Divide(hPiTransMin);
    hKTransMaxRatio->Divide(hPiTransMax);
    hKTransMinRatio->Divide(hPiTransMin);

//    hPTransMaxRatio->Rebin(2);
//    hPTransMinRatio->Rebin(2);
//    hKTransMaxRatio->Rebin(2);
//    hKTransMinRatio->Rebin(2);

    hPTransMaxOverPiTransMax_cuts.push_back(hPTransMaxRatio);
    hPTransMinOverPiTransMin_cuts.push_back(hPTransMinRatio);
    hKTransMaxOverPiTransMax_cuts.push_back(hKTransMaxRatio);
    hKTransMinOverPiTransMin_cuts.push_back(hKTransMinRatio);
  }

  // Create two canvases for P/Pi and K/Pi ratios
  TCanvas *c1 = new TCanvas("c1", "P/Pi Ratios");
  TCanvas *c2 = new TCanvas("c2", "K/Pi Ratios");

  // Draw P/Pi ratios
  c1->cd();
  TH2D *axisFrame1 = new TH2D("axisFrame1", ";p_{T} [GeV/c]; p/#pi;", 100, -0.1, 3.1, 100, 0, 0.5);
  axisFrame1->Draw("axis");

  // Draw K/Pi ratios
  c2->cd();
  TH2D *axisFrame2 = new TH2D("axisFrame2", ";p_{T} [GeV/c]; K/#pi;", 100, -0.1, 3.1, 100, 0, 0.5);
  axisFrame2->Draw("axis");

  // Colors for RT cuts
  Color_t rtColors[] = {kBlack, kRed, kBlue, kGreen, kMagenta};
  
  // Colors for different ratio types
  Color_t typeColors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+2};
  Style_t typeStyles[] = {20, 24, 21, 25, 28};

  // Draw all ratios
  for (size_t i = 0; i < RTcuts.size(); i++) {
    // Draw P/Pi ratios
    c1->cd();
    
    // Draw each type with its own style but RT-dependent color
    vector<TH1D*> pRatios = {
        hPOverPiTow_cuts[i],      // Toward
//        hPOverPiTrans_cuts[i],    // Transverse
        hPOverPiInjet_cuts[i],    // InJet
//        hPTransMaxOverPiTransMax_cuts[i], // TransMax
//        hPTransMinOverPiTransMin_cuts[i]  // TransMin
    };

    for (size_t j = 0; j < pRatios.size(); j++) {
//        pRatios[j]->SetMarkerColor(typeColors[j]);
//        pRatios[j]->SetLineColor(typeColors[j]);
        pRatios[j]->SetMarkerColor(rtColors[i]);
        pRatios[j]->SetLineColor(rtColors[i]);

        pRatios[j]->SetMarkerStyle(typeStyles[j]);
        pRatios[j]->Draw("same");
    }

    // Draw K/Pi ratios
    c2->cd();
    vector<TH1D*> kRatios = {
        hKOverPiTow_cuts[i],      // Toward
//        hKOverPiTrans_cuts[i],    // Transverse
        hKOverPiInjet_cuts[i],    // InJet
//        hKTransMaxOverPiTransMax_cuts[i], // TransMax
//        hKTransMinOverPiTransMin_cuts[i]  // TransMin
    };

    for (size_t j = 0; j < kRatios.size(); j++) {
//        kRatios[j]->SetMarkerColor(typeColors[j]);
//        kRatios[j]->SetLineColor(typeColors[j]);
        kRatios[j]->SetMarkerColor(rtColors[i]);
        kRatios[j]->SetLineColor(rtColors[i]);
        kRatios[j]->SetMarkerStyle(typeStyles[j]);
        kRatios[j]->Draw("same");
    }
  }

  // Add legends to both canvases with both RT cuts and ratio types
  for (auto c : {c1, c2}) {
    c->cd();
    TLegend *legend = new TLegend(0.65, 0.65, 0.85, 0.85);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    
    // Add ratio type entries (using first RT cut as example)
    legend->AddEntry(hPOverPiTow_cuts[0], "Toward", "p");
//    legend->AddEntry(hPOverPiTrans_cuts[0], "Transverse", "p");
    legend->AddEntry(hPOverPiInjet_cuts[0], "InJet", "p");
//    legend->AddEntry(hPTransMaxOverPiTransMax_cuts[0], "TransMax", "p");
//    legend->AddEntry(hPTransMinOverPiTransMin_cuts[0], "TransMin", "p");
    
    legend->Draw();
  }

  // Add RT cut labels if needed
  for (auto c : {c1, c2}) {
    c->cd();
    TLegend *legend = new TLegend(0.2, 0.65, 0.4, 0.85);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    for (size_t i = 0; i < RTcuts.size(); i++) {
      legend->AddEntry(hPOverPiTow_cuts[i], 
                      Form("RT %.1f-%.1f", RTcuts[i].first, RTcuts[i].second), 
                      "p");
    }
    legend->Draw();
  }

}
