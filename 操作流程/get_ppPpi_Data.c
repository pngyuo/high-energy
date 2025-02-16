TGraphErrors *convertAsymmErrToErr(TGraphAsymmErrors *gr_in){
	TGraphErrors *gr = new TGraphErrors();
	for(int i(0); i<gr_in->GetN(); i++){
		double x, y, x_err, y_err;
		gr_in->GetPoint(i, x, y);
		x_err = gr_in->GetErrorX(i);
		y_err = gr_in->GetErrorY(i);
		gr->SetPoint(i, x, y);
		gr->SetPointError(i, x_err, y_err);
	}
	return gr;
}

void convertToInvSpectra(TGraphErrors *gr_in){
	for(int i(0); i<gr_in->GetN(); i++){
		double pt, y, pt_err, y_err;
		gr_in->GetPoint(i, pt, y);
		pt_err = gr_in->GetErrorX(i);
		y_err = gr_in->GetErrorY(i);
//   double pi=3.1415926;
		gr_in->SetPoint(i, pt, y);
		gr_in->SetPointError(i, pt_err, y_err);
}
}

void get_ppPpi_Data(){


	TFile *input = new TFile("HEPData-ins2626034-v1-Table_28.root");//修改一：改为下载的名字
  //output file name convention: system_energy_particleSpecies.root
	TFile *output = new TFile("K_Tr.root","recreate");//修改二：把ppi_Tr_25_5.root改为你想要输出的名字，最好与提取的文件数据相关，避免之后混淆

	TDirectoryFile *file;


	file = (TDirectoryFile*)input->Get("Table 28");//修改三：提取的数据是table什么就改为什么，有数字也有字母的

	TGraphAsymmErrors *gr_in = (TGraphAsymmErrors*)file->Get("Graph1D_y1");//修改四：有些数据一个table里面有很多数据，第一个就是y1，第二个就是y2，一般一个table就一个数据图，因此一般不需要修改
			TGraphErrors *gr_out = convertAsymmErrToErr(gr_in);
			convertToInvSpectra(gr_out);
			output->cd();
			gr_out->Write("Pt");//修改五：可改可不改，用于之后数据转为图形的时候用的，也可以改为数据相关的名字
		}




