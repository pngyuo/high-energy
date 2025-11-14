
TF1 *fitfun = new TF1("fitfun", "x*[0]*([1]-1.)*([1]-2.)/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt(x*x+[3]*[3])-[3])/[1]/[2], -[1])", 0, 20);


void getAvgPtALICE(TGraphAsymmErrors *grIn, int type=1){

  double mass=0.134;

  if(type==2) mass=0.493;
  if(type==3) mass=0.938;
  fitfun->SetParameters(0.1,20,0.1, mass);
  for(int i=0; i<grIn->GetN(); i++){
    double x, y;
    grIn->GetPoint(i, x, y);

    TString tableName=Form("Table %i",2*type-1);
    TString grName=Form("Graph1D_y%i",i+1);
    if(i==grIn->GetN()-1) {
      tableName=Form("Table %i",2*type);
      grName=Form("Graph1D_y%i",1);
    }
    //cout<<tableName<<" "<<grName<<endl;
    TGraphAsymmErrors *grTemp = getALICE(tableName, grName);
    grTemp->Fit("fitfun");
    //get data mean within pt~(0,3)
    y=fitfun->Mean(0,3);
    grIn->SetPointY(i, y);
    delete grTemp;
  }
}

