#include <Det_MWPC.h>

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <sstream>
#include "TLine.h"
#include <string>
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "MwpcPedestalValue.h"
#include "MwpcClusterThresholdSecondApproach.h"
using namespace std;

static TH2D *histMult[12];
static const int nType=10;
static TH1D *histType[nType][12];
static string nameType[nType]={"Singlet","Doublet","Triplet","Quartet","Quintet","Sextet","Septet","OCtet","Nonet","Dectet"};
static int binNum=200;
static int binUp=500;
Long_t Det_MWPC::histos_clusterthresholdfromcharge(){
  ostringstream sstr_name, sstr_title;
  for(int i=0;i<12;i++){
    sstr_name<<"Total charge vs cluster's type for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Cluster's Type;Total charge";
    histMult[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),nType,0.5,0.5+nType,binNum,0,binUp);
    sstr_name.str("");
    sstr_title.str("");    
    
    for(int k=0;k<nType;k++){
      sstr_name<<"Total charge distribution for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<" for "<<nameType[k];
      sstr_title<<sstr_name.str()<<";Total charge;Counts";
      histType[k][i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),binNum,0,binUp);
      sstr_name.str("");
      sstr_title.str("");
    }
  }
  return 0;
}

Long_t Det_MWPC::startup_clusterthresholdfromcharge(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  return 0;
};


Long_t Det_MWPC::process_clusterthresholdfromcharge(){
  int* data[12]={treeRaw->ADC_C2X_L,treeRaw->ADC_C2X_R,treeRaw->ADC_C2Y_L,treeRaw->ADC_C2Y_R,treeRaw->ADC_C3X_L,treeRaw->ADC_C3X_R,treeRaw->ADC_C3Y_L,treeRaw->ADC_C3Y_R,treeRaw->ADC_C4X_L,treeRaw->ADC_C4X_R,treeRaw->ADC_C4Y_L,treeRaw->ADC_C4Y_R};
  const double* pedValue[12]={pedValueC2XL,pedValueC2XR,pedValueC2YL,pedValueC2YR,pedValueC3XL,pedValueC3XR,pedValueC3YL,pedValueC3YR,pedValueC4XL,pedValueC4XR,pedValueC4YL,pedValueC4YR};
  const double* pedSigma[12]={pedSigmaC2XL,pedSigmaC2XR,pedSigmaC2YL,pedSigmaC2YR,pedSigmaC3XL,pedSigmaC3XR,pedSigmaC3YL,pedSigmaC3YR,pedSigmaC4XL,pedSigmaC4XR,pedSigmaC4YL,pedSigmaC4YR};

  MwpcEvent *event[12];
  ostringstream sstr;
  vector<int> nSt[12];
  double *totCh[12];
  for(int i=0;i<12;i++){
    event[i]=new MwpcEvent(((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,!((i/4==0)&&((i%4/2))==0));
    totCh[i]=event[i]->getClusterTotCharge(event[i],data[i],pedValue[i],pedSigma[i],0);
    nSt[i]=event[i]->getClusterNStrip(event[i],data[i],pedValue[i],pedSigma[i],0);

    for(UInt_t j=0;j<nSt[i].size();j++){
      histMult[i]->Fill(nSt[i][j],totCh[i][j]);
      for(int k=0;k<nType;k++){
        if(nSt[i][j]==k+1) histType[k][i]->Fill(totCh[i][j]);
      }
    }
  }
  return 0;
}

Long_t Det_MWPC::done_clusterthresholdfromcharge(){
  ostringstream sstr;
  TCanvas *canv[nType];
  for(int k=0;k<10;k++){
    sstr<<"Total charge distribution for "<<nameType[k];
    canv[k]=new TCanvas(sstr.str().c_str(),sstr.str().c_str(),2000,900);
    sstr.str("");
    canv[k]->Divide(4,3);
  }

  double lowRange,upRange;
  double mean[nType][12],mean_err[nType][12];
  double sigma[nType][12],sigma_err[nType][12];
  for(int k=0;k<10;k++){
    for(int i=0;i<12;i++){
      canv[k]->cd(i+1);
      //gPad->SetLogx();
      lowRange=histType[k][i]->GetBinLowEdge(histType[k][i]->GetMaximumBin()-5);
      upRange=histType[k][i]->GetBinLowEdge(histType[k][i]->GetMaximumBin()+6);
      histType[k][i]->Fit("gaus","Q","",lowRange,upRange);
      mean[k][i]=histType[k][i]->GetFunction("gaus")->GetParameter(1);
      mean_err[k][i]=histType[k][i]->GetFunction("gaus")->GetParError(1);
      sigma[k][i]=histType[k][i]->GetFunction("gaus")->GetParameter(2);
      sigma_err[k][i]=histType[k][i]->GetFunction("gaus")->GetParError(2);
      histType[k][i]->Draw();
    }
  sstr<<"Total charge distribution for "<<nameType[k]<<".pdf";
  canv[k]->SaveAs(sstr.str().c_str());
  sstr.str("");
  }

  double ChargeCut[12][nType],ChargeCut_err[12][nType];
  for(int i=0;i<12;i++){
    for(int k=0;k<nType;k++){
      ChargeCut[i][k]=mean[k][i]+3*sigma[k][i];
      ChargeCut_err[i][k]=sqrt(mean_err[k][i]*mean_err[k][i]+9*sigma[k][i]*sigma[k][i]);
    }
  }

  double x[nType],x_err[nType];
  for(int k=0;k<nType;k++){
     x[k]=k+1;
     x_err[k]=0;
  }

  ostringstream sstr_title;
  TGraphErrors *gr[12];
  for(int i=0;i<12;i++){
    gr[i]=new TGraphErrors(nType,x,ChargeCut[i],x_err,ChargeCut_err[i]);
    sstr_title<<"Cut of total charge vs culster's type for MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    gr[i]->SetTitle(sstr_title.str().c_str());
    sstr_title.str("");
    gr[i]->GetXaxis()->SetTitle("Number of Strips for a clusters");
    gr[i]->GetXaxis()->SetLimits(0,nType+1);

    gr[i]->GetYaxis()->SetTitle("Cut of total charge");
    gr[i]->SetMinimum(0);
    gr[i]->SetMaximum(500);

    gr[i]->SetMarkerSize(1.5);
  }

  TF1 *func[12];
  for(int i=0;i<12;i++){
    sstr<<"func"<<i;
    func[i]=new TF1(sstr.str().c_str(),"[0]+[1]*x",0.5,nType+0.5);
    sstr.str("");
  }

  Double_t* par[12],*parErr[12];
  TCanvas *canvCut=new TCanvas("Cut of total charge vs culster's type","Cut of total charge vs culster's type",2000,900);
  canvCut->Divide(4,3);
  for(int i=0;i<12;i++){
    canvCut->cd(i+1);
    sstr<<"func"<<i;
    gr[i]->Fit(sstr.str().c_str(),"R");
    par[i]=func[i]->GetParameters();
    parErr[i]=func[i]->GetParErrors();
    sstr.str("");
    gr[i]->Draw("ap");
  }
  canvCut->SaveAs("Cut of total charge vs culster's type.pdf");

  double ChargeCutCalc[12][nType],ChargeCutCalcErr[12][nType];
  for(int i=0;i<12;i++){
    for(int j=0;j<nType;j++){
      ChargeCutCalc[i][j]=func[i]->Eval(j+1);
      ChargeCutCalcErr[i][j]=sqrt(parErr[i][0]*parErr[i][0]+(j+1)*(j+1)*parErr[i][1]*parErr[i][1]);
    }
  }   
  
  TGraphErrors *grChargeCut[12];
  for(int i=0;i<12;i++){
    grChargeCut[i]=new TGraphErrors(nType,x,ChargeCutCalc[i],x_err,ChargeCutCalcErr[i]);
    grChargeCut[i]->SetMarkerColor(kRed);
    grChargeCut[i]->SetLineColor(kRed);
  }

  TCanvas *canvMult=new TCanvas("Charge vs cluster's type","Charge vs cluster's type",2000,900);
  canvMult->Divide(4,3);
  for(int i=0;i<12;i++){
    canvMult->cd(i+1);
    gPad->SetLogz();
    histMult[i]->GetXaxis()->SetLabelOffset(99);
    histMult[i]->GetXaxis()->SetTitleOffset(0.5);
    histMult[i]->GetXaxis()->SetTitleSize(0.06);
    histMult[i]->GetYaxis()->SetLabelSize(0.06);
    histMult[i]->GetYaxis()->SetTitleOffset(0.8);
    histMult[i]->GetYaxis()->SetTitleSize(0.06);
    histMult[i]->GetZaxis()->SetLabelSize(0.045);
    histMult[i]->Draw("colz");
    grChargeCut[i]->Draw("p");
    double x, y;
    string type[10]={"S","D","T","Qa","Qi","Sx","Sp","Oc","No","Dc"};
    y = gPad->GetUymin() + 1*histMult[i]->GetYaxis()->GetBinWidth(1);
    TText t;
    t.SetTextAngle(60);
    t.SetTextSize(0.06);
    t.SetTextAlign(33);
    t.SetTextColor(kRed);
    for (int j=0;j<10;j++) {
      x = histMult[j]->GetXaxis()->GetBinCenter(j+1);
      t.DrawText(x,y,type[j].c_str());
    }
  }
  canvMult->SaveAs("Charge vs cluster's type.pdf");

  ofstream output("MwpcClusterThresholdFromTotalCharge.h");
  for(int i=0;i<12;i++){
    output<<"const double totChargeCutC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<'['<<((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8<<"]={"; 
    for(int j=0;j<((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8;j++){
      if(j<((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-1) output<<func[i]->Eval(j+1)<<", ";
      else output<<func[i]->Eval(j+1)<<"};";
    }
    output<<endl;
    output<<"const double totChargeCutErrC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<'['<<((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8<<"]={";
    for(int j=0;j<((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8;j++){
      if(j<((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-1) output<<sqrt(parErr[i][0]*parErr[i][0]+(j+1)*(j+1)*parErr[i][1]*parErr[i][1])<<", ";
      else output<<sqrt(parErr[i][0]*parErr[i][0]+(j+1)*(j+1)*parErr[i][1]*parErr[i][1])<<"};";
    }
    output<<endl;
  }

  return 0;
};
