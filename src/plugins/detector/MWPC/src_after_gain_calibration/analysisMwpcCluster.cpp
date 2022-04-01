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
#include "MwpcPedestalValue.h"
using namespace std;

static TH1D *histNClus[12];
static TH2D *histMult[12];
const static int nType=10;
static int clusType[12][nType]={{0}};
Long_t Det_MWPC::histos_cluster(){
  ostringstream sstr_name, sstr_title;
  for(int i=0;i<12;i++){
    sstr_name<<"Distribution of number of clusters for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Number of clusters; Counts";
    histNClus[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),11,-0.5,10.5);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Total charge vs cluster's type for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Cluster's Type;Total charge";
    histMult[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),nType,0.5,nType+0.5,100,0,3000);
    sstr_name.str("");
    sstr_title.str("");
  }

      
  return 0;
}

Long_t Det_MWPC::startup_cluster(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  return 0;
};


Long_t Det_MWPC::process_cluster(){
  int* data[12]={treeRaw->ADC_C2X_L,treeRaw->ADC_C2X_R,treeRaw->ADC_C2Y_L,treeRaw->ADC_C2Y_R,treeRaw->ADC_C3X_L,treeRaw->ADC_C3X_R,treeRaw->ADC_C3Y_L,treeRaw->ADC_C3Y_R,treeRaw->ADC_C4X_L,treeRaw->ADC_C4X_R,treeRaw->ADC_C4Y_L,treeRaw->ADC_C4Y_R};
  const double* pedValue[12]={pedValueC2XL,pedValueC2XR,pedValueC2YL,pedValueC2YR,pedValueC3XL,pedValueC3XR,pedValueC3YL,pedValueC3YR,pedValueC4XL,pedValueC4XR,pedValueC4YL,pedValueC4YR};
  const double* pedSigma[12]={pedSigmaC2XL,pedSigmaC2XR,pedSigmaC2YL,pedSigmaC2YR,pedSigmaC3XL,pedSigmaC3XR,pedSigmaC3YL,pedSigmaC3YR,pedSigmaC4XL,pedSigmaC4XR,pedSigmaC4YL,pedSigmaC4YR}; 

  MwpcEvent *event[12];
  ostringstream sstr;
  int *mult[12];
  vector<int> nSt[12];
  double *totCh[12];
  for(int i=0;i<12;i++){
    event[i]=new MwpcEvent(((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,!((i/4==0)&&((i%4/2))==0));
/*
    if(i==9){
      for(int j=0;j<72;j++){
        if(j==16 || j==21 || j==66 || j==67 || j==68 || j==69)  data[i][j]=0;
      }
    }
*/
    histNClus[i]->Fill(event[i]->getNCluster(event[i],data[i],pedValue[i],pedSigma[i],0));
    mult[i]=event[i]->multipletNType(data[i],pedValue[i],pedSigma[i],0,nType);
    nSt[i]=event[i]->getClusterNStrip(event[i],data[i],pedValue[i],pedSigma[i],0);
    totCh[i]=event[i]->getClusterTotCharge(event[i],data[i],pedValue[i],pedSigma[i],0);
    for(int j=0;j<nType;j++){
      clusType[i][j]+=mult[i][j];
    }   
    for(UInt_t j=0;j<nSt[i].size();j++){
      histMult[i]->Fill(nSt[i][j],totCh[i][j]);
    }
  }
  return 0;
}

Long_t Det_MWPC::done_cluster(){
  TCanvas *canvNClus=new TCanvas("Distribution of number of clusters","Distribution of number of clusters",2000,900);
  canvNClus->Divide(4,3);
  for(int i=0;i<12;i++){
    canvNClus->cd(i+1);
    histNClus[i]->Draw();
  }
  canvNClus->SaveAs("Distribution of number of clusters.pdf");
   
  int totNClus[12]={0};
  for(int i=0;i<12;i++){
    for(int j=0;j<nType;j++){
      totNClus[i]+=clusType[i][j];
    }
  }

  TH1D *histClusType[12];
  ostringstream sstr_name, sstr_title;
  for(int i=0;i<12;i++){
    sstr_name<<"Statistic of cluster's type for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Type of clusters; occupation (%)";
    histClusType[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),nType,0.5,nType+0.5);
    sstr_name.str("");
    sstr_title.str("");
    for(int j=0;j<nType;j++){
      histClusType[i]->SetBinContent(j+1,(double)clusType[i][j]/(double)totNClus[i]*100);
    }
  }

  TCanvas *canvClusType=new TCanvas("Statistic of cluster's type","Statistic of cluster's type",2000,900);
  canvClusType->Divide(4,3);
  for(int i=0;i<12;i++){
    canvClusType->cd(i+1);
    histClusType[i]->SetMinimum(0);
    histClusType[i]->SetMaximum(60);
    //histClusType[i]->GetXaxis()->SetLabelOffset(99);
    histClusType[i]->GetXaxis()->SetTitleOffset(0.8);
    histClusType[i]->GetXaxis()->SetTitleSize(0.06);
    histClusType[i]->GetYaxis()->SetLabelSize(0.06);
    histClusType[i]->GetYaxis()->SetTitleOffset(0.8);
    histClusType[i]->GetYaxis()->SetTitleSize(0.06);
    histClusType[i]->Draw();
    double x, y;
    string type[10]={"S","D","T","Qa","Qi","Sx","Sp","Oc","No","Dc"};
    y = gPad->GetUymin() + 6*histClusType[i]->GetYaxis()->GetBinWidth(1);
    TText t;
    t.SetTextAngle(60);
    t.SetTextSize(0.06);
    t.SetTextAlign(33);
    t.SetTextColor(kRed);
    for (int j=0;j<10;j++) {
      x = histClusType[j]->GetXaxis()->GetBinCenter(j+1);
      t.DrawText(x,y,type[j].c_str());
    }
  }
  canvClusType->SaveAs("Statistic of cluster's type.pdf");

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
 
  return 0;
};
