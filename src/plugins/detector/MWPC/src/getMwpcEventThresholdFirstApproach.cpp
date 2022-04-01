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
#include "TGraphErrors.h"
#include "TMarker.h"
#include "MwpcPedestalValue.h"
using namespace std;

const static int nPoint=20;
static int totNEvent=0;
static int eventType[8][12][nPoint]={{{0}}};

Long_t Det_MWPC::histos_eventthresholdFirstApproach(){
  return 0;
}

Long_t Det_MWPC::startup_eventthresholdFirstApproach(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  return 0;
};


Long_t Det_MWPC::process_eventthresholdFirstApproach(){
  int* data[12]={treeRaw->ADC_C2X_L,treeRaw->ADC_C2X_R,treeRaw->ADC_C2Y_L,treeRaw->ADC_C2Y_R,treeRaw->ADC_C3X_L,treeRaw->ADC_C3X_R,treeRaw->ADC_C3Y_L,treeRaw->ADC_C3Y_R,treeRaw->ADC_C4X_L,treeRaw->ADC_C4X_R,treeRaw->ADC_C4Y_L,treeRaw->ADC_C4Y_R};
  const double* pedValue[12]={pedValueC2XL,pedValueC2XR,pedValueC2YL,pedValueC2YR,pedValueC3XL,pedValueC3XR,pedValueC3YL,pedValueC3YR,pedValueC4XL,pedValueC4XR,pedValueC4YL,pedValueC4YR};
  const double* pedSigma[12]={pedSigmaC2XL,pedSigmaC2XR,pedSigmaC2YL,pedSigmaC2YR,pedSigmaC3XL,pedSigmaC3XR,pedSigmaC3YL,pedSigmaC3YR,pedSigmaC4XL,pedSigmaC4XR,pedSigmaC4YL,pedSigmaC4YR}; 

  MwpcEvent *event[12];
  ostringstream sstr;

  int nHit[12][nPoint],*mult[12][nPoint];
  totNEvent++;
  for(int i=0;i<12;i++){
    event[i]=new MwpcEvent(((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,!((i/4==0)&&((i%4/2))==0));
    for(int j=0;j<nPoint;j++){
      nHit[i][j]=event[i]->getNHit(data[i],pedValue[i],pedSigma[i],0.5*j);
      mult[i][j]=event[i]->multiplet3Type(data[i],pedValue[i],pedSigma[i],0.5*j);
    }
  }

  for(int i=0;i<12;i++){
    for(int j=0;j<nPoint;j++){
      if(nHit[i][j]==0){
        eventType[0][i][j]++;
      }
      else{
        if(mult[i][j][0]!=0 && mult[i][j][1]==0 && mult[i][j][2]==0) eventType[1][i][j]++;
        else if(mult[i][j][0]==0 && mult[i][j][1]!=0 && mult[i][j][2]==0) eventType[2][i][j]++;
        else if(mult[i][j][0]==0 && mult[i][j][1]==0 && mult[i][j][2]!=0) eventType[3][i][j]++; 
        else if(mult[i][j][0]!=0 && mult[i][j][1]!=0 && mult[i][j][2]==0) eventType[4][i][j]++;
        else if(mult[i][j][0]!=0 && mult[i][j][1]==0 && mult[i][j][2]!=0) eventType[5][i][j]++;
        else if(mult[i][j][0]==0 && mult[i][j][1]!=0 && mult[i][j][2]!=0) eventType[6][i][j]++;
        else eventType[7][i][j]++;
      }
    }
  }

  return 0;
}


Long_t Det_MWPC::done_eventthresholdFirstApproach(){
  double ocpt[8][12][nPoint];
  double ocpt_err[8][12][nPoint];
  for(int i=0;i<8;i++){
    for(int j=0;j<12;j++){
      for(int k=0;k<nPoint;k++){
        ocpt[i][j][k]=(double)eventType[i][j][k]/totNEvent*100;
        ocpt_err[i][j][k]=1./totNEvent*sqrt(eventType[i][j][k]*(totNEvent-eventType[i][j][k])/(double)totNEvent)*100;
      }
    }
  }
  double x[nPoint],x_err[nPoint];
  for(int k=0;k<nPoint;k++){
    x[k]=0.5*k;
    x_err[k]=0;
  }

  TGraphErrors *gr[8][12];
  for(int i=0;i<8;i++){
    for(int j=0;j<12;j++){
      gr[i][j]=new TGraphErrors(nPoint,x,ocpt[i][j],x_err,ocpt_err[i][j]);
    }
  }


  ostringstream sstr_title;
  for(int j=0;j<12;j++){
    for(int i=0;i<8;i++){
      sstr_title<<"OCPT vs factor for "<<"MwpcC"<<j/4+2<<(char)((j%4)/2+88)<<(char)((j%4%2)*6+76);
      gr[i][j]->SetTitle(sstr_title.str().c_str());
      sstr_title.str("");

      gr[i][j]->GetXaxis()->SetLimits(-0.5,0.5*nPoint);
      gr[i][j]->GetXaxis()->SetTitle("factor of #sigma");

      gr[i][j]->SetMinimum(0);
      gr[i][j]->SetMaximum(75);
      gr[i][j]->GetYaxis()->SetTitle("Occupation (%)");

      if(i==7) gr[i][j]->SetMarkerColor(9);
      else gr[i][j]->SetMarkerColor(1+i);
      gr[i][j]->SetMarkerStyle(20+i);
      gr[i][j]->SetMarkerSize(0.8);
    }
  }

  TCanvas *canv=new TCanvas("Occupation vs NSigma","Occupation vs NSigma",2000,1500);
  canv->Divide(4,3);
  for(int j=0;j<12;j++){
    canv->cd(j+1);
    for(int i=0;i<8;i++){
      if(i==0){
        gr[i][j]->Draw("ap");
      }
      else{
        gr[i][j]->Draw("p");
      }
    }
  }

  TLegend *legend=new TLegend(0.9,0.6,1,1);
  legend->AddEntry(gr[0][0],"NoH","p");
  legend->AddEntry(gr[1][0],"S","p");
  legend->AddEntry(gr[2][0],"D","p");
  legend->AddEntry(gr[3][0],"T","p");
  legend->AddEntry(gr[4][0],"SD","p");
  legend->AddEntry(gr[5][0],"ST","p");
  legend->AddEntry(gr[6][0],"DT","p");
  legend->AddEntry(gr[7][0],"SDT","p");

  canv->cd(4);
  legend->Draw();

  ostringstream sstr;
  sstr<<"Occupation vs "<<nPoint*0.5<<"Sigma.pdf";
  canv->SaveAs(sstr.str().c_str());
  sstr.str("");

  double singMiniSigma[12],singMiniOccp[12];

  int tempSingMini[12];
  for(int j=0;j<12;j++){
    for(int k=5;k<nPoint;k++){
      if(k==5) tempSingMini[j]=eventType[1][j][k];
      else{
        if(tempSingMini[j]>eventType[1][j][k]) tempSingMini[j]=eventType[1][j][k];
        else{
          singMiniSigma[j]=(k-1)*0.5;
          singMiniOccp[j]=(double)tempSingMini[j]/totNEvent*100;
          break;
        }
      }
    }
  }

  TMarker *markSingMini[12];
  for(int j=0;j<12;j++){
    markSingMini[j]=new TMarker(singMiniSigma[j],singMiniOccp[j],20);
    markSingMini[j]->SetMarkerColor(kBlack);
    markSingMini[j]->SetMarkerSize(1.5);
  }

  TCanvas *canvSingMini=new TCanvas("Occupation vs NSigma for Singlet","Occupation vs NSigma for Singlet",2000,1500);
  canvSingMini->Divide(4,3);
  for(int j=0;j<12;j++){
    canvSingMini->cd(j+1);
    gr[1][j]->SetMaximum(30);
    gr[1][j]->Draw("ap");
    markSingMini[j]->Draw(); 
  }

  sstr<<"Occupation vs "<<"Sigma for Singlet for First Approach.pdf";
  canvSingMini->SaveAs(sstr.str().c_str());
  sstr.str("");


  ofstream output("MwpcEventThresholdFirstApproach.h");
  output<<"const double singMiniSigmaEventFirstApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<singMiniSigma[j]<<", ";
    else
      output<<singMiniSigma[j]<<"};"<<endl;
  }

  output<<"const double singMiniOccpEventFirstApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<singMiniOccp[j]<<", ";
    else
      output<<singMiniOccp[j]<<"};"<<endl;
  }
  output.close();

  return 0;
}
