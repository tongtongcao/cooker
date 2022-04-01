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
#include "MwpcEventThresholdFirstApproach.h"
using namespace std;

const static int nPoint=10;
static int totNEvent=0;
static int singNEvent[12][nPoint]={{0}};
Long_t Det_MWPC::histos_eventthresholdSecondApproach(){
  return 0;
}

Long_t Det_MWPC::startup_eventthresholdSecondApproach(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  return 0;
};


Long_t Det_MWPC::process_eventthresholdSecondApproach(){
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
      nHit[i][j]=event[i]->getNHit(data[i],pedValue[i],pedSigma[i],singMiniSigmaEventFirstApproach[i]-0.5+0.1*j);
      mult[i][j]=event[i]->multiplet3Type(data[i],pedValue[i],pedSigma[i],singMiniSigmaEventFirstApproach[i]-0.5+0.1*j);
    }
  }

  for(int i=0;i<12;i++){
    for(int j=0;j<nPoint;j++){
      if(mult[i][j][0]!=0 && mult[i][j][1]==0 && mult[i][j][2]==0) singNEvent[i][j]++;
    }
  }

  return 0;
}


Long_t Det_MWPC::done_eventthresholdSecondApproach(){
  double singOcpt[12][nPoint];
  double singOcpt_err[12][nPoint];
    for(int j=0;j<12;j++){
      for(int k=0;k<nPoint;k++){
        singOcpt[j][k]=(double)singNEvent[j][k]/totNEvent*100;
        singOcpt_err[j][k]=1./totNEvent*sqrt(singNEvent[j][k]*(totNEvent-singNEvent[j][k])/(double)totNEvent)*100;
    }
  }

  double singX[12][nPoint],singX_err[12][nPoint];
  for(int j=0;j<12;j++){
    for(int k=0;k<nPoint;k++){
      singX[j][k]=singMiniSigmaEventFirstApproach[j]-0.5+0.1*k;
      singX_err[j][k]=0;
    }
  }

  TGraphErrors *grSing[12];
  for(int j=0;j<12;j++){
    grSing[j]=new TGraphErrors(nPoint,singX[j],singOcpt[j],singX_err[j],singOcpt_err[j]);
  }


  ostringstream sstr_title;
  for(int j=0;j<12;j++){
    sstr_title<<"OCPT vs factor for Singlet of "<<"MwpcC"<<j/4+2<<(char)((j%4)/2+88)<<(char)((j%4%2)*6+76);
    grSing[j]->SetTitle(sstr_title.str().c_str());
    sstr_title.str("");

    grSing[j]->GetXaxis()->SetLimits(singMiniSigmaEventFirstApproach[j]-0.5-0.1,singMiniSigmaEventFirstApproach[j]+0.5+0.1);
    grSing[j]->GetXaxis()->SetTitle("factor of #sigma");

    grSing[j]->SetMinimum(0);
    grSing[j]->SetMaximum(15);
    grSing[j]->GetYaxis()->SetTitle("Occupation (%)");

    grSing[j]->SetMarkerColor(2);
    grSing[j]->SetMarkerStyle(21);
    grSing[j]->SetMarkerSize(0.8);
  }

  double singMiniSigma[12],singMiniOccp[12];

  int tempSingMini[12];
  for(int j=0;j<12;j++){
    for(int k=0;k<nPoint;k++){
      if(k==0) tempSingMini[j]=singNEvent[j][k];
      else{
        if(tempSingMini[j]>singNEvent[j][k]) tempSingMini[j]=singNEvent[j][k];
        else{
          singMiniSigma[j]=singMiniSigmaEventFirstApproach[j]-0.5+0.1*(k-1);
          singMiniOccp[j]=(double)tempSingMini[j]/totNEvent*100;
          break;
        }
      }
    }
  }

  TMarker *markSing[12];
  for(int j=0;j<12;j++){
    markSing[j]=new TMarker(singMiniSigma[j],singMiniOccp[j],20);
    markSing[j]->SetMarkerColor(kBlack);
    markSing[j]->SetMarkerSize(1.5);
  }

  TCanvas *canvSing=new TCanvas("Occupation vs NSigma for Singlet","Occupation vs NSigma for Singlet",2000,1500);
  canvSing->Divide(4,3);
  for(int j=0;j<12;j++){
    canvSing->cd(j+1);
    grSing[j]->Draw("ap");
    markSing[j]->Draw();
  }

  ostringstream sstr;
  sstr<<"Occupation vs "<<"Sigma for Singlet for Second Approach.pdf";
  canvSing->SaveAs(sstr.str().c_str());
  sstr.str("");

  ofstream output("MwpcEventThresholdSecondApproach.h");
  output<<"const double singMiniSigmaEventSecondApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<singMiniSigma[j]<<", ";
    else
      output<<singMiniSigma[j]<<"};"<<endl;
  }

  output<<"const double singMiniOccpEventSecondApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<singMiniOccp[j]<<", ";
    else
      output<<singMiniOccp[j]<<"};"<<endl;
  }
  output.close();

  return 0;
};
