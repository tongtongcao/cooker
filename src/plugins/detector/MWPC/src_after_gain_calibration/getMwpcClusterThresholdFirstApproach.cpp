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
const static int nType=4;
static int clusType[nType][12][nPoint]={{{0}}};
Long_t Det_MWPC::histos_clusterthresholdFirstApproach(){
  return 0;
}

Long_t Det_MWPC::startup_clusterthresholdFirstApproach(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  return 0;
};


Long_t Det_MWPC::process_clusterthresholdFirstApproach(){
  int* data[12]={treeRaw->ADC_C2X_L,treeRaw->ADC_C2X_R,treeRaw->ADC_C2Y_L,treeRaw->ADC_C2Y_R,treeRaw->ADC_C3X_L,treeRaw->ADC_C3X_R,treeRaw->ADC_C3Y_L,treeRaw->ADC_C3Y_R,treeRaw->ADC_C4X_L,treeRaw->ADC_C4X_R,treeRaw->ADC_C4Y_L,treeRaw->ADC_C4Y_R};
  const double* pedValue[12]={pedValueC2XL,pedValueC2XR,pedValueC2YL,pedValueC2YR,pedValueC3XL,pedValueC3XR,pedValueC3YL,pedValueC3YR,pedValueC4XL,pedValueC4XR,pedValueC4YL,pedValueC4YR};
  const double* pedSigma[12]={pedSigmaC2XL,pedSigmaC2XR,pedSigmaC2YL,pedSigmaC2YR,pedSigmaC3XL,pedSigmaC3XR,pedSigmaC3YL,pedSigmaC3YR,pedSigmaC4XL,pedSigmaC4XR,pedSigmaC4YL,pedSigmaC4YR}; 

  MwpcEvent *event[12];
  ostringstream sstr;
  int *mult[12];
  for(int i=0;i<12;i++){
    event[i]=new MwpcEvent(((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,!((i/4==0)&&((i%4/2))==0));
    for(int j=0;j<nPoint;j++){
      mult[i]=event[i]->multipletNType(data[i],pedValue[i],pedSigma[i],0.5*j,nType);
      for(int k=0;k<nType;k++){
        clusType[k][i][j]+=mult[i][k];
      }
    }
  }
  return 0;
}


Long_t Det_MWPC::done_clusterthresholdFirstApproach(){
  int totNClus[12]={0};
  for(int i=0;i<12;i++){
    for(int k=0;k<nType;k++){
      totNClus[i]+=clusType[k][i][0];
    }
  }
  double ocpt[nType][12][nPoint];
  double ocpt_err[nType][12][nPoint];
  for(int k=0;k<nType;k++){
    for(int i=0;i<12;i++){
      for(int j=0;j<nPoint;j++){
        ocpt[k][i][j]=(double)clusType[k][i][j]/totNClus[i]*100;
        ocpt_err[k][i][j]=1./totNClus[i]*sqrt(clusType[k][i][j]*(totNClus[i]-clusType[k][i][j])/(double)totNClus[i])*100;
      }
    }
  }
  double x[nPoint],x_err[nPoint];
  for(int j=0;j<nPoint;j++){
    x[j]=0.5*j;
    x_err[j]=0;
  }

  TGraphErrors *gr[nType][12];
  for(int k=0;k<nType;k++){
    for(int i=0;i<12;i++){
      gr[k][i]=new TGraphErrors(nPoint,x,ocpt[k][i],x_err,ocpt_err[k][i]);
    }
  }


  ostringstream sstr_title;
  for(int i=0;i<12;i++){
    for(int k=0;k<nType;k++){
      sstr_title<<"OCPT vs factor for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
      gr[k][i]->SetTitle(sstr_title.str().c_str());
      sstr_title.str("");

      gr[k][i]->GetXaxis()->SetLimits(-0.5,0.5*nPoint);
      gr[k][i]->GetXaxis()->SetTitle("factor of #sigma");

      gr[k][i]->SetMinimum(0);
      gr[k][i]->SetMaximum(75);
      gr[k][i]->GetYaxis()->SetTitle("Occupation (%)");

      gr[k][i]->SetMarkerColor(1+k);
      gr[k][i]->SetMarkerStyle(20+k);
      gr[k][i]->SetMarkerSize(0.8);
    }
  }

  TCanvas *canv=new TCanvas("Occupation vs NSigma for clusters","Occupation vs NSigma for clusters",2000,1500);
  canv->Divide(4,3);
  for(int i=0;i<12;i++){
    canv->cd(i+1);
    for(int k=0;k<nType;k++){
      if(k==0){
        gr[k][i]->Draw("ap");
      }
      else{
        gr[k][i]->Draw("p");
      }
    }
  }

  TLegend *legend=new TLegend(0.9,0.6,1,1);
  legend->AddEntry(gr[0][0],"S","p");
  legend->AddEntry(gr[1][0],"D","p");
  legend->AddEntry(gr[2][0],"T","p");
  legend->AddEntry(gr[3][0],"Q","p");

  canv->cd(4);
  legend->Draw();

  ostringstream sstr;
  sstr<<"Occupation vs "<<nPoint*0.5<<"Sigma for clusters.pdf";
  canv->SaveAs(sstr.str().c_str());
  sstr.str("");

  double singMiniSigma[12],singMiniOccp[12];

  double tempSingMini[12];
  for(int i=0;i<12;i++){
    for(int j=3;j<nPoint;j++){
      if(j==3) tempSingMini[i]=ocpt[0][i][j];
      else{
        if(tempSingMini[i]>ocpt[0][i][j]+0.5) tempSingMini[i]=ocpt[0][i][j];
        else{
          singMiniSigma[i]=(j-1)*0.5;
          singMiniOccp[i]=tempSingMini[i];
          break;
        }
      }
    }
  }

  TMarker *markSingMini[12];
  for(int i=0;i<12;i++){
    markSingMini[i]=new TMarker(singMiniSigma[i],singMiniOccp[i],20);
    markSingMini[i]->SetMarkerColor(kRed);
    markSingMini[i]->SetMarkerSize(1.5);
  }

  TCanvas *canvSingMini=new TCanvas("Occupation vs NSigma for Singlet for Cluster","Occupation vs NSigma for Singlet for Cluster",2000,1500);
  canvSingMini->Divide(4,3);
  for(int i=0;i<12;i++){
    canvSingMini->cd(i+1);
    gr[0][i]->SetMaximum(75);
    gr[0][i]->Draw("ap");
    markSingMini[i]->Draw();
  }

  sstr<<"Occupation vs "<<"Sigma for Singlet for First Approach for Cluster.pdf";
  canvSingMini->SaveAs(sstr.str().c_str());
  sstr.str("");

  double doubMiniSigma[12],doubMiniOccp[12];

  int tempDoubMini[12];
  for(int i=0;i<12;i++){
    for(int j=0;j<nPoint;j++){
      if(j==0) tempDoubMini[i]=clusType[1][i][j];
      else{
        if(tempDoubMini[i]>clusType[1][i][j]) tempDoubMini[i]=clusType[1][i][j];
        else{
          doubMiniSigma[i]=(j-1)*0.5;
          doubMiniOccp[i]=(double)tempDoubMini[i]/totNClus[i]*100;
          break;
        }
      }
    }
  }

  TMarker *markDoubMini[12];
  for(int i=0;i<12;i++){
    markDoubMini[i]=new TMarker(doubMiniSigma[i],doubMiniOccp[i],20);
    markDoubMini[i]->SetMarkerColor(kBlack);
    markDoubMini[i]->SetMarkerSize(1.5);
  }

  TCanvas *canvDoubMini=new TCanvas("Occupation vs NSigma for Doublet for Cluster","Occupation vs NSigma for Doublet for Cluster",2000,1500);
  canvDoubMini->Divide(4,3);
  for(int i=0;i<12;i++){
    canvDoubMini->cd(i+1);
    gr[1][i]->SetMaximum(50);
    gr[1][i]->Draw("ap");
    markDoubMini[i]->Draw();
  }

  sstr<<"Occupation vs "<<"Sigma for Doublet for First Approach for Cluster.pdf";
  canvDoubMini->SaveAs(sstr.str().c_str());
  sstr.str("");

  ofstream output("MwpcClusterThresholdFirstApproach.h");
  output<<"const double singMiniSigmaClusterFirstApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<singMiniSigma[j]<<", ";
    else
      output<<singMiniSigma[j]<<"};"<<endl;
  }

  output<<"const double singMiniOccpClusterFirstApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<singMiniOccp[j]<<", ";
    else
      output<<singMiniOccp[j]<<"};"<<endl;
  }

  output<<"const double doubMiniSigmaClusterFirstApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<doubMiniSigma[j]<<", ";
    else
      output<<doubMiniSigma[j]<<"};"<<endl;
  }

  output<<"const double doubMiniOccpClusterFirstApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<doubMiniOccp[j]<<", ";
    else
      output<<doubMiniOccp[j]<<"};"<<endl;
  }
  output.close();

  return 0;
}
