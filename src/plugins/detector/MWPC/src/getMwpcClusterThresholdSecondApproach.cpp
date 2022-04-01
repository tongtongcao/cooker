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
#include "MwpcClusterThresholdFirstApproach.h"
using namespace std;

const static int nPoint=10;
const static int nType=4;
static int singNClus[12][nPoint]={{0}},doubNClus[12][nPoint]={{0}};
static int clusNoCut[nType][12]={{0}};
Long_t Det_MWPC::histos_clusterthresholdSecondApproach(){
  return 0;
}

Long_t Det_MWPC::startup_clusterthresholdSecondApproach(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  return 0;
};


Long_t Det_MWPC::process_clusterthresholdSecondApproach(){
  int* data[12]={treeRaw->ADC_C2X_L,treeRaw->ADC_C2X_R,treeRaw->ADC_C2Y_L,treeRaw->ADC_C2Y_R,treeRaw->ADC_C3X_L,treeRaw->ADC_C3X_R,treeRaw->ADC_C3Y_L,treeRaw->ADC_C3Y_R,treeRaw->ADC_C4X_L,treeRaw->ADC_C4X_R,treeRaw->ADC_C4Y_L,treeRaw->ADC_C4Y_R};
  const double* pedValue[12]={pedValueC2XL,pedValueC2XR,pedValueC2YL,pedValueC2YR,pedValueC3XL,pedValueC3XR,pedValueC3YL,pedValueC3YR,pedValueC4XL,pedValueC4XR,pedValueC4YL,pedValueC4YR};
  const double* pedSigma[12]={pedSigmaC2XL,pedSigmaC2XR,pedSigmaC2YL,pedSigmaC2YR,pedSigmaC3XL,pedSigmaC3XR,pedSigmaC3YL,pedSigmaC3YR,pedSigmaC4XL,pedSigmaC4XR,pedSigmaC4YL,pedSigmaC4YR}; 

  MwpcEvent *event[12];
  ostringstream sstr;
  int *multSing[12], *multDoub[12];
  int *multNoCut[12];
  for(int i=0;i<12;i++){
    event[i]=new MwpcEvent(((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,!((i/4==0)&&((i%4/2))==0));
    multNoCut[i]=event[i]->multipletNType(data[i],pedValue[i],pedSigma[i],0,nType);
    for(int k=0;k<nType;k++){
      clusNoCut[k][i]+=multNoCut[i][k];
    }
    for(int j=0;j<nPoint;j++){
      multSing[i]=event[i]->multipletNType(data[i],pedValue[i],pedSigma[i],singMiniSigmaClusterFirstApproach[i]+0.1*j-0.5,nType);
      singNClus[i][j]+=multSing[i][0];
      multDoub[i]=event[i]->multipletNType(data[i],pedValue[i],pedSigma[i],doubMiniSigmaClusterFirstApproach[i]+0.1*j-0.5,nType);
      doubNClus[i][j]+=multDoub[i][1];
    }
  }
  
  return 0;
}


Long_t Det_MWPC::done_clusterthresholdSecondApproach(){
  ostringstream sstr_title,sstr;

  int totNClus[12]={0};
  for(int i=0;i<12;i++){
    for(int k=0;k<nType;k++){
      totNClus[i]+=clusNoCut[k][i];
    }
  }

  double singOcpt[12][nPoint];
  double singOcpt_err[12][nPoint];
    for(int i=0;i<12;i++){
      for(int j=0;j<nPoint;j++){
        singOcpt[i][j]=(double)singNClus[i][j]/totNClus[i]*100;
        singOcpt_err[i][j]=1./totNClus[i]*sqrt(singNClus[i][j]*(totNClus[i]-singNClus[i][j])/(double)totNClus[i])*100;
    }
  }

  double singX[12][nPoint],singX_err[12][nPoint];
  for(int i=0;i<12;i++){
    for(int j=0;j<nPoint;j++){
      singX[i][j]=singMiniSigmaClusterFirstApproach[i]-0.5+0.1*j;
      singX_err[i][j]=0;
    }
  }

  TGraphErrors *grSing[12];
  for(int i=0;i<12;i++){
    grSing[i]=new TGraphErrors(nPoint,singX[i],singOcpt[i],singX_err[i],singOcpt_err[i]);
  }

  for(int i=0;i<12;i++){
    sstr_title<<"OCPT vs factor for Singlet of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    grSing[i]->SetTitle(sstr_title.str().c_str());
    sstr_title.str("");

    grSing[i]->GetXaxis()->SetLimits(singMiniSigmaClusterFirstApproach[i]-0.5-0.1,singMiniSigmaClusterFirstApproach[i]+0.5+0.1);
    grSing[i]->GetXaxis()->SetTitle("factor of #sigma");

    grSing[i]->SetMinimum(0);
    grSing[i]->SetMaximum(20);
    grSing[i]->GetYaxis()->SetTitle("Occupation (%)");

    grSing[i]->SetMarkerColor(1);
    grSing[i]->SetMarkerStyle(21);
    grSing[i]->SetMarkerSize(0.8);
  }

  double singMiniSigma[12],singMiniOccp[12];

  double tempSingMini[12];
  for(int i=0;i<12;i++){
    for(int j=0;j<nPoint;j++){
      if(j==0) tempSingMini[i]=singOcpt[i][j];
      else{
        if(tempSingMini[i]>singOcpt[i][j]+0.1) tempSingMini[i]=singOcpt[i][j];
        else{
          singMiniSigma[i]=singMiniSigmaClusterFirstApproach[i]-0.5+0.1*(j-1);
          singMiniOccp[i]= singOcpt[i][j-1];
          break;
        }
      }
    }
  }

  TMarker *markSing[12];
  for(int i=0;i<12;i++){
    markSing[i]=new TMarker(singMiniSigma[i],singMiniOccp[i],20);
    markSing[i]->SetMarkerColor(kRed);
    markSing[i]->SetMarkerSize(1.5);
  }

  TCanvas *canvSing=new TCanvas("Occupation vs NSigma for Singlet for Cluster","Occupation vs NSigma for Singlet for Cluster",2000,1500);
  canvSing->Divide(4,3);
  for(int i=0;i<12;i++){
    canvSing->cd(i+1);
    grSing[i]->Draw("ap");
    markSing[i]->Draw();
  }

  sstr<<"Occupation vs "<<"Sigma for Singlet for Second Approach for Cluster.pdf";
  canvSing->SaveAs(sstr.str().c_str());
  sstr.str("");

  double doubOcpt[12][nPoint];
  double doubOcpt_err[12][nPoint];
    for(int i=0;i<12;i++){
      for(int j=0;j<nPoint;j++){
        doubOcpt[i][j]=(double)doubNClus[i][j]/totNClus[i]*100;
        doubOcpt_err[i][j]=1./totNClus[i]*sqrt(doubNClus[i][j]*(totNClus[i]-doubNClus[i][j])/(double)totNClus[i])*100;
    }
  }

  double doubX[12][nPoint],doubX_err[12][nPoint];
  for(int i=0;i<12;i++){
    for(int j=0;j<nPoint;j++){
      doubX[i][j]=doubMiniSigmaClusterFirstApproach[i]-0.5+0.1*j;
      doubX_err[i][j]=0;
    }
  }

  TGraphErrors *grDoub[12];
  for(int i=0;i<12;i++){
    grDoub[i]=new TGraphErrors(nPoint,doubX[i],doubOcpt[i],doubX_err[i],doubOcpt_err[i]);
  }

  for(int i=0;i<12;i++){
    sstr_title<<"OCPT vs factor for Doublet of "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    grDoub[i]->SetTitle(sstr_title.str().c_str());
    sstr_title.str("");

    grDoub[i]->GetXaxis()->SetLimits(doubMiniSigmaClusterFirstApproach[i]-0.5-0.1,doubMiniSigmaClusterFirstApproach[i]+0.5+0.1);
    grDoub[i]->GetXaxis()->SetTitle("factor of #sigma");

    grDoub[i]->SetMinimum(0);
    grDoub[i]->SetMaximum(15);
    grDoub[i]->GetYaxis()->SetTitle("Occupation (%)");

    grDoub[i]->SetMarkerColor(2);
    grDoub[i]->SetMarkerStyle(21);
    grDoub[i]->SetMarkerSize(0.8);
  }

  double doubMiniSigma[12],doubMiniOccp[12];

  double tempDoubMini[12];
  for(int i=0;i<12;i++){
    for(int j=0;j<nPoint;j++){
      if(j==0) tempDoubMini[i]=doubOcpt[i][j];
      else{
        if(tempDoubMini[i]>doubOcpt[i][j]) tempDoubMini[i]=doubOcpt[i][j];
        else{
          doubMiniSigma[i]=doubMiniSigmaClusterFirstApproach[i]-0.5+0.1*(j-1);
          doubMiniOccp[i]= doubOcpt[i][j-1];
          break;
        }
      }
    }
  }

  TMarker *markDoub[12];
  for(int i=0;i<12;i++){
    markDoub[i]=new TMarker(doubMiniSigma[i],doubMiniOccp[i],20);
    markDoub[i]->SetMarkerColor(kBlack);
    markDoub[i]->SetMarkerSize(1.5);
  }

  TCanvas *canvDoub=new TCanvas("Occupation vs NSigma for Doublet for Cluster","Occupation vs NSigma for Doublet for Cluster",2000,1500);
  canvDoub->Divide(4,3);
  for(int i=0;i<12;i++){
    canvDoub->cd(i+1);
    grDoub[i]->Draw("ap");
    markDoub[i]->Draw();
  }

  sstr<<"Occupation vs "<<"Sigma for Doublet for Second Approach for Cluster.pdf";
  canvDoub->SaveAs(sstr.str().c_str());
  sstr.str("");
                                    
  ofstream output("MwpcClusterThresholdSecondApproach.h");
  output<<"const double singMiniSigmaClusterSecondApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<singMiniSigma[j]<<", ";
    else
      output<<singMiniSigma[j]<<"};"<<endl;
  }

  output<<"const double singMiniOccpClusterSecondApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<singMiniOccp[j]<<", ";
    else
      output<<singMiniOccp[j]<<"};"<<endl;
  }

  output<<"const double doubMiniSigmaClusterSecondApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<doubMiniSigma[j]<<", ";
    else
      output<<doubMiniSigma[j]<<"};"<<endl;
  }

  output<<"const double doubMiniOccpClusterSecondApproach[12]={";
  for(int j=0;j<12;j++){
    if(j<11)
      output<<doubMiniOccp[j]<<", ";
    else
      output<<doubMiniOccp[j]<<"};"<<endl;
  }
  output.close();

  return 0;
};
