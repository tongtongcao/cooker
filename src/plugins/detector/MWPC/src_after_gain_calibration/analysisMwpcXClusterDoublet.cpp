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
#include "MwpcPedestalValue.h"
#include "MwpcClusterThresholdSecondApproach.h"
using namespace std;

static TH2D *histXDoub[6];
static TH1D *histXDoubCutL[6];
static int clusType[6][2]={{0}};
Long_t Det_MWPC::histos_xclusterdoub(){
  ostringstream sstr_name, sstr_title;
  for(int i=0;i<6;i++){
    sstr_name<<"Total charge for L and R of Cut for "<<"MwpcC"<<i/2+2<<'X'<<(char)((i%2%2)*6+76)<<" for XDoublet";
    sstr_title<<sstr_name.str()<<";L and R of Cut;Total charge";
    histXDoub[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),2,0.5,2.5,100,0,3000);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Total charge distribution for L of Cut for "<<"MwpcC"<<i/2+2<<'X'<<(char)((i%2%2)*6+76)<<" for XDoublet";
    sstr_title<<sstr_name.str()<<";Total charge;Counts";
    histXDoubCutL[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,0,3000);
    sstr_name.str("");
    sstr_title.str("");
  }

  return 0;
}

Long_t Det_MWPC::startup_xclusterdoub(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  return 0;
};


Long_t Det_MWPC::process_xclusterdoub(){
  int* data[6]={treeRaw->ADC_C2X_L,treeRaw->ADC_C2X_R,treeRaw->ADC_C3X_L,treeRaw->ADC_C3X_R,treeRaw->ADC_C4X_L,treeRaw->ADC_C4X_R};
  const double* pedValue[6]={pedValueC2XL,pedValueC2XR,pedValueC3XL,pedValueC3XR,pedValueC4XL,pedValueC4XR};
  const double* pedSigma[6]={pedSigmaC2XL,pedSigmaC2XR,pedSigmaC3XL,pedSigmaC3XR,pedSigmaC4XL,pedSigmaC4XR}; 

  MwpcEvent *event[6];
  ostringstream sstr;
  vector<vector<int> > clusStNumNoCut,clusStNumCut;
  double *totChNoCut[6],*totChCut[6];
  for(int i=0;i<6;i++){
    event[i]=new MwpcEvent((i/2+7)*8,(i/2)!=0);

    totChCut[i]=event[i]->getClusterTotCharge(event[i],data[i],pedValue[i],pedSigma[i],doubMiniSigmaXClusterSecondApproach[2*i-i%2]);
    clusStNumCut=event[i]->clusterStripNumber(data[i],pedValue[i],pedSigma[i],doubMiniSigmaXClusterSecondApproach[2*i-i%2]);
 
    for(int j=0;j<(int)clusStNumCut.size();j++){
      if((int)clusStNumCut[j].size()==2) histXDoub[i]->Fill(2,totChCut[i][j]);
    }

    totChNoCut[i]=event[i]->getClusterTotCharge(event[i],data[i],pedValue[i],pedSigma[i],0);
    clusStNumNoCut=event[i]->clusterStripNumber(data[i],pedValue[i],pedSigma[i],0);
    for(int j=0;j<(int)clusStNumNoCut.size();j++){
      if((int)clusStNumNoCut[j].size()==2){
        for(int k=0;k<(int)clusStNumCut.size();k++){
          if(clusStNumCut[k].size()==2){
            if(clusStNumNoCut[j][0]!=clusStNumCut[k][0]){
              histXDoub[i]->Fill(1,totChNoCut[i][j]);
              histXDoubCutL[i]->Fill(totChNoCut[i][j]);
            }
          }
        }
      }
    }
  }

  return 0;
}

Long_t Det_MWPC::done_xclusterdoub(){

  TCanvas *canvXDoub=new TCanvas("Total charge for L and R of Cut for XDoublet","Total charge for L and R of Cut for XDoublet",1000,900);
  canvXDoub->Divide(2,3);
  for(int i=0;i<6;i++){
    canvXDoub->cd(i+1);
    gPad->SetLogz();
    histXDoub[i]->GetXaxis()->SetLabelOffset(99);
    histXDoub[i]->GetXaxis()->SetTitleOffset(0.5);
    histXDoub[i]->GetXaxis()->SetTitleSize(0.06);
    histXDoub[i]->GetYaxis()->SetLabelSize(0.06);
    histXDoub[i]->GetYaxis()->SetTitleOffset(0.8);
    histXDoub[i]->GetYaxis()->SetTitleSize(0.06);
    histXDoub[i]->GetZaxis()->SetLabelSize(0.045);
    histXDoub[i]->Draw("colz");
    double x, y;
    string type[2]={"L","R"};
    y = gPad->GetUymin() + 1*histXDoub[i]->GetYaxis()->GetBinWidth(1);
    TText t;
    t.SetTextAngle(60);
    t.SetTextSize(0.06);
    t.SetTextAlign(33);
    t.SetTextColor(kRed);
    for (int j=0;j<2;j++) {
      x = histXDoub[i]->GetXaxis()->GetBinCenter(j+1);
      t.DrawText(x,y,type[j].c_str());
    }
  }

  canvXDoub->SaveAs("total charge for L and R of Cut for XDoublet.pdf");  
  
  int temp=0;
  int xDoubChargeCut[6];
  for(int i=0;i<6;i++){
    for(int j=5;j<3000;j++){
      if(j==5) temp=histXDoubCutL[i]->GetBinContent(j+1);
      else{
        if(temp>histXDoubCutL[i]->GetBinContent(j+1)) temp=histXDoubCutL[i]->GetBinContent(j+1);
        else{
          xDoubChargeCut[i]=j;
          break;
       }
      }
    }
  } 
  
  TMarker *marker[6];
  for(int i=0;i<6;i++){
    marker[i]=new TMarker(xDoubChargeCut[i]*3000/100,histXDoubCutL[i]->GetBinContent(xDoubChargeCut[i]),20);
    marker[i]->SetMarkerSize(1.5);
    marker[i]->SetMarkerColor(kRed);
  }

  TCanvas *canvXDoubCutL=new TCanvas("Total charge distribution for L of Cut for XDoublet","Total charge distribution for L and R of Cut for XDoublet",1000,900);
  canvXDoubCutL->Divide(2,3);
  for(int i=0;i<6;i++){
    canvXDoubCutL->cd(i+1);
    gPad->SetLogy();
    histXDoubCutL[i]->Draw();
    marker[i]->Draw();
  }
  canvXDoubCutL->SaveAs("total charge distribution for L Cut for XDoublet.pdf");

  ofstream output("MwpcXClusterTotalChargeCut.h");
  output<<"const double doubMiniTotalChargeXCluster[6]={";
  for(int j=0;j<6;j++){
    if(j<5)
      output<<xDoubChargeCut[j]*3000/100<<", ";
    else
      output<<xDoubChargeCut[j]*3000/100<<"};"<<endl;
  }
  output.close();

  return 0;
};
