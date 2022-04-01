#include <Det_MWPC.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <sstream>
#include "TLine.h"
#include <string>
#include "TH1D.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TTree.h"
#include "TF1.h"
using namespace std;

static TH1D *histNClus[12];
static TH1D *histNSt[12];
static TH1D *histDoubChargeRatio[12];
static TH1D *histTripChargeRatio[12];

static int numClusterOnlyX[6],numClusterOnlyY[6]={0},numClusterBoth[6]={0};
Long_t Det_MWPC::histos_analysisExtractedCluster3(){
  ostringstream sstr_name,sstr_title;
  for(int i=0;i<12;i++){
    sstr_name<<"Number of clusters for C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Number of clusters; Counts";
    histNClus[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),7,-0.5,6.5);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Number of strips for C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Number of strips; Counts";
    histNSt[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),10,0.5,10.5);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Charge ratio for doublet of C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Charge ratio; Counts";
    histDoubChargeRatio[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),19,0,20);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Charge ratio for Triplet of C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Charge ratio; Counts";
    histTripChargeRatio[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),19,0,20);
    sstr_name.str("");
    sstr_title.str("");
  }
  return 0;
}

Long_t Det_MWPC::startup_analysisExtractedCluster3(){
  getBranchObject("CaliMwpcInfo",(TObject **) &treeCali);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);

  return 0;
};

Long_t Det_MWPC::process_analysisExtractedCluster3(){
  for(int i=0;i<6;i++){
    if(treeCali->thTCCutNumCluster[i/2*4+i%2]!=0 && treeCali->thTCCutNumCluster[i/2*4+i%2+2]==0) numClusterOnlyX[i]++;
    if(treeCali->thTCCutNumCluster[i/2*4+i%2]==0 && treeCali->thTCCutNumCluster[i/2*4+i%2+2]!=0) numClusterOnlyY[i]++;
    if(treeCali->thTCCutNumCluster[i/2*4+i%2]==0 && treeCali->thTCCutNumCluster[i/2*4+i%2+2]==0) numClusterBoth[i]++;
  }

  for(int i=0;i<12;i++){
    histNClus[i]->Fill(treeCali->thTCCutNumCluster[i]);
    for(UInt_t j=0;j<treeCali->thTCCutClusterStNum[i].size();j++){
      histNSt[i]->Fill(treeCali->thTCCutClusterStNum[i][j].size());
      if(treeCali->thTCCutClusterStNum[i][j].size()==2){
        if(treeCali->thTCCutClusterStCharge[i][j][0]>treeCali->thTCCutClusterStCharge[i][j][1]) histDoubChargeRatio[i]->Fill(treeCali->thTCCutClusterStCharge[i][j][0]/treeCali->thTCCutClusterStCharge[i][j][1]);
        else histDoubChargeRatio[i]->Fill(treeCali->thTCCutClusterStCharge[i][j][1]/treeCali->thTCCutClusterStCharge[i][j][0]);
      }
      if(treeCali->thTCCutClusterStNum[i][j].size()==3){
        if(treeCali->thTCCutClusterStCharge[i][j][0]>treeCali->thTCCutClusterStCharge[i][j][1] && treeCali->thTCCutClusterStCharge[i][j][0]>treeCali->thTCCutClusterStCharge[i][j][2]){
          if(treeCali->thTCCutClusterStCharge[i][j][1]>treeCali->thTCCutClusterStCharge[i][j][2]) histTripChargeRatio[i]->Fill(treeCali->thTCCutClusterStCharge[i][j][0]/treeCali->thTCCutClusterStCharge[i][j][1]);
          else histTripChargeRatio[i]->Fill(treeCali->thTCCutClusterStCharge[i][j][0]/treeCali->thTCCutClusterStCharge[i][j][2]);
        }
        if(treeCali->thTCCutClusterStCharge[i][j][1]>treeCali->thTCCutClusterStCharge[i][j][0] && treeCali->thTCCutClusterStCharge[i][j][1]>treeCali->thTCCutClusterStCharge[i][j][2]){
          if(treeCali->thTCCutClusterStCharge[i][j][0]>treeCali->thTCCutClusterStCharge[i][j][2]) histTripChargeRatio[i]->Fill(treeCali->thTCCutClusterStCharge[i][j][1]/treeCali->thTCCutClusterStCharge[i][j][0]);
          else histTripChargeRatio[i]->Fill(treeCali->thTCCutClusterStCharge[i][j][1]/treeCali->thTCCutClusterStCharge[i][j][2]);
        }
        if(treeCali->thTCCutClusterStCharge[i][j][2]>treeCali->thTCCutClusterStCharge[i][j][0] && treeCali->thTCCutClusterStCharge[i][j][2]>treeCali->thTCCutClusterStCharge[i][j][1]){
          if(treeCali->thTCCutClusterStCharge[i][j][0]>treeCali->thTCCutClusterStCharge[i][j][1]) histTripChargeRatio[i]->Fill(treeCali->thTCCutClusterStCharge[i][j][2]/treeCali->thTCCutClusterStCharge[i][j][0]);
          else histTripChargeRatio[i]->Fill(treeCali->thTCCutClusterStCharge[i][j][2]/treeCali->thTCCutClusterStCharge[i][j][1]);
        }
      }
    }
  }  
  return 0;
}

Long_t Det_MWPC::done_analysisExtractedCluster3(){
  ostringstream sstr_name,sstr_title;
  TH1D *histCompNClusXY[6];
  for(int i=0;i<6;i++){
    sstr_name<<"Comparison for number of clusters between X and Y for C"<<i/2+2<<(char)((i%2)*6+76)<<"after Total-Charge and Threshold Cuts";
    sstr_title<<sstr_name.str()<<";X/Y/Both;Counts";
    histCompNClusXY[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),3,0,3);
    sstr_name.str("");
    sstr_title.str("");
  }

  for(int i=0;i<6;i++){ 
    histCompNClusXY[i]->SetBinContent(1,numClusterOnlyX[i]);
    histCompNClusXY[i]->SetBinContent(2,numClusterOnlyY[i]);
    histCompNClusXY[i]->SetBinContent(3,numClusterBoth[i]);
  }

  TCanvas *canvCompNClusXY=new TCanvas("Comparison for number of clusters between X and Y after Total-Charge and Threshold Cuts","Comparison for number of clusters between X and Y after Total-Charge and Threshold Cuts",1000,900);
  canvCompNClusXY->Divide(2,3);
  for(int i=0;i<6;i++){
    canvCompNClusXY->cd(i+1);
    histCompNClusXY[i]->SetMinimum(0);
    histCompNClusXY[i]->SetMaximum(60000);
    histCompNClusXY[i]->Draw();
    histCompNClusXY[i]->GetXaxis()->SetLabelOffset(99);
    double x, y;
    string type[3]={"Only in X","Only in Y","Two Sides"};
    y = gPad->GetUymin() + 5000*histCompNClusXY[i]->GetYaxis()->GetBinWidth(1);
    TText t;
    t.SetTextAngle(30);
    t.SetTextSize(0.06);
    t.SetTextAlign(33);
    t.SetTextColor(kRed);
    for (int j=0;j<3;j++) {
      x = histCompNClusXY[i]->GetXaxis()->GetBinCenter(j+1)+0.3;
      t.DrawText(x,y,type[j].c_str());
  }
  }
  canvCompNClusXY->SaveAs("Comparison for number of clusters between X and Y after Total-Charge and Threshold Cuts.pdf");

  TCanvas *canvNClus=new TCanvas("Distribution of number of clusters","Distribution of number of clusters",2000,900);
  canvNClus->Divide(4,3);
  for(int i=0;i<12;i++){
    canvNClus->cd(i+1);
    histNClus[i]->Draw();
  }
  canvNClus->SaveAs("Distribution of number of clusters.pdf");

  TCanvas *canvNSt=new TCanvas("Distribution of number of strips","Distribution of number of strips",2000,900);
  canvNSt->Divide(4,3);
  for(int i=0;i<12;i++){
    canvNSt->cd(i+1);
    histNSt[i]->Draw();
  }
  canvNSt->SaveAs("Distribution of number of strips.pdf");

  TCanvas *canvDoubChargeRatio=new TCanvas("Distribution of charge ratio for Doublet","Distribution of charge ratio for Doublet",2000,900);
  canvDoubChargeRatio->Divide(4,3);
  for(int i=0;i<12;i++){
    canvDoubChargeRatio->cd(i+1);
    histDoubChargeRatio[i]->Draw();
  }
  canvDoubChargeRatio->SaveAs("Distribution of charge ratio for Doublet.pdf");

  TCanvas *canvTripChargeRatio=new TCanvas("Distribution of charge ratio for Triplet","Distribution of charge ratio for Triplet",2000,900);
  canvTripChargeRatio->Divide(4,3);
  for(int i=0;i<12;i++){
    canvTripChargeRatio->cd(i+1);
    histTripChargeRatio[i]->Draw();
  }
  canvTripChargeRatio->SaveAs("Distribution of charge ratio for Triplet.pdf");
  return 0;
};

