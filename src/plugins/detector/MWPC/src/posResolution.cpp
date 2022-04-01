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

static TH2D *histPosVsTCSing;
static TH2D *histPosVsTCDoub;
static TH2D *histPosVsTCTrip;
static TH1D *histPosCentroid;
static TH1D *histPosChargeRatio;
static TH1D *histDiffPos;
static TH1D *histPosLocal[12];
static TH1D *histPosLocalGap[12][6];
Long_t Det_MWPC::histos_posResolution(){
  ostringstream sstr_name,sstr_title;
  histPosVsTCSing=new TH2D("Position vs total-charge for singlet","strip 30 of C3XR;Total-charge;Position",50,0,3000,50,29,30);
  histPosVsTCDoub=new TH2D("Position vs total-charge for doublet","strip 30 of C3XR;Total-charge;Position",50,0,3000,50,29,30);
  histPosVsTCTrip=new TH2D("Position vs total-charge for triplet","strip 30 of C3XR;Total-charge;Position",50,0,3000,50,29,30);

  histPosCentroid=new TH1D("Position by Centroid","Position;Position;Counts",1120,0,56);
  histPosChargeRatio=new TH1D("Position by Charge Ratio","Position;Position;Counts",1120,0,56);
  histDiffPos=new TH1D("Difference of Position","Difference of Position;Diffrencen of Position;Counts",500,-0.1,0.1);

  for(int i=0;i<12;i++){
    sstr_name<<"Position for C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Position (mm); Counts";
    histPosLocal[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*40,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*40);
    sstr_name.str("");
    sstr_title.str("");

    for(int j=0;j<6;j++){
      sstr_name<<"Position for gap "<<(!(i%2))*6+j+1<<" of C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
      sstr_title<<"Gap "<<(!(i%2))*6+j+1<<";Position (mm); Counts";
      histPosLocalGap[i][j]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*40,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*40);
      sstr_name.str("");
      sstr_title.str("");
    }
  }




  return 0;
}

Long_t Det_MWPC::startup_posResolution(){
  getBranchObject("CookMwpcInfoV1",(TObject **) &treeCookV1);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);


  return 0;
}

Long_t Det_MWPC::process_posResolution(){
  if(treeCookV1->flagLeftRight!=0){
    if(treeCookV1->numSt[1]==1) histPosVsTCSing->Fill(treeCookV1->totalCharge[1],treeCookV1->posChargeRatio[1]);
    if(treeCookV1->numSt[1]==2) histPosVsTCDoub->Fill(treeCookV1->totalCharge[1],treeCookV1->posChargeRatio[1]);
    if(treeCookV1->numSt[1]==3) histPosVsTCTrip->Fill(treeCookV1->totalCharge[1],treeCookV1->posChargeRatio[1]);

    if(treeCookV1->numSt[1]>1 && treeCookV1->gapNum==2){
      histPosCentroid->Fill(treeCookV1->posCentroid[1]);
      histPosChargeRatio->Fill(treeCookV1->posChargeRatio[1]);
      histDiffPos->Fill(treeCookV1->posChargeRatio[1]-treeCookV1->posCentroid[1]);
    }
    for(int i=0;i<12;i++){
      if(i==9){
        if(treeCookV1->posChargeRatioLocal[i]!=-1000 && treeCookV1->flagDeadStripCluster==0){
          histPosLocal[i]->Fill(treeCookV1->posChargeRatioLocal[i]);
          for(int j=0;j<6;j++){
            if(treeCookV1->gapNum==(!(i%2))*6+j+1) histPosLocalGap[i][j]->Fill(treeCookV1->posChargeRatioLocal[i]);
          }
        }
      }
      else{
        if(treeCookV1->posChargeRatioLocal[i]!=-1000){
          histPosLocal[i]->Fill(treeCookV1->posChargeRatioLocal[i]);
          for(int j=0;j<6;j++){
            if(treeCookV1->gapNum==(!(i%2))*6+j+1) histPosLocalGap[i][j]->Fill(treeCookV1->posChargeRatioLocal[i]);
          }
        }
      }
    }
  }

  return 0;
}

Long_t Det_MWPC::done_posResolution(){
  ostringstream sstr;

  TCanvas *canvSing=new TCanvas("Position vs total-charge for singlet","Position vs total-charge for singlet",1000,1000);
  histPosVsTCSing->Draw("colz");
  canvSing->SaveAs("Position vs total-charge for singlet.pdf");

  TCanvas *canvDoub=new TCanvas("Position vs total-charge for doublet","Position vs total-charge for doublet",1000,1000);
  histPosVsTCDoub->Draw("colz");
  canvDoub->SaveAs("Position vs total-charge for doublet.pdf");

  TCanvas *canvTrip=new TCanvas("Position vs total-charge for triplet","Position vs total-charge for triplet",1000,1000);
  histPosVsTCTrip->Draw("colz");
  canvTrip->SaveAs("Position vs total-charge for triplet.pdf");

  TCanvas *canvPos=new TCanvas("Position","Position",1000,300);
  histPosCentroid->Draw();
  histPosChargeRatio->SetLineColor(kRed);
  histPosChargeRatio->Draw("same");
  TLegend *legend=new TLegend(0.75,0.7,0.9,0.9);
  legend->AddEntry(histPosCentroid,"Centroid","l");
  legend->AddEntry(histPosChargeRatio,"Charge Ratio","l");
  legend->Draw();
  canvPos->SaveAs("Position.pdf");

  TCanvas *canvDiffPos=new TCanvas("Difference of Position","Difference of Position",1000,600);
  histDiffPos->Draw("same");
  canvDiffPos->SaveAs("Difference of Position.pdf");

  TFile *file=new TFile("hist.root","recreate");

  TCanvas *canvPosLocal=new TCanvas("Local position","Local position",2000,900);
  canvPosLocal->Divide(4,3);
  for(int i=0;i<12;i++){
    canvPosLocal->cd(i+1);
    histPosLocal[i]->Draw();
    histPosLocal[i]->Write();
  }
  canvPosLocal->SaveAs("Local position.pdf");

  TCanvas *canvPosLocalGap[12];
  for(int i=0;i<12;i++){
    sstr<<"Position for C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    canvPosLocalGap[i]=new TCanvas(sstr.str().c_str(),sstr.str().c_str(),1000,900);
    sstr.str("");
    canvPosLocalGap[i]->Divide(2,3);
    for(int j=0;j<6;j++){
      canvPosLocalGap[i]->cd(j+1);
      histPosLocalGap[i][j]->Draw();
histPosLocalGap[i][j]->Write();
    }
    sstr<<"Position for C"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<".pdf";
    canvPosLocalGap[i]->SaveAs(sstr.str().c_str());
    sstr.str("");
  }

file->Write();
  return 0;
}

