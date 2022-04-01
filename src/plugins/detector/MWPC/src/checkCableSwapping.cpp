#include <Det_MWPC.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
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

static string inputFilePath(getenv("mwpc_inputFiles"));
static string inputFile=inputFilePath+"/gapNumTable/run3994";
static ifstream input(inputFile.c_str());
static int runNum, eventNum, gapNumTOF1, gapNumTOF2;

static TH2D *histCompPosC3YC4Y[12];
static TH2D *histCompPosC3XC4X[12];
static TH2D *histCompPosC3XC4XAfterCableSwapping[12];
static TH2D *histCompPosC3YC4YWithTOFInequalEventsBasedOnTOF2[12];
static TH2D *histCompPosC3XC4XWithTOFInequalEventsBasedOnTOF2[12];
static TH2D *histCompPosC3YC4YWithTOFInequalEventsBasedOnTOF1[12];
static TH2D *histCompPosC3XC4XWithTOFInequalEventsBasedOnTOF1[12];
Long_t Det_MWPC::histos_checkCableSwapping(){
  ostringstream sstr_name,sstr_title;
  for(int i=0;i<12;i++){
    sstr_name<<"Comparion of position between C3Y and C4Y for gap "<<i+1;
    sstr_title<<sstr_name.str()<<";Position of C3Y; Position of C4Y";
    histCompPosC3YC4Y[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,0,17,100,0,17);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Comparion of position between C3X and C4X for gap "<<i+1;
    sstr_title<<sstr_name.str()<<";Position of C3X; Position of C4X";
    histCompPosC3XC4X[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,0,65,100,0,73);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Comparion of position between C3X and C4X for gap "<<i+1<<" after cable swapping";
    sstr_title<<sstr_name.str()<<";Position of C3X; Position of C4X";
    histCompPosC3XC4XAfterCableSwapping[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,0,65,100,0,73);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Comparion of position between C3Y and C4Y for gap "<<i+1<<"with TOF-inequal events based on TOF2";
    sstr_title<<sstr_name.str()<<";Position of C3Y; Position of C4Y";
    histCompPosC3YC4YWithTOFInequalEventsBasedOnTOF2[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,0,17,100,0,17);
    sstr_name.str("");
    sstr_title.str("");
    
    sstr_name<<"Comparion of position between C3X and C4X for gap "<<i+1<<"with TOF-inequal events based on TOF2";
    sstr_title<<sstr_name.str()<<";Position of C3X; Position of C4X";
    histCompPosC3XC4XWithTOFInequalEventsBasedOnTOF2[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,0,65,100,0,73);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Comparion of position between C3Y and C4Y for gap "<<i+1<<"with TOF-inequal events based on TOF1";
    sstr_title<<sstr_name.str()<<";Position of C3Y; Position of C4Y";
    histCompPosC3YC4YWithTOFInequalEventsBasedOnTOF1[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,0,17,100,0,17);
    sstr_name.str("");
    sstr_title.str("");
    
    sstr_name<<"Comparion of position between C3X and C4X for gap "<<i+1<<"with TOF-inequal events based on TOF1";
    sstr_title<<sstr_name.str()<<";Position of C3X; Position of C4X";
    histCompPosC3XC4XWithTOFInequalEventsBasedOnTOF1[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,0,65,100,0,73);
    sstr_name.str("");
    sstr_title.str("");
  }
    
  return 0;
}

Long_t Det_MWPC::startup_checkCableSwapping(){
  getBranchObject("CaliMwpcInfo",(TObject **) &treeCali);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);

  return 0;
};

static int numAcptEvent[7]={0};
static bool condTwoSide, condOneSide, condBad, condGood, condUnmatch, condMatch, condInequal, condEqual;

Long_t Det_MWPC::process_checkCableSwapping(){
  input>>runNum>>eventNum>>gapNumTOF1>>gapNumTOF2;
  if(runNum!=(int)treeCali->run || eventNum!=(int)treeCali->event){
    cout<<"Error: Run number or event number are not consistent between two files!"<<endl;
  }

  condTwoSide=treeCali->thTCCutNumCluster[4]==1 && treeCali->thTCCutNumCluster[6]==1 && treeCali->thTCCutNumCluster[8]==1 && treeCali->thTCCutNumCluster[10]==1 && treeCali->thTCCutNumCluster[5]==1 && treeCali->thTCCutNumCluster[7]==1 && treeCali->thTCCutNumCluster[9]==1 && treeCali->thTCCutNumCluster[11]==1;
  condOneSide=((treeCali->thTCCutNumCluster[4]==1 && treeCali->thTCCutNumCluster[6]==1 && treeCali->thTCCutNumCluster[8]==1 && treeCali->thTCCutNumCluster[10]==1) || (treeCali->thTCCutNumCluster[5]==1 && treeCali->thTCCutNumCluster[7]==1 && treeCali->thTCCutNumCluster[9]==1 && treeCali->thTCCutNumCluster[11]==1)) && !condTwoSide;

  condBad=(gapNumTOF2==-1);
  condGood=(gapNumTOF2!=-1);

  condUnmatch=treeCali->thTCCutNumCluster[4]==1 && treeCali->thTCCutNumCluster[6]==1 && treeCali->thTCCutNumCluster[8]==1 && treeCali->thTCCutNumCluster[10]==1 && (gapNumTOF2==1 || gapNumTOF2==2 || gapNumTOF2==3 || gapNumTOF2==4 || gapNumTOF2==5 || gapNumTOF2==6) && (gapNumTOF1==1 || gapNumTOF1==2 || gapNumTOF1==3 || gapNumTOF1==4 || gapNumTOF1==5 || gapNumTOF1==6) || treeCali->thTCCutNumCluster[5]==1 && treeCali->thTCCutNumCluster[7]==1 && treeCali->thTCCutNumCluster[9]==1 && treeCali->thTCCutNumCluster[11]==1 && (gapNumTOF2==7 || gapNumTOF2==8 || gapNumTOF2==9 || gapNumTOF2==10 || gapNumTOF2==11 || gapNumTOF2==12) && (gapNumTOF1==7 || gapNumTOF1==8 || gapNumTOF1==9 || gapNumTOF1==10 || gapNumTOF1==11 || gapNumTOF1==12);
  condMatch=treeCali->thTCCutNumCluster[4]==1 && treeCali->thTCCutNumCluster[6]==1 && treeCali->thTCCutNumCluster[8]==1 && treeCali->thTCCutNumCluster[10]==1 && ((gapNumTOF2==7 || gapNumTOF2==8 || gapNumTOF2==9 || gapNumTOF2==10 || gapNumTOF2==11 || gapNumTOF2==12) || (gapNumTOF1==7 || gapNumTOF1==8 || gapNumTOF1==9 || gapNumTOF1==10 || gapNumTOF1==11 || gapNumTOF1==12)) || treeCali->thTCCutNumCluster[5]==1 && treeCali->thTCCutNumCluster[7]==1 && treeCali->thTCCutNumCluster[9]==1 && treeCali->thTCCutNumCluster[11]==1 && ((gapNumTOF2==1 || gapNumTOF2==2 || gapNumTOF2==3 || gapNumTOF2==4 || gapNumTOF2==5 || gapNumTOF2==6) || (gapNumTOF1==1 || gapNumTOF1==2 || gapNumTOF1==3 || gapNumTOF1==4 || gapNumTOF1==5 || gapNumTOF1==6));

  condInequal=(gapNumTOF1!=gapNumTOF2);
  condEqual=(gapNumTOF1==gapNumTOF2);

  if(condTwoSide && condBad) numAcptEvent[0]++;

  if(condTwoSide && condGood && condInequal) numAcptEvent[1]++;

  if(condTwoSide && condGood && condEqual) numAcptEvent[2]++;

  if(condOneSide && condBad) numAcptEvent[3]++;
  
  if(condOneSide && condGood && condUnmatch) numAcptEvent[4]++;

  if(condOneSide && condGood && condMatch && condInequal) numAcptEvent[5]++;

  if(condOneSide && condGood && condMatch && condEqual) numAcptEvent[6]++;


  if(condOneSide && condGood && condMatch && condEqual && treeCali->thTCCutNumCluster[5]==1 && treeCali->thTCCutNumCluster[7]==1 && treeCali->thTCCutNumCluster[9]==1 && treeCali->thTCCutNumCluster[11]==1 && (gapNumTOF2==1 || gapNumTOF2==2 || gapNumTOF2==3 || gapNumTOF2==4 || gapNumTOF2==5 || gapNumTOF2==6)){
    for(int i=1;i<=6;i++){
      if(gapNumTOF2==i){
        histCompPosC3YC4Y[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[7][0],treeCali->thTCCutClusterPosCentroid[11][0]);
        histCompPosC3XC4X[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[5][0],treeCali->thTCCutClusterPosCentroid[9][0]);
        histCompPosC3XC4XAfterCableSwapping[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[5][0],treeCali->thTCCutClusterPosCentroid[9][0]);
      }
    }
  }
  if(condOneSide && condGood && condMatch && condEqual && treeCali->thTCCutNumCluster[4]==1 && treeCali->thTCCutNumCluster[6]==1 && treeCali->thTCCutNumCluster[8]==1 && treeCali->thTCCutNumCluster[10]==1 && (gapNumTOF2==7 || gapNumTOF2==8 || gapNumTOF2==9 || gapNumTOF2==10 || gapNumTOF2==11 || gapNumTOF2==12)){
    for(int i=7;i<=12;i++){
      if(gapNumTOF2==i){
        histCompPosC3YC4Y[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[6][0],treeCali->thTCCutClusterPosCentroid[10][0]);
        histCompPosC3XC4X[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[4][0],treeCali->thTCCutClusterPosCentroid[8][0]);
        if(gapNumTOF2!=7 && gapNumTOF2!=11) histCompPosC3XC4XAfterCableSwapping[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[4][0],treeCali->thTCCutClusterPosCentroid[8][0]);
        else if(gapNumTOF2==7) histCompPosC3XC4XAfterCableSwapping[6]->Fill(treeCali->thTCCutClusterPosCentroidC3XLGap7Swap[0],treeCali->thTCCutClusterPosCentroidC4XLGap7Swap[0]);
        else histCompPosC3XC4XAfterCableSwapping[10]->Fill(treeCali->thTCCutClusterPosCentroid[4][0],treeCali->thTCCutClusterPosCentroidC4XLGap11Swap[0]);
      }
    }
  }   

  if(condOneSide && condGood && condMatch && condInequal && treeCali->thTCCutNumCluster[5]==1 && treeCali->thTCCutNumCluster[7]==1 && treeCali->thTCCutNumCluster[9]==1 && treeCali->thTCCutNumCluster[11]==1 && (gapNumTOF2==1 || gapNumTOF2==2 || gapNumTOF2==3 || gapNumTOF2==4 || gapNumTOF2==5 || gapNumTOF2==6)){
    for(int i=1;i<=6;i++){
      if(gapNumTOF2==i){
        histCompPosC3YC4YWithTOFInequalEventsBasedOnTOF2[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[7][0],treeCali->thTCCutClusterPosCentroid[11][0]);
        histCompPosC3XC4XWithTOFInequalEventsBasedOnTOF2[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[5][0],treeCali->thTCCutClusterPosCentroid[9][0]);
      }
    }
  }
  if(condOneSide && condGood && condMatch && condInequal && treeCali->thTCCutNumCluster[4]==1 && treeCali->thTCCutNumCluster[6]==1 && treeCali->thTCCutNumCluster[8]==1 && treeCali->thTCCutNumCluster[10]==1 && (gapNumTOF2==7 || gapNumTOF2==8 || gapNumTOF2==9 || gapNumTOF2==10 || gapNumTOF2==11 || gapNumTOF2==12)){
    for(int i=7;i<=12;i++){
      if(gapNumTOF2==i){
        histCompPosC3YC4YWithTOFInequalEventsBasedOnTOF2[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[6][0],treeCali->thTCCutClusterPosCentroid[10][0]);
        histCompPosC3XC4XWithTOFInequalEventsBasedOnTOF2[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[4][0],treeCali->thTCCutClusterPosCentroid[8][0]);
      }
    }
  }

  if(condOneSide && condGood && condMatch && condInequal && treeCali->thTCCutNumCluster[5]==1 && treeCali->thTCCutNumCluster[7]==1 && treeCali->thTCCutNumCluster[9]==1 && treeCali->thTCCutNumCluster[11]==1 && (gapNumTOF1==1 || gapNumTOF1==2 || gapNumTOF1==3 || gapNumTOF1==4 || gapNumTOF1==5 || gapNumTOF1==6)){
    for(int i=1;i<=6;i++){
      if(gapNumTOF1==i){
        histCompPosC3YC4YWithTOFInequalEventsBasedOnTOF1[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[7][0],treeCali->thTCCutClusterPosCentroid[11][0]);
        histCompPosC3XC4XWithTOFInequalEventsBasedOnTOF1[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[5][0],treeCali->thTCCutClusterPosCentroid[9][0]);
      }
    }
  }
  if(condOneSide && condGood && condMatch && condInequal && treeCali->thTCCutNumCluster[4]==1 && treeCali->thTCCutNumCluster[6]==1 && treeCali->thTCCutNumCluster[8]==1 && treeCali->thTCCutNumCluster[10]==1 && (gapNumTOF1==7 || gapNumTOF1==8 || gapNumTOF1==9 || gapNumTOF1==10 || gapNumTOF1==11 || gapNumTOF1==12)){
    for(int i=7;i<=12;i++){
      if(gapNumTOF1==i){
        histCompPosC3YC4YWithTOFInequalEventsBasedOnTOF1[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[6][0],treeCali->thTCCutClusterPosCentroid[10][0]);
        histCompPosC3XC4XWithTOFInequalEventsBasedOnTOF1[i-1]->Fill(treeCali->thTCCutClusterPosCentroid[4][0],treeCali->thTCCutClusterPosCentroid[8][0]);
      }
    }
  }
  return 0;
}

Long_t Det_MWPC::done_checkCableSwapping(){
  TH1D *histCompTwoFiles=new TH1D("Comparison between two files","Comparison between two files;;Counts",7,0,7);
  for(int i=0;i<7;i++){
    histCompTwoFiles->SetBinContent(i+1,numAcptEvent[i]);
  }
  TCanvas *canvCompTwoFiles=new TCanvas("Comparison between two files","Comparison between two files",1400,600);
  histCompTwoFiles->Draw();
  histCompTwoFiles->GetXaxis()->SetLabelOffset(99);
  double x, y;
  string type[7]={"F1-T;F2-B","F1-T;F2-GI","F1-T;F2-GE","F1-O;F2-B","F1-O;F2-GU","F1-O;F2-GMI","F1-O;F2-GME"};
  y = gPad->GetUymin() + 15000*histCompTwoFiles->GetYaxis()->GetBinWidth(1);
  TText t;
  t.SetTextAngle(30);
  t.SetTextSize(0.06);
  t.SetTextAlign(33);
  t.SetTextColor(kRed);
  for (int j=0;j<7;j++) {
    x = histCompTwoFiles->GetXaxis()->GetBinCenter(j+1)+0.5;
    t.DrawText(x,y,type[j].c_str());
  }
  canvCompTwoFiles->SaveAs("Comparison between two files.pdf");
 
  TCanvas *canvCompPosC3YC4Y=new TCanvas("Comparison of position between C3Y and C4Y","Comparison of position between C3Y and C4Y",1200,900);
  canvCompPosC3YC4Y->Divide(4,3);
  for(int i=0;i<12;i++){
    canvCompPosC3YC4Y->cd(i+1);
    gPad->SetLogz();
    histCompPosC3YC4Y[i]->Draw("colz");
  }
  canvCompPosC3YC4Y->SaveAs("Comparison of position between C3Y and C4Y.pdf");

  TCanvas *canvCompPosC3XC4X=new TCanvas("Comparison of position between C3X and C4X","Comparison of position between C3X and C4X",1200,900/8.*9);
  canvCompPosC3XC4X->Divide(4,3);
  for(int i=0;i<12;i++){
    canvCompPosC3XC4X->cd(i+1);
    gPad->SetLogz();
    histCompPosC3XC4X[i]->Draw("colz");
  }
  canvCompPosC3XC4X->SaveAs("Comparison of position between C3X and C4X.pdf");

  TCanvas *canvCompPosC3XC4XAfterCableSwapping=new TCanvas("Comparison of position between C3X and C4X after cable swapping","Comparison of position between C3X and C4X after cable swapping",1200,900/8.*9);
  canvCompPosC3XC4XAfterCableSwapping->Divide(4,3);
  for(int i=0;i<12;i++){
    canvCompPosC3XC4XAfterCableSwapping->cd(i+1);
    gPad->SetLogz();
    histCompPosC3XC4XAfterCableSwapping[i]->Draw("colz");
  }
  canvCompPosC3XC4XAfterCableSwapping->SaveAs("Comparison of position between C3X and C4X after cable swapping.pdf");
   
  TCanvas *canvCompPosC3YC4YWithTOFInequalEventsBasedOnTOF2=new TCanvas("Comparison of position between C3Y and C4Y with TOF-inequal events based on TOF2","Comparison of position between C3Y and C4Y with TOF-inequal events based on TOF2",1200,900);
  canvCompPosC3YC4YWithTOFInequalEventsBasedOnTOF2->Divide(4,3);
  for(int i=0;i<12;i++){
    canvCompPosC3YC4YWithTOFInequalEventsBasedOnTOF2->cd(i+1);
    gPad->SetLogz();
    histCompPosC3YC4YWithTOFInequalEventsBasedOnTOF2[i]->Draw("colz");
  }
  canvCompPosC3YC4YWithTOFInequalEventsBasedOnTOF2->SaveAs("Comparison of position between C3Y and C4Y with TOF-inequal events based on TOF2.pdf");

  TCanvas *canvCompPosC3XC4XWithTOFInequalEventsBasedOnTOF2=new TCanvas("Comparison of position between C3X and C4X with TOF-inequal events based on TOF2","Comparison of position between C3X and C4X with TOF-inequal events based on TOF2",1200,900/8.*9);
  canvCompPosC3XC4XWithTOFInequalEventsBasedOnTOF2->Divide(4,3);
  for(int i=0;i<12;i++){
    canvCompPosC3XC4XWithTOFInequalEventsBasedOnTOF2->cd(i+1);
    gPad->SetLogz();
    histCompPosC3XC4XWithTOFInequalEventsBasedOnTOF2[i]->Draw("colz");
  }
  canvCompPosC3XC4XWithTOFInequalEventsBasedOnTOF2->SaveAs("Comparison of position between C3X and C4X with TOF-inequal events based on TOF2.pdf");

  TCanvas *canvCompPosC3YC4YWithTOFInequalEventsBasedOnTOF1=new TCanvas("Comparison of position between C3Y and C4Y with TOF-inequal events based on TOF1","Comparison of position between C3Y and C4Y with TOF-inequal events based on TOF1",1200,900);
  canvCompPosC3YC4YWithTOFInequalEventsBasedOnTOF1->Divide(4,3);
  for(int i=0;i<12;i++){
    canvCompPosC3YC4YWithTOFInequalEventsBasedOnTOF1->cd(i+1);
    gPad->SetLogz();
    histCompPosC3YC4YWithTOFInequalEventsBasedOnTOF1[i]->Draw("colz");
  }
  canvCompPosC3YC4YWithTOFInequalEventsBasedOnTOF1->SaveAs("Comparison of position between C3Y and C4Y with TOF-inequal events based on TOF1.pdf");

  TCanvas *canvCompPosC3XC4XWithTOFInequalEventsBasedOnTOF1=new TCanvas("Comparison of position between C3X and C4X with TOF-inequal events based on TOF1","Comparison of position between C3X and C4X with TOF-inequal events based on TOF1",1200,900/8.*9);
  canvCompPosC3XC4XWithTOFInequalEventsBasedOnTOF1->Divide(4,3);
  for(int i=0;i<12;i++){
    canvCompPosC3XC4XWithTOFInequalEventsBasedOnTOF1->cd(i+1);
    gPad->SetLogz();
    histCompPosC3XC4XWithTOFInequalEventsBasedOnTOF1[i]->Draw("colz");
  }
  canvCompPosC3XC4XWithTOFInequalEventsBasedOnTOF1->SaveAs("Comparison of position between C3X and C4X with TOF-inequal events based on TOF1.pdf");

  return 0;
};

