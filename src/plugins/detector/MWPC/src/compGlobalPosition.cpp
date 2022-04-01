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
#include "mwpcLocalGlobal.h"
using namespace std;

static string inputFilePath(getenv("mwpc_inputFiles"));
static string inputFile=inputFilePath+"/pos_run3994_TRIUMF";
static ifstream input(inputFile.c_str());
static int eventNum, gapNum;
static double mom,c2Y,c2Z,c3X,c3Y,c4X,c4Y;

static TH1D *histc2YDiff;
static TH1D *histc2ZDiff;
static TH1D *histc3XDiff;
static TH1D *histc3YDiff;
static TH1D *histc4XDiff;
static TH1D *histc4YDiff;

static TH1D *histc2YDiffGap[12];
static TH1D *histc2ZDiffGap[12];

Long_t Det_MWPC::histos_compGlobalPosition(){
  ostringstream sstr_name,sstr_title;
  histc2YDiff=new TH1D("Diff. of C2Y","Diff. of C2Y;Diff (mm); Counts",100,-10,10);
  histc2ZDiff=new TH1D("Diff. of C2Z","Diff. of C2Z;Diff (mm); Counts",100,-10,10);
  histc3XDiff=new TH1D("Diff. of C3X","Diff. of C3X;Diff (mm); Counts",100,-10,10);
  histc3YDiff=new TH1D("Diff. of C3Y","Diff. of C3Y;Diff (mm); Counts",100,-10,10);
  histc4XDiff=new TH1D("Diff. of C4X","Diff. of C4X;Diff (mm); Counts",100,-10,10);
  histc4YDiff=new TH1D("Diff. of C4Y","Diff. of C4Y;Diff (mm); Counts",100,-10,10);

  for(int i=0;i<12;i++){
    sstr_name<<"Diff. of C2Y for gap "<<i+1;
    sstr_title<<"Gap "<<i+1<<";Diff (mm); Counts";
    histc2YDiffGap[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,-10,10);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Diff. of C2Z for gap "<<i+1;
    sstr_title<<"Gap "<<i+1<<";Diff (mm); Counts";
    histc2ZDiffGap[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),100,-10,10);
    sstr_name.str("");
    sstr_title.str("");
  }

  return 0;
}

Long_t Det_MWPC::startup_compGlobalPosition(){
  getBranchObject("CookMwpcInfoV1",(TObject **) &treeCookV1);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);


  return 0;
}

Long_t Det_MWPC::process_compGlobalPosition(){
  input>>eventNum>>gapNum>>mom>>c2Y>>c2Z>>c3X>>c3Y>>c4X>>c4Y;
  if((UInt_t)eventNum!=treeCookV1->event){
    cout<<"Error: Event number are not consistent between two files!"<<endl;
    cout<<"For MWPC file: Event number - "<<treeCookV1->event<<endl;
    cout<<"For TRIUMF file: Event number - "<<eventNum<<endl;
  }

  if(c2Y!=-10000 && treeCookV1->posC2Y!=-10000 && treeCookV1->posC3Y!=-10000 && treeCookV1->posC4Y!=-10000){
    histc2YDiff->Fill(treeCookV1->posC2Y-c2Y*10);
    histc2ZDiff->Fill(treeCookV1->posC2Z-c2Z*10);
    histc3XDiff->Fill(treeCookV1->posC3X-c3X*10);
    histc3YDiff->Fill(treeCookV1->posC3Y-c3Y*10);
    histc4XDiff->Fill(treeCookV1->posC4X-c4X*10);
    histc4YDiff->Fill(treeCookV1->posC4Y-c4Y*10);
    
    for(int i=0;i<12;i++){
      if(gapNum==i+1){
       histc2YDiffGap[i]->Fill(treeCookV1->posC2Y-c2Y*10);
       histc2ZDiffGap[i]->Fill(treeCookV1->posC2Z-c2Z*10);
      }
    }
  }

  return 0;
}

Long_t Det_MWPC::done_compGlobalPosition(){
  TCanvas *canv=new TCanvas("Difference of global position","Difference of global position",1000,900);
  canv->Divide(2,3);
  canv->cd(1);
  histc2ZDiff->Draw();
  canv->cd(2);
  histc2YDiff->Draw();
  canv->cd(3);
  histc3XDiff->Draw();
  canv->cd(4);
  histc3YDiff->Draw();
  canv->cd(5);
  histc4XDiff->Draw();
  canv->cd(6);
  histc4YDiff->Draw();
  canv->SaveAs("Difference of global position.pdf");

  TLine *line;

  TCanvas *canvC2Y=new TCanvas("Difference of global position for C2Y","Difference of global position for C2Y",1500,1200);
  canvC2Y->Divide(3,4);
  for(int i=0;i<12;i++){
    canvC2Y->cd(i+1);
    histc2YDiffGap[i]->Draw();
    line=new TLine(c2YCenter[i],0,c2YCenter[i],10000);
    line->SetLineColor(kRed);
    line->Draw();
  }
  canvC2Y->SaveAs("Difference of global position for C2Y.pdf");

  TCanvas *canvC2Z=new TCanvas("Difference of global position for C2Z","Difference of global position for C2Z",1500,1200);
  canvC2Z->Divide(3,4);
  for(int i=0;i<12;i++){
    canvC2Z->cd(i+1);
    histc2ZDiffGap[i]->Draw();
    line=new TLine(c2ZCenter[i],0,c2ZCenter[i],10000);
    line->SetLineColor(kRed);
    line->Draw();
  }
  canvC2Z->SaveAs("Difference of global position for C2Z.pdf");

  return 0;
}

