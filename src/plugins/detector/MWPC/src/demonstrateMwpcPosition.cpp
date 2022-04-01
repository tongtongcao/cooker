////// To use cut functions of diff vs. sum of XY-total-charge, a full directory should be given since the functions are saved in a root file. Here, the full directory is /home/caot/trek-e36/cooker/src/plugins/detector/MWPC/src/MwpcDiffVsSumTotalChargeCutFunction.root 
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
static MwpcPosition *positionX2Y2C2L[6];
static int ievent=0;
Long_t Det_MWPC::histos_demonstratePosition(){
  return 0;
}

Long_t Det_MWPC::startup_demonstratePosition(){
  getBranchObject("CaliMwpcInfo",(TObject **) &treeCali);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  gStyle->SetNdivisions(72,"X");
  gStyle->SetNdivisions(16,"Y");  
  gStyle->SetTickLength(0,"X");
  gStyle->SetTickLength(0,"Y");

  return 0;
};


Long_t Det_MWPC::process_demonstratePosition(){
  MwpcPosition *position[6];
  MwpcPosition *positionCentroid[6];
  vector<vector<int> > stNumX[6],stNumY[6];
  vector<vector<double> > stChargeX[6],stChargeY[6];
  vector<double> pairTotalChargeX[6],pairTotalChargeY[6];
  vector<double> posCentroidX[6],posCentroidY[6], posErrCentroidX[6],posErrCentroidY[6];
  for(int i=0;i<6;i++){
    stNumX[i].clear();
    stChargeX[i].clear();
    pairTotalChargeX[i].clear();
    stNumY[i].clear();
    stChargeY[i].clear();
    pairTotalChargeY[i].clear();
  }
  ostringstream sstr;
  vector<TMarker*> markPosCentroid[6];
  vector<TLine*> markPosErrCentroidLeft[6],markPosErrCentroidRight[6],markPosErrCentroidTop[6],markPosErrCentroidBottom[6];

  if(treeCali->event==140 || treeCali->event==463 || treeCali->event==596 || treeCali->event==813 || treeCali->event==1995 || treeCali->event==2473){
    sstr<<"MwpcC2LPosition";
    positionX2Y2C2L[ievent]=new MwpcPosition(treeCali->event,sstr.str().c_str(),56);
    sstr.str("");
    positionX2Y2C2L[ievent]->defineHist();
    positionX2Y2C2L[ievent]->setBinContent(treeCali->thTCPairCutPairClusterStNum[0],treeCali->thTCPairCutPairClusterStCharge[0],treeCali->thTCPairCutPairClusterTC[0],treeCali->thTCPairCutPairClusterStNum[2],treeCali->thTCPairCutPairClusterStCharge[2],treeCali->thTCPairCutPairClusterTC[2]); 
    ievent++;
  }
  if(treeCali->event<10){
    for(int i=0;i<6;i++){
      sstr<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76)<<"Position";
      position[i]=new MwpcPosition(treeCali->event,sstr.str().c_str(),(i/2+7)*8);
      sstr.str("");
      stNumX[i]=treeCali->thTCPairCutPairClusterStNum[i/2*4+i%2];
      stChargeX[i]=treeCali->thTCPairCutPairClusterStCharge[i/2*4+i%2];
      pairTotalChargeX[i]=treeCali->thTCPairCutPairClusterTC[i/2*4+i%2];
      stNumY[i]=treeCali->thTCPairCutPairClusterStNum[i/2*4+i%2+2];
      stChargeY[i]=treeCali->thTCPairCutPairClusterStCharge[i/2*4+i%2+2];
      pairTotalChargeY[i]=treeCali->thTCPairCutPairClusterTC[i/2*4+i%2+2];
      position[i]->defineHist();
      position[i]->setBinContent(stNumX[i],stChargeX[i],pairTotalChargeX[i],stNumY[i],stChargeY[i],pairTotalChargeY[i]);
   
      sstr<<"MwpcC"<<i/2+2<<(char)((i%2)*6+76)<<"PositionCentroid";
      positionCentroid[i]=new MwpcPosition(treeCali->event,sstr.str().c_str(),(i/2+7)*8);
      sstr.str("");
      posCentroidX[i]=treeCali->thTCPairCutPairPosCentroid[i/2*4+i%2];
      posErrCentroidX[i]=treeCali->thTCPairCutPairPosErrCentroid[i/2*4+i%2];   
      posCentroidY[i]=treeCali->thTCPairCutPairPosCentroid[i/2*4+i%2+2];
      posErrCentroidY[i]=treeCali->thTCPairCutPairPosErrCentroid[i/2*4+i%2+2]; 
      positionCentroid[i]->defineHist();
      markPosCentroid[i]= positionCentroid[i]->markPosition(posCentroidX[i],posErrCentroidX[i],posCentroidY[i],posErrCentroidY[i]);
      markPosErrCentroidLeft[i]=positionCentroid[i]->markPositionErrorLeft(posCentroidX[i],posErrCentroidX[i],posCentroidY[i],posErrCentroidY[i]);
      markPosErrCentroidRight[i]=positionCentroid[i]->markPositionErrorRight(posCentroidX[i],posErrCentroidX[i],posCentroidY[i],posErrCentroidY[i]);
      markPosErrCentroidTop[i]=positionCentroid[i]->markPositionErrorTop(posCentroidX[i],posErrCentroidX[i],posCentroidY[i],posErrCentroidY[i]);
      markPosErrCentroidBottom[i]=positionCentroid[i]->markPositionErrorBottom(posCentroidX[i],posErrCentroidX[i],posCentroidY[i],posErrCentroidY[i]);
    }
    sstr<<"MwpcPosition"<<treeCali->event;
    TCanvas *canv=new TCanvas(sstr.str().c_str(),sstr.str().c_str(),2500./3*2/16*72,2500);
    sstr.str("");
    TPad *pad[6];
    pad[0]=new TPad("pad0","pad0",0,2./3,0.5*56./72,1);
    pad[1]=new TPad("pad1","pad1",0.5,2./3,0.5+0.5*56./72,1);
    pad[2]=new TPad("pad2","pad2",0,1./3,0.5*64/72,2./3); 
    pad[3]=new TPad("pad3","pad3",0.5,1./3,0.5+0.5*64/72,2./3); 
    pad[4]=new TPad("pad4","pad4",0,0,0.5,1./3); 
    pad[5]=new TPad("pad5","pad5",0.5,0,1,1./3);            
    
    for(int i=0;i<6;i++){
      if(i==0 || i==1){
        pad[i]->SetLeftMargin(0.02*72/56);
        pad[i]->SetRightMargin(0.1-0.02*72/56);
        pad[i]->SetTopMargin(0.05);
        pad[i]->SetBottomMargin(0.05);
        pad[i]->Draw();
      }
      if(i==2 || i==3){
        pad[i]->SetLeftMargin(0.02*72/64);
        pad[i]->SetRightMargin(0.1-0.02*72/64);
        pad[i]->SetTopMargin(0.05);
        pad[i]->SetBottomMargin(0.05);
        pad[i]->Draw();
      }
      if(i==4 || i==5){
        pad[i]->SetLeftMargin(0.02);
        pad[i]->SetRightMargin(0.08);
        pad[i]->SetTopMargin(0.05);
        pad[i]->SetBottomMargin(0.05);
        pad[i]->Draw();
      }
    }

    TLine *line;
    for(int i=0;i<6;i++){
      pad[i]->cd();
      position[i]->getHist()->Draw("colz");
      if(i==0 || i==1){
        for(int m=1;m<56;m++){
          line=new TLine(m,0,m,16);
          line->SetLineWidth(0.01);
          line->Draw();
        }  
        for(int m=1;m<16;m++){
          line=new TLine(0,m,56,m);
          line->SetLineWidth(0.01);
          line->Draw();
        }
      }

      if(i==2 || i==3){
        for(int m=1;m<64;m++){
          line=new TLine(m,0,m,16);
          line->SetLineWidth(0.01);
          line->Draw();
        }
        for(int m=1;m<16;m++){
          line=new TLine(0,m,64,m);
          line->SetLineWidth(0.01);
          line->Draw();
        }
      }

      if(i==4 || i==5){
        for(int m=1;m<72;m++){
          line=new TLine(m,0,m,16);
          line->SetLineWidth(0.01);
          line->Draw();
        }
        for(int m=1;m<16;m++){
          line=new TLine(0,m,72,m);
          line->SetLineWidth(0.01);
          line->Draw();
        }
      }
    }
     
   sstr<<"MwpcPosition"<<treeCali->event<<".pdf";
   canv->SaveAs(sstr.str().c_str());
   sstr.str("");

   sstr<<"MwpcPositionCentroid"<<treeCali->event;
   TCanvas *canvCentroid=new TCanvas(sstr.str().c_str(),sstr.str().c_str(),2500./3*2/16*72,2500);
   sstr.str("");

    for(int i=0;i<6;i++){
      pad[i]->cd();
      positionCentroid[i]->getHist()->Draw();
      if(i==0 || i==1){
        for(int m=1;m<56;m++){
          line=new TLine(m,0,m,16);
          line->SetLineWidth(0.01);
          line->Draw();
        }
        for(int m=1;m<16;m++){
          line=new TLine(0,m,56,m);
          line->SetLineWidth(0.01);
          line->Draw();
        }
      }

      if(i==2 || i==3){
        for(int m=1;m<64;m++){
          line=new TLine(m,0,m,16);
          line->SetLineWidth(0.01);
          line->Draw();
        }
        for(int m=1;m<16;m++){
          line=new TLine(0,m,64,m);
          line->SetLineWidth(0.01);
          line->Draw();
        }
      }

      if(i==4 || i==5){
        for(int m=1;m<72;m++){
          line=new TLine(m,0,m,16);
          line->SetLineWidth(0.01);
          line->Draw();
        }
        for(int m=1;m<16;m++){
          line=new TLine(0,m,72,m);
          line->SetLineWidth(0.01);
          line->Draw();
        }
      }
      for(UInt_t j=0;j<markPosCentroid[i].size();j++){
        markPosCentroid[i][j]->SetMarkerColor(kRed);
        markPosCentroid[i][j]->SetMarkerSize(3);
        markPosCentroid[i][j]->Draw();
        markPosErrCentroidLeft[i][j]->SetLineColor(kRed);
        markPosErrCentroidLeft[i][j]->Draw();
        markPosErrCentroidRight[i][j]->SetLineColor(kRed);
        markPosErrCentroidRight[i][j]->Draw();
        markPosErrCentroidTop[i][j]->SetLineColor(kRed);
        markPosErrCentroidTop[i][j]->Draw();
        markPosErrCentroidBottom[i][j]->SetLineColor(kRed);
        markPosErrCentroidBottom[i][j]->Draw();
      }
    }
   sstr<<"MwpcPositionCentroid"<<treeCali->event<<".pdf";
   canv->SaveAs(sstr.str().c_str());
   sstr.str("");
  }

  return 0;
}

Long_t Det_MWPC::done_demonstratePosition(){
  TLine *line;
  TCanvas *canvX2Y2=new TCanvas("demonstration of X2Y2","demonstration of X2Y2",1800./3*2/16*56,1800);
  canvX2Y2->Divide(2,3);
  for(int i=0;i<6;i++){
    canvX2Y2->cd(i+1);
    positionX2Y2C2L[i]->getHist()->Draw("colz");
    for(int m=1;m<56;m++){
      line=new TLine(m,0,m,16);
      line->SetLineWidth(0.01);
      line->Draw();
    }
    for(int m=1;m<16;m++){
      line=new TLine(0,m,56,m);
      line->SetLineWidth(0.01);
      line->Draw();
    }
  }
  canvX2Y2->SaveAs("demonstration of X2Y2.pdf");

  return 0;
};

