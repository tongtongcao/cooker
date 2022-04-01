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
#include "MwpcPedestalValue.h"
using namespace std;

static TH1D *histNHit[12];
static TH2D *histMult[12];
static int totNEvent=0;
static int hit0[12]={0},sing[12]={0},doub[12]={0},trip[12]={0},singdoub[12]={0},singtrip[12]={0},doubtrip[12]={0},singdoubtrip[12]={0};
TH2D *histAdcVsStrip[12];

Long_t Det_MWPC::histos_event(){
  ostringstream sstr_name, sstr_title;
  for(int i=0;i<12;i++){
    sstr_name<<"ADC Vs Strip for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Strip; ADC";
    histAdcVsStrip[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,150,0,5000);
    sstr_name.str("");
    sstr_title.str("");

    sstr_name<<"Distribution of number of hits for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Number of hits; Counts";
    histNHit[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),11,-0.5,10.5);
    sstr_name.str("");
    sstr_title.str("");
    
    sstr_name<<"Distribution of multiplets for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Singlet/Doublet/>=Triplet; Number of S/D/>=T";
    histMult[i]=new TH2D(sstr_name.str().c_str(),sstr_title.str().c_str(),3,0.5,3.5,10,0.5,10.5);
    sstr_name.str("");
    sstr_title.str("");

  }

      
  return 0;
}

Long_t Det_MWPC::startup_event(){
  getBranchObject("RawMwpcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  return 0;
};


Long_t Det_MWPC::process_event(){
  int* data[12]={treeRaw->ADC_C2X_L,treeRaw->ADC_C2X_R,treeRaw->ADC_C2Y_L,treeRaw->ADC_C2Y_R,treeRaw->ADC_C3X_L,treeRaw->ADC_C3X_R,treeRaw->ADC_C3Y_L,treeRaw->ADC_C3Y_R,treeRaw->ADC_C4X_L,treeRaw->ADC_C4X_R,treeRaw->ADC_C4Y_L,treeRaw->ADC_C4Y_R};
  const double* pedValue[12]={pedValueC2XL,pedValueC2XR,pedValueC2YL,pedValueC2YR,pedValueC3XL,pedValueC3XR,pedValueC3YL,pedValueC3YR,pedValueC4XL,pedValueC4XR,pedValueC4YL,pedValueC4YR};
  const double* pedSigma[12]={pedSigmaC2XL,pedSigmaC2XR,pedSigmaC2YL,pedSigmaC2YR,pedSigmaC3XL,pedSigmaC3XR,pedSigmaC3YL,pedSigmaC3YR,pedSigmaC4XL,pedSigmaC4XR,pedSigmaC4YL,pedSigmaC4YR}; 

  MwpcEvent *event[12], *event_pedSub[12],*event_pedSub_sigmaCut[12];
  ostringstream sstr;
  if(treeRaw->event<200){
    for(int i=0;i<12;i++){
      sstr<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);  
      event[i]=new MwpcEvent(treeRaw->event,sstr.str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,!((i/4==0)&&((i%4/2))==0));
      sstr.str("");
      event[i]->setHist();
      event[i]->addHistData(data[i]);
      sstr<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<"_pedSub";           
      event_pedSub[i]=new MwpcEvent(treeRaw->event,sstr.str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,!((i/4==0)&&((i%4/2))==0));
      sstr.str("");
      event_pedSub[i]->setHist();
      event_pedSub[i]->addHistData_pedSub(data[i],pedValue[i]);
      sstr<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76)<<"_pedSub_sigmaCut";
      event_pedSub_sigmaCut[i]=new MwpcEvent(treeRaw->event,sstr.str(),((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,-0.5,((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8-0.5,!((i/4==0)&&((i%4/2))==0));
      sstr.str("");
      event_pedSub_sigmaCut[i]->setHist();
      event_pedSub_sigmaCut[i]->addHistData_pedSub_sigmaCut(data[i],pedValue[i],pedSigma[i],3);
    }

    sstr<<"MwpcEvent"<<treeRaw->event;
    TCanvas *canv=new TCanvas(sstr.str().c_str(),sstr.str().c_str(),2000,1200);
    sstr.str("");
    canv->Divide(4,3);
    for(int i=0;i<12;i++){
        canv->cd(i+1);
        event[i]->getHist()->SetMinimum(0);
        event[i]->getHist()->SetLineColor(kRed);
        event[i]->getHist()->SetLineWidth(2);
        event[i]->getHist()->Draw();
        event_pedSub[i]->getHist()->SetLineColor(kBlue);
        event_pedSub[i]->getHist()->SetLineWidth(2);
        event_pedSub[i]->getHist()->Draw("same");
        event_pedSub_sigmaCut[i]->getHist()->SetLineColor(kGreen);
        event_pedSub_sigmaCut[i]->getHist()->SetLineWidth(2);
        //event_pedSub_sigmaCut[i]->getHist()->Draw("same");
    } 
    TLegend *legend=new TLegend(0.7,0.8,1,1);
    legend->AddEntry(event[0]->getHist(),"Original Data","l");
    legend->AddEntry(event_pedSub[0]->getHist(),"After Ped Sub","l");
    //legend->AddEntry(event_pedSub_sigmaCut[0]->getHist(),"+ 3#sigma Cut","l");
    canv->cd(1);
    legend->Draw();
    sstr<<"MwpcEvent"<<treeRaw->event<<".png";
    canv->SaveAs(sstr.str().c_str());
    canv->Close();
  }  

  int nHit[12],*mult[12];
  totNEvent++;
  double *adcValue[12];
  for(int i=0;i<12;i++){
    event[i]=new MwpcEvent(((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8,!((i/4==0)&&((i%4/2))==0));
    adcValue[i]=event[i]->getAdcValue(data[i],pedValue[i],pedSigma[i],3);
    for(int j=0;j<((i%4)/2*(-5)+(i/4)*(!((i%4)/2))+7)*8;j++){
      histAdcVsStrip[i]->Fill(j,adcValue[i][j]);
    }
    nHit[i]=event[i]->getNHit(data[i],pedValue[i],pedSigma[i],3);
    mult[i]=event[i]->multiplet3Type(data[i],pedValue[i],pedSigma[i],3);
    histNHit[i]->Fill(nHit[i]);
    for(int j=0;j<3;j++){
      histMult[i]->Fill(j+1,mult[i][j]);
    }
  }

  for(int i=0;i<12;i++){
    if(nHit[i]==0){
      hit0[i]++;
      }
    else{
      if(mult[i][0]!=0 && mult[i][1]==0 && mult[i][2]==0) sing[i]++;
      else if(mult[i][0]==0 && mult[i][1]!=0 && mult[i][2]==0) doub[i]++;
      else if(mult[i][0]==0 && mult[i][1]==0 && mult[i][2]!=0) trip[i]++; 
      else if(mult[i][0]!=0 && mult[i][1]!=0 && mult[i][2]==0) singdoub[i]++;
      else if(mult[i][0]!=0 && mult[i][1]==0 && mult[i][2]!=0) singtrip[i]++;
      else if(mult[i][0]==0 && mult[i][1]!=0 && mult[i][2]!=0) doubtrip[i]++;
      else singdoubtrip[i]++;
    }
  }

  return 0;
}


Long_t Det_MWPC::done_event(){
  TCanvas *canvNHit=new TCanvas("Distribution of number of hits","Distribution of number of hits",2000,900);
  canvNHit->Divide(4,3);
  for(int i=0;i<12;i++){
    canvNHit->cd(i+1);
    histNHit[i]->Draw();
  }
  canvNHit->SaveAs("Distribution of number of hits.pdf");
    
  TCanvas *canvMult=new TCanvas("Distribution of multiplets","Distribution of multiplets",2000,900);
  canvMult->Divide(4,3);
  for(int i=0;i<12;i++){
    canvMult->cd(i+1);
    gPad->SetLogz();
    histMult[i]->GetXaxis()->SetLabelOffset(99);
    histMult[i]->GetXaxis()->SetTitleOffset(0.5);
    histMult[i]->GetXaxis()->SetTitleSize(0.06);
    histMult[i]->GetYaxis()->SetLabelSize(0.06);
    histMult[i]->GetYaxis()->SetTitleOffset(0.8);
    histMult[i]->GetYaxis()->SetTitleSize(0.06);
    histMult[i]->GetZaxis()->SetLabelSize(0.045);
    histMult[i]->Draw("colz");
    double x, y;
    string type[3]={"Singlet","Doublet",">=Triplet"};
    y = gPad->GetUymin() + 4*histMult[i]->GetYaxis()->GetBinWidth(1);
    TText t;
    t.SetTextAngle(60);
    t.SetTextSize(0.06);
    t.SetTextAlign(33);
    for (int j=0;j<3;j++) {
      x = histMult[j]->GetXaxis()->GetBinCenter(j+1);
      t.DrawText(x,y,type[j].c_str());
    }
  }
  canvMult->SaveAs("Distribution of multiplets.pdf");

  TH1D *histMultStat[12];
  ostringstream sstr_name, sstr_title;
  for(int i=0;i<12;i++){
    sstr_name<<"Statistic of multiplets for "<<"MwpcC"<<i/4+2<<(char)((i%4)/2+88)<<(char)((i%4%2)*6+76);
    sstr_title<<sstr_name.str()<<";Type of multiplets; occupation (%)";
    histMultStat[i]=new TH1D(sstr_name.str().c_str(),sstr_title.str().c_str(),8,-0.5,7.5);
    sstr_name.str("");
    sstr_title.str("");
    histMultStat[i]->SetBinContent(1,hit0[i]/(double)totNEvent*100);
    histMultStat[i]->SetBinContent(2,sing[i]/(double)totNEvent*100);
    histMultStat[i]->SetBinContent(3,doub[i]/(double)totNEvent*100);
    histMultStat[i]->SetBinContent(4,trip[i]/(double)totNEvent*100);
    histMultStat[i]->SetBinContent(5,singdoub[i]/(double)totNEvent*100);
    histMultStat[i]->SetBinContent(6,singtrip[i]/(double)totNEvent*100);
    histMultStat[i]->SetBinContent(7,doubtrip[i]/(double)totNEvent*100);
    histMultStat[i]->SetBinContent(8,singdoubtrip[i]/(double)totNEvent*100);    
  }

  TCanvas *canvMultStat=new TCanvas("Statistic of multiplets","Statistic of multiplets",2000,900);
  canvMultStat->Divide(4,3);
  for(int i=0;i<12;i++){
    canvMultStat->cd(i+1);
    histMultStat[i]->SetMinimum(0);
    histMultStat[i]->SetMaximum(50); 
    histMultStat[i]->GetXaxis()->SetLabelOffset(99);
    histMultStat[i]->GetXaxis()->SetTitleOffset(0.5);
    histMultStat[i]->GetXaxis()->SetTitleSize(0.06);
    histMultStat[i]->GetYaxis()->SetLabelSize(0.06);
    histMultStat[i]->GetYaxis()->SetTitleOffset(0.8);
    histMultStat[i]->GetYaxis()->SetTitleSize(0.06);
    histMultStat[i]->Draw();
    double x, y;
    string type[8]={"NoH","S","D","T","SD","ST","DT","SDT"};
    y = gPad->GetUymin() + 6*histMultStat[i]->GetYaxis()->GetBinWidth(1);
    TText t;
    t.SetTextAngle(60);
    t.SetTextSize(0.06);
    t.SetTextAlign(33);
    t.SetTextColor(kRed);
    for (int j=0;j<8;j++) {
      x = histMultStat[j]->GetXaxis()->GetBinCenter(j+1);
      t.DrawText(x,y,type[j].c_str());
    }
  }
  canvMultStat->SaveAs("Statistic of multiplets.pdf");

  TCanvas *canvAdcVsStrip=new TCanvas("Adc Vs Strip","Adc Vs Strip",2000,1500);
  canvAdcVsStrip->Divide(4,3);
  for(int i=0;i<12;i++){
    canvAdcVsStrip->cd(i+1);
    gPad->SetLogz();
    histAdcVsStrip[i]->Draw("colz");
  }
  canvAdcVsStrip->SaveAs("Adc Vs Strip.pdf");

  return 0;
};
