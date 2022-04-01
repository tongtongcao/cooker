#include <Det_ToF1.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
using namespace std;
#define NSPERTDC 0.025


Long_t Det_ToF1::histos_lar()
{
  for(int iTof=0;iTof<12;iTof++){
    string name("TDiff_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    string title("Half Time Diff between Tof1U and D Tdc;#frac{T_{u}-T_{d}}{2} (ns)");
    h1Diff[iTof]=new TH1D(name.c_str(),title.c_str(),200,-10.,10.);
  }
  string name("GapVSTDiff");
  h2GapVSTDiff=new TH2D("GapVSTDiff","Gap vs. half time diff;#frac{T_{u}-T_{d}}{2} (ns);Gap",200,-10.,10.,12,0.5,12.5);
  h2GapVSTDiffCorr=new TH2D("GapVSTDiffCorr","Gap vs. corrected half time diff;#frac{T_{u}-T_{d}}{2} (ns);Gap",200,-10.,10.,12,0.5,12.5);

  return 0;

}

Long_t Det_ToF1::startup_lar()
{
  getBranchObject("RawTof1Info",(TObject **) &treeRaw);
  gStyle->SetOptStat(0);
  shiftLeftRight[0]= -0.5;
  shiftLeftRight[1]= -1.15;
  shiftLeftRight[2]= 0.45;
  shiftLeftRight[3]= -0.55;
  shiftLeftRight[4]= -2.2;
  shiftLeftRight[5]= -1.55;
  shiftLeftRight[6]= -1.7;
  shiftLeftRight[7]= -2.75;
  shiftLeftRight[8]= -1.75;
  shiftLeftRight[9]= -1.75;
  shiftLeftRight[10]= -2.4;
  shiftLeftRight[11]= -2.4;

  return 0;
};


Long_t Det_ToF1::process_lar()
{
  for(int iTof=0;iTof<12;iTof++){
    if(treeRaw->TDC_TOF1U[iTof]>0)
      if(treeRaw->TDC_TOF1D[iTof]>0){
	double diff=(treeRaw->TDC_TOF1U[iTof]-treeRaw->TDC_TOF1D[iTof])/2.0*NSPERTDC;
	h1Diff[iTof]->Fill(diff);
	h2GapVSTDiff->Fill(diff,iTof+1);
	double corrLeft=treeRaw->TDC_TOF1U[iTof]*NSPERTDC-shiftLeftRight[iTof];
	double corrRight=treeRaw->TDC_TOF1D[iTof]*NSPERTDC+shiftLeftRight[iTof];
	h2GapVSTDiffCorr->Fill((corrLeft-corrRight)/2.0,iTof+1);
      }  
  }
  //  cout<<treeBeam->TDC_Trig[0][0]<<endl;
  return 0;
};


Long_t Det_ToF1::done_lar()
{
  for(int iTof=0;iTof<12;iTof++){
    Int_t iBin=h1Diff[iTof]->GetMaximumBin();
    Double_t max=h1Diff[iTof]->GetBinContent(iBin);
    Double_t threshold=max/5.0;
    Double_t low,high,middle;
    for(Int_t i=iBin-1;i>0;i--){
      if(h1Diff[iTof]->GetBinContent(i)<threshold){
	low=h1Diff[iTof]->GetBinCenter(i);
	break;
      }
    }
    for(Int_t i=iBin+1,n=h1Diff[iTof]->GetNbinsX();i<n;i++){
      if(h1Diff[iTof]->GetBinContent(i)<threshold){
	high=h1Diff[iTof]->GetBinCenter(i);
	break;      
      }
    }
    middle=low+high;
    middle/=2.0;
    cout<<"shiftLeftRight["<<iTof<<"]= "<<middle<<";"<<endl;
  }
  return 0;
};
