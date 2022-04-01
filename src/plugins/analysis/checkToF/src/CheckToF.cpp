#include <CheckToF.h>

#include<cmath>

#include <iostream>
#include <TStyle.h>
#include <TProfile.h>
using namespace std;

CheckToF::CheckToF(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_):Plugin(in_,out_,inf_,outf_,p_)
{
  // Set defaults for various options
};


CheckToF::~CheckToF()
{
};
Long_t CheckToF::set_tofDiff(int id, double offset)
{
  cout<<"set_tofDiff:"<<id<<" "<<offset<<endl;
  tofDiff[id]=offset;
  return 0;
};

Long_t CheckToF::histos()
{
  h2TimeVSAngle=new TH2D("TimeVSAngle","Time difference between Tof2 and 1 vs. #theta;#theta;time",100,0.0,180.0,100,0.0,100.0);
  h1TimeDiff=new TH1D("TimeDiff","Time difference between Tof2 and Tof1;Tof2-Tof1",100,0.0,100.0);
  h2TimeDiffVSGap=new TH2D("TimeDiffVSGap","Time difference between Tof VS. gap;gap;Tof2-Tof1",12,0.5,12.5,100,200.0,1000.0);
  h2MomentumVSCz=new TH2D("MomentumVSCz","Momentum VS cos(#theta);cos(#theta);Momentum (MeV/c)",100,-1.0,1.0,100,100,300);
  h2DeltaTVSCzMu=new TH2D("DeltaTVSCz_Mu","TOF time difference VS cos(#theta);cos(#theta);Tof2-Tof1 (ns)",100,-1,1,100,5,15);
  h2DeltaTVSCzPi=new TH2D("DeltaTVSCz_Pi","TOF time difference VS cos(#theta);cos(#theta);Tof2-Tof1 (ns)",100,-1,1,100,5,15);
  h2DeltaTVSCzPositron=new TH2D("DeltaTVSCz_Positron","TOF time difference VS cos(#theta);cos(#theta);Tof2-Tof1 (ns)",100,-1,1,100,5,15);


  for(int i=0;i<12;i++){
    for(int j=0;j<2;j++){
      TString name;
      name.Form("TOF1VSTOF2_%c_Gap%d",'A'+j,i+1);
      TString title;
      title.Form("Tof 1 vs tof2 %c for gap %d;tof2;tof1",'a'+j,i+1);
      h2TOF1VSTOF2[j][i]=new TH2D(name,title,100,0,3000,100,0,3000);
      name.Form("DeltaTVSCz_%c_Gap%d",'A'+j,i+1);
      title.Form("Tof2-Tof1 VS cos(#theta) for momentum between [230,240];cos(#theta);flight time (ns)");
      h2DeltaTVSCz[j][i]=new TH2D(name,title,100,-1.0,1.0,100,5,25);
    }
  }
  nameTrack="/home/hlu/data/track/track-run2800.dat";
  fileTrack.open(nameTrack);
  UInt_t iEvent,iGap;
  Double_t momentum;
  Double_t x,y,z;
  Double_t cx,cy,cz;
  while(fileTrack.good()){
    fileTrack>>iEvent>>momentum>>x>>y>>z>>cx>>cy>>cz>>iGap;
    if(fileTrack.good()){
      mapMomentum[iEvent]=momentum;
      mapGap[iEvent]=iGap;
      mapCz[iEvent]=cz;
      h2MomentumVSCz->Fill(cz,momentum);
    }
  }
  fileTrack.close();
  return 0;

}

Long_t CheckToF::startup()
{

  getBranchObject("RawTof1Info",(TObject **) &treeTof1);
  getBranchObject("RawTof2Info",(TObject **) &treeTof2);
  gStyle->SetOptStat(0);
  tofDiff[0]=13.6354;
  tofDiff[12]=15.6867;
  tofDiff[1]=16.0202;
  tofDiff[13]=15.94;
  tofDiff[2]=10.7359;
  tofDiff[14]=12.1655;
  tofDiff[3]=13.3166;
  tofDiff[15]=14.8267;
  tofDiff[4]=12.1622;
  tofDiff[16]=15.8662;
  tofDiff[5]=10.6512;
  tofDiff[17]=14.6548;
  tofDiff[6]=13.5342;
  tofDiff[18]=16.1861;
  tofDiff[7]=12.6835;
  tofDiff[19]=14.1794;
  tofDiff[8]=12.3781;
  tofDiff[20]=18.131;
  tofDiff[9]=0;
  tofDiff[21]=15.8414;
  tofDiff[10]=12.118;
  tofDiff[22]=14.335;
  tofDiff[11]=10.7256;
  tofDiff[23]=12.7888;
  //listEvent
  return 0;
};


Long_t CheckToF::process()
{
  const float calibratedTD=9.0568;
  int goodGap=-1;
  UInt_t event=treeTof1->event;
  if(mapGap.find(event)!=mapGap.end()){
    goodGap=mapGap[treeTof1->event];
  }

  float timeTof1=0,timeTof2=0;
  for(int i=0;i<12;i++){
    timeTof1=-1;
    timeTof2=-1;
    if(treeTof1->TDC_TOF1U[i]>0 && treeTof1->TDC_TOF1D[i]>0){
      bool findA=false;
      bool findB=false;
      timeTof1=treeTof1->TDC_TOF1U[i]+treeTof1->TDC_TOF1D[i];
      timeTof1/=2.0;
      //      cout<<"tof1:"<<timeTof1<<endl;
      if(treeTof2->TDC_TOF2AO[i]>0 && treeTof2->TDC_TOF2AI[i]>0){
	timeTof2=treeTof2->TDC_TOF2AO[i]+treeTof2->TDC_TOF2AI[i];
	timeTof2/=2.0;
	h2TOF1VSTOF2[0][i]->Fill(timeTof2,timeTof1);
	findA=true;
	//	cout<<"tof2:"<<timeTof2<<endl;
      }
      if(treeTof2->TDC_TOF2BO[i]>0 && treeTof2->TDC_TOF2BI[i]>0){
	timeTof2=treeTof2->TDC_TOF2BO[i]+treeTof2->TDC_TOF2BI[i];
	timeTof2/=2.0;
	h2TOF1VSTOF2[1][i]->Fill(timeTof2,timeTof1);
	findB=true;
	//	cout<<"tof2:"<<timeTof2<<endl;
      }
      if(timeTof1>0 && timeTof2>0){

	h1TimeDiff->Fill(timeTof2-timeTof1);
	UInt_t myEffGap=0;
	if(goodGap==i+1){
	  if(findB){
	    if(mapMomentum[event]>=230 && mapMomentum[event]<=240){
	      h2DeltaTVSCz[1][i]->Fill(mapCz[event],(timeTof2-timeTof1)*0.025);
	    }
	    myEffGap=i+12;
	  }
	  if(findA){
	    if(mapMomentum[event]>=230 && mapMomentum[event]<=240){
	      h2DeltaTVSCz[0][i]->Fill(mapCz[event],(timeTof2-timeTof1)*0.025);
	    }
	    myEffGap=i;
	  }
	  h2TimeDiffVSGap->Fill(i+1,timeTof2-timeTof1-tofDiff[i]+calibratedTD);
	  if(mapMomentum[event]>=230 && mapMomentum[event]<=240){
	      h2DeltaTVSCzMu->Fill(mapCz[event],(timeTof2-timeTof1)*0.025-tofDiff[myEffGap]+calibratedTD);
	  }
	  if(mapMomentum[event]>=195 && mapMomentum[event]<=210){
	    h2DeltaTVSCzPi->Fill(mapCz[event],(timeTof2-timeTof1)*0.025-tofDiff[myEffGap]+calibratedTD);
	  }
	  if(mapMomentum[event]>240){
	    h2DeltaTVSCzPositron->Fill(mapCz[event],(timeTof2-timeTof1)*0.025-tofDiff[myEffGap]+calibratedTD);
	  }

	}
      }
    }
  }
  //  cout<<treeBeam->TDC_Trig[0][0]<<endl;
  return 0;
};

Double_t fitCz(TH2D* h2,const TString& name){
  if(h2->GetEntries()<100) return 0;
  TProfile* h1=(TProfile*)h2->ProfileX();
  h1->SetName(name);
  
  h1->Fit("pol1","Q","",-0.15,0.15);
  return h1->GetFunction("pol1")->GetParameter(0);
}
Long_t CheckToF::done()
{
  TString name;
  for(UInt_t i=0;i<12;i++){
    cout<<"tofDiff["<<i<<"]="<<fitCz(h2DeltaTVSCz[0][i],name.Format("Cz_A_Gap%d",i))<<";"<<endl;
    cout<<"tofDiff["<<i+12<<"]="<<fitCz(h2DeltaTVSCz[1][i],name.Format("Cz_B_Gap%d",i))<<";"<<endl;

    /*
    cout<<"<tofDiff id=\""<<i<<"\">"<<fitCz(h2DeltaTVSCz[0][i],name.Format("Cz_A_Gap%d",i))<<"</tofDiff>"<<endl;
    cout<<"<tofDiff id=\""<<i+12<<"\">"<<fitCz(h2DeltaTVSCz[1][i],name.Format("Cz_A_Gap%d",i))<<"</tofDiff>"<<endl;
    */

  }
  return 0;
};

Long_t CheckToF::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};

/////////////////////////////////////////////////////////////////////////////////////////////////// 

extern "C"{
  Plugin *factory(TTree *in_, TTree *out_, TFile *inf_, TFile * outf_,TObject *p_)
  {
      return (Plugin *) new CheckToF(in_,out_,inf_,outf_,p_);
  }
}


ClassImp(CheckToF);
