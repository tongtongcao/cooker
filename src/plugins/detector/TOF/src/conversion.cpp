#include <Det_TOF.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
using namespace std;



Long_t Det_TOF::histos_conversion()
{
  h1Multi=new TH1D("Mulplicity_TOF","Number of gaps coincidence between TOF1 and TOF2;Number of gaps",13,-0.5,12.5);
  h1MultiLow=new TH1D("MulplicityLow_TOF","Number of gaps with gap of TOF1+1=gap of TOF2;Number of gaps",13,-0.5,12.5);
  h1MultiUp=new TH1D("MulplicityUp_TOF","Number of gaps with gap of TOF1-1=gap of TOF2;Number of gaps",13,-0.5,12.5);
  h1MultiTotal=new TH1D("MulplicityTotal_TOF","Number of total good gaps;Number of gaps",13,-0.5,12.5);

  //  h1Diff=new TH1D("Time_Diff","Time difference between TOF1 and TOF2;TOF2-TOF1 (ns)",100,0,100);
  return 0;

}

Long_t Det_TOF::startup_conversion()
{
  getBranchObject("hitTof1",(TObject **) &treeTof1);
  getBranchObject("hitTof2",(TObject **) &treeTof2);

  // Make new branch on output tree                                                            
  //  makeBranch("hitTof1",(TObject **)&treeCali);

  gStyle->SetOptStat(0);

  return 0;
};


Long_t Det_TOF::process_conversion()
{
  static UInt_t nEvent=0;
  nEvent++;
  //  cout<<nEvent<<endl;
  UInt_t nGap=0;
  UInt_t nGapTof1=treeTof1->gap.size();
  UInt_t nGapTof2=treeTof2->gap.size();
  UInt_t iGapTof1=0;
  UInt_t iGapTof2=0;
  UInt_t nMatch=0;
  UInt_t nTotal=0;
  while(iGapTof1<nGapTof1 && iGapTof2<nGapTof2){
    UInt_t indexGapTof1=treeTof1->gap[iGapTof1];
    UInt_t indexGapTof2=treeTof2->gap[iGapTof2];
    
    if(indexGapTof1==indexGapTof2){
      nMatch++;
      iGapTof1++;
      iGapTof2++;
    }
    else if(indexGapTof1<indexGapTof2){
      iGapTof1++;
    }
    else{
      iGapTof2++;
    }
  }
  h1Multi->Fill(nMatch);
  nTotal+=nMatch; 
  if(nMatch==0){
    iGapTof1=0;
    iGapTof2=0;
    UInt_t nLow=0;
    while(iGapTof1<nGapTof1 && iGapTof2<nGapTof2){
      UInt_t indexGapTof1=treeTof1->gap[iGapTof1];
      UInt_t indexGapTof2=treeTof2->gap[iGapTof2];
      if((indexGapTof1+1==indexGapTof2)||(indexGapTof1==12 && indexGapTof2==1)){
	nLow++;
	iGapTof1++;
	iGapTof2++;
      }
      else if(indexGapTof1+1<indexGapTof2){
	iGapTof1++;
      }
      else{
	iGapTof2++;
      }
    }
    h1MultiLow->Fill(nLow);
    nTotal+=nLow;
    if(nLow==0){
      iGapTof1=0;
      iGapTof2=0;
      UInt_t nUp=0;
      while(iGapTof1<nGapTof1 && iGapTof2<nGapTof2){
	UInt_t indexGapTof1=treeTof1->gap[iGapTof1];
	UInt_t indexGapTof2=treeTof2->gap[iGapTof2];
	if(indexGapTof1-1==indexGapTof2 || (indexGapTof1==1 && indexGapTof2==12)){
	  nUp++;
	  iGapTof1++;
	  iGapTof2++;
	}
	else if(indexGapTof1-1<indexGapTof2){
	  iGapTof1++;
	}
	else{
	  iGapTof2++;
	}
      }
      h1MultiUp->Fill(nUp);
      
      nTotal+=nUp;
    }
  }
  h1MultiTotal->Fill(nTotal);
  if(nTotal==0){
    cout<<endl;
    cout<<" tof1:"<<nGapTof1<<endl;

    for(int iGapTof1=0;iGapTof1<nGapTof1;iGapTof1++){
      cout<<" "<<treeTof1->gap[iGapTof1];
    }
    cout<<endl;
    cout<<" tof2:"<<nGapTof2<<endl;
    for(int iGapTof2=0;iGapTof2<nGapTof2;iGapTof2++){
      cout<<" "<<treeTof2->gap[iGapTof2];
    }

  }
  //  cout<<treeBeam->TDC_Trig[0][0]<<endl;
  return 0;
};


Long_t Det_TOF::done_conversion()
{

  return 0;
};
