#include <Det_ToF1.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <sstream>
using namespace std;
#define NSPERTDC 0.025


Long_t Det_ToF1::histos_timewalk()
{
  double range[12][2]={
    {1500,4000},
    {2500,3000},
    {1500,4000},
    {3000,4000},
    {3500,3000},
    {1000,3000},//gap 6
    {3000,2500},
    {2000,2000},
    {3000,2500},
    {2500,2000},
    {1000,2500},
    {1000,2500}
  };

  for(int iTof=0;iTof<12;iTof++){
    for(int iEnd=0;iEnd<2;iEnd++){
      string name("Timewalk_");
      name+=iTof/10+'0';
      name+=iTof%10+'0';
      name+='_';
      name+=iEnd+'0';
      tof1Timewalk[iTof][iEnd]=new TrekTimewalk(name);
      tof1Timewalk[iTof][iEnd]->setUpperX(range[iTof][iEnd]);
      tof1Timewalk[iTof][iEnd]->setNBinsX(100);
    }
  }

  return 0;

}

Long_t Det_ToF1::startup_timewalk()
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
  offset[0]=2.62555;
  veff[0]=9.51727;
  offset[1]=2.9544;
veff[1]=8.58775;
 offset[2]=2.24266;
veff[2]=8.5935;
 offset[3]=3.41194;
veff[3]= 7.43159;
 offset[4]=4.34046;veff[4]= 9.73211;
 offset[5]=3.95966;veff[5]= 7.71025;
 offset[6]=3.81489;veff[6]= 9.35991;
 offset[7]=3.7243;veff[7]= 9.88155;
 offset[8]=4.66487;veff[8]= 9.45376;
 offset[9]=4.77853;veff[9]= 8.89083;
 offset[10]=3.27032;veff[10]= 9.60284;
 offset[11]=2.86695;veff[11]= 9.77261;

 listPedestal[0]=45.0186;
 listPedestal[12]=60.0226;
 listPedestal[1]=67.484;
 listPedestal[13]=48.4007;
 listPedestal[2]=56.1659;
 listPedestal[14]=79.1373;
 listPedestal[3]=37.9318;
 listPedestal[15]=44.1452;
 listPedestal[4]=59.794;
 listPedestal[16]=62.797;
 listPedestal[5]=60.4463;
 listPedestal[17]=72.4111;
 listPedestal[6]=52.0727;
 listPedestal[18]=82.3057;
 listPedestal[7]=66.2435;
 listPedestal[19]=65.7156;
 listPedestal[8]=63.2392;
 listPedestal[20]=69.8056;
 listPedestal[9]=51.6298;
 listPedestal[21]=64.618;
 listPedestal[10]=67.5943;
 listPedestal[22]=56.4271;
 listPedestal[11]=49.2037;
 listPedestal[23]=63.7336;

  memset(indexGap,0,sizeof(indexGap));
  static ifstream fileTrack("/home/hlu/data/track/sebasft-3994c3-all.txt");
  string nameLine;  
  unsigned currentEvent=0;
  unsigned currentGap=0;
  for(unsigned i=0;i<12;i++)
    countGap[i]=0;
  for(unsigned i=0;i<500000;i++)
    zGap[i]=1000000;
  while(fileTrack.good()){
    getline(fileTrack,nameLine);
    if(fileTrack.good()){
      istringstream buffer(nameLine);
      unsigned tmpEvent,tag;
      buffer>>tmpEvent>>tag;
      if(tag==0){
	currentEvent=tmpEvent;
	unsigned dummy,gap;
	buffer>>dummy>>dummy>>dummy>>dummy>>dummy>>gap;
	currentGap=gap-1;
	indexGap[tmpEvent]=gap;
      }
      else if(tag==1){

	float dummy;
	float x,y,z;
	for(unsigned i=0;i<8;i++)
	  buffer>>dummy;
	buffer>>x>>y>>z;
	zGap[tmpEvent]=z;
	countGap[currentGap]++;
      }
    }

  }
  for(unsigned i=0;i<12;i++){
    cout<<i+1<<" "<<countGap[i]<<endl;
  }
  return 0;
};


Long_t Det_ToF1::process_timewalk()
{
  if(indexGap[treeRaw->event]==0) return 0;
  unsigned iTof=indexGap[treeRaw->event]-1;
  double corrLeft=treeRaw->TDC_TOF1U[iTof]*NSPERTDC-shiftLeftRight[iTof];
  double corrRight=treeRaw->TDC_TOF1D[iTof]*NSPERTDC+shiftLeftRight[iTof];
  double diff=corrLeft-corrRight;
  diff/=2.0;

  diff*=veff[iTof];
  diff+=offset[iTof];
  Double_t position=zGap[treeRaw->event];

  if(position<100){
    diff-=zGap[treeRaw->event];
    tof1Timewalk[iTof][0]->addData(diff,treeRaw->ADC_TOF1U[iTof]-listPedestal[iTof]);
    tof1Timewalk[iTof][1]->addData(-diff,treeRaw->ADC_TOF1D[iTof]-listPedestal[iTof+12]);
  }
  return 0;
};


Long_t Det_ToF1::done_timewalk()
{
  return 0;
};
