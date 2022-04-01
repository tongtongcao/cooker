#include <Det_ToF1.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <sstream>
using namespace std;
#define NSPERTDC 0.025


Long_t Det_ToF1::histos_veff()
{
  for(int iTof=0;iTof<12;iTof++){
    string name("Veff_");
    name+=iTof/10+'0';
    name+=iTof%10+'0';
    string title("Hit Position VS Half Time Difference;#frac{1}{2}(T_{U}-T{D}) (ns);Hit Position (cm)");
    h2HitVSDiff[iTof]=new TH2D(name.c_str(),title.c_str(),100,-2.,2.,100,-20.0,20.0);

  }

  return 0;

}

Long_t Det_ToF1::startup_veff()
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


Long_t Det_ToF1::process_veff()
{
  if(indexGap[treeRaw->event]==0) return 0;
  unsigned iTof=indexGap[treeRaw->event]-1;
  double corrLeft=treeRaw->TDC_TOF1U[iTof]*NSPERTDC-shiftLeftRight[iTof];
  double corrRight=treeRaw->TDC_TOF1D[iTof]*NSPERTDC+shiftLeftRight[iTof];
  double diff=corrLeft-corrRight;
  diff/=2.0;
  Double_t position=zGap[treeRaw->event];
  if(position<100)
    h2HitVSDiff[iTof]->Fill(diff,zGap[treeRaw->event]);

return 0;
};


Long_t Det_ToF1::done_veff()
{
  return 0;
};
