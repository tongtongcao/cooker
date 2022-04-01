#include <Target.h>

#include <iostream>
#include <cmath>
#include <fstream>
#include <TStyle.h>
using namespace std;

// Constructor/destructor
Target::Target(TTree *in, TTree *out,TFile *inf_, TFile *outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

Target::~Target()
{
};

Long_t Target::histos()
{

  char Title_TARGET[256][100];	
  char Name_TARGET[256][100];		

  for(int i=0;i<256;i++){
    sprintf(Title_TARGET[i],"Raw ADC High Gain (Ch. %d)  --  TARGET",i); 
    sprintf(Name_TARGET[i],"ADC_High (Ch. %d) - TARGET",i);	
    hADC_HG[i] = new TH1D(Name_TARGET[i],Title_TARGET[i],512,0,4096);

    sprintf(Title_TARGET[i],"Raw ADC Low Gain (Ch. %d)  --  TARGET",i); 
    sprintf(Name_TARGET[i],"ADC_Low (Ch. %d) - TARGET",i);	
    hADC_LG[i] = new TH1D(Name_TARGET[i],Title_TARGET[i],512,0,4096);
 
    sprintf(Title_TARGET[i],"Raw TDC (LE) (Ch. %d)  --  TARGET",i); 
    sprintf(Name_TARGET[i],"TDC_LE (Ch. %d) - TARGET",i);	
    hTDC_LE[i] = new TH1D(Name_TARGET[i],Title_TARGET[i],512,0,4096);

    sprintf(Title_TARGET[i],"Raw TDC (TE) (Ch. %d)  --  TARGET",i); 
    sprintf(Name_TARGET[i],"TDC_TE (Ch. %d) - TARGET",i);	
    hTDC_TE[i] = new TH1D(Name_TARGET[i],Title_TARGET[i],512,0,4096);
  }


  return 0;
};

Long_t Target::startup()
{

  //getBranchObject("Tree",(TObject **) &treeRaw);  
	
  getBranchObject("RawTargetInfo",(TObject **) &treeRaw);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
  gStyle->SetOptStat(0);

  return 0;
};

Long_t Target::process()
{
  for(int i=0;i<256;i++){
    hADC_HG[i]->Fill(treeRaw->ADC_High_TARGET[i]);
    hADC_LG[i]->Fill(treeRaw->ADC_Low_TARGET[i]);

    if(treeRaw->TDC_LE_TARGET[i].size()>0) hTDC_LE[i]->Fill(treeRaw->TDC_LE_TARGET[i][0]);
    if(treeRaw->TDC_TE_TARGET[i].size()>0) hTDC_TE[i]->Fill(treeRaw->TDC_TE_TARGET[i][0]);

  }
  return 0;
};

Long_t Target::done()
{
  return 0;
};

Long_t Target::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};




// Standard recipe functions
//Long_t Target::startup()
//{
//   return 0;
//};

//Long_t Target::process()
//{
//   return 0;
//};

//Long_t Target::finalize()
//{
//   return 0;
//};

// Cooker command line option function
//Long_t Target::cmdline(char *cmd)
//{
   // Add command line handling flags here

//   return 0;
//};


extern "C"
{
   Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
   {
      return (Plugin *) new Target(in,out,inf_,outf_,p);
   }
}

ClassImp(Target);
