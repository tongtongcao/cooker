#include <Det_SFT.h>

#include <iostream>
#include <cmath>
#include <fstream>
#include <TStyle.h>
using namespace std;

// Constructor/destructor
Det_SFT::Det_SFT(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

Det_SFT::~Det_SFT()
{
};


// Standard recipe functions

Long_t Det_SFT::histos()
{

  char Title_SFT[128][100];	
  char Name_SFT[128][100];		

  for(int i=0;i<128;i++){
    sprintf(Title_SFT[i],"Raw ADC High Gain (Ch. %d)  --  SFT",i); 
    sprintf(Name_SFT[i],"ADC_High (Ch. %d) - Name_SFT",i);	
    hADC_HG_SFT[i] = new TH1D(Name_SFT[i],Title_SFT[i],512,0,4096);

    sprintf(Title_SFT[i],"Raw ADC Low Gain (Ch. %d)  --  SFT",i); 
    sprintf(Name_SFT[i],"ADC_Low (Ch. %d) - SFT",i);	
    hADC_LG_SFT[i] = new TH1D(Name_SFT[i],Title_SFT[i],512,0,4096);
 
    sprintf(Title_SFT[i],"Raw TDC (LE) (Ch. %d)  --  SFT",i); 
    sprintf(Name_SFT[i],"TDC_LE (Ch. %d) - SFT",i);	
    hTDC_LE_SFT[i] = new TH1D(Name_SFT[i],Title_SFT[i],512,0,4096);

    sprintf(Title_SFT[i],"Raw TDC (TE) (Ch. %d)  --  SFT",i); 
    sprintf(Name_SFT[i],"TDC_TE (Ch. %d) - SFT",i);	
    hTDC_TE_SFT[i] = new TH1D(Name_SFT[i],Title_SFT[i],512,0,4096);
  }


  return 0;
};


Long_t Det_SFT::startup()
{
  getBranchObject("RawSftInfo",(TObject **) &treeRaw);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
  //gStyle->SetOptStat(0);
  
  return 0;
};

Long_t Det_SFT::process()
{
  for(int i=0;i<128;i++){
    hADC_HG_SFT[i]->Fill(treeRaw->ADC_High_SFT[i]);
    hADC_LG_SFT[i]->Fill(treeRaw->ADC_Low_SFT[i]);

    if(treeRaw->TDC_LE_SFT[i].size()>0) hTDC_LE_SFT[i]->Fill(treeRaw->TDC_LE_SFT[i][0]);
    if(treeRaw->TDC_TE_SFT[i].size()>0) hTDC_TE_SFT[i]->Fill(treeRaw->TDC_TE_SFT[i][0]);

  }
  
  return 0;
};

Long_t Det_SFT::done()
{
  return 0;
};

//Long_t Det_SFT::finalize()
//{
//   return 0;
//};

// Cooker command line option function
Long_t Det_SFT::cmdline(char *cmd)
{
   // Add command line handling flags here

   return 0;
};


extern "C"
{
   Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
   {
      return (Plugin *) new Det_SFT(in,out,inf_,outf_,p);
   }
}

ClassImp(Det_SFT);
