#include <Det_ToF2.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
using namespace std;

Det_ToF2::Det_ToF2(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_):Plugin(in_,out_,inf_,outf_,p_)
{
  // Set defaults for various options
  treeCali=0;
};


Det_ToF2::~Det_ToF2()
{
};

Long_t Det_ToF2::histos()
{
  return 0;

}

Long_t Det_ToF2::startup()
{
  getBranchObject("RawTof2Info",(TObject **) &treeRaw);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
  gStyle->SetOptStat(0);

  return 0;
};


Long_t Det_ToF2::process()
{
  return 0;
};


Long_t Det_ToF2::done()
{
  return 0;
};

Long_t Det_ToF2::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};

/////////////////////////////////////////////////////////////////////////////////////////////////// 

extern "C"{
  Plugin *factory(TTree *in_, TTree *out_, TFile *inf_, TFile * outf_,TObject *p_)
  {
      return (Plugin *) new Det_ToF2(in_,out_,inf_,outf_,p_);
  }
}


ClassImp(Det_ToF2);
