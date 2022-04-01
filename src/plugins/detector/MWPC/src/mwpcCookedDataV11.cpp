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

Long_t Det_MWPC::histos_mwpcCookedDataV11(){
  ostringstream sstr_name,sstr_title;

  return 0;
}

static vector<double> c2X,c2Y,c2Z,c3X,c3Y,c3Z,c4X,c4Y,c4Z;
static vector<double> c2XErr,c2YErr,c2ZErr,c3XErr,c3YErr,c3ZErr,c4XErr,c4YErr,c4ZErr;
static UInt_t run,event;
static int gapNum;
Long_t Det_MWPC::startup_mwpcCookedDataV11(){
  getBranchObject("CookMwpcInfoV1",(TObject **) &treeCookV1);
  out->Branch("run",&run,"run/i");
  out->Branch("event",&event,"event/i");
  out->Branch("gapNum",&gapNum,"gapNum/I");
  out->Branch("c2X",&c2X);
  out->Branch("c2Y",&c2Y);
  out->Branch("c2Z",&c2Z);
  out->Branch("c3X",&c3X);
  out->Branch("c3Y",&c3Y);
  out->Branch("c3Z",&c3Z);
  out->Branch("c4X",&c4X);
  out->Branch("c4Y",&c4Y);
  out->Branch("c4Z",&c4Z);

  out->Branch("c2XErr",&c2XErr);
  out->Branch("c2YErr",&c2YErr);
  out->Branch("c2ZErr",&c2ZErr);
  out->Branch("c3XErr",&c3XErr);
  out->Branch("c3YErr",&c3YErr);
  out->Branch("c3ZErr",&c3ZErr);
  out->Branch("c4XErr",&c4XErr);
  out->Branch("c4YErr",&c4YErr);
  out->Branch("c4ZErr",&c4ZErr);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);


  return 0;
}

Long_t Det_MWPC::process_mwpcCookedDataV11(){
  run=treeCookV1->run;
  event=treeCookV1->event;
  gapNum=treeCookV1->gapNum;
  c2X.clear();
  c2Y.clear();
  c2Z.clear();
  c3X.clear();
  c3Y.clear();
  c3Z.clear();
  c4X.clear();
  c4Y.clear();
  c4Z.clear();

  c2XErr.clear();
  c2YErr.clear();
  c2ZErr.clear();
  c3XErr.clear();
  c3YErr.clear();
  c3ZErr.clear();
  c4XErr.clear();
  c4YErr.clear();
  c4ZErr.clear();

  if(treeCookV1->posC2X!=-10000){
    c2X.push_back(treeCookV1->posC2X);
    c2XErr.push_back(1);
  }
  if(treeCookV1->posC2Y!=-10000){
    c2Y.push_back(treeCookV1->posC2Y);
    c2YErr.push_back(1);
  }
  if(treeCookV1->posC2Z!=-10000){
    c2Z.push_back(treeCookV1->posC2Z);
    c2ZErr.push_back(1);
  }

  if(treeCookV1->posC3X!=-10000){ 
    c3X.push_back(treeCookV1->posC3X);
    c3XErr.push_back(1);
  }
  if(treeCookV1->posC3Y!=-10000){ 
    c3Y.push_back(treeCookV1->posC3Y);
    c3YErr.push_back(1);
  }
  if(treeCookV1->posC3Z!=-10000){
    c3Z.push_back(treeCookV1->posC3Z);
    c3ZErr.push_back(1);
  }

  if(treeCookV1->posC4X!=-10000){ 
    c4X.push_back(treeCookV1->posC4X);
    c4XErr.push_back(1);
  }
  if(treeCookV1->posC4Y!=-10000){ 
    c4Y.push_back(treeCookV1->posC4Y);
    c4YErr.push_back(1);
  }
  if(treeCookV1->posC4Z!=-10000){
    c4Z.push_back(treeCookV1->posC4Z);
    c4ZErr.push_back(1);
  }

  return 0;
}

Long_t Det_MWPC::done_mwpcCookedDataV11(){


  return 0;
}

