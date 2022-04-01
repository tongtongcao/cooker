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

Long_t Det_MWPC::histos_mwpcCookedDataV21(){
  ostringstream sstr_name,sstr_title;

  return 0;
}

static UInt_t run,event;
static UInt_t nHitsC2X[12], nHitsC2Y[12], nHitsC2Z[12], nHitsC3X[12], nHitsC3Y[12], nHitsC3Z[12], nHitsC4X[12], nHitsC4Y[12], nHitsC4Z[12];
static double posC2X[12][28], posC2Y[12][8], posC2Z[12][28], posC3X[12][32], posC3Y[12][8], posC3Z[12][1], posC4X[12][36], posC4Y[12][8], posC4Z[12][1];

Long_t Det_MWPC::startup_mwpcCookedDataV21(){
  getBranchObject("CookMwpcInfoV2",(TObject **) &treeCookV2);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);

  out->Branch("run",&run,"run/i");
  out->Branch("event",&event,"event/i");

  out->Branch("nHitsC2X",nHitsC2X,"nHitsC2X[12]/i");
  out->Branch("nHitsC2Y",nHitsC2Y,"nHitsC2Y[12]/i");
  out->Branch("nHitsC2Z",nHitsC2Z,"nHitsC2Z[12]/i");
  out->Branch("nHitsC3X",nHitsC3X,"nHitsC3X[12]/i");
  out->Branch("nHitsC3Y",nHitsC3Y,"nHitsC3Y[12]/i");
  out->Branch("nHitsC3Z",nHitsC3Z,"nHitsC3Z[12]/i");
  out->Branch("nHitsC4X",nHitsC4X,"nHitsC4X[12]/i");
  out->Branch("nHitsC4Y",nHitsC4Y,"nHitsC4Y[12]/i");
  out->Branch("nHitsC4Z",nHitsC4Z,"nHitsC4Z[12]/i");

  out->Branch("posC2X",posC2X,"posC2X[12][28]/D");
  out->Branch("posC2Y",posC2Y,"posC2Y[12][8]/D");
  out->Branch("posC2Z",posC2Z,"posC2Z[12][28]/D");
  out->Branch("posC3X",posC3X,"posC3X[12][32]/D");
  out->Branch("posC3Y",posC3Y,"posC3Y[12][8]/D");
  out->Branch("posC3Z",posC3Z,"posC3Z[12][1]/D");
  out->Branch("posC4X",posC4X,"posC4X[12][36]/D");
  out->Branch("posC4Y",posC4Y,"posC4Y[12][8]/D");
  out->Branch("posC4Z",posC4Z,"posC4Z[12][1]/D");




  return 0;
}

Long_t Det_MWPC::process_mwpcCookedDataV21(){
  run=treeCookV2->run;
  event=treeCookV2->event;


  for(int i=0;i<12;i++){
    nHitsC2X[i]=treeCookV2->posC2X[i].size();
    nHitsC2Y[i]=treeCookV2->posC2Y[i].size();
    nHitsC2Z[i]=treeCookV2->posC2Z[i].size();
    nHitsC3X[i]=treeCookV2->posC3X[i].size();
    nHitsC3Y[i]=treeCookV2->posC3Y[i].size();
    nHitsC3Z[i]=treeCookV2->posC3Z[i].size();
    nHitsC4X[i]=treeCookV2->posC4X[i].size();
    nHitsC4Y[i]=treeCookV2->posC4Y[i].size();
    nHitsC4Z[i]=treeCookV2->posC4Z[i].size();

    for(int j=0;j<nHitsC2X[i];j++){
      posC2X[i][j]=treeCookV2->posC2X[i][j];
    }
    for(int j=0;j<nHitsC2Y[i];j++){
      posC2Y[i][j]=treeCookV2->posC2Y[i][j];
    }
    for(int j=0;j<nHitsC2Z[i];j++){
      posC2Z[i][j]=treeCookV2->posC2Z[i][j];
    }

    for(int j=0;j<nHitsC3X[i];j++){
      posC3X[i][j]=treeCookV2->posC3X[i][j];
    }
    for(int j=0;j<nHitsC3Y[i];j++){
      posC3Y[i][j]=treeCookV2->posC3Y[i][j];
    }
    for(int j=0;j<nHitsC3Z[i];j++){
      posC3Z[i][j]=treeCookV2->posC3Z[i][j];
    }

    for(int j=0;j<nHitsC4X[i];j++){
      posC4X[i][j]=treeCookV2->posC4X[i][j];
    }
    for(int j=0;j<nHitsC4Y[i];j++){
      posC4Y[i][j]=treeCookV2->posC4Y[i][j];
    }
    for(int j=0;j<nHitsC4Z[i];j++){
      posC4Z[i][j]=treeCookV2->posC4Z[i][j];
    }
  }

  return 0;
}

Long_t Det_MWPC::done_mwpcCookedDataV21(){


  return 0;
}

