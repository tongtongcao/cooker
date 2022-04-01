#include <Det_tracking.h>
#include "TClonesArray.h"
#include <iostream>
#include <cmath>

// Constructor/destructor
Det_tracking::Det_tracking(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p){
}

Det_tracking::~Det_tracking(){}

// Standard recipe functions
Long_t Det_tracking::histos(){
  return 0;
}

Long_t Det_tracking::startup(){
  getBranchObject("caliTargetTofMwpc",(TObject **) &caliTargetTofMwpc);
  trackArr=0;
  out->Branch("tracks",&trackArr,160000,0);
  uSft=new sft(kSft, kPlaneDirX);
  uSft->init(sftX, sftY, sftZ, sftXErr, sftYErr, sftZErr);
  uAc=new ac(kAc, kPlaneDirX);
  

  return 0;
}

Long_t Det_tracking::process(){
 
  mwpc *mwpc2=new mwpc(kMwpc2,kPlaneDirX);
  mwpc *mwpc3=new mwpc(kMwpc3,kPlaneDirZ);
  mwpc *mwpc4=new mwpc(kMwpc4,kPlaneDirZ);
  mwpc2->init(caliTargetTofMwpc->nHitsC2X, caliTargetTofMwpc->nHitsC2Y, caliTargetTofMwpc->nHitsC2Z, caliTargetTofMwpc->posC2X, caliTargetTofMwpc->posC2Y, caliTargetTofMwpc->posC2Z);
  mwpc3->init(caliTargetTofMwpc->nHitsC3X, caliTargetTofMwpc->nHitsC3Y, caliTargetTofMwpc->nHitsC3Z, caliTargetTofMwpc->posC3X, caliTargetTofMwpc->posC3Y, caliTargetTofMwpc->posC3Z);
  mwpc4->init(caliTargetTofMwpc->nHitsC4X, caliTargetTofMwpc->nHitsC4Y, caliTargetTofMwpc->nHitsC4Z, caliTargetTofMwpc->posC4X, caliTargetTofMwpc->posC4Y, caliTargetTofMwpc->posC4Z);

  tof2 *uTof2=new tof2(kTof2,kPlaneDirZ, caliTargetTofMwpc->gapNumTof2, caliTargetTofMwpc->timeAITof2, caliTargetTofMwpc->timeAOTof2, caliTargetTofMwpc->timeBITof2, caliTargetTofMwpc->timeBOTof2);

  tof1 *uTof1=new tof1(kTof1,kPlaneDirX, caliTargetTofMwpc->gapNumTof1, caliTargetTofMwpc->timeUTof1, caliTargetTofMwpc->timeDTof1);

  target *uTarget=new target(kTarget, kPlaneDirX, caliTargetTofMwpc->tag, caliTargetTofMwpc->vertX, caliTargetTofMwpc->vertY, caliTargetTofMwpc->vertPhi, caliTargetTofMwpc->vertXErr, caliTargetTofMwpc->vertYErr, caliTargetTofMwpc->vertPhiErr);
  uTarget->setNHits(0);

//  partId=kPionPlusId;
  partId=kMuonId;
  // partId=kPositronId;
  trackFinder *trackFind=new trackFinder(partId, kParticleMass[partId], kParticleCharge[partId], uTarget, uSft, uTof1, uAc, mwpc2, mwpc3, mwpc4, uTof2);

  TClonesArray *theTracks=new TClonesArray("track", 1000, kTRUE);
  theTracks->SetOwner(kTRUE);
  theTracks->Clear();
  trackFind->processHits(theTracks);

  trackArr->vectTrack.clear(); 

  int nHits[100];
  double chi2PerNDF[100];
  for(int j=0; j<theTracks->GetEntries(); j++){
    trackArr->vectTrack.push_back(*((track*)theTracks->At(j)));
  }

  theTracks->SetOwner(kTRUE);
  theTracks->Delete();
  delete trackFind;
  trackFind=NULL;


  delete mwpc2;
  mwpc2=NULL;
  delete mwpc3;
  mwpc3=NULL;
  delete mwpc4;
  mwpc4=NULL;
  delete uTof2;
  uTof2=NULL;
  delete uTof1;
  uTof1=NULL;

  delete uTarget;
  uTarget=NULL;

  return 0;
}

Long_t Det_tracking::done(){
  delete uSft;
  delete uAc;

  return 0;
}

// Cooker command line option function
Long_t Det_tracking::cmdline(char *cmd){
   // Add command line handling flags here

  return 0;
};


extern "C"
{
  Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p){
    return (Plugin *) new Det_tracking(in,out,inf_,outf_,p);
  }
}

ClassImp(Det_tracking);
