#include <Det_tracking.h>
#include <iostream>
using namespace std;

Long_t Det_tracking::histos_basicDataType(){
  return 0;
}

static UInt_t run,event;

static Int_t nTracks=MAXNTRACKS;
static Int_t nHits[MAXNTRACKS];
static Double_t chi2PerNDF[MAXNTRACKS];
static Int_t flagTargetToAC[MAXNTRACKS];
static Bool_t flagTofMatch[MAXNTRACKS];

static Double_t sMwpc4[MAXNTRACKS];
static Double_t tMwpc4[MAXNTRACKS];
static Double_t eLossMwpc4[MAXNTRACKS];
static Double_t mXMwpc4[MAXNTRACKS];
static Double_t mYMwpc4[MAXNTRACKS];
static Double_t mZMwpc4[MAXNTRACKS];
static Double_t mXErrMwpc4[MAXNTRACKS];
static Double_t mYErrMwpc4[MAXNTRACKS];
static Double_t mZErrMwpc4[MAXNTRACKS];
static Double_t xMwpc4[MAXNTRACKS];
static Double_t yMwpc4[MAXNTRACKS];
static Double_t zMwpc4[MAXNTRACKS];
static Double_t nXMwpc4[MAXNTRACKS];
static Double_t nYMwpc4[MAXNTRACKS];
static Double_t nZMwpc4[MAXNTRACKS];
static Double_t pMwpc4[MAXNTRACKS];

Long_t Det_tracking::startup_basicDataType(){
  getBranchObject("tracks",(TObject **) &trackArr);


  out->Branch("run",&run,"run/i");
  out->Branch("event",&event,"event/i");

  out->Branch("nTracks",&nTracks,"nTracks/I");
  out->Branch("nHits",nHits,"nHits[nTracks]/I");
  out->Branch("chi2PerNDF",chi2PerNDF,"chi2PerNDF[nTracks]/D");
  out->Branch("flagTargetToAC",flagTargetToAC,"flagTargetToAC[nTracks]/I");
  out->Branch("flagTofMatch",flagTofMatch,"flagTofMatch[nTracks]/O");

  out->Branch("sMwpc4",sMwpc4,"sMwpc4[nTracks]/D");
  out->Branch("tMwpc4",tMwpc4,"tMwpc4[nTracks]/D");
  out->Branch("eLossMwpc4",eLossMwpc4,"eLossMwpc4[nTracks]/D");
  out->Branch("mXMwpc4",mXMwpc4,"mXMwpc4[nTracks]/D");
  out->Branch("mYMwpc4",mYMwpc4,"mYMwpc4[nTracks]/D");
  out->Branch("mZMwpc4",mZMwpc4,"mZMwpc4[nTracks]/D");
  out->Branch("mXErrMwpc4",mXErrMwpc4,"mXErrMwpc4[nTracks]/D");
  out->Branch("mYErrMwpc4",mYErrMwpc4,"mYErrMwpc4[nTracks]/D");
  out->Branch("mZErrMwpc4",mZErrMwpc4,"mZErrMwpc4[nTracks]/D");
  out->Branch("xMwpc4",xMwpc4,"xMwpc4[nTracks]/D");
  out->Branch("yMwpc4",yMwpc4,"yMwpc4[nTracks]/D");
  out->Branch("zMwpc4",zMwpc4,"zMwpc4[nTracks]/D");
  out->Branch("nXMwpc4",nXMwpc4,"nXMwpc4[nTracks]/D");
  out->Branch("nYMwpc4",nYMwpc4,"nYMwpc4[nTracks]/D");
  out->Branch("nZMwpc4",nZMwpc4,"nZMwpc4[nTracks]/D");
  out->Branch("pMwpc4",pMwpc4,"pMwpc4[nTracks]/D");


  return 0;
}


Long_t Det_tracking::process_basicDataType(){

cout<<  trackArr->vectTrack.size()<<endl;
  nTracks=trackArr->vectTrack.size();
  for(int j=0; j<nTracks; j++){
    //-----------------------------------Muon---------------------------------------//
    nHits[j]=trackArr->vectTrack[j].getNHits();
    chi2PerNDF[j]=trackArr->vectTrack[j].getChi2PerNDF();
    flagTargetToAC[j]=trackArr->vectTrack[j].getFlagTargetToAC();
    flagTofMatch[j]=trackArr->vectTrack[j].getFlagToFMatch();

    sMwpc4[j]=trackArr->vectTrack[j].getMwpc4PathLocation();
    tMwpc4[j]=trackArr->vectTrack[j].getMwpc4PathTime();
    eLossMwpc4[j]=trackArr->vectTrack[j].getMwpc4EnergyLoss();
    mXMwpc4[j]=trackArr->vectTrack[j].getMwpc4MX();
    mYMwpc4[j]=trackArr->vectTrack[j].getMwpc4MY();
      mZMwpc4[j]=trackArr->vectTrack[j].getMwpc4MZ();
      mXErrMwpc4[j]=trackArr->vectTrack[j].getMwpc4MXErr();
      mYErrMwpc4[j]=trackArr->vectTrack[j].getMwpc4MYErr();
      mZErrMwpc4[j]=trackArr->vectTrack[j].getMwpc4MZErr();
      xMwpc4[j]=trackArr->vectTrack[j].getMwpc4SX();
      yMwpc4[j]=trackArr->vectTrack[j].getMwpc4SY();
      zMwpc4[j]=trackArr->vectTrack[j].getMwpc4SZ();
      nXMwpc4[j]=trackArr->vectTrack[j].getMwpc4SNx();
      nYMwpc4[j]=trackArr->vectTrack[j].getMwpc4SNy();
      nZMwpc4[j]=trackArr->vectTrack[j].getMwpc4SNz();
      pMwpc4[j]=trackArr->vectTrack[j].getMwpc4SP();

cout<<pMwpc4[j]<<endl;


  }
  return 0;
}


Long_t Det_tracking::done_basicDataType(){


  return 0;
}



