#include <Det_tracking.h>
#include <iostream>
using namespace std;

Long_t Det_tracking::histos_basicDataType(){
  return 0;
}

  static UInt_t run,event;
static Int_t nTracks;
static Double_t fTof2SX[nTracks];
static Double_t fTof2SY[nTracks];
static Double_t fTof2SZ[nTracks];
static Double_t fTof2SNx[nTracks];
static Double_t fTof2SNy[nTracks];
static Double_t fTof2SNz[nTracks];
static Double_t fTof2SP[nTracks];

Long_t Det_tracking::startup_basicDataType(){
  getBranchObject("tracks",(TObject **) &trackArr);

  out->Branch("run",&run,"run/i");
  out->Branch("event",&event,"event/i");
  out->Branch("nTracks",&nTracks,"nTracks/I");
  out->Branch("nHits",nHits,"nHits[nTracks]/I");
  out->Branch("fTof2SX",fTof2SX,"fTof2SX[nTracks]/D");
  out->Branch("fTof2SY",fTof2SY,"fTof2SY[nTracks]/D");
  out->Branch("fTof2SZ",fTof2SZ,"fTof2SZ[nTracks]/D");
  out->Branch("fTof2SNx",fTof2SNx,"fTof2SNx[nTracks]/D");
  out->Branch("fTof2SNy",fTof2SNy,"fTof2SNy[nTracks]/D");

  out->Branch("fTof2SNz",fTof2SNz,"fTof2SNz[nTracks]/D");
  out->Branch("fTof2SP",fTof2SP,"fTof2SP[nTracks]/D");
  return 0;
}


Long_t Det_tracking::process_basicDataType(){

      // cout<<  trackArr->vectTrack.size()<<endl;
nTracks=trackArr->vectTrack.size(); 
for(int j=0; j<nTracks; j++){
      fTof2SX[j] = trackArr->vectTrack[j].getTof2SX();
      fTof2SY[j] = trackArr->vectTrack[j].getTof2SY();
      fTof2SZ[j] = trackArr->vectTrack[j].getTof2SZ();
      fTof2SNx[j] = trackArr->vectTrack[j].getTof2SNx();
      fTof2SNy[j] = trackArr->vectTrack[j].getTof2SNy();
      fTof2SNz[j] = trackArr->vectTrack[j].getTof2SNz();
      fTof2SP[j] = trackArr->vectTrack[j].getTof2SP();

}
  
  return 0;
}


Long_t Det_tracking::done_basicDataType(){


  return 0;
}



