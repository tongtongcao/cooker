#include <Det_tracking.h>
#include <iostream>
using namespace std;

Long_t Det_tracking::histos_basicDataType(){
  return 0;
}

  static UInt_t run,event;
static Int_t nTracks=MAXNTRACKS;
static Int_t nHits[MAXNTRACKS];
Int_t fgapNumTof2[MAXNTRACKS];
Int_t fgapNumTof1[MAXNTRACKS];
static Double_t fTof1SX[MAXNTRACKS];
static Double_t fTof1SY[MAXNTRACKS];
static Double_t fTof1SZ[MAXNTRACKS];
static Double_t fTof1SNx[MAXNTRACKS];
static Double_t fTof1SNy[MAXNTRACKS];
static Double_t fTof1SNz[MAXNTRACKS];


static Double_t fTof2SX[MAXNTRACKS];
static Double_t fTof2SY[MAXNTRACKS];
static Double_t fTof2SZ[MAXNTRACKS];
static Double_t fTof2SNx[MAXNTRACKS];
static Double_t fTof2SNy[MAXNTRACKS];
static Double_t fTof2SNz[MAXNTRACKS];
static Double_t fTof2SP[MAXNTRACKS];
static Double_t fVertMX[MAXNTRACKS];
static Double_t fVertMY[MAXNTRACKS];
static Double_t fVertMZ[MAXNTRACKS];
static Double_t fVertMPhi[MAXNTRACKS];
static Double_t fSftSX[MAXNTRACKS];
static Double_t fSftSY[MAXNTRACKS];
static Double_t fSftSZ[MAXNTRACKS];
static Double_t fVertSX[MAXNTRACKS];
static Double_t fVertSY[MAXNTRACKS];
static Double_t fVertSZ[MAXNTRACKS];

static Double_t fVertSPhi[MAXNTRACKS];

static Double_t fSftSNx[MAXNTRACKS];
static Double_t fSftSNy[MAXNTRACKS];
static Double_t fSftSNz[MAXNTRACKS];

static Double_t fTof1TimeMean[MAXNTRACKS];
static Double_t fTof2TimeMean[MAXNTRACKS];
static Double_t fTof1PathLocation[MAXNTRACKS];
static Double_t fTof2PathLocation[MAXNTRACKS];
static Double_t fTof1PathTime[MAXNTRACKS];
static Double_t fTof2PathTime[MAXNTRACKS];
static Double_t fTof1SP[MAXNTRACKS];
static Double_t fChi2PerNDF[MAXNTRACKS];
static Int_t fFlagTargetToAC[MAXNTRACKS];
static Bool_t fFlagToFMatch[MAXNTRACKS];

static Double_t fTof1TimeD[MAXNTRACKS];
static Double_t fTof1TimeU[MAXNTRACKS];
static Double_t fTof2TimeI[MAXNTRACKS];
static Double_t fTof2TimeO[MAXNTRACKS];
//c4
static Double_t fMwpc4SX[MAXNTRACKS];
static Double_t fMwpc4SY[MAXNTRACKS];
static Double_t fMwpc4SZ[MAXNTRACKS];
static Double_t fMwpc4SNx[MAXNTRACKS];
static Double_t fMwpc4SNy[MAXNTRACKS];
static Double_t fMwpc4SNz[MAXNTRACKS];
static Double_t fMwpc4SP[MAXNTRACKS];
static Double_t fMwpc4PathLocation[MAXNTRACKS];
static Double_t fMwpc4PathTime[MAXNTRACKS];
//c2
static Double_t fMwpc2SP[MAXNTRACKS];
static Double_t fMwpc2SX[MAXNTRACKS];
static Double_t fMwpc2SY[MAXNTRACKS];
static Double_t fMwpc2SZ[MAXNTRACKS];
static Double_t fMwpc2SNx[MAXNTRACKS];
static Double_t fMwpc2SNy[MAXNTRACKS];
static Double_t fMwpc2SNz[MAXNTRACKS];
static Double_t fMwpc2PathLocation[MAXNTRACKS];
static Double_t fMwpc2PathTime[MAXNTRACKS];

//c3
static Double_t fMwpc3SP[MAXNTRACKS];
static Double_t fMwpc3SX[MAXNTRACKS];
static Double_t fMwpc3SY[MAXNTRACKS];
static Double_t fMwpc3SZ[MAXNTRACKS];
static Double_t fMwpc3SNx[MAXNTRACKS];
static Double_t fMwpc3SNy[MAXNTRACKS];
static Double_t fMwpc3SNz[MAXNTRACKS];
static Double_t fMwpc3PathLocation[MAXNTRACKS];
static Double_t fMwpc3PathTime[MAXNTRACKS];

static Double_t fVertSP[MAXNTRACKS];
static Double_t fVertSPPionPlus[MAXNTRACKS];
static Double_t fVertSPPositron[MAXNTRACKS];
static Double_t fTof2SPPionPlus[MAXNTRACKS];
static Double_t fTof2SPPositron[MAXNTRACKS];

static Double_t fMwpc2EnergyLoss[MAXNTRACKS];
static Double_t fMwpc3EnergyLoss[MAXNTRACKS];
static Double_t fMwpc4EnergyLoss[MAXNTRACKS];


Long_t Det_tracking::startup_basicDataType(){
  getBranchObject("tracks",(TObject **) &trackArr);
  getBranchObject("RawAcInfo",(TObject **) &treeRaw);

  out->Branch("run",&run,"run/i");
  out->Branch("event",&event,"event/i");
  out->Branch("nTracks",&nTracks,"nTracks/I");
  out->Branch("fgapNumTof2",fgapNumTof2,"fgapNumTof2[nTracks]/I");
  out->Branch("fgapNumTof1",fgapNumTof1,"fgapNumTof1[nTracks]/I");
  out->Branch("nHits",nHits,"nHits[nTracks]/I");
  out->Branch("fTof1TimeD",fTof1TimeD,"fTof1TimeD[nTracks]/D");
  out->Branch("fTof1TimeU",fTof1TimeU,"fTof1TimeU[nTracks]/D");
  out->Branch("fTof2TimeI",fTof2TimeI,"fTof2TimeI[nTracks]/D");
  out->Branch("fTof2TimeO",fTof2TimeO,"fTof2TimeO[nTracks]/D");
  out->Branch("fTof1SX",fTof1SX,"fTof1SX[nTracks]/D");
  out->Branch("fTof1SY",fTof1SY,"fTof1SY[nTracks]/D");
  out->Branch("fTof1SZ",fTof1SZ,"fTof1SZ[nTracks]/D");
  out->Branch("fTof1SNx",fTof1SNx,"fTof1SNx[nTracks]/D");
  out->Branch("fTof1SNy",fTof1SNy,"fTof1SNy[nTracks]/D");
  out->Branch("fTof1SNz",fTof1SNz,"fTof1SNz[nTracks]/D");


  out->Branch("fTof2SX",fTof2SX,"fTof2SX[nTracks]/D");
  out->Branch("fTof2SY",fTof2SY,"fTof2SY[nTracks]/D");
  out->Branch("fTof2SZ",fTof2SZ,"fTof2SZ[nTracks]/D");
  out->Branch("fTof2SNx",fTof2SNx,"fTof2SNx[nTracks]/D");
  out->Branch("fTof2SNy",fTof2SNy,"fTof2SNy[nTracks]/D");
  out->Branch("fTof2SNz",fTof2SNz,"fTof2SNz[nTracks]/D");
  out->Branch("fTof2SP",fTof2SP,"fTof2SP[nTracks]/D");
  out->Branch("fVertSX",fVertSX,"fVertSX[nTracks]/D");
  out->Branch("fVertSZ",fVertSZ,"fVertSZ[nTracks]/D");
  out->Branch("fVertSY",fVertSY,"fVertSY[nTracks]/D");
  out->Branch("fVertSP",fVertSP,"fVertSP[nTracks]/D");
  out->Branch("fVertSPPositron",fVertSPPositron,"fVertSPPositron[nTracks]/D");
  out->Branch("fVertSPPionPlus",fVertSPPionPlus,"fVertSPPionPlus[nTracks]/D");
  out->Branch("fVertSPhi",fVertMPhi,"fVertMPhi[nTracks]/D");
  out->Branch("fVertMX",fVertMX,"fVertMX[nTracks]/D");
  out->Branch("fVertMY",fVertMY,"fVertMY[nTracks]/D");
  out->Branch("fVertMPhi",fVertMPhi,"fVertMPhi[nTracks]/D");
  out->Branch("fSftSX",fSftSX,"fSftSX[nTracks]/D");
  out->Branch("fSftSZ",fSftSZ,"fSftSZ[nTracks]/D");
  out->Branch("fSftSY",fSftSY,"fSftSY[nTracks]/D");
  out->Branch("fSftSNx",fSftSNx,"fSftSNx[nTracks]/D");
  out->Branch("fSftSNz",fSftSNz,"fSftSNz[nTracks]/D");
  out->Branch("fSftSNy",fSftSNy,"fSftSNy[nTracks]/D");

  out->Branch("fTof1TimeMean",fTof1TimeMean,"fTof1TimeMean[nTracks]/D");
  out->Branch("fTof2TimeMean",fTof2TimeMean,"fTof2TimeMean[nTracks]/D");
  
  out->Branch("fTof1SP",fTof1SP,"fTof1SP[nTracks]/D");
  out->Branch("fTof1PathLocation",fTof1PathLocation,"fTof1PathLocation[nTracks]/D");
  out->Branch("fTof2PathLocation",fTof2PathLocation,"fTof2PathLocation[nTracks]/D");

  out->Branch("fTof1PathTime",fTof1PathTime,"fTof1PathTime[nTracks]/D");
  out->Branch("fTof2PathTime",fTof2PathTime,"fTof2PathTime[nTracks]/D");
  out->Branch("fChi2PerNDF",fChi2PerNDF,"fChi2PerNDF[nTracks]/D");
  out->Branch("fFlagTargetToAC",fFlagTargetToAC,"fFlagTargetToAC[nTracks]/D");
  out->Branch("fFlagToFMatch",fFlagToFMatch,"fFlagToFMatch[nTracks]/D");

  out->Branch("fMwpc3SP",fMwpc3SP,"fMwpc3SP[nTracks]/D");
  out->Branch("fMwpc3SX",fMwpc3SX,"fMwpc3SX[nTracks]/D");
  out->Branch("fMwpc3SY",fMwpc3SY,"fMwpc3SY[nTracks]/D");
  out->Branch("fMwpc3SZ",fMwpc3SZ,"fMwpc3SZ[nTracks]/D");
  out->Branch("fMwpc3SNx",fMwpc3SNx,"fMwpc3SNx[nTracks]/D");
  out->Branch("fMwpc3SNy",fMwpc3SNy,"fMwpc3SNy[nTracks]/D");
  out->Branch("fMwpc3SNz",fMwpc3SNz,"fMwpc3SNz[nTracks]/D");
  out->Branch("fMwpc3PathLocation",fMwpc3PathLocation,"fMwpc3PathLocation[nTracks]/D");
  out->Branch("fMwpc3PathTime",fMwpc3PathTime,"fMwpc3PathTime[nTracks]/D");

  out->Branch("fMwpc4SP",fMwpc4SP,"fMwpc4SP[nTracks]/D");
  out->Branch("fMwpc4SX",fMwpc4SX,"fMwpc4SX[nTracks]/D");
  out->Branch("fMwpc4SY",fMwpc4SY,"fMwpc4SY[nTracks]/D");
  out->Branch("fMwpc4SZ",fMwpc4SZ,"fMwpc4SZ[nTracks]/D");
  out->Branch("fMwpc4SNx",fMwpc4SNx,"fMwpc4SNx[nTracks]/D");
  out->Branch("fMwpc4SNy",fMwpc4SNy,"fMwpc4SNy[nTracks]/D");
  out->Branch("fMwpc4SNz",fMwpc4SNz,"fMwpc4SNz[nTracks]/D");
  out->Branch("fMwpc4PathLocation",fMwpc4PathLocation,"fMwpc4PathLocation[nTracks]/D");
  out->Branch("fMwpc4PathTime",fMwpc4PathTime,"fMwpc4PathTime[nTracks]/D");
  
  out->Branch("fMwpc2SP",fMwpc2SP,"fMwpc2SP[nTracks]/D");
  out->Branch("fMwpc2SX",fMwpc2SX,"fMwpc2SX[nTracks]/D");
  out->Branch("fMwpc2SY",fMwpc2SY,"fMwpc2SY[nTracks]/D");
  out->Branch("fMwpc2SZ",fMwpc2SZ,"fMwpc2SZ[nTracks]/D");
  out->Branch("fMwpc2SNx",fMwpc2SNx,"fMwpc2SNx[nTracks]/D");
  out->Branch("fMwpc2SNy",fMwpc2SNy,"fMwpc2SNy[nTracks]/D");
  out->Branch("fMwpc2SNz",fMwpc2SNz,"fMwpc2SNz[nTracks]/D");
  out->Branch("fMwpc2PathLocation",fMwpc2PathLocation,"fMwpc2PathLocation[nTracks]/D");
  out->Branch("fMwpc2PathTime",fMwpc2PathTime,"fMwpc2PathTime[nTracks]/D");
  
  out->Branch("fMwpc2EnergyLoss",fMwpc2EnergyLoss,"fMwpc2EnergyLoss[nTracks]/D");
  out->Branch("fMwpc3EnergyLoss",fMwpc3EnergyLoss,"fMwpc3EnergyLoss[nTracks]/D");
  out->Branch("fMwpc4EnergyLoss",fMwpc4EnergyLoss,"fMwpc4EnergyLoss[nTracks]/D");
  return 0;
}


Long_t Det_tracking::process_basicDataType(){

      // cout<<  trackArr->vectTrack.size()<<endl;
nTracks=trackArr->vectTrack.size(); 
run=treeRaw->run; 
event=treeRaw->event; 
std::cout<<"Run:  "<<treeRaw->run<<"   event:  "<<treeRaw->event<<std::endl;
for(int j=0; j<nTracks; j++){
      nHits[j]=trackArr->vectTrack[j].getNHits();
      fChi2PerNDF[j] = trackArr->vectTrack[j].getChi2PerNDF();
      fFlagTargetToAC[j] = trackArr->vectTrack[j].getFlagTargetToAC();
      fFlagToFMatch[j] = trackArr->vectTrack[j].getFlagToFMatch();
      fTof1SX[j] = trackArr->vectTrack[j].getTof1SX();
      fTof1SY[j] = trackArr->vectTrack[j].getTof1SY();
      fTof1SZ[j] = trackArr->vectTrack[j].getTof1SZ();
      fTof1SNx[j] = trackArr->vectTrack[j].getTof1SNx();
      fTof1SNy[j] = trackArr->vectTrack[j].getTof1SNy();
      fTof1SNz[j] = trackArr->vectTrack[j].getTof1SNz();
      fTof1SP[j] = trackArr->vectTrack[j].getTof1SP();
      fTof1TimeD[j] = trackArr->vectTrack[j].getTof1TimeD();
      fTof1TimeU[j] = trackArr->vectTrack[j].getTof1TimeU();


      fgapNumTof1[j] = trackArr->vectTrack[j].getTof1GapNum();
      fgapNumTof2[j] = trackArr->vectTrack[j].getTof2GapNum();
      fTof2SX[j] = trackArr->vectTrack[j].getTof2SX();
      fTof2SY[j] = trackArr->vectTrack[j].getTof2SY();
      fTof2SZ[j] = trackArr->vectTrack[j].getTof2SZ();
      fTof2SNx[j] = trackArr->vectTrack[j].getTof2SNx();
      fTof2SNy[j] = trackArr->vectTrack[j].getTof2SNy();
      fTof2SNz[j] = trackArr->vectTrack[j].getTof2SNz();
      fTof2SP[j] = trackArr->vectTrack[j].getTof2SP();
      fTof2TimeI[j] = trackArr->vectTrack[j].getTof2TimeI();
      fTof2TimeO[j] = trackArr->vectTrack[j].getTof2TimeO();

      fVertMX[j] = trackArr->vectTrack[j].getVertMX();
      fVertMY[j] = trackArr->vectTrack[j].getVertMY();
      fVertMPhi[j] = trackArr->vectTrack[j].getVertMPhi();
      
      fSftSY[j] = trackArr->vectTrack[j].getSftSY();
      fSftSZ[j] = trackArr->vectTrack[j].getSftSZ();
      fSftSX[j] = trackArr->vectTrack[j].getSftSX();
     
      fVertSX[j] = trackArr->vectTrack[j].getVertSX();
      fVertSY[j] = trackArr->vectTrack[j].getVertSY();
      fVertSZ[j] = trackArr->vectTrack[j].getVertSZ();
      fVertSP[j] = trackArr->vectTrack[j].getVertSP();
      fVertSPPositron[j] = trackArr->vectTrack[j].getVertSPPositron();
      fVertSPPionPlus[j] = trackArr->vectTrack[j].getVertSPPionPlus();
      fVertSPhi[j] = trackArr->vectTrack[j].getVertSPhi();

      fSftSNy[j] = trackArr->vectTrack[j].getSftSNy();
      fSftSNz[j] = trackArr->vectTrack[j].getSftSNz();
      fSftSNx[j] = trackArr->vectTrack[j].getSftSNx();
    
      fTof1TimeMean[j] = trackArr->vectTrack[j].getTof1TimeMean();
      fTof2TimeMean[j] = trackArr->vectTrack[j].getTof2TimeMean();
      fTof1PathTime[j] = trackArr->vectTrack[j].getTof1PathTime();
      fTof1PathLocation[j] = trackArr->vectTrack[j].getTof1PathLocation();
      fTof2PathTime[j] = trackArr->vectTrack[j].getTof2PathTime();
      fTof2PathLocation[j] = trackArr->vectTrack[j].getTof2PathLocation();
      
      fMwpc4SX[j] = trackArr->vectTrack[j].getMwpc4SX();
      fMwpc4SY[j] = trackArr->vectTrack[j].getMwpc4SY();
      fMwpc4SZ[j] = trackArr->vectTrack[j].getMwpc4SZ();
      fMwpc4SNx[j] = trackArr->vectTrack[j].getMwpc4SNx();
      fMwpc4SNy[j] = trackArr->vectTrack[j].getMwpc4SNy();
      fMwpc4SNz[j] = trackArr->vectTrack[j].getMwpc4SNz();
      fMwpc4SP[j] = trackArr->vectTrack[j].getMwpc4SP();
      fMwpc4PathTime[j] = trackArr->vectTrack[j].getMwpc4PathTime();
      fMwpc4PathLocation[j] = trackArr->vectTrack[j].getMwpc4PathLocation();

      fMwpc2SP[j] = trackArr->vectTrack[j].getMwpc2SP();
      fMwpc2SX[j] = trackArr->vectTrack[j].getMwpc2SX();
      fMwpc2SY[j] = trackArr->vectTrack[j].getMwpc2SY();
      fMwpc2SZ[j] = trackArr->vectTrack[j].getMwpc2SZ();
      fMwpc2SNx[j] = trackArr->vectTrack[j].getMwpc2SNx();
      fMwpc2SNy[j] = trackArr->vectTrack[j].getMwpc2SNy();
      fMwpc2SNz[j] = trackArr->vectTrack[j].getMwpc2SNz();
      fMwpc2PathTime[j] = trackArr->vectTrack[j].getMwpc2PathTime();
      fMwpc2PathLocation[j] = trackArr->vectTrack[j].getMwpc2PathLocation();

      fMwpc3SP[j] = trackArr->vectTrack[j].getMwpc3SP();
      fMwpc3SX[j] = trackArr->vectTrack[j].getMwpc3SX();
      fMwpc3SY[j] = trackArr->vectTrack[j].getMwpc3SY();
      fMwpc3SZ[j] = trackArr->vectTrack[j].getMwpc3SZ();
      fMwpc3SNx[j] = trackArr->vectTrack[j].getMwpc3SNx();
      fMwpc3SNy[j] = trackArr->vectTrack[j].getMwpc3SNy();
      fMwpc3SNz[j] = trackArr->vectTrack[j].getMwpc3SNz();
      fMwpc3PathTime[j] = trackArr->vectTrack[j].getMwpc3PathTime();
      fMwpc3PathLocation[j] = trackArr->vectTrack[j].getMwpc3PathLocation();

      fMwpc2EnergyLoss[j] = trackArr->vectTrack[j].getMwpc2EnergyLoss();
      fMwpc3EnergyLoss[j] = trackArr->vectTrack[j].getMwpc3EnergyLoss();
      fMwpc4EnergyLoss[j] = trackArr->vectTrack[j].getMwpc4EnergyLoss();



}
  
  return 0;
}


Long_t Det_tracking::done_basicDataType(){


  return 0;
}



