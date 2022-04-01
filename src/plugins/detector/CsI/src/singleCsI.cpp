#include <singleCsI.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TList.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TFrame.h>
#include <TSpectrum.h>
Double_t WaveformCsI::waveformSingleSignal(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return 0;
  Double_t termFirst=par[0]/(1+TMath::Exp(-par[2]*(x[0]-par[6])));
  Double_t termSecond=(x[0]-par[1])/(par[3]*par[3]);
  Double_t termThird=TMath::Exp(-(x[0]-par[1])/par[3])+par[4]*TMath::Exp(-(x[0]-par[1])/par[5]);
  Double_t value=termFirst*termSecond*termThird;
  return value;
}
Double_t WaveformCsI::waveformSignal(Double_t *x,Double_t* par){
  Double_t value=0;
  for(UInt_t iWave=0;iWave<mNWave;iWave++){
    value+=waveformSingleSignal(x,&par[iWave*7]);
  }
  return value;
}
Double_t WaveformCsI::waveform(Double_t *x,Double_t* par){
  Double_t value=waveformSignal(x,par);
  UInt_t iPar=mNWave*7+1;
  value+=par[iPar-1];
  return value;
}
Double_t WaveformCsI::waveformCut(Double_t *x,Double_t* par){
  Double_t value=waveform(x,par);
  if(value>1023.0) value=1023.0;
  return value;
}
void SingleCsI::setData(const vector<UShort_t>& listD){
  for(UInt_t iS=0,nS=listD.size();iS<nS;iS++){
    mListData.push_back(listD[iS]);
  }
}
static Double_t waveformPure(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return 0;
  Double_t termFirst=par[0]/(1+TMath::Exp(-par[2]*(x[0]-par[6])));
  Double_t termSecond=(x[0]-par[1])/(par[3]*par[3]);
  Double_t termThird=TMath::Exp(-(x[0]-par[1])/par[3])+par[4]*TMath::Exp(-(x[0]-par[1])/par[5]);
  Double_t value=termFirst*termSecond*termThird;
  return value;
}

static Double_t derivativeSingle(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return 0;
  Double_t firstExp=par[4]*TMath::Exp(-(x[0]-par[1])/par[5]);
  Double_t secondExp=TMath::Exp(-(x[0]-par[1])/par[3]);
  Double_t thirdExp=TMath::Exp(par[2]*par[6]-par[2]*x[0]);
  Double_t termFirst=par[0]*(x[0]-par[1])*(-firstExp/par[5]-secondExp/par[3])/(par[3]*par[3]*(thirdExp+1));
  Double_t termSecond=par[0]*(firstExp+secondExp)/(par[3]*par[3]*(thirdExp+1));
  Double_t termThird=par[0]*par[2]*(x[0]-par[1])*(firstExp+secondExp)*thirdExp/(par[3]*par[3]*(thirdExp+1)*(thirdExp+1));
  Double_t value=termFirst+termSecond+termThird;
  return value;
}
static Double_t findMaximumX(TF1* f1,Double_t lowX,Double_t highX){

  Double_t range=highX-lowX;
  Double_t interval=sqrt(range);
  int NPx=range/interval;


  Double_t valueMax=f1->Eval(lowX),xMax=lowX;

  for(int iX=1;iX<=NPx;iX++){
    Double_t x=lowX+interval*iX;
    Double_t value=f1->Eval(x);

    if(value>valueMax){
      valueMax=value;
      xMax=x;
    }
    else
      break;
  }
  interval=1;
  Double_t valueTry=f1->Eval(xMax-interval);
  if(valueTry>valueMax){
    interval=-1;
    xMax+=interval;
    valueMax=valueTry;
  }
  while(true){
    Double_t x=xMax+interval;
    valueTry=f1->Eval(x);

    if(valueTry>valueMax){
      valueMax=valueTry;
      xMax=x;
    }
    else{
      return (x+xMax)/2.0;
    }
  }

}
static Double_t findPeak(Double_t *par){
  if(waveformPure(&par[6],par)<0) return 0.0;
  TF1* fPure=new TF1("pure",waveformPure,7,1,500);
  fPure->SetParameters(par);
  double xPeak=findMaximumX(fPure,par[1],500);
  delete fPure;
  return xPeak;
}
static Double_t findPeak(Double_t *par,Double_t low,Double_t high){
  if(waveformPure(&par[6],par)<0) return 0.0;
  TF1* fPure=new TF1("pure",waveformPure,7,low,high);
  fPure->SetParameters(par);
  double xPeak=findMaximumX(fPure,low,high);
  delete fPure;
  return xPeak;
}


static Double_t findTime(Double_t *par){
  if(waveformPure(&par[6],par)<0) return 0.0;
  TF1* fDeriv=new TF1("deriv",derivativeSingle,7,1,250);
  fDeriv->SetParameters(par);
  double xTime=findMaximumX(fDeriv,par[1],250);
  delete fDeriv;
  return xTime;
}
static Double_t findEnergy(Double_t *par){
  if(waveformPure(&par[6],par)<0) return 0.0;
  Double_t xPeak=findPeak(par);
  if(xPeak>249.9999||xPeak<0.00001) return 0;
  Double_t value=waveformPure(&xPeak,par);
  return value;
}
static Double_t quickFindTime(TH1D* h){
  UInt_t iBin=h->GetMaximumBin();
  UInt_t step=sqrt(iBin);
  UInt_t nStep=iBin/step+1;
  Double_t valueZero=h->GetBinContent(1);

  nStep*=step;
  if(nStep>250) nStep=250;

  for(UInt_t iStep=nStep;iStep>0;iStep-=step){
    Double_t myValue=h->GetBinContent(iStep);

    if(myValue<valueZero+5){

      for(int iMini=iStep,nMini=iStep+step;iMini<nMini;iMini++){

	if(h->GetBinContent(iMini)>valueZero+5){
	  return h->GetBinCenter(iMini);
	}
      }
      return h->GetBinCenter(iStep+step);
    }
  }

  for(int iMini=1,nMini=step;iMini<nMini;iMini++){

    if(h->GetBinContent(iMini)>valueZero+5){
      return h->GetBinCenter(iMini);
    }
  }
  return 0;
}
char* SingleCsI::nameCsI(UInt_t index){
  if(index<1){
    sprintf(mName,"0000");
  }
  UInt_t iClock=(index-1)/(2*2*16)+1;
  UInt_t iInside=index-1-(iClock-1)*2*2*16;
  UInt_t iFB=iInside/(2*16);
  UInt_t iUD=(iInside-iFB*2*16)/16;
  UInt_t iPaddle=iInside%16;
  char NameFB[3]="fb";
  char NameUD[3]="ud";
  sprintf(mName,"%i%c%c%i",iClock,NameFB[iFB],NameUD[iUD],iPaddle+1);
  return mName;
}
Double_t formReference(Double_t*x, Double_t*par){
  if(x[0]<=par[0]) return par[9];
  if(x[0]-par[0]<=5) return par[9]+(x[0]-par[0])*(par[10]-par[9])/5.0;
  if(x[0]<=par[1]) return par[10];
  if(x[0]-par[1]<=5) return par[10]-(x[0]-par[1])*(par[10]-par[9])/5.0;
  if(x[0]<=par[2]) return par[9];
  if(x[0]-par[2]<=5) return par[9]+(x[0]-par[2])*(par[11]-par[9])/5.0;
  if(x[0]<=par[3]) return par[11];
  if(x[0]-par[3]<=5) return par[11]+(x[0]-par[3])*(par[10]-par[11])/5.0;
  if(x[0]<=par[4]) return par[10];
  if(x[0]-par[4]<=5) return par[10]-(x[0]-par[4])*(par[10]-par[12])/5.0;
  if(x[0]<=par[5]) return par[12];
  if(x[0]-par[5]<=5) return par[12]+(x[0]-par[5])*(par[13]-par[12])/5.0;
  if(x[0]<=par[6]) return par[13];
  if(x[0]-par[6]<=5) return par[13]+(x[0]-par[6])*(par[9]-par[13])/5.0;
  if(x[0]<=par[7]) return par[9];
  if(x[0]-par[7]<=5) return par[9]+(x[0]-par[7])*(par[10]-par[9])/5.0;
  if(x[0]<=par[8]) return par[10];
  return par[9];

}
Double_t tryReferenceForm(TH1D* h1){
  TF1* f1=new TF1("waveReference",formReference,1,250,14);
  double par[]={68,82,94,107,132,144,157,169,182,112,876,220,440,650};
  f1->SetParameters(par);
  h1->Fit(f1,"QR");
  UInt_t iLow=1;
  UInt_t iUp=h1->GetNbinsX();
  TF1* funFit=h1->GetFunction("waveReference");
  Double_t chi2=0;
  for(UInt_t i=iLow;i<=iUp;i++){
    Double_t x=h1->GetBinCenter(i);
    Double_t y=h1->GetBinContent(i);
    Double_t error=h1->GetBinError(i);
    Double_t value=funFit->Eval(x);
    Double_t thisChi2=(value-y)*(value-y)/error/error;
    chi2+=thisChi2;
  }
  chi2/=(iUp-1-14);
  double max=h1->GetBinContent(h1->GetMaximumBin());
  delete f1;
  delete funFit;
  return chi2/max/max*400;
}
void SingleCsI::quickDrawWaves(TH1D* h1){
  char name[256];
  sprintf(name,"canvas_run%d_%dCsI_%s",mRunNo,mEventNo,nameCsI(mIndexCsI));
  TCanvas* c1=new TCanvas(name,name,1200,900);
  if(h1!=0){
    h1->SetMarkerStyle(2);
    h1->SetMarkerSize(1.2);
    h1->SetMinimum(0);
    h1->DrawCopy();
    Double_t xTime=mListTime[0];
    TLine* lineTime=new TLine(xTime,0,xTime,h1->GetBinContent(h1->FindBin(xTime)));
    lineTime->SetLineColor(kRed);
    Double_t peak=mListPeak[0];
    TLine* linePeak=new TLine(peak,0,peak,mListEnergy[0]+mBaseline);
    linePeak->SetLineColor(kBlue);
    lineTime->Draw("same");
    linePeak->Draw("same");
  }
  TLine* lineBase=new TLine(1,mBaseline,250,mBaseline);
  lineBase->Draw("same");
  c1->Write();
  //  c1->Delete();
  //  delete c1;
  //  sprintf(name,"wave_run%d_%dCsI_%s.png",mRunNo,mEventNo,nameCsI(mIndexCsI));
  //  c1->SaveAs(name);
}
static Double_t quickFindEnergy(Double_t width,Double_t height){
  Double_t RatioWidth[88]={0,0.000532561,0.000662612,0.000785248,0.000881907,0.000974722,0.00107656,0.00116987,0.00126172,0.00134836,0.00144428,0.0015228,0.00160406,0.0016911,0.00178789,0.00188698,0.00198058,0.00205644,0.00214122,0.00223334,0.00232779,0.00243305,0.00254952,0.00265186,0.00274375,0.00285149,0.00296217,0.00307013,0.00318701,0.0032979,0.00343094,0.00355336,0.00367467,0.00380957,0.00394861,0.0040957,0.00423333,0.00438593,0.00453805,0.00470083,0.00486912,0.00505448,0.00523489,0.00541428,0.00560778,0.00580845,0.00601038,0.00623907,0.00647674,0.00672392,0.00695933,0.00724047,0.00750577,0.00782522,0.00810785,0.0084258,0.00875854,0.00911275,0.00949442,0.00987883,0.0103332,0.0107683,0.0112071,0.0117181,0.0122396,0.0128687,0.0134567,0.0141223,0.0147963,0.015518,0.0163686,0.0172566,0.0181466,0.0192216,0.0203826,0.0215858,0.0229179,0.0243755,0.0260653,0.0279032,0.0299566,0.0322744,0.0348206,0.0378138,0.0412497,0.0451933,0.049826,0.055305};
  Double_t RatioHeight[88]={1,1.02041,1.03093,1.04167,1.05263,1.06383,1.07527,1.08696,1.0989,1.11111,1.1236,1.13636,1.14943,1.16279,1.17647,1.19048,1.20482,1.21951,1.23457,1.25,1.26582,1.28205,1.2987,1.31579,1.33333,1.35135,1.36986,1.38889,1.40845,1.42857,1.44928,1.47059,1.49254,1.51515,1.53846,1.5625,1.5873,1.6129,1.63934,1.66667,1.69492,1.72414,1.75439,1.78571,1.81818,1.85185,1.88679,1.92308,1.96078,2,2.04082,2.08333,2.12766,2.17391,2.22222,2.27273,2.32558,2.38095,2.43902,2.5,2.5641,2.63158,2.7027,2.77778,2.85714,2.94118,3.0303,3.125,3.22581,3.33333,3.44828,3.57143,3.7037,3.84615,4,4.16667,4.34783,4.54545,4.7619,5,5.26316,5.55556,5.88235,6.25,6.66667,7.14286,7.69231,8.33333};
  Double_t myRatio=width/height;
  Double_t lowRatioWidth,highRatioWidth;
  Double_t lowRatioHeight,highRatioHeight;
  if(myRatio>RatioWidth[87]){
    lowRatioWidth=RatioWidth[86];
    highRatioWidth=RatioWidth[87];
    lowRatioHeight=RatioHeight[86];
    highRatioHeight=RatioHeight[87];
  }
  else{
    UInt_t indexTry=UInt_t(myRatio/RatioWidth[87]*87);
    if(RatioWidth[indexTry]<myRatio){
      while(RatioWidth[indexTry]<myRatio){
	indexTry++;
      }
    }
    else{
      while(RatioWidth[indexTry]>myRatio){
	indexTry--;
      }
      indexTry+=1;//the first larger one
    }
    lowRatioWidth=RatioWidth[indexTry-1];
    lowRatioHeight=RatioHeight[indexTry-1];
    highRatioWidth=RatioWidth[indexTry];
    highRatioHeight=RatioHeight[indexTry];
  }
  return (lowRatioHeight+(myRatio-lowRatioWidth)/(highRatioWidth-lowRatioWidth)*(highRatioHeight-lowRatioHeight))*height;


}
bool SingleCsI::quick_fit(bool draw){
  static UInt_t countCsI[768];

  static UInt_t nEvent=0;
  char name[256];
  sprintf(name,"wave_run%d_%dCsI_%s",mRunNo,mEventNo,nameCsI(mIndexCsI));
  char title[256];
  sprintf(title,"fAdc waveform");
  TH1D* h1=new TH1D(name,title,250,0.5,250.5);
  for(UInt_t iData=0,nData=mListData.size();iData<nData;iData++){
    h1->SetBinContent(iData+1,mListData[iData]);
    h1->SetBinError(iData+1,1.1);
  }
  mNWave=1;
  UInt_t indexMax=h1->GetMaximumBin();
  Double_t valueMax=h1->GetBinContent(indexMax);
  Double_t xMax=h1->GetBinCenter(indexMax);
  if(valueMax<1023){
    mBaseline=h1->GetBinContent(1);
    mListPeak.push_back(xMax);
    mListEnergy.push_back(valueMax-mBaseline);
    mListTime.push_back(quickFindTime(h1));
  }
  else{
    mBaseline=h1->GetBinContent(1);
    UInt_t iLow,iHigh;
    for(iLow=indexMax-1;iLow>0;iLow--){
      if(h1->GetBinContent(iLow)<1023) break;
    }
    UInt_t nBin=h1->GetNbinsX();
    for(iHigh=indexMax+1;iHigh<=nBin;iHigh++){
      if(h1->GetBinContent(iHigh)<1023) break;
    }
    mListPeak.push_back((iLow+iHigh)*0.5);
    mListEnergy.push_back(quickFindEnergy(iHigh-iLow+1,1023-mBaseline));
    mListTime.push_back(quickFindTime(h1));
  }
  if(draw && countCsI[mIndexCsI-1]<5){
    quickDrawWaves(h1);
    countCsI[mIndexCsI-1]++;
  }
  nEvent++;
  delete(h1);
  return true;
}
bool SingleCsI::fit(){
  static UInt_t nEvent=0;
  char name[256];
  sprintf(name,"wave_run%d_%dCsI_%s",mRunNo,mEventNo,nameCsI(mIndexCsI));
  char title[256];
  sprintf(title,"fAdc waveform");
  TH1D* h1=new TH1D(name,title,250,0.5,250.5);
  for(UInt_t iData=0,nData=mListData.size();iData<nData;iData++){
    h1->SetBinContent(iData+1,mListData[iData]);
    h1->SetBinError(iData+1,1.1);
  }
  mNWave=1;
  tryFit(h1);
  Double_t chi2Reduced=findChi2(h1);
  if(chi2Reduced>80){
    mNWave=2;
    tryFit(h1);
    Double_t chi2Reduced=findChi2(h1);
  }
  calculateTdcAdc(h1);
  mBaseline=h1->GetFunction("waveCut")->GetParameter(mNWave*7);

  if(nEvent<20)
    drawWaves(h1);

  nEvent++;
  delete(h1);
  return true;
}
void SingleCsI::findLocalMax(TH1D* h1){
  TF1* f1=h1->GetFunction("waveCut");
  int nBin=h1->GetNbinsX();
  double maxValue=0;
  double newPeak;
  for(int iBin=1;iBin<=nBin;iBin++){
    float vBin=h1->GetBinContent(iBin);
    float value=f1->Eval(h1->GetBinCenter(iBin));
    if(vBin-value>maxValue){
      maxValue=vBin-value;
      newPeak=h1->GetBinCenter(iBin);
    }
  }
  mListLocalMax.push_back(newPeak);
}
void SingleCsI::tryFit(TH1D* h1){
  const UInt_t NPar=mNWave*7+1;
  WaveformCsI* myWave=new WaveformCsI(mNWave);
  TF1* f1=new TF1("waveCut",myWave,&WaveformCsI::waveformCut,1,250,NPar,"WaveformCsI","waveformCut");
  Double_t par[NPar];
  if(mNWave==1){
    f1->SetParameters(5035,26.69,0.1782,12.280,32.91,12.28,h1->GetBinCenter(h1->GetMaximumBin()),72.79);    
    f1->SetLineStyle(6);
    f1->SetLineColor(1);
    f1->SetLineWidth(3);
    h1->Fit(f1,"Q");
    for(UInt_t i=0;i<NPar;i++){
      mPar[i]=h1->GetFunction("waveCut")->GetParameter(i);
    }
  }
  else{
    findLocalMax(h1);
    f1->SetParameter(0,5035);
    f1->SetParameter(1,26.69);
    f1->SetParameter(2,0.1782);
    f1->SetParameter(3,12.280);
    f1->SetParameter(4,32.91);
    f1->SetParameter(5,12.28);
    f1->SetParameter(6,mListLocalMax[0]);    
    
    for(int i=0,n=h1->GetFunction("waveCut")->GetNpar();i<n;i++){
      f1->SetParameter(7+i,h1->GetFunction("waveCut")->GetParameter(i));
    }

    f1->SetLineStyle(6);
    f1->SetLineColor(1);
    f1->SetLineWidth(3);
    h1->Fit(f1,"Q");
    for(UInt_t i=0;i<NPar;i++){
      mPar[i]=h1->GetFunction("waveCut")->GetParameter(i);
    }
  }
  delete f1;

}
void SingleCsI::quickCalculateTdcAdc(TH1D* h1){
  Double_t par[8];
  for(UInt_t iPar=0;iPar<8;iPar++){
    par[iPar]=h1->GetFunction("waveCut")->GetParameter(iPar);
  }
  UInt_t iMax=h1->GetMaximumBin();
  UInt_t iLow,iHigh;
  for(iLow=iMax-1;iLow>0;iLow--){
    if(h1->GetBinContent(iLow)<1023) break;
  }
  UInt_t nBin=h1->GetNbinsX();
  for(iHigh=iMax+1;iHigh<=nBin;iHigh++){
    if(h1->GetBinContent(iHigh)<1023) break;
  }

  Double_t peak=findPeak(par,iLow,iHigh);
  mListPeak.push_back(peak);

  Double_t energy=waveformPure(&peak,par);;
  mListEnergy.push_back(energy);
  Double_t time=findTime(par);
  mListTime.push_back(time);

}

void SingleCsI::calculateTdcAdc(TH1D* h1){
  Double_t par[mNWave*7+1];
  par[mNWave*7]=h1->GetFunction("waveCut")->GetParameter(mNWave*7);
  for(UInt_t iPar=0;iPar<7*mNWave;iPar++){
    par[iPar]=h1->GetFunction("waveCut")->GetParameter(iPar);
  }

  for(UInt_t i=0;i<mNWave;i++){
    Double_t peak=findPeak(&par[i*7]);
    mListPeak.push_back(peak);
    Double_t energy=waveformPure(&peak,par);;
    mListEnergy.push_back(energy);

    Double_t time=findTime(&par[i*7]);
    mListTime.push_back(time);

  }

}
void SingleCsI::drawWaves(TH1D* h1){
  char name[256];
  sprintf(name,"canvas_run%d_%dCsI_%s",mRunNo,mEventNo,nameCsI(mIndexCsI));
  TCanvas* c1=new TCanvas(name,name,1200,900);
  if(h1!=0){
    h1->SetMarkerStyle(2);
    h1->SetMarkerSize(1.2);
    h1->SetMinimum(0);
    h1->DrawCopy();

    WaveformCsI* myWave=new WaveformCsI(mNWave);
    TF1* f2=new TF1("waveform",myWave,&WaveformCsI::waveform,1,250,mNWave*7+1,"WaveformCsI","waveform");
    Double_t par[mNWave*7+1];
    par[mNWave*7]=h1->GetFunction("waveCut")->GetParameter(mNWave*7);
    for(UInt_t iPar=0;iPar<7*mNWave;iPar++){
      par[iPar]=h1->GetFunction("waveCut")->GetParameter(iPar);
    }

    for(UInt_t i=0;i<mNWave;i++){
      char name[128];
      sprintf(name,"wave_%i",i+1);
      TF1* fWavelet=new TF1(name,waveformPure,1,250,7);
      for(int iPar=0;iPar<7;iPar++){
	f2->FixParameter(i*7+iPar,par[i*7+iPar]);
	fWavelet->FixParameter(iPar,par[i*7+iPar]);
      }
      Double_t xTime=mListTime[i];
      TLine* lineTime=new TLine(xTime,0,xTime,h1->GetBinContent(h1->FindBin(xTime)));
      lineTime->SetLineColor(kRed);
      Double_t peak=mListPeak[i];
      TLine* linePeak=new TLine(peak,0,peak,mListEnergy[i]);
      linePeak->SetLineColor(kBlue);
      fWavelet->DrawCopy("same");
      lineTime->Draw("same");
      linePeak->Draw("same");

    }
    f2->FixParameter(mNWave*7,par[mNWave*7]);
    f2->SetLineStyle(9);
    f2->SetLineColor(kBlue);
    h1->GetListOfFunctions()->Add(f2);
  }

  c1->Write();
  //  c1->Delete();
  //  delete c1;
  //  sprintf(name,"wave_run%d_%dCsI_%s.png",mRunNo,mEventNo,nameCsI(mIndexCsI));
  //  c1->SaveAs(name);
}
Double_t SingleCsI::findChi2(TH1D* h1){
  UInt_t iLow=1;
  UInt_t iUp=h1->GetNbinsX();
  TF1* funFit=h1->GetFunction("waveCut");
  Double_t chi2=0;
  for(UInt_t i=iLow;i<=iUp;i++){
    Double_t x=h1->GetBinCenter(i);
    Double_t y=h1->GetBinContent(i);
    Double_t error=h1->GetBinError(i);
    Double_t value=funFit->Eval(x);
    Double_t thisChi2=(value-y)*(value-y)/error/error;
    chi2+=thisChi2;

  }
  chi2/=(iUp-1-h1->GetFunction("waveCut")->GetNpar());
  return chi2;
}
