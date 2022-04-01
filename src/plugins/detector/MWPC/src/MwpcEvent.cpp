#include "MwpcEvent.h"
#include <fstream>
#include <string>
#include "TF1.h"
#include <vector>
#include "TCanvas.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TLatex.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include "C2X_Strip_Trans.h"
#include "RPArray.h"
using namespace std;

void MwpcEvent::setHist(){
 ostringstream name,title;
 name<<eventName<<"_Event"<<eventNumber;
 title<<eventName<<"_Event"<<eventNumber<<";Strip Channel;ADC";
 eventHist=new TH1D(name.str().c_str(),title.str().c_str(),numBin,lowBound,upBound);
 name.str("");
 title.str("");
}

void MwpcEvent::addHistData(int* value){
  for(int i=0;i<numBin;i++){
    if(stripTransStat==1)
      eventHist->SetBinContent(i+1,value[i]);
    else
      eventHist->SetBinContent(i+1,value[C2XStTrans[i]]);
  }
}

void MwpcEvent::addHistData_pedSub(int* value, const double* pdd){
  for(int i=0;i<numBin;i++){
    if(stripTransStat==1){
      if(value[i]-pdd[i]>0)
        eventHist->SetBinContent(i+1,value[i]-pdd[i]);
      else
        eventHist->SetBinContent(i+1,0);
    }
    else{
      if(value[C2XStTrans[i]]-pdd[i]>0)
        eventHist->SetBinContent(i+1,value[C2XStTrans[i]]-pdd[i]);
      else
        eventHist->SetBinContent(i+1,0);
    }
  }
}

void MwpcEvent::addHistData_pedSub_sigmaCut(int* value, const double* pdd, const double *sigma, double nSigma){
  for(int i=0;i<numBin;i++){
    if(stripTransStat==1){
      if(value[i]-pdd[i]>nSigma*sigma[i])
        eventHist->SetBinContent(i+1,value[i]-pdd[i]);
      else
        eventHist->SetBinContent(i+1,0);
    }
    else{
      if(value[C2XStTrans[i]]-pdd[i]>nSigma*sigma[i])
        eventHist->SetBinContent(i+1,value[C2XStTrans[i]]-pdd[i]);
      else
        eventHist->SetBinContent(i+1,0);
    }
  }
}

void MwpcEvent::drawHist(){
  eventHist->Draw();
}

void MwpcEvent::saveHist(){
  ostringstream name;
  TCanvas *canv;
  name<<eventName<<"_Event"<<eventNumber;
  canv=new TCanvas(name.str().c_str(),name.str().c_str(),1000,300);
  eventHist->Draw();
  name.str("");
  name<<eventName<<"_Event"<<eventNumber<<".png";
  canv->SaveAs(name.str().c_str());
  name.str("");
}

double* MwpcEvent::getAdcValue(int* value, const double* pdd){
  double *adcValue=new double[numBin];
  for(int i=0;i<numBin;i++){
    if(stripTransStat==1){
      if(value[i]-pdd[i]>0)
        adcValue[i]=value[i]-pdd[i];
      else
        adcValue[i]=0;
    }
    else{
      if(value[C2XStTrans[i]]-pdd[i]>0)
        adcValue[i]=value[C2XStTrans[i]]-pdd[i];
      else
        adcValue[i]=0;
    }
  }
  return adcValue;
}

double* MwpcEvent::getAdcValue(int* value, const double* pdd, const double *sigma, double nSigma){
  double *adcValue=new double[numBin];
  for(int i=0;i<numBin;i++){
    if(stripTransStat==1){
      if(value[i]-pdd[i]>nSigma*sigma[i])
        adcValue[i]=value[i]-pdd[i];
      else
        adcValue[i]=0;
    }
    else{
      if(value[C2XStTrans[i]]-pdd[i]>nSigma*sigma[i])
        adcValue[i]=value[C2XStTrans[i]]-pdd[i];
      else
        adcValue[i]=0;
    }
  }
  return adcValue;
}

int MwpcEvent::getNHit(int* value, const double* pdd, const double *sigma, double nSigma){
  int nhit=0;
  for(int i=0;i<numBin;i++){
    if(stripTransStat==1){
      if(value[i]-pdd[i]>nSigma*sigma[i])
        nhit++;
    }
    else{
      if(value[C2XStTrans[i]]-pdd[i]>nSigma*sigma[i])
        nhit++;
    }
  }
  return nhit;
}

int* MwpcEvent::multiplet3Type(int* value, const double* pdd, const double* sigma, double nSigma){
  int hitStat[numBin+1];
  
  for(int i=0;i<numBin;i++){
    if(stripTransStat==1){
      if(value[i]-pdd[i]>nSigma*sigma[i])
        hitStat[i]=1;
      else
        hitStat[i]=0;
    }
    else{
      if(value[C2XStTrans[i]]-pdd[i]>nSigma*sigma[i])
        hitStat[i]=1;
      else
        hitStat[i]=0;
    }
  }
  hitStat[numBin]=0;

  int temp=0,stat=0;
  vector<int> vect_temp;
  for(int i=0;i<numBin+1;i++){
    if(hitStat[i]==0){
      if(stat==1) vect_temp.push_back(temp);
      temp=0;
      stat=0;
    }
    else{
      temp++;
      stat=1;
    }
  }

  int *mult=new int[3]{0,0,0};
  for(int i=0;i<(int)vect_temp.size();i++){
    if(vect_temp[i]==1) mult[0]++;
    else if(vect_temp[i]==2) mult[1]++;
    else mult[2]++;
  }
  return mult;
}

int* MwpcEvent::multiplet4Type(int* value, const double* pdd, const double* sigma, double nSigma){
  int hitStat[numBin+1];

  for(int i=0;i<numBin;i++){
    if(stripTransStat==1){
      if(value[i]-pdd[i]>nSigma*sigma[i])
        hitStat[i]=1;
      else
        hitStat[i]=0;
    }
    else{
      if(value[C2XStTrans[i]]-pdd[i]>nSigma*sigma[i])
        hitStat[i]=1;
      else
        hitStat[i]=0;
    }
  }
  hitStat[numBin]=0;

  int temp=0,stat=0;
  vector<int> vect_temp;
  for(int i=0;i<numBin+1;i++){
    if(hitStat[i]==0){
      if(stat==1) vect_temp.push_back(temp);
      temp=0;
      stat=0;
    }
    else{
      temp++;
      stat=1;
    }
  }

  int *mult=new int[4]{0,0,0,0};
  for(int i=0;i<(int)vect_temp.size();i++){
    if(vect_temp[i]==1) mult[0]++;
    else if(vect_temp[i]==2) mult[1]++;
    else if(vect_temp[i]==3) mult[2]++;
    else mult[3]++;
  }
  return mult;
}

int* MwpcEvent::multipletNType(int* value, const double* pdd, const double* sigma, double nSigma,int nType){
  int hitStat[numBin+1];

  for(int i=0;i<numBin;i++){
    if(stripTransStat==1){
      if(value[i]-pdd[i]>nSigma*sigma[i])
        hitStat[i]=1;
      else
        hitStat[i]=0;
    }
    else{
      if(value[C2XStTrans[i]]-pdd[i]>nSigma*sigma[i])
        hitStat[i]=1;
      else
        hitStat[i]=0;
    }
  }
  hitStat[numBin]=0;

  int temp=0,stat=0;
  vector<int> vect_temp;
  for(int i=0;i<numBin+1;i++){
    if(hitStat[i]==0){
      if(stat==1) vect_temp.push_back(temp);
      temp=0;
      stat=0;
    }
    else{
      temp++;
      stat=1;
    }
  }

  int *mult=new int[nType];
  for(int j=0;j<nType;j++){
    mult[j]=0;
  }

  for(int i=0;i<(int)vect_temp.size();i++){
    for(int j=0;j<nType;j++){
      // if(vect_temp[i]<nType && vect_temp[i]==j+1) mult[j]++;
      // else if(vect_temp[i]>=nType) mult[nType-1]++;
      if(vect_temp[i]==j+1) mult[j]++;
    }
  }

  return mult;
}

vector<vector<int> > MwpcEvent::clusterStripNumber(int* value, const double* pdd, const double* sigma, double nSigma){
  int hitStat[numBin+1];
  for(int i=0;i<numBin;i++){
    if(stripTransStat==1){
      if(value[i]-pdd[i]>nSigma*sigma[i])
        hitStat[i]=1;
      else
        hitStat[i]=0;
    }
    else{
      if(value[C2XStTrans[i]]-pdd[i]>nSigma*sigma[i])
        hitStat[i]=1;
      else
        hitStat[i]=0;
    }
  }
  hitStat[numBin]=0;

  int stat=0;
  vector<int> vect_stNum;
  vector<vector<int> > vect_clus;
  for(int i=0;i<numBin+1;i++){
    if(hitStat[i]==0){
      if(stat==1) vect_clus.push_back(vect_stNum);
      vect_stNum.clear();
      stat=0;
    }
    else{
      vect_stNum.push_back(i);
      stat=1;
    }
  }
 return vect_clus;
}

int MwpcEvent::getNCluster(MwpcEvent* event,int* value, const double* pdd, const double* sigma,double nSigma){
  vector<vector<int> > clusterStNum=event->clusterStripNumber(value, pdd, sigma, nSigma);
  return clusterStNum.size();
} 

vector<int> MwpcEvent::getClusterNStrip(MwpcEvent* event,int* value, const double* pdd, const double* sigma,double nSigma){
  vector<vector<int> > clusterStNum=event->clusterStripNumber(value, pdd, sigma, nSigma);
  vector<int> nStrip;
  for(int i=0;i<(int)clusterStNum.size();i++){
    nStrip.push_back(clusterStNum[i].size());
  }
  return nStrip;
}

double* MwpcEvent::getClusterTotCharge(MwpcEvent* event,int* value, const double* pdd, const double* sigma, double nSigma){
  double* AdcValue=event->getAdcValue(value, pdd, sigma, nSigma);
  vector<vector<int> > clusterStNum=event->clusterStripNumber(value, pdd, sigma, nSigma);
  double *totCharge=new double[clusterStNum.size()];
  double temp;
  for(int i=0;i<(int)clusterStNum.size();i++){
    temp=0;
    for(int j=0;j<(int)clusterStNum[i].size();j++){
      temp+=AdcValue[clusterStNum[i][j]];
    }
    totCharge[i]=temp;
  }

  return totCharge;
}

//Logic of this function: Firstly, set status for each strip after pedestal subtraction; Secondly, extract clusters and save strip numbers of each cluster into a vector -- Nest vector is used here, the first layer is to save strip numbers of each cluster, and the second layer is to save all clusters; Thirdly, thrshold cut for singlet and doublet is applied to filter clusters so as to obtain a first-filtered vector; Finally, total charge cut is applied to further filter clusters so as to obtain a second-filtered vector.
vector<vector<int> > MwpcEvent::clusterStripNumber(MwpcEvent* event,int* value, const double* pdd, const double* sigma, bool statusTotalChargeCut, const double* totalChargeCut, bool statusSingletThreshold, double singletThreshold, bool statusDoubletThresholdTotalChargeCut, double doubletThreshold, double doubletTotalChargeCut){
  // Set status for each strip after pedestal subtraction
  int hitStat[numBin+1];
  for(int i=0;i<numBin;i++){
    if(stripTransStat==1){
      if(value[i]-pdd[i]>0)
        hitStat[i]=1;
      else
        hitStat[i]=0;
    }
    else{
      if(value[C2XStTrans[i]]-pdd[i]>0)
        hitStat[i]=1;
      else
        hitStat[i]=0;
    }
  }
  hitStat[numBin]=0;

  // Extract clusters and save strip numbers of each cluster into a vector -- Nest vector is used here, the first layer is to save strip numbers of each cluster, and the second layer is to save all clusters; Thirdly, thrshold cut for singlet and doublet is applied to filter clusters so as to obtain a first-filtered vector
  int stat=0;
  vector<int> vect_stNum;
  vector<vector<int> > vect_clus;
  vect_stNum.clear();
  vect_clus.clear();
  for(int i=0;i<numBin+1;i++){
    if(hitStat[i]==0){
      if(stat==1) vect_clus.push_back(vect_stNum);
      vect_stNum.clear();
      stat=0;
    }
    else{
      vect_stNum.push_back(i);
      stat=1;
    }
  }

 // Thrshold cut for singlet and doublet is applied to filter clusters so as to obtain a first-filtered vector 
  double* AdcValue=event->getAdcValue(value, pdd, sigma, 0); // Get Adc values after pedestal subtraction
  vector<vector<int> > vect_clus_after_thresholdcut;
  vect_clus_after_thresholdcut.clear();
  for(int i=0;i<(int)vect_clus.size();i++){
    //Singlet threshold cut
    if(vect_clus[i].size()==1){
      if(statusSingletThreshold==1){ // Status == 1 -- Apply threshold cut for singlet; Status == 0 -- Dont apply
        if(AdcValue[vect_clus[i][0]]>sigma[vect_clus[i][0]]*singletThreshold){ 
          vect_clus_after_thresholdcut.push_back(vect_clus[i]); // Keep if Adc values is large than threshold
          continue;
        }
        else continue;
      }
      else{
        vect_clus_after_thresholdcut.push_back(vect_clus[i]);
        continue;
      }
    }
    
    // Doublet threshold cut
    if(vect_clus[i].size()==2){
      if(statusDoubletThresholdTotalChargeCut==1){ // Status == 1 -- Apply threshold-total-charge cut for doublet; Status == 0 -- Dont apply
        if((AdcValue[vect_clus[i][0]]<sigma[vect_clus[i][0]]*doubletThreshold || AdcValue[vect_clus[i][1]]<sigma[vect_clus[i][1]]*doubletThreshold) && AdcValue[vect_clus[i][0]]+AdcValue[vect_clus[i][1]]<doubletTotalChargeCut) continue; // Cut off if one of two Adc values is less than threshold while total charge is less than total-charge cut value 
        else{
          vect_clus_after_thresholdcut.push_back(vect_clus[i]);
          continue;
        }
      }
      else{
        vect_clus_after_thresholdcut.push_back(vect_clus[i]);
        continue;
      }
    }

    // Other multiplets
    else{
      vect_clus_after_thresholdcut.push_back(vect_clus[i]);
      continue;    
    }
  }

  //Total charge cut is applied to further filter clusters so as to obtain a second-filtered vector. 
  double totalcharge;
  vector<vector<int> > vect_clus_after_thresholdcut_after_totalchargecut;
  vect_clus_after_thresholdcut_after_totalchargecut.clear();
  for(int i=0;i<(int)vect_clus_after_thresholdcut.size();i++){
    if(statusTotalChargeCut==1){ // Status == 1 -- Apply total charge cut; Status == 0 -- Dont apply
      totalcharge=0;
      for(int j=0;j<(int)vect_clus_after_thresholdcut[i].size();j++){
        totalcharge+=AdcValue[vect_clus_after_thresholdcut[i][j]];
      }
      if(totalcharge>totalChargeCut[vect_clus_after_thresholdcut[i].size()]){// Keep if total charge is larger than total-charge cut value
        vect_clus_after_thresholdcut_after_totalchargecut.push_back(vect_clus_after_thresholdcut[i]); 
        continue;
      }
      else continue;
    }
    else{
      vect_clus_after_thresholdcut_after_totalchargecut.push_back(vect_clus_after_thresholdcut[i]);
      continue;
    }
  }
 return vect_clus_after_thresholdcut_after_totalchargecut;
}
      
int MwpcEvent::getNCluster(MwpcEvent* event, int* value, const double* pdd, const double* sigma, bool statusTotalChargeCut, const double* totalChargeCut, bool statusSingletThreshold, double singletThreshold, bool statusDoubletThresholdTotalChargeCut, double doubletThreshold, double doubletTotalChargeCut){
  vector<vector<int> > clusterStNum=event->clusterStripNumber(event, value, pdd, sigma, statusTotalChargeCut, totalChargeCut, statusSingletThreshold, singletThreshold, statusDoubletThresholdTotalChargeCut, doubletThreshold, doubletTotalChargeCut);
  return clusterStNum.size();
} 

vector<int> MwpcEvent::getClusterNStrip(MwpcEvent* event, int* value, const double* pdd, const double* sigma, bool statusTotalChargeCut, const double* totalChargeCut, bool statusSingletThreshold, double singletThreshold, bool statusDoubletThresholdTotalChargeCut, double doubletThreshold, double doubletTotalChargeCut){
  vector<vector<int> > clusterStNum=event->clusterStripNumber(event, value, pdd, sigma, statusTotalChargeCut, totalChargeCut, statusSingletThreshold, singletThreshold, statusDoubletThresholdTotalChargeCut, doubletThreshold, doubletTotalChargeCut);
  vector<int> nStrip;
  nStrip.clear();
  for(int i=0;i<(int)clusterStNum.size();i++){
    nStrip.push_back(clusterStNum[i].size());
  }
  return nStrip;
}

double* MwpcEvent::getClusterTotCharge(MwpcEvent* event, int* value, const double* pdd, const double* sigma, bool statusTotalChargeCut, const double* totalChargeCut, bool statusSingletThreshold, double singletThreshold, bool statusDoubletThresholdTotalChargeCut, double doubletThreshold, double doubletTotalChargeCut){
  double* AdcValue=event->getAdcValue(value, pdd, sigma, 0);
  vector<vector<int> > clusterStNum=event->clusterStripNumber(event, value, pdd, sigma, statusTotalChargeCut, totalChargeCut, statusSingletThreshold, singletThreshold, statusDoubletThresholdTotalChargeCut, doubletThreshold, doubletTotalChargeCut);
  double *totCharge=new double[clusterStNum.size()];
  double temp;
  for(int i=0;i<(int)clusterStNum.size();i++){
    temp=0;
    for(int j=0;j<(int)clusterStNum[i].size();j++){
      temp+=AdcValue[clusterStNum[i][j]];
    }
    totCharge[i]=temp;
  }
  return totCharge;
}

vector<vector<double> > MwpcEvent::clusterStripCharge(MwpcEvent* event, int* value, const double* pdd, const double* sigma, bool statusTotalChargeCut, const double* totalChargeCut, bool statusSingletThreshold, double singletThreshold, bool statusDoubletThresholdTotalChargeCut, double doubletThreshold, double doubletTotalChargeCut){
  double* AdcValue=event->getAdcValue(value, pdd, sigma, 0);
  vector<vector<int> > clusterStNum=event->clusterStripNumber(event, value, pdd, sigma, statusTotalChargeCut, totalChargeCut, statusSingletThreshold, singletThreshold, statusDoubletThresholdTotalChargeCut, doubletThreshold, doubletTotalChargeCut);
  vector<vector<double> >clusterStCharge; 
  clusterStCharge.clear();
  vector<double> stCharge;
  for(int i=0;i<(int)clusterStNum.size();i++){
    stCharge.clear();
    for(int j=0;j<(int)clusterStNum[i].size();j++){
      stCharge.push_back(AdcValue[clusterStNum[i][j]]);
    }
    clusterStCharge.push_back(stCharge);
  }
  return clusterStCharge;
}







vector<double> MwpcEvent::getEachClusterStCharge(MwpcEvent* event,int* value, const double* pdd,vector<int> stNum){
  double* AdcValue=event->getAdcValue(value, pdd);
  vector<double> stCharge;
  stCharge.clear();
  for(UInt_t i=0;i<stNum.size();i++){
    stCharge.push_back(AdcValue[stNum[i]]);
  }
  return stCharge;
}

double MwpcEvent::getEachClusterPositionCentroid(vector<int> stNum,vector<double> stCharge){
  double sumCharge=0,sumChargeMultStNum=0;
  for(UInt_t i=0;i<stNum.size();i++){
    sumCharge+=stCharge[i];
    sumChargeMultStNum+=stCharge[i]*stNum[i];
  }
  if(sumCharge==0){
    cout<<"There is 0 strip for this cluster; Error takes place."<<endl;
    exit;
  }
  else
    return sumChargeMultStNum/sumCharge+0.5;
}

double MwpcEvent::getEachClusterPositionCentroid(MwpcEvent* event,int* value, const double* pdd,vector<int> stNum){
  double* AdcValue=event->getAdcValue(value, pdd);
  double sumCharge=0,sumChargeMultStNum=0;
  for(UInt_t i=0;i<stNum.size();i++){
    sumCharge+=AdcValue[stNum[i]];
    sumChargeMultStNum+=AdcValue[stNum[i]]*stNum[i];
  }
  if(sumCharge==0){
    cout<<"There is 0 strip for this cluster; Error takes place."<<endl;
    exit;
  }
  else
    return sumChargeMultStNum/sumCharge+0.5;
}

double MwpcEvent::getEachClusterPositionErrorCentroid(vector<double> stCharge){
  double stError=0.5;
  double sumCharge=0,sumChargeSq=0;
  for(UInt_t i=0;i<stCharge.size();i++){
    sumCharge+=stCharge[i];
    sumChargeSq+=stCharge[i]*stCharge[i];
  }

  if(sumCharge==0){
    cout<<"There is 0 strip for this cluster; Error takes place."<<endl;
    exit;
  }
  else
    return sqrt(sumChargeSq)*stError/sumCharge;
}

double MwpcEvent::getEachClusterPositionErrorCentroid(MwpcEvent* event,int* value, const double* pdd,vector<int> stNum){
  double stError=0.5;
  double* AdcValue=event->getAdcValue(value, pdd);
  double sumCharge=0,sumChargeSq=0;
  for(UInt_t i=0;i<stNum.size();i++){
    sumCharge+=AdcValue[stNum[i]];
    sumChargeSq+=AdcValue[stNum[i]]*AdcValue[stNum[i]];
  }
  if(sumCharge==0){
    cout<<"There is 0 strip for this cluster; Error takes place."<<endl;
    exit;
  }
  else
    return sqrt(sumChargeSq)*stError/sumCharge;
}

int MwpcEvent::getEachClusterStNumWithHighestADC(vector<int> stNum,vector<double> stCharge){
  double temp=0;
  int num;
  for(UInt_t i=0;i<stNum.size();i++){
    if(temp<stCharge[i]){
      temp=stCharge[i];
      num=stNum[i];
    }
  }
  return num;
}

double MwpcEvent::getEachClusterPositionCentroidUpto3St(vector<int> stNum,vector<double> stCharge){
  double sumCharge=0,sumChargeMultStNum=0;
  if(stNum.size()<3){
    for(UInt_t i=0;i<stNum.size();i++){
      sumCharge+=stCharge[i];
      sumChargeMultStNum+=stCharge[i]*stNum[i];
    }
  }
  else{
    double temp=0;
    UInt_t index;
    for(UInt_t i=0;i<stNum.size();i++){
      if(temp<stCharge[i]){
        temp=stCharge[i];
        index=i;
      }
    }
    if(index==0){
      for(UInt_t i=0;i<=2;i++){
        sumCharge+=stCharge[i];
        sumChargeMultStNum+=stCharge[i]*stNum[i];
      }
    }
    else if(index==stNum.size()-1){
      for(UInt_t i=stNum.size()-3;i<=stNum.size()-1;i++){
        sumCharge+=stCharge[i];
        sumChargeMultStNum+=stCharge[i]*stNum[i];
      }
    }
    else{
      for(UInt_t i=index-1;i<=index+1;i++){
        sumCharge+=stCharge[i];
        sumChargeMultStNum+=stCharge[i]*stNum[i];
      }
    }
  }

  if(sumCharge==0){
    cout<<"There is 0 strip for this cluster; Error takes place."<<endl;
    exit;
  }
  else
    return sumChargeMultStNum/sumCharge+0.5;
}

double MwpcEvent::getEachClusterPositionErrorCentroidUpto3St(vector<double> stCharge){
  double stError=0.5;
  double sumCharge=0,sumChargeSq=0;
  if(stCharge.size()<3){
    for(UInt_t i=0;i<stCharge.size();i++){
      sumCharge+=stCharge[i];
      sumChargeSq+=stCharge[i]*stCharge[i];
    }
  }
  else{
    double temp=0;
    UInt_t index;
    for(UInt_t i=0;i<stCharge.size();i++){
      if(temp<stCharge[i]){
        temp=stCharge[i];
        index=i;
      }
    }
    if(index==0){
      for(UInt_t i=0;i<=2;i++){
        sumCharge+=stCharge[i];
        sumChargeSq+=stCharge[i]*stCharge[i];
      }
    }
    else if(index==stCharge.size()-1){
      for(UInt_t i=stCharge.size()-3;i<=stCharge.size()-1;i++){
        sumCharge+=stCharge[i];
        sumChargeSq+=stCharge[i]*stCharge[i];
      }
    }
    else{
      for(UInt_t i=index-1;i<=index+1;i++){
        sumCharge+=stCharge[i];
        sumChargeSq+=stCharge[i]*stCharge[i];
      }
    }
  }
  if(sumCharge==0){
    cout<<"There is 0 strip for this cluster; Error takes place."<<endl;
    exit;
  }
  else
    return sqrt(sumChargeSq)*stError/sumCharge;
}

double MwpcEvent::getEachClusterPositionChargeRatio(vector<int> stNum,vector<double> stCharge){
  if(stNum.size()==0){
    cout<<"There is 0 strip for this cluster; Error takes place."<<endl;
    exit;
  }

  double ratio;
  if(stNum.size()==1) return stNum[0]+0.5;
  else{
    double temp=0;
    UInt_t index;
    for(UInt_t i=0;i<stNum.size();i++){
      if(temp<stCharge[i]){
        temp=stCharge[i];
        index=i;
      }
    }
    if(index==0){
      ratio=(stCharge[0]-stCharge[1])/stCharge[0];
    }   
    else if(index==stNum.size()-1){
      ratio=stCharge[stNum.size()-1]/(stCharge[stNum.size()-1]-stCharge[stNum.size()-2]);
    }
    else{
      ratio=(stCharge[index]-stCharge[index+1])/(stCharge[index]-stCharge[index-1]);
    }
 
    int digit1,digit2,digit3,digit4;
    if(ratio>=ratioArray[0])
      return stNum[index];
    else{
      for(int i=0;i<5;i++){
        if(i==4){
          if(ratio<ratioArray[1000*i] && ratio>=ratioArray[1000*(i+1)-1]){
            digit1=i;
            break;
          }
        }
        else{
          if(ratio<ratioArray[1000*i] && ratio>=ratioArray[1000*(i+1)]){
            digit1=i;
            break;
          }
        }
      }
      for(int i=0;i<10;i++){
        if(digit1==4 && i==9){
          if(ratio<ratioArray[digit1*1000+100*i] && ratio>=ratioArray[digit1*1000+100*(i+1)-1]){
            digit2=i;
            break;
          }
        }
        else{
          if(ratio<ratioArray[digit1*1000+100*i] && ratio>=ratioArray[digit1*1000+100*(i+1)]){
            digit2=i;
            break;
          }
        }
      }
      for(int i=0;i<10;i++){
        if(digit1==4 && digit2==9 && i==9){
          if(ratio<ratioArray[digit1*1000+digit2*100+10*i] && ratio>=ratioArray[digit1*1000+digit2*100+10*(i+1)-1]){
            digit3=i;
            break;
          }
        }
        else{
          if(ratio<ratioArray[digit1*1000+digit2*100+10*i] && ratio>=ratioArray[digit1*1000+digit2*100+10*(i+1)]){
            digit3=i;
            break;
          }
        }
      }
      for(int i=0;i<10;i++){
        if(digit1==4 && digit2==9 && digit3==9 && i==9){
          if(ratio<ratioArray[digit1*1000+digit2*100+digit3*10+i]){
            cout<<"ratio < 0, it is impossible. Error takes place."<<endl;
            exit;
          } 
        }
        else{
          if(ratio<ratioArray[digit1*1000+digit2*100+digit3*10+i] && ratio>=ratioArray[digit1*1000+digit2*100+digit3*10+i+1]){
            digit4=i;
            break;
          }
        }
      }
      return stNum[index]+posArray[digit1*1000+digit2*100+digit3*10+digit4];
    }
  }
}

MwpcEvent::~MwpcEvent(){}
ClassImp(MwpcEvent);
