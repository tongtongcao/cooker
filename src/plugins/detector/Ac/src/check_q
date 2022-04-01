#include <Det_Ac.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
using namespace std;

Long_t Det_Ac::histos_ped(){
  h1Pedestal=new TH1D("Pedestal","Pedestal for AC;AC Paddle [U1-12,D1-12]",24,0.5,24.5);
  for(int iAc=0;iAc<12;iAc++){
    string name("AcU_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcU[iAc]=new TrekAcAdcPedestal(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcD_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcD[iAc]=new TrekAcAdcPedestal(150,0,1500,name);
  }
  return 0;

}

Long_t Det_Ac::startup_ped(){
  getBranchObject("RawAcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0);
  return 0;
};

Long_t Det_Ac::process_ped(){
  for(int iAc = 0; iAc < 12; iAc++){
    if(treeRaw->ADC_ACU[iAc]>0)
      adcAcU[iAc]->addData(treeRaw->ADC_ACU[iAc]);
    if(treeRaw->ADC_ACD[iAc]>0)
      adcAcD[iAc]->addData(treeRaw->ADC_ACD[iAc]);
  }
  return 0;
};

Long_t Det_Ac::done_ped(){
  double pedestal[24];
  double pedSigma[24];
  ofstream file1, file2;
  int runNo = treeRaw->run;
  int evtNo = tree->GetEntries();
  int n1 = 120000, n2 = 140000;
  cout << " --- CHECKING THE No. OF ENTRIES " << tree->GetEntries() << endl;
  if(evtNo > 100000){
    file1.open("avePedNo"+std::to_string(runNo)+"D.txt", ios::app);
    file2.open("avePedNo"+std::to_string(runNo)+"U.txt", ios::app);
    for(int i = 0; i < 5; i+=1){
      for(int iAc=0;iAc<12;iAc++){
        double ped1U=adcAcU[iAc]->getPedestal();
        double err1U=adcAcU[iAc]->getError();
        double sigmU=adcAcU[iAc]->getSigma();
        double errSU=adcAcU[iAc]->getSigmaError();
        h1Pedestal->SetBinContent(1+iAc,ped1U);
        h1Pedestal->SetBinError(1+iAc,err1U);
        h1Pedestal->SetBinError(1+iAc,sigmU);
        h1Pedestal->SetBinError(1+iAc,errSU);
        pedestal[iAc]=ped1U;
        pedSigma[iAc]=sigmU;
        cout<<ped1U<<"&";
       
        double ped1D=adcAcD[iAc]->getPedestal();
        double err1D=adcAcD[iAc]->getError();
        double sigmD=adcAcD[iAc]->getSigma();
        double errSD=adcAcD[iAc]->getSigmaError();
        h1Pedestal->SetBinContent(13+iAc,ped1D);
        h1Pedestal->SetBinError(13+iAc,err1D);
        h1Pedestal->SetBinError(1+iAc,sigmD);
        h1Pedestal->SetBinError(1+iAc,errSD);
        cout<<ped1D<<"\\"<<endl;
        pedestal[iAc+12]=ped1D;
        pedSigma[iAc+12]=sigmD;
        file1<< pedestal[iAc+12] <<"\t";
        file2<< pedestal[iAc] <<"\t";
        /*file3<< pedSigma[iAc] <<"\t";
        file4<< pedSigma[iAc+12] <<"\t";
        cout<<"listPedestal["<<iAc<<"]="<<pedestal[iAc]<<";"<<endl;
        cout<<"listPedestal["<<iAc+12<<"]="<<pedestal[iAc+12]<<";"<<endl;*/
      } 
      file1 << "\n";
      file2 << "\n";
    }
    if(evtNo >= n2){
      for(int i = 0; i < 2; i+=1){
        for(int iAc=0;iAc<12;iAc++){
          double ped1U=adcAcU[iAc]->getPedestal();
          double err1U=adcAcU[iAc]->getError();
          double sigmU=adcAcU[iAc]->getSigma();
          double errSU=adcAcU[iAc]->getSigmaError();
          h1Pedestal->SetBinContent(1+iAc,ped1U);
          h1Pedestal->SetBinError(1+iAc,err1U);
          h1Pedestal->SetBinError(1+iAc,sigmU);
          h1Pedestal->SetBinError(1+iAc,errSU);
          pedestal[iAc]=ped1U;
          pedSigma[iAc]=sigmU;
          cout<<ped1U<<"&";
         
          double ped1D=adcAcD[iAc]->getPedestal();
          double err1D=adcAcD[iAc]->getError();
          double sigmD=adcAcD[iAc]->getSigma();
          double errSD=adcAcD[iAc]->getSigmaError();
          h1Pedestal->SetBinContent(13+iAc,ped1D);
          h1Pedestal->SetBinError(13+iAc,err1D);
          h1Pedestal->SetBinError(1+iAc,sigmD);
          h1Pedestal->SetBinError(1+iAc,errSD);
          cout<<ped1D<<"\\"<<endl;
          pedestal[iAc+12]=ped1D;
          pedSigma[iAc+12]=sigmD;
          file1<< pedestal[iAc+12] <<"\t";
          file2<< pedestal[iAc] <<"\t";
          /*file3<< pedSigma[iAc] <<"\t";
          file4<< pedSigma[iAc+12] <<"\t";
          cout<<"listPedestal["<<iAc<<"]="<<pedestal[iAc]<<";"<<endl;
          cout<<"listPedestal["<<iAc+12<<"]="<<pedestal[iAc+12]<<";"<<endl;*/
        } 
        file1 << "\n";
        file2 << "\n";
      }
    }else if(evtNo >= n1){
      for(int i = 0; i < 1; i+=1){
        for(int iAc=0;iAc<12;iAc++){
          double ped1U=adcAcU[iAc]->getPedestal();
          double err1U=adcAcU[iAc]->getError();
          double sigmU=adcAcU[iAc]->getSigma();
          double errSU=adcAcU[iAc]->getSigmaError();
          h1Pedestal->SetBinContent(1+iAc,ped1U);
          h1Pedestal->SetBinError(1+iAc,err1U);
          h1Pedestal->SetBinError(1+iAc,sigmU);
          h1Pedestal->SetBinError(1+iAc,errSU);
          pedestal[iAc]=ped1U;
          pedSigma[iAc]=sigmU;
          cout<<ped1U<<"&";
         
          double ped1D=adcAcD[iAc]->getPedestal();
          double err1D=adcAcD[iAc]->getError();
          double sigmD=adcAcD[iAc]->getSigma();
          double errSD=adcAcD[iAc]->getSigmaError();
          h1Pedestal->SetBinContent(13+iAc,ped1D);
          h1Pedestal->SetBinError(13+iAc,err1D);
          h1Pedestal->SetBinError(1+iAc,sigmD);
          h1Pedestal->SetBinError(1+iAc,errSD);
          cout<<ped1D<<"\\"<<endl;
          pedestal[iAc+12]=ped1D;
          pedSigma[iAc+12]=sigmD;
          file1<< pedestal[iAc+12] <<"\t";
          file2<< pedestal[iAc] <<"\t";
          /*file3<< pedSigma[iAc] <<"\t";
          file4<< pedSigma[iAc+12] <<"\t";
          cout<<"listPedestal["<<iAc<<"]="<<pedestal[iAc]<<";"<<endl;
          cout<<"listPedestal["<<iAc+12<<"]="<<pedestal[iAc+12]<<";"<<endl;*/
        } 
        file1 << "\n";
        file2 << "\n";
      }
    }
  }
		  
  file1.close();
  file2.close();
  /*
  file1.close();
  file2 << "\n";
  file2.close();
  file3 << "\n";
  file3.close();
  file4 << "\n";
  file4.close();
  */
  for(int i = 0; i < treeRaw->event; i+=1){
    if(i%10000==1)	cout<<"**** "<<i<<" events done"<<endl;
  }
  for(int iAc=0;iAc<12;iAc++){
    cout <<" " << pedSigma[iAc] << " \t " << pedSigma[iAc+12] <<endl;
  }

  cout << "pedU = [";
  for(int iAc=0;iAc<12;iAc++){
    cout<< pedestal[iAc] <<",";
  }
  cout <<"] \n";
  cout << "pedD = [";
  for(int iAc=0;iAc<12;iAc++){
    cout<< pedestal[iAc+12] <<",";
  }
  cout <<"] \n";

  cout << "sigU = [";
  for(int iAc=0;iAc<12;iAc++){
    cout<< pedSigma[iAc] <<",";
  }
  cout <<"] \n";
  cout << "sigD = [";
  for(int iAc=0;iAc<12;iAc++){
    cout<< pedSigma[iAc+12] <<",";
  }
  cout <<"] \n";

  cout<<"  " <<endl;
  for(int iAc=0;iAc<12;iAc++){
    cout<<"listPedsigma["<<iAc<<"]="   <<pedSigma[iAc]<<";"<<endl;
    cout<<"listPedsigma["<<iAc+12<<"]="<<pedSigma[iAc+12]<<";"<<endl;
  }
  cout<<"  " <<endl;
  for(int iAc=0;iAc<12;iAc++){
    cout <<" " << pedSigma[iAc] << " \t " << pedSigma[iAc+12] <<endl;
  }
  cout<<"  " <<endl;
  for(int iAc=0;iAc<12;iAc++){
    cout<<pedestal[iAc]<<"    "<<pedestal[iAc+12]<<endl;
  }
  return 0;
};












