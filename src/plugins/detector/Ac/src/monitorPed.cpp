/*
 * Simple program to monitor the ped shits during a run
 */
#include <Det_Ac.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <TStyle.h>
#include <string>
#include <TSpectrum.h>
#include <TH1D.h>
#include <TDirectory.h>
#include <TList.h>


using namespace std;

Long_t Det_Ac::monitor_ped(){
  h1Pedestal=new TH1D("Pedestal","Pedestal for AC;AC Paddle [U1-12,D1-12]",24,0.5,24.5);
  for(int iAc=0;iAc<12;iAc++){
    string name("AcU1_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcU1[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcD1_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcD1[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcU2_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcU2[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcD2_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcD2[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcU3_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcU3[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcD3_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcD3[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcU4_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcU4[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcD4_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcD4[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcU5_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcU5[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcD5_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcD5[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcU6_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcU6[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcD6_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcD6[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcU7_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcU7[iAc]=new genericPedMon(150,0,1500,name);
  }
  for(int iAc=0;iAc<12;iAc++){
    string name("AcD7_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    adcAcD7[iAc]=new genericPedMon(150,0,1500,name);
  }
  return 0;
}

Long_t Det_Ac::monitor_begin(){
  getBranchObject("RawAcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0);
  return 0;
};
Long_t Det_Ac::monitor_done(){
  double pedestal1[24], pedestal2[24],pedestal3[24],pedestal4[24],pedestal5[24];
  double pedestal6[24], pedestal7[24];
  double pedSigma1[24], pedSigma2[24], pedSigma3[24], pedSigma4[24], pedSigma5[24], pedSigma6[24];
  double pedSigma7[24];
  double fWHM1[24], fWHM2[24], fWHM3[24], fWHM4[24], fWHM5[24], fWHM6[24], fWHM7[24];
  ofstream file1, file2, file3, file4;
  int runNo = treeRaw->run;
  int n1 = 120000, n2 = 140000;
  int evtNo = tree->GetEntries();
  if(evtNo < 100000){
    cout << " Error! \n Not enough events for an evaluation " << endl;
    //abort();
  }else if (evtNo > 100000){
    file1.open("pedRun"+std::to_string(runNo)+"D.txt", ios::app);
    file2.open("pedRun"+std::to_string(runNo)+"U.txt", ios::app);
    file3.open("fwhmRn"+std::to_string(runNo)+"D.txt", ios::app);
    file4.open("fwhmRn"+std::to_string(runNo)+"U.txt", ios::app);
    for(int n = 0; n < 20000; n+=1){
      tree->GetEntry(n);
      for(int iAc=0;iAc<12;iAc++){
        if(treeRaw->ADC_ACU[iAc]>0)
          adcAcU1[iAc]->addData(treeRaw->ADC_ACU[iAc]);
        if(treeRaw->ADC_ACD[iAc]>0)
          adcAcD1[iAc]->addData(treeRaw->ADC_ACD[iAc]);
      }
    }
    for(int n = 20000; n < 40000; n+=1){
      tree->GetEntry(n);
      for(int iAc=0;iAc<12;iAc++){
        if(treeRaw->ADC_ACU[iAc]>0)
          adcAcU2[iAc]->addData(treeRaw->ADC_ACU[iAc]);
        if(treeRaw->ADC_ACD[iAc]>0)
          adcAcD2[iAc]->addData(treeRaw->ADC_ACD[iAc]);
      }
    }
    for(int n = 40000; n < 60000; n+=1){
      tree->GetEntry(n);
      for(int iAc=0;iAc<12;iAc++){
        if(treeRaw->ADC_ACU[iAc]>0)
          adcAcU3[iAc]->addData(treeRaw->ADC_ACU[iAc]);
        if(treeRaw->ADC_ACD[iAc]>0)
          adcAcD3[iAc]->addData(treeRaw->ADC_ACD[iAc]);
      }
    }
    for(int n = 60000; n < 80000; n+=1){
      tree->GetEntry(n);
      for(int iAc=0;iAc<12;iAc++){
        if(treeRaw->ADC_ACU[iAc]>0)
          adcAcU4[iAc]->addData(treeRaw->ADC_ACU[iAc]);
        if(treeRaw->ADC_ACD[iAc]>0)
          adcAcD4[iAc]->addData(treeRaw->ADC_ACD[iAc]);
      }
    }
    for(int n = 80000; n < 100000; n+=1){
      tree->GetEntry(n);
      for(int iAc=0;iAc<12;iAc++){
        if(treeRaw->ADC_ACU[iAc]>0)
          adcAcU5[iAc]->addData(treeRaw->ADC_ACU[iAc]);
        if(treeRaw->ADC_ACD[iAc]>0)
          adcAcD5[iAc]->addData(treeRaw->ADC_ACD[iAc]);
      }
    }
    // write peds to file
    for(int iAc = 0; iAc < 12; iAc+=1){
      TSpectrum *s1 = new TSpectrum(4);
      double ped1U=adcAcU1[iAc]->getPedestal(s1);
      double err1U=adcAcU1[iAc]->getError();
      double sigmU=adcAcU1[iAc]->getSigma();
      double errSU=adcAcU1[iAc]->getSigmaError();
      double fwhmU=adcAcU1[iAc]->getFWHM();
      h1Pedestal->SetBinContent(1+iAc,ped1U);
      h1Pedestal->SetBinError(1+iAc,err1U);
      h1Pedestal->SetBinError(1+iAc,sigmU);
      h1Pedestal->SetBinError(1+iAc,errSU);
      pedestal1[iAc]=ped1U;
      pedSigma1[iAc]=sigmU;
      fWHM1[iAc]=fwhmU;
      
      double ped1D=adcAcD1[iAc]->getPedestal(s1);
      double err1D=adcAcD1[iAc]->getError();
      double sigmD=adcAcD1[iAc]->getSigma();
      double errSD=adcAcD1[iAc]->getSigmaError();
      double fwhmD=adcAcD1[iAc]->getFWHM();
      h1Pedestal->SetBinContent(13+iAc,ped1D);
      h1Pedestal->SetBinError(13+iAc,err1D);
      h1Pedestal->SetBinError(1+iAc,sigmD);
      h1Pedestal->SetBinError(1+iAc,errSD);
      pedestal1[iAc+12]=ped1D;
      pedSigma1[iAc+12]=sigmD;
      fWHM1[iAc+12]=fwhmD;
      file1<< pedestal1[iAc+12] <<"\t";
      file2<< pedestal1[iAc] <<"\t";
      file3<< fWHM1[iAc+12] <<"\t";
      file4<< fWHM1[iAc] <<"\t";
    }
    file1 << "\n";
    file2 << "\n";
    file3 << "\n";
    file4 << "\n";
    for(int iAc = 0; iAc < 12; iAc+=1){
      TSpectrum *s2 = new TSpectrum(4);
      double ped1U=adcAcU2[iAc]->getPedestal(s2);
      double err1U=adcAcU2[iAc]->getError();
      double sigmU=adcAcU2[iAc]->getSigma();
      double errSU=adcAcU2[iAc]->getSigmaError();
      double fwhmU=adcAcU2[iAc]->getFWHM();
      h1Pedestal->SetBinContent(1+iAc,ped1U);
      h1Pedestal->SetBinError(1+iAc,err1U);
      h1Pedestal->SetBinError(1+iAc,sigmU);
      h1Pedestal->SetBinError(1+iAc,errSU);
      pedestal2[iAc]=ped1U;
      pedSigma2[iAc]=sigmU;
      fWHM2[iAc]=fwhmU;

      
      double ped1D=adcAcD2[iAc]->getPedestal(s2);
      double err1D=adcAcD2[iAc]->getError();
      double sigmD=adcAcD2[iAc]->getSigma();
      double errSD=adcAcD2[iAc]->getSigmaError();
      double fwhmD=adcAcD2[iAc]->getFWHM();
      h1Pedestal->SetBinContent(13+iAc,ped1D);
      h1Pedestal->SetBinError(13+iAc,err1D);
      h1Pedestal->SetBinError(1+iAc,sigmD);
      h1Pedestal->SetBinError(1+iAc,errSD);
      pedestal2[iAc+12]=ped1D;
      pedSigma2[iAc+12]=sigmD;
      fWHM2[iAc+12]=fwhmD;
      file1<< pedestal2[iAc+12] <<"\t";
      file2<< pedestal2[iAc] <<"\t";
      file3<< fWHM2[iAc+12] <<"\t";
      file4<< fWHM2[iAc] <<"\t";

    }
    file1 << "\n";
    file2 << "\n";
    file3 << "\n";
    file4 << "\n";
    for(int iAc = 0; iAc < 12; iAc+=1){
      TSpectrum *s3 = new TSpectrum(4);
      double ped1U=adcAcU3[iAc]->getPedestal(s3);
      double err1U=adcAcU3[iAc]->getError();
      double sigmU=adcAcU3[iAc]->getSigma();
      double errSU=adcAcU3[iAc]->getSigmaError();
      double fwhmU=adcAcU3[iAc]->getFWHM();
      h1Pedestal->SetBinContent(1+iAc,ped1U);
      h1Pedestal->SetBinError(1+iAc,err1U);
      h1Pedestal->SetBinError(1+iAc,sigmU);
      h1Pedestal->SetBinError(1+iAc,errSU);
      pedestal3[iAc]=ped1U;
      pedSigma3[iAc]=sigmU;
      fWHM3[iAc]=fwhmU;
      
      double ped1D=adcAcD3[iAc]->getPedestal(s3);
      double err1D=adcAcD3[iAc]->getError();
      double sigmD=adcAcD3[iAc]->getSigma();
      double errSD=adcAcD3[iAc]->getSigmaError();
      double fwhmD=adcAcD3[iAc]->getFWHM();
      h1Pedestal->SetBinContent(13+iAc,ped1D);
      h1Pedestal->SetBinError(13+iAc,err1D);
      h1Pedestal->SetBinError(1+iAc,sigmD);
      h1Pedestal->SetBinError(1+iAc,errSD);
      pedestal3[iAc+12]=ped1D;
      pedSigma3[iAc+12]=sigmD;
      fWHM3[iAc+12]=fwhmD;
      file1<< pedestal1[iAc+12] <<"\t";
      file2<< pedestal1[iAc] <<"\t";
      file3<< fWHM3[iAc+12] <<"\t";
      file4<< fWHM3[iAc] <<"\t";
    }
    file1 << "\n";
    file2 << "\n";
    file3 << "\n";
    file4 << "\n";
    for(int iAc = 0; iAc < 12; iAc+=1){
      TSpectrum *s4 = new TSpectrum(4);
      double ped1U=adcAcU4[iAc]->getPedestal(s4);
      double err1U=adcAcU4[iAc]->getError();
      double sigmU=adcAcU4[iAc]->getSigma();
      double errSU=adcAcU4[iAc]->getSigmaError();
      double fwhmU=adcAcU4[iAc]->getFWHM();
      h1Pedestal->SetBinContent(1+iAc,ped1U);
      h1Pedestal->SetBinError(1+iAc,err1U);
      h1Pedestal->SetBinError(1+iAc,sigmU);
      h1Pedestal->SetBinError(1+iAc,errSU);
      pedestal4[iAc]=ped1U;
      pedSigma4[iAc]=sigmU;
      fWHM4[iAc]=fwhmU;
      
      double ped1D=adcAcD4[iAc]->getPedestal(s4);
      double err1D=adcAcD4[iAc]->getError();
      double sigmD=adcAcD4[iAc]->getSigma();
      double errSD=adcAcD4[iAc]->getSigmaError();
      double fwhmD=adcAcD4[iAc]->getFWHM();
      h1Pedestal->SetBinContent(13+iAc,ped1D);
      h1Pedestal->SetBinError(13+iAc,err1D);
      h1Pedestal->SetBinError(1+iAc,sigmD);
      h1Pedestal->SetBinError(1+iAc,errSD);
      pedestal4[iAc+12]=ped1D;
      pedSigma4[iAc+12]=sigmD;
      fWHM4[iAc+12]=fwhmD;
      file1<< pedestal4[iAc+12] <<"\t";
      file2<< pedestal4[iAc] <<"\t";
      file3<< fWHM4[iAc+12] <<"\t";
      file4<< fWHM4[iAc] <<"\t";
    }
    file1 << "\n";
    file2 << "\n";
    file3 << "\n";
    file4 << "\n";
    for(int iAc = 0; iAc < 12; iAc+=1){
      TSpectrum *s5 = new TSpectrum(4);
      double ped1U=adcAcU5[iAc]->getPedestal(s5);
      double err1U=adcAcU5[iAc]->getError();
      double sigmU=adcAcU5[iAc]->getSigma();
      double errSU=adcAcU5[iAc]->getSigmaError();
      double fwhmU=adcAcU5[iAc]->getFWHM();
      h1Pedestal->SetBinContent(1+iAc,ped1U);
      h1Pedestal->SetBinError(1+iAc,err1U);
      h1Pedestal->SetBinError(1+iAc,sigmU);
      h1Pedestal->SetBinError(1+iAc,errSU);
      pedestal5[iAc]=ped1U;
      pedSigma5[iAc]=sigmU;
      fWHM5[iAc]=fwhmU;
      
      double ped1D=adcAcD5[iAc]->getPedestal(s5);
      double err1D=adcAcD5[iAc]->getError();
      double sigmD=adcAcD5[iAc]->getSigma();
      double errSD=adcAcD5[iAc]->getSigmaError();
      double fwhmD=adcAcD5[iAc]->getFWHM();
      h1Pedestal->SetBinContent(13+iAc,ped1D);
      h1Pedestal->SetBinError(13+iAc,err1D);
      h1Pedestal->SetBinError(1+iAc,sigmD);
      h1Pedestal->SetBinError(1+iAc,errSD);
      pedestal5[iAc+12]=ped1D;
      pedSigma5[iAc+12]=sigmD;
      fWHM5[iAc+12]=fwhmD;
      file1<< pedestal5[iAc+12] <<"\t";
      file2<< pedestal5[iAc] <<"\t";
      file3<< fWHM5[iAc+12] <<"\t";
      file4<< fWHM5[iAc] <<"\t";
    }
    file1 << "\n";
    file2 << "\n";
    file3 << "\n";
    file4 << "\n";
    if(evtNo >= n2){
      for(int n = 100000; n < 120000; n+=1){
        tree->GetEntry(n);
        for(int iAc=0;iAc<12;iAc++){
          if(treeRaw->ADC_ACU[iAc]>0)
            adcAcU6[iAc]->addData(treeRaw->ADC_ACU[iAc]);
          if(treeRaw->ADC_ACD[iAc]>0)
            adcAcD6[iAc]->addData(treeRaw->ADC_ACD[iAc]);
        }
      }
      for(int iAc = 0; iAc < 12; iAc+=1){
        TSpectrum *s6 = new TSpectrum(4);
        double ped1U=adcAcU6[iAc]->getPedestal(s6);
        double err1U=adcAcU6[iAc]->getError();
        double sigmU=adcAcU6[iAc]->getSigma();
        double errSU=adcAcU6[iAc]->getSigmaError();
        double fwhmU=adcAcU6[iAc]->getFWHM();
        h1Pedestal->SetBinContent(1+iAc,ped1U);
        h1Pedestal->SetBinError(1+iAc,err1U);
        h1Pedestal->SetBinError(1+iAc,sigmU);
        h1Pedestal->SetBinError(1+iAc,errSU);
        pedestal6[iAc]=ped1U;
        pedSigma6[iAc]=sigmU;
        fWHM6[iAc]=fwhmU;
        
        double ped1D=adcAcD6[iAc]->getPedestal(s6);
        double err1D=adcAcD6[iAc]->getError();
        double sigmD=adcAcD6[iAc]->getSigma();
        double errSD=adcAcD6[iAc]->getSigmaError();
        double fwhmD=adcAcD6[iAc]->getFWHM();
        h1Pedestal->SetBinContent(13+iAc,ped1D);
        h1Pedestal->SetBinError(13+iAc,err1D);
        h1Pedestal->SetBinError(1+iAc,sigmD);
        h1Pedestal->SetBinError(1+iAc,errSD);
        pedestal6[iAc+12]=ped1D;
        pedSigma6[iAc+12]=sigmD;
        fWHM6[iAc+12]=fwhmD;
        file1<< pedestal6[iAc+12] <<"\t";
        file2<< pedestal6[iAc] <<"\t";
        file3<< fWHM6[iAc+12] <<"\t";
        file4<< fWHM6[iAc] <<"\t";
      }
      file1 << "\n";
      file2 << "\n";
      file3 << "\n";
      file4 << "\n";
      for(int n = 120000; n < 140000; n+=1){
        tree->GetEntry(n);
        for(int iAc=0;iAc<12;iAc++){
          if(treeRaw->ADC_ACU[iAc]>0)
            adcAcU7[iAc]->addData(treeRaw->ADC_ACU[iAc]);
          if(treeRaw->ADC_ACD[iAc]>0)
            adcAcD7[iAc]->addData(treeRaw->ADC_ACD[iAc]);
        }
      }
      for(int iAc = 0; iAc < 12; iAc+=1){
        TSpectrum *s7 = new TSpectrum(4);
        double ped1U=adcAcU7[iAc]->getPedestal(s7);
        double err1U=adcAcU7[iAc]->getError();
        double sigmU=adcAcU7[iAc]->getSigma();
        double errSU=adcAcU7[iAc]->getSigmaError();
        double fwhmU=adcAcU7[iAc]->getFWHM();
        h1Pedestal->SetBinContent(1+iAc,ped1U);
        h1Pedestal->SetBinError(1+iAc,err1U);
        h1Pedestal->SetBinError(1+iAc,sigmU);
        h1Pedestal->SetBinError(1+iAc,errSU);
        pedestal7[iAc]=ped1U;
        pedSigma7[iAc]=sigmU;
        fWHM7[iAc]=fwhmU;
        
        double ped1D=adcAcD7[iAc]->getPedestal(s7);
        double err1D=adcAcD7[iAc]->getError();
        double sigmD=adcAcD7[iAc]->getSigma();
        double errSD=adcAcD7[iAc]->getSigmaError();
        double fwhmD=adcAcD7[iAc]->getFWHM();
        h1Pedestal->SetBinContent(13+iAc,ped1D);
        h1Pedestal->SetBinError(13+iAc,err1D);
        h1Pedestal->SetBinError(1+iAc,sigmD);
        h1Pedestal->SetBinError(1+iAc,errSD);
        pedestal7[iAc+12]=ped1D;
        pedSigma7[iAc+12]=sigmD;
        fWHM7[iAc+12]=fwhmD;
        file1<< pedestal7[iAc+12] <<"\t";
        file2<< pedestal7[iAc] <<"\t";
        file3<< fWHM7[iAc+12] <<"\t";
        file4<< fWHM7[iAc] <<"\t";
      }
      file1 << "\n";
      file2 << "\n";
      file3 << "\n";
      file4 << "\n";
    }else if(evtNo < n2){
      for(int n = 100000; n < 120000; n+=1){
        tree->GetEntry(n);
        for(int iAc=0;iAc<12;iAc++){
          if(treeRaw->ADC_ACU[iAc]>0)
            adcAcU6[iAc]->addData(treeRaw->ADC_ACU[iAc]);
          if(treeRaw->ADC_ACD[iAc]>0)
            adcAcD6[iAc]->addData(treeRaw->ADC_ACD[iAc]);
        }
      }
      for(int iAc = 0; iAc < 12; iAc+=1){
        TSpectrum *s6 = new TSpectrum(4);
        double ped1U=adcAcU6[iAc]->getPedestal(s6);
        double err1U=adcAcU6[iAc]->getError();
        double sigmU=adcAcU6[iAc]->getSigma();
        double errSU=adcAcU6[iAc]->getSigmaError();
        double fwhmU=adcAcU6[iAc]->getFWHM();
        h1Pedestal->SetBinContent(1+iAc,ped1U);
        h1Pedestal->SetBinError(1+iAc,err1U);
        h1Pedestal->SetBinError(1+iAc,sigmU);
        h1Pedestal->SetBinError(1+iAc,errSU);
        pedestal6[iAc]=ped1U;
        pedSigma6[iAc]=sigmU;
        fWHM6[iAc]=fwhmU;
        
        double ped1D=adcAcD6[iAc]->getPedestal(s6);
        double err1D=adcAcD6[iAc]->getError();
        double sigmD=adcAcD6[iAc]->getSigma();
        double errSD=adcAcD6[iAc]->getSigmaError();
        double fwhmD=adcAcD6[iAc]->getFWHM();
        h1Pedestal->SetBinContent(13+iAc,ped1D);
        h1Pedestal->SetBinError(13+iAc,err1D);
        h1Pedestal->SetBinError(1+iAc,sigmD);
        h1Pedestal->SetBinError(1+iAc,errSD);
        pedestal6[iAc+12]=ped1D;
        pedSigma6[iAc+12]=sigmD;
        fWHM6[iAc+12]=fwhmD;
        file1<< pedestal6[iAc+12] <<"\t";
        file2<< pedestal6[iAc] <<"\t";
        file3<< fWHM6[iAc+12] <<"\t";
        file4<< fWHM6[iAc] <<"\t";
      }
      file1 << "\n";
      file2 << "\n";
      file3 << "\n";
      file4 << "\n";
    }
  } 
  file1.close();
  file2.close();
  file3.close();
  file4.close();
  return 0;
};

