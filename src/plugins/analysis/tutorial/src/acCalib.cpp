#include <tutorial.h>

#include<cmath>
#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <vector>
using namespace std;

Long_t tutorial::hist_calibS(){
  h1Pedestal=new TH1D("Peak","Calibration Peak for AC;AC Paddle [U1-12,D1-12]",24,0.5,24.5);
  for(int iAc=0;iAc<12;iAc++){
    string name("AcU_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    pedADC_AcU[iAc]=new peakSelect(100,0,1500,name);

  }
  
  for(int iAc=0;iAc<12;iAc++){
    string name("AcD_");
    name+=(iAc+1)/10+'0';
    name+=(iAc+1)%10+'0';
    pedADC_AcD[iAc]=new peakSelect(100,0,1500,name);
  }

  return 0;

}

Long_t tutorial::start_calibS(){
  getBranchObject("RawAcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0);
  return 0;
};

Long_t tutorial::calib_peakS(){
  for(int iAc = 0; iAc < 12; iAc++){
    if(treeRaw->ADC_ACU[iAc] > 0)
      pedADC_AcU[iAc]->addData(treeRaw->ADC_ACU[iAc]);
    if(treeRaw->ADC_ACD[iAc] > 0)
      pedADC_AcD[iAc]->addData(treeRaw->ADC_ACD[iAc]);
  }
  return 0;
};

Long_t tutorial::done_calibS(){
  
  ofstream file1, file2, file3, file4;
  int runNo = treeRaw->run;
  //file1.open("/data/trek/E36/gautam/crnRoot/cooker/datfiles/pedestal_run/PedRun"+std::to_string(runNo)+"U.txt", ios::app);
  //file2.open("/data/trek/E36/gautam/crnRoot/cooker/datfiles/pedestal_run/PedRun"+std::to_string(runNo)+"D.txt", ios::app);
  file1.open("/data/trek/E36/gautam/crnRoot/cooker/datfiles/Pro_run/ProRun"+std::to_string(runNo)+"U.txt", ios::app);
  file2.open("/data/trek/E36/gautam/crnRoot/cooker/datfiles/Pro_run/ProRun"+std::to_string(runNo)+"D.txt", ios::app);
  file3.open("/data/trek/E36/gautam/crnRoot/cooker/datfiles/Pro_run/FWHM/ProRun_sigma"+std::to_string(runNo)+"U.txt", ios::app);
  file4.open("/data/trek/E36/gautam/crnRoot/cooker/datfiles/Pro_run/FWHM/ProRun_sigma"+std::to_string(runNo)+"D.txt", ios::app);
//  file3.open("datfiles/sigma"+std::to_string(runNo)+"U.txt", ios::app);
//  file4.open("datfiles/sigma"+std::to_string(runNo)+"D.txt", ios::app);

  double pedestalU[24],pedestalD[24],sigmaU[24],sigmaD[24];
  for(int iAc = 0; iAc < 12; iAc++){
    double ped1U=pedADC_AcU[iAc]->getPedestal();
    double err1U=pedADC_AcU[iAc]->getError();
    double sigmU=pedADC_AcU[iAc]->getSigma();
    double errSU=pedADC_AcU[iAc]->getSigmaError();
    h1Pedestal->SetBinContent(1+iAc,ped1U);
    h1Pedestal->SetBinError(1+iAc,err1U);
    h1Pedestal->SetBinError(1+iAc,sigmU);
    h1Pedestal->SetBinError(1+iAc,errSU);
    sigmU = 2*sqrt(2*log(2))*sigmU;
    pedestalU[iAc]=ped1U;
    sigmaU[iAc] = sigmU;
    double ped1D=pedADC_AcD[iAc]->getPedestal();
    double err1D=pedADC_AcD[iAc]->getError();
    double sigmD=pedADC_AcD[iAc]->getSigma();
    double errSD=pedADC_AcD[iAc]->getSigmaError();
    h1Pedestal->SetBinContent(13+iAc,ped1D);
    h1Pedestal->SetBinError(13+iAc,err1D);
    h1Pedestal->SetBinError(1+iAc,sigmD);
    h1Pedestal->SetBinError(1+iAc,errSD);
    pedestalD[iAc]=ped1D;
    sigmD = 2*sqrt(2*log(2))*sigmD;
    sigmaD[iAc] = sigmD;
  }
  int size_of_array = sizeof(pedestalU)/sizeof(pedestalU[0]);
  std::vector<double> pedestal_U(pedestalU, pedestalU+size_of_array);
  std::vector<double> sigma_U(sigmaU, sigmaU+size_of_array);
  pedestal_U.insert(pedestal_U.begin(),runNo);
  sigma_U.insert(sigma_U.begin(),runNo);
  

  std::vector<double> pedestal_D(pedestalD, pedestalD+size_of_array);
  std::vector<double> sigma_D(sigmaD, sigmaD+size_of_array);
  pedestal_D.insert(pedestal_D.begin(),runNo);
  sigma_D.insert(sigma_D.begin(),runNo);
  

  for(Int_t j=0; j<13 ;j++){
    cout <<"pedestal_U : "  <<  pedestal_U[j] << " sigma_U :  "<<sigma_U[j]<<"    pedestal_D : "  <<  pedestal_D[j] << " sigma_D :  "<<sigma_D[j]<<endl;
  
    file1 << pedestal_U[j] << "\t" ;
    file2 << pedestal_D[j] << "\t" ;
    file3 << sigma_U[j] << "\t" ;
    file4 << sigma_D[j] << "\t" ;
  }
  file1 << "\n";
  file2 << "\n";
  file3 << "\n";
  file4 << "\n";

  return 0;
}
