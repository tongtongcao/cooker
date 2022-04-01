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

Long_t Det_Ac::conversionAc_histos()
{
  h_sum=new TH1D("h_sum","adc",100,-1000,4500);
  return 0;

}


Long_t Det_Ac::conversionAc_begin(){
	getBranchObject("RawAcInfo",(TObject **) &treeRaw);
	treeCali=new CRTCaliAc();
	makeBranch("hitAc",(TObject **)&treeCali);
gStyle->SetOptStat(0);
	return 0;
};
Long_t Det_Ac::conversionAc_calib(){
	int runNo = treeRaw->run;
	Int_t Run_U,Run_D;
	Double_t U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U11, U12;
	Double_t D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12;
	ifstream infile;
	infile.open("/data/trek/E36/gautam/crnRoot/cooker/datfiles/Pro_run/ProRun"+std::to_string(runNo)+"U.txt");
//	infile.open("/data/trek/E36/gautam/crnRoot/cooker/datfiles/pedestal_run/ProRun3399U.txt");
	infile>>fixed>>Run_U>>U1>>U2>>U3>>U4>>U5>>U6>>U7>>U8>>U9>>U10>>U11>>U12;
	infile.close();
	ifstream infile1;
	infile1.open("/data/trek/E36/gautam/crnRoot/cooker/datfiles/Pro_run/ProRun"+std::to_string(runNo)+"D.txt");
//	infile1.open("/data/trek/E36/gautam/crnRoot/cooker/datfiles/pedestal_run/ProRun3399D.txt");
	infile1>>fixed>>Run_D>>D1>>D2>>D3>>D4>>D5>>D6>>D7>>D8>>D9>>D10>>D11>>D12
     
;
	infile1.close();
	Double_t array_U[12] = {U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U11, U12};
	Double_t array_D[12] = {D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12}; 

  treeCali->run=treeRaw->run;
  treeCali->event=treeRaw->event;
  treeCali->adcACU.clear();
  treeCali->adcACD.clear();
//  treeCali->cer_adc_sum.clear();
//  double pedU[12] = {602.439,643.602,626.907,514.168,937.094,610.519,666.728,719.734,569.338,719.425,666.703,712.182};
 // double pedD[12] = {388.124,465.105,370.873,481.46,382.436,551.753,461.122,315.574,563.614,475.83,433.043,571.538};

double thr[24] = {324.1, 331.9, 318.6, 540.7, 381.8, 398.3,
                    589.1, 368.7, 361.3, 359.7, 322.5, 324.6,
                    318.8, 461.8, 392.8, 329.0, 388.2, 407.4,
                    339.9, 409.1, 400.5, 371.0, 308.0, 379.1};
       Float_t my_cer_adc_sum = 0;
  for(int iAc=0;iAc<12;iAc++){
//    cout<<" gap:  "<<iAc<<"  Up:  "<<array_U[iAc]<<"  Down  "<<array_D[iAc]<<endl;
     myAdcAcU[iAc] = treeRaw->ADC_ACU[iAc] - array_U[iAc];
     myAdcAcD[iAc]  = treeRaw->ADC_ACD[iAc] - array_D[iAc];
      
        h_sum->Fill(myAdcAcU[iAc]);
        h_sum->Fill(myAdcAcD[iAc]);
      

      // my_cer_adc_sum += (myAdcAcU + myAdcAcD);

     
      treeCali->adcACU.push_back(myAdcAcU[iAc]);
      treeCali->adcACD.push_back(myAdcAcD[iAc]);
     

  //   treeCali->cer_adc_sum.push_back(myAdcAcU);
//     treeCali->cer_adc_sum.push_back(myAdcAcD);



    }
 
 return 0;
};

