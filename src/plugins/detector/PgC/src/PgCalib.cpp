#include <Det_PgC.h>

#include<cmath>
#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <vector>
using namespace std;

Long_t Det_PgC::hist_calibS(){
  for(int iAc=0;iAc<12;iAc++){
    for(int ii=0;ii<7;ii++){
      string name("Adc_PGC");
      name+=(iAc+1)/10+'0';
      name+=(iAc+1)%10+'0';
      name+= "_";
      name+=(ii+1)+'0';
      name+="_pgc";
      pedADC_PgC[iAc][ii]=new peakSelect(200,0,1600,name);

    }
  }  
  return 0;

}

Long_t Det_PgC::start_calibS(){
  getBranchObject("RawPgcInfo",(TObject **) &treeRaw);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
  gStyle->SetOptStat(0);
  return 0;
};

Long_t Det_PgC::calib_peakS(){
  for(int iAc = 0; iAc < 12; iAc++){
    for(int ii = 0;ii < 7; ii++){
      if(treeRaw->ADC_PGC_Gap[iAc][ii] > 20) pedADC_PgC[iAc][ii]->addData(treeRaw->ADC_PGC_Gap[iAc][ii]);
    }
  }
  return 0;
};

Long_t Det_PgC::done_calibS(){
  
  int runNo = treeRaw->run;
          TString fname_data[12];

	 for(int z=0; z<12; z++){

                  fname_data[z] = Form("/data/trek/E36/gautam/crnRoot/cooker/datfiles/PgC_kmu2/PgC_ped%d_%d.dat",runNo,z+1);
                    }

double pedestal_m[12][7],sigma_m[12][7];
  double pedestal[7],sigma[7];
  for(int iAc = 0; iAc < 12; iAc++){
       const char* fname_tmp=fname_data[iAc];
      FILE *parameter = fopen(fname_tmp,"w");
    for(int ii = 0;ii < 7; ii++){
      pedestal_m[iAc][ii]=pedADC_PgC[iAc][ii]->getPedestal();
      pedestal[ii] = pedestal_m[iAc][ii];

    fprintf(parameter,"%8.4f\n", pedestal[ii]);
 
    }
    fclose(parameter);
  } 
  
  return 0;
}
