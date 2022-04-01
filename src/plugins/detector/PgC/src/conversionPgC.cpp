/*
 * Simple program to monitor the ped shits during a run
 */
#include <Det_PgC.h>
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
#include <sstream>

using namespace std;

Long_t Det_PgC::conversionPgC_histos()
{
  h_sum=new TH1D("h_sum","adc",100,-1000,4500);
  return 0;

}

static Int_t run,event;
float adcPGC[12][7];
Long_t Det_PgC::conversionPgC_begin(){
	getBranchObject("RawPgcInfo",(TObject **) &treeRaw);
	treeCali=new CRTCaliPgC();
	makeBranch("hitPgC",(TObject **)&treeCali);
 // out->Branch("run",&run,"run/I");
 // out->Branch("event",&event,"event/I");
 // out->Branch("adcPGC",adcPGC,"adcPGC[12][7]/F");




         gStyle->SetOptStat(0);
	return 0;
};

Long_t Det_PgC::conversionPgC_calib(){
	
    int runNo = treeRaw->run;
          
for (int i=0; i<12; i++)
     {
       int j=0;
       ifstream file(Form("/data/trek/E36/gautam/crnRoot/cooker/datfiles/PgC_kmu2/PgC_ped%d_%d.dat",runNo, i+1), ios::in);
 
       string line;
       if ( file.is_open() )
     {
       while ( ! file.eof() )
         {
           getline(file,line);
           if(line.size()>0)
         {
           istringstream vars(line);
 
           vars>>pgc_ped[i][j];

            j++;
 
         }
         }
     }
       file.close();
 
     }
  //run=treeRaw->run;
//  event=treeRaw->event;
  treeCali->run=treeRaw->run;
  treeCali->event=treeRaw->event;
 //  adcPGC.clear();
treeCali->adcPgC[12][7]=0;  
// adcPGC[12][7]=0;  

for(int iPgC=0;iPgC<12;iPgC++){

       for(int jPgC=0;jPgC<7;jPgC++){
     myAdcPgC[iPgC][jPgC] = treeRaw->ADC_PGC_Gap[iPgC][jPgC]- pgc_ped[iPgC][jPgC];
      
        h_sum->Fill(myAdcPgC[iPgC][jPgC]- pgc_ped[iPgC][jPgC]);
      
  //      cout<<iPgC<<"   "<<jPgC<<"  PgC: "<<pgc_ped[iPgC][jPgC]<<endl;
//      adcPGC.push_back(myAdcPgC[iPgC][jPgC]);
//      treeCali->adcPgC.push_back(myAdcPgC[iPgC][jPgC]);
      treeCali->adcPgC[iPgC][jPgC] = myAdcPgC[iPgC][jPgC];
//     adcPGC[iPgC][jPgC] = myAdcPgC[iPgC][jPgC];
     
    }
} 
 return 0;
};

