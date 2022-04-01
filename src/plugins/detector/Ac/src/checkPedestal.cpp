#include <Det_Ac.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPad.h>
using namespace std;

Long_t Det_Ac::histos_checkPed(){
  for(int iAc = 0; iAc < 12; iAc++){
    string name("Ac1U_");
    name+=iAc/10+'0';
    name+=iAc%10+'0';
    h1Adc[iAc]=new TH1D(name.c_str(),"ADC Upstr Spectrum",150,0,4100);
  }
  for(int iAc = 0; iAc < 12; iAc++){
    string name("Ac1D_");
    name+=iAc/10+'0';
    name+=iAc%10+'0';
    h1Adc[iAc+12]=new TH1D(name.c_str(),"ADC Dnstr Spectrum",150,0,4100);
  }
  return 0;

}

Long_t Det_Ac::startup_checkPed(){
  getBranchObject("RawAcInfo",(TObject **) &treeRaw);
  gStyle->SetOptStat(0);
  return 0;
};


Long_t Det_Ac::process_checkPed(){
  for(int iAc = 0; iAc < 12; iAc++){
    if(treeRaw->ADC_ACU[iAc] > 0)
      h1Adc[iAc]->Fill(treeRaw->ADC_ACU[iAc]);
    if(treeRaw->ADC_ACD[iAc] > 0)
      h1Adc[iAc+12]->Fill(treeRaw->ADC_ACD[iAc]);
  }

  return 0;
};

Long_t Det_Ac::done_checkPed(){
  vector<double> listPedestal(24,0);
  listPedestal[0]=602.439;
  listPedestal[12]=388.124;
  listPedestal[1]=643.602;
  listPedestal[13]=465.105;
  listPedestal[2]=626.907;
  listPedestal[14]=370.873;
  listPedestal[3]=514.168;
  listPedestal[15]=481.46;
  listPedestal[4]=937.094;
  listPedestal[16]=382.436;
  listPedestal[5]=610.519;
  listPedestal[17]=551.753;
  listPedestal[6]=666.728;
  listPedestal[18]=461.122;
  listPedestal[7]=719.734;
  listPedestal[19]=315.574;
  listPedestal[8]=569.338;
  listPedestal[20]=563.614;
  listPedestal[9]=719.425;
  listPedestal[21]=475.83;
  listPedestal[10]=666.703;
  listPedestal[22]=433.043;
  listPedestal[11]=712.182;
  listPedestal[23]=571.538;

  double thr[24] = {424.1, 431.9, 418.6, 640.7, 481.8, 498.3,
                    689.1, 468.7, 461.3, 459.7, 422.5, 424.6,
                    418.8, 561.8, 492.8, 429.0, 488.2, 507.4,
                    439.9, 509.1, 500.5, 471.0, 408.0, 479.1};
  //Values to be gain adjusted
  double gAdj[24] = {1650., 1700., 1850., 1650., 2100., 1650.,
                     1600., 1750., 1750., 1850., 1780., 1600.,
                     1600., 1550., 1710., 1650., 1500., 1700.,
                     1750., 1650., 1770., 1710., 1650., 1750.};
  
  TCanvas* c1=new TCanvas("ACU","ACU",3508,2480);
  TCanvas* c2=new TCanvas("ACD","ACD",3508,2480);
  c1->Divide(3,4);
  c2->Divide(3,4);
  for(int iAc = 0; iAc < 12; iAc++){
    c1->cd(iAc+1);
    h1Adc[iAc]->Draw();
    double ymax=h1Adc[iAc]->GetMaximum();
    TLine* line=new TLine(listPedestal[iAc],0,listPedestal[iAc],ymax);
    line->SetLineColor(kRed);
    line->Draw();
    gPad->SetLogy();
    TLine* lin2=new TLine(thr[iAc]+listPedestal[iAc],0,thr[iAc]+listPedestal[iAc],ymax);
    lin2->SetLineColor(kMagenta);
    //lin2->Draw();
    TLine* lin3=new TLine(gAdj[iAc],0,gAdj[iAc],ymax);
    lin3->SetLineColor(kSpring);
    //lin3->Draw();
    c2->cd(iAc+1);
    h1Adc[iAc+12]->Draw();
    ymax=h1Adc[iAc+12]->GetMaximum();
    line=new TLine(listPedestal[iAc+12],0,listPedestal[iAc+12],ymax);
    line->SetLineColor(kRed);
    line->Draw();
    lin2=new TLine(thr[iAc+12]+listPedestal[iAc+12],0,thr[iAc+12]+listPedestal[iAc+12],ymax);
    lin2->SetLineColor(kMagenta);
    //lin2->Draw();
    gPad->SetLogy();
    
  }
  c1->Write();
  c2->Write();
  c1->SaveAs("checkPed1U.png");
  c2->SaveAs("checkPed1D.png");

  return 0;
};
