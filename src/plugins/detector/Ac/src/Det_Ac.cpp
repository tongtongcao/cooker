//#include <bits/stdc++.h> //conflicts with preprocessor directives when compiling cooker.
                           // --> Investigation needed

#include <Det_Ac.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <TStyle.h>

using namespace std;

Det_Ac::Det_Ac(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_):Plugin(in_,out_,inf_,outf_,p_){
  // Set defaults for various options
  treeCali=0;
  tree = new TTree();
  tree = in_;
    

};


Det_Ac::~Det_Ac(){
};

Long_t Det_Ac::histos(){

trackMom =  new TH1D("trackMom","trackMom",100, 0,0.4);
  for(int i = 0; i < 12; i++){
    string name("ACu_ch");
    name+=(i+1)/10+'0';
    name+=(i+1)%10+'0';
    name+="_ac";
    string title("adc Ch for ac u");
    title+=(i+1)/10+'0';
    title+=(i+1)%10+'0';
    title+="Ch_ADC";
    h1adc_ACu[i] = dH1(name.c_str(),title.c_str(),100,0,4000);

    name=string("ACd_ch");
    name+=(i+1)/10+'0';
    name+=(i+1)%10+'0';
    name+="_ac";
    title=string("adc Ch for ac d");
    title+=(i+1)/10+'0';
    title+=(i+1)%10+'0';
    title+="Ch_ADC";
    h1adc_ACd[i] = dH1(name.c_str(),title.c_str(),100,0,4000);

    name=string("pedSubtr_U");
    name+=(i+1)/10+'0';
    name+=(i+1)%10+'0';
    name+="_ac";
    title=string("adc Ch pedestal subtr U");
    title+=(i+1)/10+'0';
    title+=(i+1)%10+'0';
    title+="Ch_ADC";
    calib_ACu[i] = new TH1D(name.c_str(),title.c_str(),100,0,4000);

    name=string("pedSubtr_D");
    name+=(i+1)/10+'0';
    name+=(i+1)%10+'0';
    name+="_ac";
    title=string("adc Ch pedestal subtr D");
    title+=(i+1)/10+'0';
    title+=(i+1)%10+'0';
    title+="Ch_ADC";
    calib_ACd[i] = new TH1D(name.c_str(),title.c_str(),100,0,4000);
  gStyle->SetOptStat(0);
  }
  return 0;
}

Long_t Det_Ac::startup(){
  getBranchObject("RawAcInfo",(TObject **) &treeRaw);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
//  getBranchObject("testBranch",(TObject **) &ptr);
  gStyle->SetOptStat(0);

  return 0;
};


Long_t Det_Ac::process(){
  //cout << " **** TESTING AC TDC VALUES : " << treeRaw->ADC_ACU[2] << "  ***** AC U TDC *** \n";
  // fill raw ac histograms for each gap
//trackMom->Fill(ptr->pVertpi0);
  for(int i = 0; i < 12; i++){
    h1adc_ACu[i]->Fill(treeRaw->ADC_ACU[i]);
    h1adc_ACd[i]->Fill(treeRaw->ADC_ACD[i]);
  gStyle->SetOptStat(0);
  }
  return 0;
};

Long_t Det_Ac::calib(){
  // average pedastals for each ac PMT
  // from ped runs
  double pedU[12] = {618.044, 665.372, 640.553, 446.387, 951.85, 622.559, 659.741, 729.028, 587.575, 727.765, 678.153, 719.818}; 
  double pedD[12] = {405.804, 497.161, 375.981, 495.873, 402.038, 575.475, 485.304, 337.158, 584.047, 492.762, 456.359, 590.692};
  double sigU[12] = {79.9772,77.692,90.0308,244.858,82.3735,92.813,183.242,89.4695,83.3293,99.9577,92.183,97.0597};
  double sigD[12] = {77.8834,219.51,99.0386,82.2189,68.9708,69.3355,66.893,68.5849,71.2415,74.3353,63.8502,75.6117};
/*
  //ped from calibration run: Run
  double ped[24] = {791.482, 789.921, 688.526, 751.938, 1043.54, 766.328, 791.82, 781.094,
                    805.547, 986.17, 849.183, 872.944, 595.823, 659.219, 578.115, 691.537,
                    581.34, 660.842, 656.739, 526.276, 698.491, 632.228, 644.929, 688.741};
*/
/*  double thr[24] = {324.1, 331.9, 318.6, 540.7, 381.8, 398.3,
                    589.1, 368.7, 361.3, 359.7, 322.5, 324.6,
                    318.8, 461.8, 392.8, 329.0, 388.2, 407.4,
                    339.9, 409.1, 400.5, 371.0, 308.0, 379.1};

*/
double thr[24] = {219.191, 236.297,245.148, 608.795,222.209,260.744,395.383,228.076,236.629,257.769,253.047,249.773,
                 214.198,594.404,238.741,224.515,201.825,214.829,202.7,210.59,213.162,207.583,194.146,218.672};	  


 
for(int i = 0; i < 12; i++){
    fgainU[i] = (treeRaw->ADC_ACU[i]); // - ped[i]);
    fgainD[i] = (treeRaw->ADC_ACD[i]); // - ped[i+12]);
    
// if(ptr->pVertpi0>0.230 && ptr->pVertpi0<0.240){
// if(treeRaw->ADC_ACU[i]>(pedU[i]+5.*sigU[i])) 
{  
 calib_ACu[i]->Fill(treeRaw->ADC_ACU[i]);
 }

//if(treeRaw->ADC_ACD[i]>(pedD[i]+5.*sigD[i])) 
{
  calib_ACd[i]->Fill(treeRaw->ADC_ACD[i]);
  }
//}
}
  return 0;
};

Long_t Det_Ac::done(){

  TCanvas* c1 = new TCanvas("c1", " AC ADC Ped Sub. Upstream ", 100, 100, 800, 650);
  c1->Divide(3,4);
  TCanvas* c2 = new TCanvas("c2", " AC ADC Ped Sub. Downstream ", 100, 100, 800, 650);
  c2->Divide(3,4);
  TCanvas* c3 = new TCanvas("c3", " AC ADC Upstream ", 100, 100, 800, 650);
  c3->Divide(3,4);
  TCanvas* c4 = new TCanvas("c4", " AC ADC Downstream ", 100, 100, 800, 650);
  c4->Divide(3,4);

  for(int i = 0; i < 12; i++){
    c1->cd(i+1);
    string htitle("AC ADC gap: ");
    htitle+=(i+1)/10+'0';
    htitle+=(i+1)%10+'0';
    htitle+="Gap";
    calib_ACu[i]->SetTitle(htitle.c_str());
    //calib_ACu[i]->SetLineColor(kMagenta-7);
    calib_ACu[i]->GetXaxis()->SetTitle("ADC Ch");
    //gPad->SetLogy();
    calib_ACu[i]->Draw();
    c2->cd(i+1);
    calib_ACd[i]->SetTitle(htitle.c_str());
   // calib_ACd[i]->SetLineColor(kMagenta-7);
    calib_ACd[i]->GetXaxis()->SetTitle("ADC Ch");
  //  gPad->SetLogy();
    calib_ACd[i]->Draw();

    c3->cd(i+1);
    h1adc_ACu[i]->SetTitle(htitle.c_str());
    h1adc_ACu[i]->SetLineColor(kMagenta-7);
    h1adc_ACu[i]->GetXaxis()->SetTitle("ADC Ch");
    gPad->SetLogy();
    h1adc_ACu[i]->Draw();
    c4->cd(i+1);
    h1adc_ACd[i]->SetTitle(htitle.c_str());
    h1adc_ACd[i]->SetLineColor(kMagenta-7);
    h1adc_ACd[i]->GetXaxis()->SetTitle("ADC Ch");
   // gPad->SetLogy();
    h1adc_ACd[i]->Draw();
//    gPad->SetLogy();
    //

}
TCanvas* c5 = new TCanvas("c5", "", 100, 100, 800, 650);
TH1D* pedSub = (TH1D*) calib_ACu[0]->Clone("pedSub");
for(int j = 1; j < 12; j++){
        pedSub->Add(calib_ACu[j]);
        pedSub->Add(calib_ACd[j]);
        pedSub->Draw("HIST");


  }
TCanvas* c6 = new TCanvas("c6", "", 100, 100, 800, 650);
 //trackMom->Draw();
  c1->Write("AC_adcUsub");
  c2->Write("AC_adcDsub");
  c3->Write("AC_adcUped");
  c4->Write("AC_adcDped");
  //c1->SaveAs("calidAC1.pdf");
  c5->Write("pedSub");
 // c6->Write("trackMome");
  return 0;
};

Long_t Det_Ac::cmdline(char *cmd){
  //add cmdline hanling here

  return 0; // 0 = all ok
};

/////////////////////////////////////////////////////////////////////////////////////////////////// 

extern "C"{
  Plugin *factory(TTree *in_, TTree *out_, TFile *inf_, TFile * outf_,TObject *p_)
  {
      return (Plugin *) new Det_Ac(in_,out_,inf_,outf_,p_);
  }
}


ClassImp(Det_Ac);
