#include <Det_TOF.h>

#include<cmath>

#include <fstream>
#include <iostream>
#include <TStyle.h>
using namespace std;

Det_TOF::Det_TOF(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_):Plugin(in_,out_,inf_,outf_,p_)
{
  // Set defaults for various options

};


Det_TOF::~Det_TOF()
{
};

Long_t Det_TOF::histos()
{
//trackM =  new TH1D("trackM","trackM",50,-10000, 20000);
trackBeta =  new TH1D("trackBeta","trackBeta",100,-3, 3);

for(int i=0;i<12;i++){
    string name1E("TOF Gap_");
    name1E+=(i+1)/10+'0';
    name1E+=(i+1)%10+'0';
    name1E+="Ene";
    string title1E("TOF_Gap_");
    title1E+=(i+1)/10+'0';
    title1E+=(i+1)%10+'0';
    title1E+="Ene";
  //  h1E_TOF[i] = dH1(Form("TOF_Mass%d",i+1),Form("TOF_Mass%d",i+1),200,-10000,20000);
    h1E_TOF[i] = dH1(Form("TOF_Mass%d",i+1),Form("TOF_Mass%d",i+1),50,-2,3);
    h1E_TOF[i]->SetTitleFont(62);
    h1E_TOF[i]->SetTitleSize(1.2);
    h1E_TOF[i]->GetXaxis()->SetTitle("N_{TOF}^{2}");
    h1E_TOF[i]->GetXaxis()->CenterTitle();
    h1E_TOF[i]->GetXaxis()->SetTitleSize(0.055);
    h1E_TOF[i]->SetTitleOffset(0.85,"X");
    h1E_TOF[i]->SetTitleFont(62,"X");
    h1E_TOF[i]->SetLabelFont(62,"X");
    h1E_TOF[i]->SetLabelFont(62,"Y");
    h1E_TOF[i]->SetTitleFont(62,"X");
    h1E_TOF[i]->SetTitleFont(62,"Y");
    h1E_TOF[i]->SetLabelSize(0.055,"X");
    h1E_TOF[i]->SetLabelSize(0.055,"Y");
}

  return 0;

}

Long_t Det_TOF::startup()
{
getBranchObject("testBranch",(TObject **) &ptr);

  return 0;
};

Double_t c = 3E10;
Double_t Mass_total =0.;
Long_t Det_TOF::process()
{
Int_t  GapNum =ptr->fgapNumTof2;

/*path_length = ptr->fTof2PathLocation - ptr->fTof1PathLocation;
path_time = ptr->fTof2PathTime - ptr->fTof1PathTime;
Beta= (path_length/path_time)*(1./3e2); 
Mass = pow(ptr->fTof2SP,2)*(1./pow(Beta,2) -1.);
*/
path_length[GapNum] = ptr->fTof2PathLocation - TMath::Abs(ptr->fTof1PathLocation);
//std::cout<<path_length[GapNum]<<"    "<< ptr->fTof2PathLocation <<"    " << ptr->fTof1PathLocation<<std::endl;

path_time[GapNum] = ptr->fTof2PathTime - TMath::Abs(ptr->fTof1PathTime);
std::cout<<" diff_Length:  "<<path_length[GapNum]<<"  TOF2_PL  "<< ptr->fTof2PathLocation <<"  TOF1_PL:   " << ptr->fTof1PathLocation<<"   diff_Time:  "<<path_time[GapNum]<<"  TOF2_time:  "<< ptr->fTof2PathTime <<"  TOF1_time:  " << ptr->fTof1PathTime<<std::endl;
Beta[GapNum] = TMath::Abs(path_length[GapNum]/(path_time[GapNum]*3e2)); 
Mass[GapNum] = pow(ptr->fTof2SP,2)/pow(c,2)*(1./pow(Beta[GapNum],2) -1.);
h1E_TOF[GapNum]->Fill(Beta[GapNum]);

//trackM->Fill(Mass[GapNum]);
trackBeta->Fill(TMath::Abs((ptr->fTof2PathLocation - TMath::Abs(ptr->fTof1PathLocation))/((ptr->fTof2PathTime-TMath::Abs(ptr->fTof1PathTime))*3e2)));

  return 0;
};


Long_t Det_TOF::done()
{

TCanvas* c7 = new TCanvas("c7", "", 100, 100, 800, 650);
c7->Divide(3,4);
for(int i =0;i<12;i++){
c7->cd(i+1);
h1E_TOF[i]->Draw("HIST");
}
TCanvas* c6 = new TCanvas("c6", "", 100, 100, 800, 650);
    trackBeta->Draw("HIST");
//c7->Write("Mass^{2}");
c6->Write("Beta");
  return 0;
};

Long_t Det_TOF::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};

/////////////////////////////////////////////////////////////////////////////////////////////////// 

extern "C"{
  Plugin *factory(TTree *in_, TTree *out_, TFile *inf_, TFile * outf_,TObject *p_)
  {
      return (Plugin *) new Det_TOF(in_,out_,inf_,outf_,p_);
  }
}


ClassImp(Det_TOF);
