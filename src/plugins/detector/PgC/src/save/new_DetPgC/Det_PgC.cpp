#include <Det_PgC.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <TStyle.h>
#include <TH1D.h>
#include <TMath.h>
using namespace std;


Det_PgC::Det_PgC(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_):Plugin(in_,out_,inf_,outf_,p_)
{
  treeCali=0;
};


Det_PgC::~Det_PgC()
{

};

Long_t Det_PgC::histos()
{
  for(int i=0;i<12;i++){
    string name1E("PgC Gap_");
    name1E+=(i+1)/10+'0';
    name1E+=(i+1)%10+'0';
    name1E+="Ene";
    string title1E("PgC_Gap_");
    title1E+=(i+1)/10+'0';
    title1E+=(i+1)%10+'0';
    title1E+="Ene";
    h1E_PgC[i] = dH1(Form("PgC_Gap_%d",i+1),Form("PgC_Gap_%d",i+1),50,0,400);
    h1E_PgC[i]->SetTitleFont(62);
    h1E_PgC[i]->SetTitleSize(1.2);
    h1E_PgC[i]->GetXaxis()->SetTitle("Energy (MeV)");
    h1E_PgC[i]->GetXaxis()->CenterTitle();
    h1E_PgC[i]->GetXaxis()->SetTitleSize(0.055);
    h1E_PgC[i]->SetTitleOffset(0.85,"X");
    h1E_PgC[i]->SetTitleFont(62,"X");
    h1E_PgC[i]->SetLabelFont(62,"X");
    h1E_PgC[i]->SetLabelFont(62,"Y");
    h1E_PgC[i]->SetTitleFont(62,"X");
    h1E_PgC[i]->SetTitleFont(62,"Y");
    h1E_PgC[i]->SetLabelSize(0.055,"X");
    h1E_PgC[i]->SetLabelSize(0.055,"Y");


    h1Mom_PgC[i] = dH1(Form("PgC_Mom_Gap_%d ",i+1),Form("PgC_Mom_Gap_%d ",i+1),50,0,400);
    h1Mom_PgC[i]->SetTitleFont(62);
    h1Mom_PgC[i]->SetTitleSize(1.2);
    h1Mom_PgC[i]->GetXaxis()->SetTitle("Track Momentum (MeV)");
    h1Mom_PgC[i]->GetXaxis()->CenterTitle();
    h1Mom_PgC[i]->GetXaxis()->SetTitleSize(0.075);
    h1Mom_PgC[i]->SetTitleOffset(0.6,"X");
    h1Mom_PgC[i]->SetTitleFont(62,"X");
    h1Mom_PgC[i]->SetLabelFont(62,"X");
    h1Mom_PgC[i]->SetLabelFont(62,"Y");
    h1Mom_PgC[i]->SetTitleFont(62,"X");
    h1Mom_PgC[i]->SetTitleFont(62,"Y");
    h1Mom_PgC[i]->SetLabelSize(0.055,"X");
    h1Mom_PgC[i]->SetLabelSize(0.055,"Y");



    string name("PgC_ch");
    name+=(i+1)/10+'0';
    name+=(i+1)%10+'0';
    name+="_PgC";
    string title("ADC");
    title+=(i+1)/10+'0';
    title+=(i+1)%10+'0';
    h1adc_PgC[i] = dH1(name.c_str(),title.c_str(),100,0,1500);
    h1adc_PgC[i]->SetTitleFont(62);
    h1adc_PgC[i]->SetTitleSize(1.2);
    h1adc_PgC[i]->GetXaxis()->SetTitle("ADC Channel");
    h1adc_PgC[i]->GetXaxis()->CenterTitle();
    h1adc_PgC[i]->GetXaxis()->SetTitleSize(0.075);
    h1adc_PgC[i]->SetTitleOffset(0.6,"X");
    h1adc_PgC[i]->SetTitleFont(62,"X");
    h1adc_PgC[i]->SetLabelFont(62,"X");
    h1adc_PgC[i]->SetLabelFont(62,"Y");
    h1adc_PgC[i]->SetTitleFont(62,"X");
    h1adc_PgC[i]->SetTitleFont(62,"Y");
    h1adc_PgC[i]->SetLabelSize(0.055,"X");
    h1adc_PgC[i]->SetLabelSize(0.055,"Y");

    string n1("PgC Gap_");
    n1+=(i+1)/10+'0';
    n1+=(i+1)%10+'0';
    n1+="mom";
    string t1("PgC Gap_");
    t1+=(i+1)/10+'0';
    t1+=(i+1)%10+'0';
    t1+="Mom";
    h1Mom_ADC[i] = dH2(n1.c_str(),t1.c_str(),100,0,1500,50,0,400);
    h1Mom_ADC_TP[i] = new TProfile(n1.c_str(),t1.c_str(),100,0,1500,0,400);
    h1Mom_ADC_TP[i]->SetLineColor(kRed);
    h1Mom_ADC_TP[i]->SetMarkerColor(kRed);
    h1Mom_ADC[i]->SetTitleFont(62);
    h1Mom_ADC[i]->SetTitleSize(1.2);
    h1Mom_ADC[i]->GetYaxis()->SetTitle("Momentum (MeV)");
    h1Mom_ADC[i]->GetYaxis()->CenterTitle();
    h1Mom_ADC[i]->GetYaxis()->SetTitleSize(0.05);
    h1Mom_ADC[i]->GetXaxis()->SetTitle("ADC channel");
    h1Mom_ADC[i]->GetXaxis()->CenterTitle();
    h1Mom_ADC[i]->GetXaxis()->SetTitleSize(0.05);
    h1Mom_ADC[i]->SetTitleOffset(0.85,"X");
    h1Mom_ADC[i]->SetTitleOffset(0.85,"Y");
    h1Mom_ADC[i]->SetTitleFont(62,"X");
    h1Mom_ADC[i]->SetLabelFont(62,"X");
    h1Mom_ADC[i]->SetLabelFont(62,"Y");
    h1Mom_ADC[i]->SetTitleFont(62,"X");
    h1Mom_ADC[i]->SetTitleFont(62,"Y");
    h1Mom_ADC[i]->SetLabelSize(0.04,"X");
    h1Mom_ADC[i]->SetLabelSize(0.04,"Y");


    h1Mom_ADC_TP[i] = new TProfile(n1.c_str(),t1.c_str(),100,0,1500,0,400);
    h1Mom_ADC_TP[i]->SetLineColor(kRed);
    h1Mom_ADC_TP[i]->SetMarkerColor(kRed);
    h1Mom_ADC_TP[i]->SetTitleFont(62);
    h1Mom_ADC_TP[i]->SetTitleSize(1.2);
    h1Mom_ADC_TP[i]->GetYaxis()->SetTitle("Momentum (MeV)");
    h1Mom_ADC_TP[i]->GetYaxis()->CenterTitle();
    h1Mom_ADC_TP[i]->GetYaxis()->SetTitleSize(0.05);
    h1Mom_ADC_TP[i]->GetXaxis()->SetTitle("ADC channel");
    h1Mom_ADC_TP[i]->GetXaxis()->CenterTitle();
    h1Mom_ADC_TP[i]->GetXaxis()->SetTitleSize(0.05);
    h1Mom_ADC_TP[i]->SetTitleOffset(0.85,"X");
    h1Mom_ADC_TP[i]->SetTitleOffset(0.85,"Y");
    h1Mom_ADC_TP[i]->SetTitleFont(62,"X");
    h1Mom_ADC_TP[i]->SetLabelFont(62,"X");
    h1Mom_ADC_TP[i]->SetLabelFont(62,"Y");
    h1Mom_ADC_TP[i]->SetTitleFont(62,"X");
    h1Mom_ADC_TP[i]->SetTitleFont(62,"Y");
    h1Mom_ADC_TP[i]->SetLabelSize(0.04,"X");
    h1Mom_ADC_TP[i]->SetLabelSize(0.04,"Y");




    for(int ii=0;ii<7;ii++){
      string namep("Adc_PGC");
      namep+=(i+1)/10+'0';
      namep+=(i+1)%10+'0';
      namep+= "_";
      namep+=(ii+1)+'0';
      namep+="_pgc";

      string titlep("PgC Gap_");
      titlep+=(i+1)/10+'0';
      titlep+=(i+1)%10+'0';
      titlep+= "_";
      titlep+=(ii+1)+'0';
      titlep+="";

      ADC_PgC[i][ii]= dH1(namep.c_str(),titlep.c_str(),100,0,1500);
      ADC_PgC[i][ii]->SetTitleFont(62);
      ADC_PgC[i][ii]->SetTitleSize(1.2);
      ADC_PgC[i][ii]->GetXaxis()->SetTitle("ADC Channel");
      ADC_PgC[i][ii]->GetXaxis()->CenterTitle();
      ADC_PgC[i][ii]->GetXaxis()->SetTitleSize(0.075);
      ADC_PgC[i][ii]->SetTitleOffset(0.6,"X");
      ADC_PgC[i][ii]->SetTitleFont(62,"X");
      ADC_PgC[i][ii]->SetLabelFont(62,"X");
      ADC_PgC[i][ii]->SetLabelFont(62,"Y");
      ADC_PgC[i][ii]->SetTitleFont(62,"X");
      ADC_PgC[i][ii]->SetTitleFont(62,"Y");
      ADC_PgC[i][ii]->SetLabelSize(0.055,"X");
      ADC_PgC[i][ii]->SetLabelSize(0.055,"Y");

    }

  }


  return 0;

}

Long_t Det_PgC::startup()
{
  getBranchObject("RawPgcInfo",(TObject **) &treeRaw);
  getBranchObject("RawBeamInfo",(TObject **) &treeBeam);
  getBranchObject("testBranch",(TObject **) &ptr);
  gStyle->SetOptStat(0);
  return 0;
};


//Double_t M = 0.511;
//PgC Calibraiton
Double_t pedestal[12][7] = { {479.2008, 569.6043, 186.0088, 314.5900, 480.9254, 489.9421, 165.9059},
			     {400.0266, 520.2240, 543.4152, 514.8561, 442.4048, 462.0687, 618.9130},
			     {219.8563, 152.1412, 208.5225, 383.0596, 440.1627, 352.4358, 505.0797},
			     {486.3342, 439.8117, 449.0930, 385.7619, 433.0558, 460.7424, 420.1256},
			     {519.7247, 577.0609, 405.0830, 489.3602, 709.1587, 471.5588, 568.0758},
			     {721.1133, 846.2933, 790.2884, 709.4451, 514.9213, 590.0196, 703.7817},
			     {615.7067, 564.9588, 264.9764, 584.1545, 694.1780, 419.1175, 276.1318},
			     {661.0245, 516.5154, 503.4639, 465.1396, 479.2645, 454.9846, 487.2282},
			     {295.9655, 243.4812, 339.5547, 315.1526, 365.1006, 198.7629, 250.5476},
			     {355.9769, 328.0053, 318.4344, 286.9747, 237.4743, 284.1670, 215.4730},
			     {427.6373, 263.0033, 329.6864, 381.3515, 270.2102, 343.0392, 359.3404},
			     {498.0114, 570.0067, 368.9883, 447.8834, 432.3074, 452.4020, 307.2901} };

//For tof2
double z_tofc[12] = {160.0, 161.3, 162.1, 162.0, 153.6, 153.9, 153.5, 161.1, 160.8, 160.0, 159.8, 160.1}; //z of each tof
double x_tofc[12] = {139.0,139.0,139.0,139.0,131.6,139.0,131.6,139.0,139.0,139.0,139.0,139.0}; //x of each tof
double x_tof[12], y_tof[12];
double t_tof[12], z_tofb[12];   
double pl_btof[12]; 
double x_btof[12], y_btof[12];  
double x_tof_max[12], x_tof_min[12];
Double_t xv_tof = 80.; Double_t yv_tof = 30.; Double_t zv_tof = 2.; //All Tof has same x and y size
Double_t MEE_tof = 68.7; Double_t Z_A_tof = 0.53768;
Double_t ploss_tof_b[12]; 
Double_t rho_tof = 1.;

//For Deg

Double_t z_degc[12] = {171.5, 172.2, 172.8, 172.8, 172.9, 173.0, 172.6, 172.2, 171.8, 171.5, 171.2, 171.4}; //z of each tof
Double_t x_degc[12] = {139.0, 139.0, 139.0, 139.0, 131.6, 139.0, 131.6, 139.0, 139.0, 139.0, 139.0, 139.0}; //x of each tof
Double_t x_deg[12],y_deg[12],z_degf[12],z_degb[12],pl_fdeg[12],pl_bdeg[12];
Double_t x_bdeg[12],y_bdeg[12],x_fdeg[12],y_fdeg[12];
Double_t x_deg_max[12], x_deg_min[12];
Double_t xv_deg = 93.; Double_t yv_deg = 41.; Double_t zv_deg =10.; //All Tof has same x and y size
Double_t MEE_deg = 57.4; Double_t Z_A_deg = 0.57034;
Double_t ploss_deg_b[12]; 
Double_t ploss_deg_f[12]; 
 Double_t rho_deg =0.89;


Double_t E_PgC[12];
Double_t P_PgC[12];
Double_t M = 105.66;



Double_t Det_PgC::pathlength(Double_t z0, Double_t z, Double_t n_z)
{
  Double_t pl = (z -z0)/n_z;
  return pl;
}

Double_t Det_PgC::x_p(Double_t x0, Double_t n_x, Double_t t_p)
{
  Double_t pos_x = x0 + n_x*t_p;
  return pos_x;
}

Double_t Det_PgC::y_p(Double_t y0, Double_t n_y, Double_t t_p)
{
  Double_t pos_y = y0 + n_y*t_p;
  return pos_y;
}

Double_t Det_PgC::dp(Double_t p0, Double_t t_p, Double_t Z_A, Double_t I, Double_t rho)
{
 Double_t E = sqrt(p0*p0 + M*M);
  Double_t K = 0.307075;
  Double_t dE = - rho*t_p*(K/2*Z_A)*(2.*log(2.*M*pow(10,6)/I) + 4.*log(p0/E) - 2.);
  Double_t dp_m = E*dE*pow((E*E - M*M),(-0.5));

  return dp_m;
}

Double_t Det_PgC::E_loss(Double_t p0, Double_t p_loss){
Double_t dE_loss = p0/sqrt(p0*p0+M*M) * p_loss;
}

Double_t E_loss_tof2[12], E_loss_degf[12], E_loss_degb[12];

Long_t Det_PgC::process() 
{
  Int_t  GapNum =ptr->fgapNumTof2; 
  //Tof2 
  x_tof_max[GapNum] = x_tofc[GapNum] + xv_tof/2.;
  x_tof_min[GapNum] = x_tofc[GapNum] - xv_tof/2.;
  z_tofb[GapNum] = z_tofc[GapNum] + zv_tof/2.; //Back z

  pl_btof[GapNum] =  pathlength(ptr->fTof2SZ/10., z_tofb[GapNum], ptr->fTof2SNz);//Pathlength back
  x_btof[GapNum] = x_p(ptr->fTof2SX/10.,ptr->fTof2SNx, pl_btof[GapNum]);//x back 
  y_btof[GapNum] = y_p(ptr->fTof2SY/10.,ptr->fTof2SNy, pl_btof[GapNum]);//y back
  //cout<<" pl_btof:  "<<pl_btof[GapNum]<<endl;
 if(x_tof_min[GapNum]<x_btof[GapNum] && x_btof[GapNum]<x_tof_max[GapNum] && -yv_tof/2<y_btof[GapNum] && y_btof[GapNum]<yv_deg/2 )
    {
      ploss_tof_b[GapNum] = dp(ptr->fTof2SP, pl_btof[GapNum], Z_A_tof, MEE_tof, rho_tof);// front plane
    }

  //deg
  x_deg_max[GapNum] = x_degc[GapNum] + xv_deg/2.;
  x_deg_min[GapNum] = x_degc[GapNum] - xv_deg/2.;

  z_degb[GapNum] = z_degc[GapNum] + zv_deg/2.; //Back z
  z_degf[GapNum] = z_degc[GapNum] - zv_deg/2.; //Back z

  pl_bdeg[GapNum] =  pathlength(ptr->fTof2SZ/10. + z_degc[GapNum] - z_tofc[GapNum], z_degb[GapNum], ptr->fTof2SNz);//Pathlength back
  pl_fdeg[GapNum] =  pathlength(ptr->fTof2SZ/10. + z_degc[GapNum] - z_tofc[GapNum], z_degf[GapNum], ptr->fTof2SNz);//Pathlength back
//cout<<" pl_bdeg:  "<<pl_bdeg[GapNum] << "  pl_fdeg:  "<<pl_fdeg[GapNum]<<endl;
  x_bdeg[GapNum] = x_p(ptr->fTof2SX/10., ptr->fTof2SNx, pl_bdeg[GapNum]);//x back 
  y_bdeg[GapNum] = y_p(ptr->fTof2SY/10., ptr->fTof2SNy, pl_bdeg[GapNum]);//y back

  x_fdeg[GapNum] = x_p(ptr->fTof2SX/10., ptr->fTof2SNx, pl_fdeg[GapNum]);//x back 
  y_fdeg[GapNum] = y_p(ptr->fTof2SY/10., ptr->fTof2SNy, pl_fdeg[GapNum]);//y back

//cout<<"  x_fdeg:   "<<x_fdeg[GapNum]<<"   x_btof:   "<<ptr->fTof2SX/10.<< "   y_fdeg:  "<<y_fdeg[GapNum]<< " y_tof:  "<<  ptr->fTof2SY/10.<<endl;

  if(x_deg_min[GapNum]<x_fdeg[GapNum]  && x_fdeg[GapNum]<x_deg_max[GapNum] && -yv_deg/2.<y_fdeg[GapNum] && y_fdeg[GapNum]<yv_deg/2.)
    {
      ploss_deg_f[GapNum] = dp(ptr->fTof2SP - TMath::Abs(ploss_tof_b[GapNum]), pl_fdeg[GapNum], Z_A_deg, MEE_deg,rho_deg);// front plane
    }
  if(x_deg_min[GapNum]<x_bdeg[GapNum] && x_bdeg[GapNum]<x_deg_max[GapNum] && -yv_deg/2<y_bdeg[GapNum] && y_bdeg[GapNum]<yv_deg/2.)
    {
      ploss_deg_b[GapNum] = dp(ptr->fTof2SP - TMath::Abs(ploss_tof_b[GapNum]) - TMath::Abs(ploss_deg_f[GapNum]), pl_bdeg[GapNum], Z_A_deg, MEE_deg, rho_deg);// front plane
    }
ploss_tof_b[GapNum] = fabs(ploss_tof_b[GapNum]);
ploss_deg_b[GapNum]= fabs(ploss_deg_b[GapNum]);
ploss_deg_f[GapNum] = fabs(ploss_deg_f[GapNum]);
  //deg

  P_PgC[GapNum] = ptr->fTof2SP - TMath::Abs(ploss_tof_b[GapNum]) - 2.*TMath::Abs(ploss_deg_f[GapNum]);
  E_PgC[GapNum] = sqrt(P_PgC[GapNum]*P_PgC[GapNum] + M*M);

/*E_loss_tof2[GapNum] = E_loss(ptr->fTof2SP, ploss_tof_b[GapNum]);
E_loss_degf[GapNum] = E_loss(ptr->fTof2SP - ploss_tof_b[GapNum], ploss_deg_f[GapNum]);
E_loss_degb[GapNum] = E_loss(ptr->fTof2SP - ploss_tof_b[GapNum] - ploss_deg_f[GapNum], ploss_deg_b[GapNum]);
  */

E_loss_tof2[GapNum] = ptr->fTof2SP/sqrt(ptr->fTof2SP*ptr->fTof2SP + M*M)*ploss_tof_b[GapNum];
E_loss_degf[GapNum] = (ptr->fTof2SP - ploss_tof_b[GapNum])/sqrt((ptr->fTof2SP - ploss_tof_b[GapNum])*(ptr->fTof2SP -ploss_tof_b[GapNum]) + M*M)*ploss_deg_f[GapNum];
E_loss_degb[GapNum] = (ptr->fTof2SP - ploss_tof_b[GapNum] - ploss_deg_f[GapNum])/sqrt((ptr->fTof2SP - ploss_tof_b[GapNum] - ploss_deg_f[GapNum])*(ptr->fTof2SP -ploss_tof_b[GapNum]- ploss_deg_f[GapNum]) + M*M)*ploss_deg_b[GapNum];



cout<<"tof_b  "<<E_loss_tof2[GapNum]<<"  deg_f: "<<E_loss_degf[GapNum]<<" deg_b   "<<E_loss_degb[GapNum]<<endl;

  Double_t E_cal[12];  
//  cout<<"tof_b  "<<ploss_tof_b[GapNum]<<"  deg_f: "<<ploss_deg_f[GapNum]<<" deg_b   "<<ploss_deg_b[GapNum]<<endl;
  h1Mom_PgC[GapNum]->Fill(P_PgC[GapNum]);
  h1E_PgC[GapNum]->Fill(E_PgC[GapNum]);
  for(int i=0;i<12;i++){
    for(int ii=0;ii<7;ii++){
      ADC_PgC[i][ii]->Fill(treeRaw->ADC_PGC_Gap[i][ii] - pedestal[i][ii]); 
      if((treeRaw->ADC_PGC_Gap[i][ii] - pedestal[i][ii])>0.)
	{
	  h1adc_PgC[i]->Add(ADC_PgC[i][ii]);
         E_cal[i] += treeRaw->ADC_PGC_Gap[i][ii] - pedestal[i][ii];
	}
    }
      h1Mom_ADC[i]->Fill(E_cal[i],P_PgC[i]);
      h1Mom_ADC_TP[i]->Fill(E_cal[i],P_PgC[i]);
  }
  return 0;
};


Long_t Det_PgC::done()
{

  TCanvas* c2 = new TCanvas("c2", " PGC Mom", 100, 100, 800, 650);
  c2->Divide(3,4);
  TCanvas* c3 = new TCanvas("c3", " PGC Ene ", 100, 100, 800, 650);
  c3->Divide(3,4);
  TCanvas* c4 = new TCanvas("c4", " Ene Adc sum", 100, 100, 800, 650);
  c4->Divide(3,4);
  TCanvas* c5 = new TCanvas("c5", " Ene Adc sum", 100, 100, 800, 650);
  c5->Divide(3,4);
  TCanvas* c6 = new TCanvas("c6", " Ene Adc sum", 100, 100, 800, 650);
  c6->Divide(3,4);
  for(int i = 0; i < 12; i++){
    c1[i] = new TCanvas(Form("c1[%d]",i));
    c1[i]->Divide(4,2);
    for(int ii = 0; ii < 7; ii++){
      c1[i]->cd(ii+1);
      ADC_PgC[i][ii]->Draw("HIST");
      gPad->SetLogy();
      c2->cd(i+1);
      h1Mom_PgC[i]->Draw("HIST");
   
      c3->cd(i+1);
      h1E_PgC[i]->Draw("HIST");
   
      c4->cd(i+1);
      h1adc_PgC[i]->Draw("HIST");
  
      gPad->SetLogy(); 
      c5->cd(i+1);
      h1Mom_ADC[i]->Draw("cloz");
      h1Mom_ADC_TP[i]->Draw("same"); 
      c6->cd(i+1);
      h1Mom_ADC_TP[i]->Draw(); 
    }
   // c1[i]->Write("PGC");
    //c1[i]->SaveAs(Form("output_root/PGC_plots/pcg_Gap_%d.pdf",i));
  }
  c2->Write("Mom");
  //c3->Write("Ene");
// c5->Write("ADC_Mom");
// c6->Write("ADC_Mom");
  /*c2->Print("TOF_MWPC/PGC_3994/MomentumE_3994.pdf[");
  c2->Print("TOF_MWPC/PGC_3994/MomentumE_3994.pdf");
  c3->Print("TOF_MWPC/PGC_3994/MomentumE_3994.pdf");
  c4->Print("TOF_MWPC/PGC_3994/MomentumE_3994.pdf");
  c5->Print("TOF_MWPC/PGC_3994/MomentumE_3994.pdf");
  c5->Print("TOF_MWPC/PGC_3994/MomentumE_3994.pdf]");
  c3->Print("TOF_MWPC/PGC_3994/Energy_3994.pdf");
  */return 0;
};

Long_t Det_PgC::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};

/////////////////////////////////////////////////////////////////////////////////////////////////// 

extern "C"{
  Plugin *factory(TTree *in_, TTree *out_, TFile *inf_, TFile * outf_,TObject *p_)
  {
    return (Plugin *) new Det_PgC(in_,out_,inf_,outf_,p_);
  }
}


ClassImp(Det_PgC);
