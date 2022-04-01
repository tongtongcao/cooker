//#include <bits/stdc++.h> //conflicts with preprocessor directives when compiling cooker.
// --> Investigation needed

#include <Det_PgC.h>
//#include </data/trek/E36/gautam/crnRoot/cooker/src/plugins/tracking/src/include/Det_tracking.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <TStyle.h>

using namespace std;


Det_PgC::Det_PgC(TTree *in_,TTree *out_,TFile *inf_, TFile * outf_,TObject *p_):Plugin(in_,out_,inf_,outf_,p_)
{
	// Set defaults for various options
	treeCali=0;
};


Det_PgC::~Det_PgC()
{

};

Long_t Det_PgC::histos()
{
	for(int i=0;i<12;i++){

		    string name1("mom");
		      name1+=(i+1)/10+'0';
		      name1+=(i+1)%10+'0';
		      name1+="_Tof2";
		      string title1("mom");
		      title1+=(i+1)/10+'0';
		      title1+=(i+1)%10+'0';
		      title1+="Tof2_Mom";
		      h1adc_tof2[i] = dH1(name1.c_str(),title1.c_str(),100,0,400);
		      

/*		for(int ii=0;ii<7;ii++){
			string name("Adc_PGC");
			name+=(i+1)/10+'0';
			name+=(i+1)%10+'0';
			name+= "_";
			name+=(ii+1)+'0';
			name+="_pgc";
			string title("Adc for pgc ");
			title+=(i+1)/10+'0';
			title+=(i+1)%10+'0';
			title+= "_";
			title+=(ii+1)+'0';
			title+="ADC";
			h1Adc_PGC[i][ii] = dH1(name.c_str() ,title.c_str() ,100,0,1500);  
		}

*/	}
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

Double_t Det_PgC::pathlength(Double_t z0, Double_t z, Double_t n_z)
  {
    Double_t pl = (z -z0)/n_z;
    return pl;
  }

Double_t Det_PgC::x_p(Double_t x0, Double_t n_x, Double_t t)
  {
    Double_t pos_x = x0 + n_x*t;
    return pos_x;
  }

Double_t Det_PgC::y_p(Double_t y0, Double_t n_y, Double_t t)
  {
    Double_t pos_y = y0 + n_y*t;
    return pos_y;
  }

Double_t Det_PgC::dp(Double_t p0, Double_t t, Double_t Z_A, Double_t I)
  {
    Double_t M = 0.511; Double_t E = sqrt(p0*p0 + M*M);
    Double_t K = 0.307075;
    Double_t dE = - t*(K/2*Z_A)*(2*log(2*M*1e6/I) + 4.*log(p0/E) - 2.);
    Double_t dp_m = E*dE*pow((E*E - M*M),(-0.5));
    return dp_m;
    }
//For tof2
double z_tofc[12] = {160.0, 161.3, 162.1, 162.0, 153.6, 153.9, 153.5, 161.1, 160.8, 160.0, 159.8, 160.1}; //z of each tof
double x_tofc[12] = {139.0,139.0,139.0,139.0,131.6,139.0,131.6,139.0,139.0,139.0,139.0,139.0}; //x of each tof
double x_tof[12], y_tof[12];
double t_tof[12], z_tofb[12];   
double pl_btof[12]; 
double x_btof[12], y_btof[12];  
double x_tof_max[12], x_tof_min[12];
Double_t xv_tof = 80.; Double_t yv_tof = 30.; Double_t zv_tof = 2.; //All Tof has same x and y size
Double_t MEE_tof = 6.87*1e-5; Double_t Z_A_tof = 0.53768;
Double_t ploss_tof_b[12]; 

  //For Deg

  Double_t z_degc[12] = {171.5, 172.2, 172.8, 172.8, 172.9, 173.0, 172.6, 172.2, 171.8, 171.5, 171.2, 171.4}; //z of each tof
  Double_t x_degc[12] = {139.0, 139.0, 139.0, 139.0, 131.6, 139.0, 131.6, 139.0, 139.0, 139.0, 139.0, 139.0}; //x of each tof
  Double_t x_deg[12],y_deg[12],t_deg[12],z_degf[12],z_degb[12],pl_fdeg[12],pl_bdeg[12];
  Double_t x_bdeg[12],y_bdeg[12],x_fdeg[12],y_fdeg[12];

  Double_t x_deg_maxf[12], x_deg_minf[12];
  Double_t x_deg_maxb[12], x_deg_minb[12];
  Double_t xv_deg = 93.; Double_t yv_deg = 41.; Double_t zv_deg =10; //All Tof has same x and y size
  Double_t MEE_deg = 5.74*1e-5; Double_t Z_A_deg = 0.57034;
  Double_t ploss_deg_b[12]; 
  Double_t ploss_deg_f[12]; 

Long_t Det_PgC::process() 
{
     Int_t  GapNum =ptr->fgapNumTof2; 
	//Tof2 
	x_tof_max[GapNum] = x_tofc[GapNum] + xv_tof/2.;
	x_tof_min[GapNum] = x_tofc[GapNum];

	z_tofb[GapNum] = z_tofc[GapNum] + zv_tof/2.; //Back z

        pl_btof[GapNum] =  pathlength(ptr->fTof2SZ/10., z_tofb[GapNum], ptr->fTof2SNz);//Pathlength back
 	x_btof[GapNum] = x_p(ptr->fTof2SX/10.,ptr->fTof2SNx, pl_btof[GapNum]);//x back 
	y_btof[GapNum] = x_p(ptr->fTof2SY/10.,ptr->fTof2SNy, pl_btof[GapNum]);//y back
//if(0.<y_btof[GapNum]<15.)cout<<Form(" y_btof[%d]: ",GapNum) <<y_btof[GapNum]<<" x0  "<<ptr->fTof2SY/10.<<endl;

//       if(x_tof_min[GapNum]<x_btof[GapNum]<x_tof_max[GapNum])cout<<"  x_ftof:  "<<x_btof[GapNum]<<"  y_ftof "<<y_btof[GapNum]<<endl;
//cout<<Form(" x_btof[%d]: ",GapNum) <<x_btof[GapNum]<<" x0  "<<ptr->fTof2SX/10.<<endl;
       if(x_tof_min[GapNum]<x_btof[GapNum] && x_btof[GapNum]<x_tof_max[GapNum] && 0.<y_btof[GapNum] && y_btof[GapNum]<15. ){

//cout<<"  GapNum:  "<<GapNum<<"  x_btof    "<<x_btof[GapNum]<<"  y_btof    "<<y_btof[GapNum]<<endl;
	ploss_tof_b[GapNum] = dp(ptr->fTof2SP, pl_btof[GapNum], Z_A_tof, MEE_tof);// front plane
cout<<"  GapNum:  "<<GapNum<<"    ploss_tof_b:  "<<ploss_tof_b[GapNum]<<"  x_btof    "<<x_btof[GapNum]<<"  y_btof    "<<y_btof[GapNum]<<endl;
	h1adc_tof2[GapNum]->Fill(ptr->fTof2SP - ploss_tof_b[GapNum]);
	}


/*	//Tof2 
	x_deg_maxf[GapNum] = x_degc[GapNum];
	x_deg_minf[GapNum] = x_degc[GapNum] - xv_deg/2.;

	x_deg_maxb[GapNum] = x_degc[GapNum] + xv_deg/2;
	x_deg_minb[GapNum] = x_degc[GapNum];


	z_degb[GapNum] = z_degc[GapNum] + zv_deg/2.; //Back z
	z_degf[GapNum] = z_degc[GapNum] - zv_deg/2.; //Back z

        pl_bdeg[GapNum] =  pathlength(ptr->fTof2SZ/10. + z_degc[GapNum] - z_tofc[GapNum], z_degb[GapNum], ptr->fTof2SNz);//Pathlength back
        pl_fdeg[GapNum] =  pathlength(ptr->fTof2SZ/10. + z_degc[GapNum] - z_tofc[GapNum], z_degf[GapNum], ptr->fTof2SNz);//Pathlength back
 
	x_bdeg[GapNum] = x_p(ptr->fTof2SX/10. + x_degc[GapNum] - x_tofc[GapNum], pl_bdeg[GapNum]);//x back 
	y_bdeg[GapNum] = x_p(ptr->fTof2SY/10., ptr->fTof2SNy,pl_bdeg[GapNum]);//y back

	x_fdeg[GapNum] = x_p(ptr->fTof2SX/10. + x_degc[GapNum] - x_tofc[GapNum], ptr->fTof2SNx,  pl_fdeg[GapNum]);//x back 
	y_fdeg[GapNum] = x_p(ptr->fTof2SY/10., ptr->fTof2SNy, pl_fdeg[GapNum]);//y back
	
        if(x_deg_minb[GapNum]<x_bdeg[GapNum]<x_deg_maxb[GapNum] && 0.<y_bdeg[GapNum]<yv_deg/2.){
	ploss_deg_b[GapNum] = dp(ptr->fTof2SP - ploss_tof_b[GapNum], pl_bdeg[GapNum], Z_A_deg, MEE_deg);// front plane
//cout<<"   xmin  "<<x_deg_minb[GapNum]<<  " xmax     "<<x_deg_maxb[GapNum]<<"  x_fdeg  "<<x_bdeg[GapNum]<<"   ymin   "<<0.<<"   ymax   "<<yv_deg/2<<"  y_fdeg   "<<y_bdeg[GapNum]<<endl;
	}

        if(x_deg_minf[GapNum]<x_fdeg[GapNum]<x_deg_maxf[GapNum] && -yv_deg/2.<y_fdeg[GapNum]<0.){
	ploss_deg_f[GapNum] = dp(ptr->fTof2SP - ploss_tof_b[GapNum] - ploss_deg_b[GapNum], pl_fdeg[GapNum], Z_A_deg, MEE_deg);// front plane
	}
//cout<<"   xmin  "<<x_deg_minf[GapNum]<<  " xmax     "<<x_deg_maxf[GapNum]<<"  x_fdeg  "<<x_fdeg[GapNum]<<"   ymin   "<<-yv_deg/2.<<"   ymax   "<<0.<<"  y_fdeg   "<<y_fdeg[GapNum]<<endl;
// cout<<"  GapNum:  "<<GapNum<<"    ploss_deg_b:  "<<ploss_deg_b[GapNum]<<"  ploss_deg_f    "<<ploss_deg_f[GapNum]<<endl;

*/	

	return 0;
};


Long_t Det_PgC::done()
{
	/* gStyle->SetOptFit(1);
	   Double_t x_max[12][7];
	   for(int i = 0; i < 12; i++){
	   c1[i] = new TCanvas(Form("c1[%d]",i));
	   c1[i]->Divide(4,2);
	   for(int ii = 0; ii < 7; ii++){
	   c1[i]->cd(ii+1);
	   string htitle("PGC ADC gap: ");
	   htitle+=(i+1)/10+'0';
	   htitle+=(i+1)%10+'0';
	   htitle+= "_";
	   htitle+=(ii+1)+'0';
	   htitle+="Gap";
	   h1Adc_PGC[i][ii]->SetTitle(htitle.c_str());
	   int binmax = h1Adc_PGC[i][ii]->GetMaximumBin();  
	   x_max[i][ii] = h1Adc_PGC[i][ii]->GetXaxis()->GetBinCenter(binmax);
	//    h1Adc_PGC[i][ii]->SetAxisRange(x_max[i][ii]-100, x_max[i][ii]+100,"X");
	if(i==3 && ii==7)
	{
	h1Adc_PGC[i][ii]->SetAxisRange(x_max[i][ii]-60, x_max[i][ii]+100,"X");
	}
	h1Adc_PGC[i][ii]->GetXaxis()->SetTitle("ADC Ch");
	h1Adc_PGC[i][ii]->Draw("HIST");
	//  h1Adc_PGC[i][ii]->Fit("gaus","","",x_max[i][ii]-30,x_max[i][ii]+30);
	gPad->SetLogy();
	}
	c1[i]->Write("PGC");
	c1[i]->SaveAs(Form("output_root/PGC_plots/pcg_Gap_%d.pdf",i));

	}*/
	return 0;
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
