#ifndef __TOF_M2__
#define __TOF_M2__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include <iostream>
#include "anaE36.h"
#include "cookedTof.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TVector.h"
#include <cmath>
#include <iomanip>
#include <utility>
#include <fstream>
#include "TString.h"
#include "actree.h"
#include "pgctree.h"
#include "TProfile.h"
//#include "greaterThan.h"
using namespace std;
class Tof_M2:public Plugin
{
   private:
CRTCaliPgC *treePgC;
//PgcInfo *treeRaw; 
CRTCaliAc *treeAc;
//AcInfo *treeRaw; 
//CRTCaliAc *treeAc;
  //AcInfo *treeRaw;                      /// Input tree with AC raw data
  BeamInfo* treeBeam;
  anaE36 *ptr;
  cookedTof* ptof = new cookedTof;
   public:
      Tof_M2(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
      virtual ~Tof_M2();
      // add funtions with return value Long_t here:
      Long_t histos();
      Long_t startup();
      Long_t process();
      Long_t done();
   int gap_Num;
   //double adc_ACu[12];
   //double adc_ACu[12];
   double ACu[12];
   double ACd[12];
   double Ac_U[12];
   double Ac_D[12];
 
   TH1F* h1acu_tof1[12][36];
   TH1F* h1acd_tof1[12][36];
   TH2D *h1adc_acuGap; 
   TH2D *h1adc_acGapuGap;  
   TH2D *h1adc_acGapdGap;  
   TH2D *h1adc_acGapuGap2;  
   TH2D *h1adc_acGapdGap2;  
   TH1D* h1tdc_ACu[12];
   TH1D* h1tdc_ACu_ntdc[12];
   TH1D* h1tdc_ACd_ntdc[12];
   TH1D* h1adc_ACu_ntdc[12];
   //TH1D* h1tdc_ACd_ntdc_ntt[12];
   //TH1D* h1adc_ACu_ntdc_ntt[12];
   TH1D* h1tdc_ACd[12];
   TH1D* h1adc_ACu[12];
   TH1D* h1adc_ACd[12];
   TH2D* h1adc_ACuGap[12];
   TH2D* h1adc_ACdGap[12];
   TH2D* h2adc_ACUD[12];
TH1D* h1AcAv[12];
double   AcAv[12];
    double p_max, p_min;
    double u, d, AI, AO, BI, BO;
    double u0, d0, AIo, AOo, BIo, BOo;
    double off;
    double ADC_PgC[12];
    double ADC_PgCf;
    double Ep_PgCf;
    double pgc_ped[12][7];
   TH1D *h1adc_PgC[12];
   TH1D *h1Mom_PgC[12];
   TH1D *ADC_PgCh[12][7];
   TH2D *h1Mom_ADC[12];
   TH1D *h1EP_PgC[12];
   TH2D *hBeta_gap[12];
   TH2D *hBetaGap;
   TH2D *hAdcGap;
   TH2D *hAdcPGap;
   TH2D *theta_pzx;
   TH2D *phi_pzx;
   TH2D *p_pzx;
  TH1F *h_GapNum; 
   TH1D *hEloss; 
   TH1D *EP_PgCh[12][7];
   TProfile *h1Mom_ADC_TP[12];
    TH1F *vert_Mom;
  double chi_min=0.,chi_max = 10;
  TH1F *final_PgCAdc;
  TH1F *final_PgCEp;
  TH1F *final_beta;
  TH1F *final_M2;
  TH1F *final_acAdc;
  TH1F *final_acAdcU;
  TH1F *final_acAdcD;
  TH1F *final_acAdcD_ntdc;
  TH1F *final_acAdcU_ntdc;
  TH1F *final_acAdc_ntdc;
  
 double  acAdcf_ntdc;
 double  acAdcu_ntdc;
double  acAdcd_ntdc;

   
  double mean_a[12][36];
   double mean_b[12][36];
   double peak_val[12][7];
  double maX_val[12][7];
   double  time_tof1A[20], time_tof1B[20],path_tof1[20];
    double TOF_cal_A[30][20], TOF_cal_B[30][20],fTof2PathLocationA,fTof2PathLocationB;
  double delta_corr_A[30][30];
  double delta_corr_B[30][30];
  double beta_A[30][30];
  double beta_B[30][30];
  double betaGap[12];
  
  int ppndex = 0;
  int pndex = 0;
  int index = 0;
  Int_t c = 3e8;
  double fTof2TimeMeanA,fTof2TimeMeanB;
  Double_t fTof2TimeMeanA1[20],fTof1TimeMean;
  Double_t fTof2TimeMeanB1[20];
  Double_t fTof2PathTimeB,fTof2PathTimeA;
  Double_t fTof2PathTimeA1[20],fTof2PathTimeB1[20],fTof1PathTime1A[20],fTof1PathTime1B[20];
  vector<double> tof2_p, c2_p; 
  double acAdcf=0.;
  double acAdcu=0.;
  double acAdcd=0.;
  double acAdcf1;
  double acAdcfU[12];
  double acAdcfD[12];
  double fTofbeta;
  double fTofM2;
  
  double sum_beta_A[12];
  double sum_beta_B[12];
  double sum_M2_A[12];
  double sum_M2_B[12];
  double F_beta;
  double F_M2;

  //PgC
Double_t pathlength(Double_t z0, Double_t z, Double_t n_z);
Double_t x_p(Double_t x0, Double_t n_x, Double_t t_p);
Double_t y_p(Double_t y0, Double_t n_y, Double_t t_p);
Double_t dp(Double_t p0, Double_t t_p, Double_t Z_A, Double_t I, Double_t rho);
Double_t E_loss(Double_t p0,Double_t p_loss);
double Mean_B[12][7];
double Max_B[12][7];
double peak_12[12];
double array_ADC[12][7];
double Array_ADC[12][7];
double tof2_mom[12];
double sf;
double ac_dsf[12];
double ac_usf[12];

 TString par;

int indexk; double maxk;
void findMaximum(double *num_thir, int  size, int *indexk, double *maxk);
//void findMaximum(vector<double>num, int  size, int *indexk, double *maxk);   
 virtual Long_t cmdline(char * cmd);

      ClassDef(Tof_M2,1);
};

#endif
