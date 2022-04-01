#include <Det_MWPC.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <TStyle.h>
#include <sstream>
#include "TLine.h"
#include <string>
#include "TH1D.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TFormula.h"
#include "TAxis.h"
#include "TGaxis.h"
using namespace std;

static TH1D *hist;
static TH2D *histPVsZ;
static TH2D *histPVsTheta;
Long_t Det_MWPC::histos_trackV2(){
  ostringstream sstr_name,sstr_title;
  histPVsTheta=new TH2D("Momentum vs theta","Momentum vs theta;#theta in zx plane;p (GeV/c)",100,0,3.14,100,0,300);
  histPVsZ=new TH2D("Momentum vs z","Momentum vs z;z (mm);p (GeV/c)",100,-100,100,100,0,300);
  hist=new TH1D("Momentum Distribution","Momentum Distribution;p (GeV/c);Counts",100,0,300);
  return 0;
}

Long_t Det_MWPC::startup_trackV2(){
  getBranchObject("CookMwpcInfoV1",(TObject **) &treeCookV1);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);

  return 0;
}

static const double p1Z=380;
static const double p2X=725.01;
static double p1X;
static double a1,b1,a1P,b1P;
static double c1,c2,c3,c4;
static double a,b,c,d;
static double alpha,beta,delta;
static double p0Z,p0X,p2Z,p0Z1,p0Z2,p0Z3,p0X1,p0X2,p0X3,p2Z1,p2Z2,p2Z3;
static double r,a2,b2,r1,a21,b21,r2,a22,b22,r3,a23,b23;
static double aCond,bCond,aCond1,bCond1,aCond2,bCond2,aCond3,bCond3; // Center of circle is located at the right of the strainge line which passes through p1 and p2
static vector<double> par;
static vector<vector<double> > vect_par;
static vector<double> point;
static vector<vector<double> > vect_point;
static int nTrack;
static vector<int> vect_nTrack;

static double s1,s0,s2;
static double pL,pT;

static int numTrack=0;
Long_t Det_MWPC::process_trackV2(){
//cout<<treeCookV1->event<<"    "<<treeCookV1->posC2Y<<"    "<<treeCookV1->posC2Z<<"    "<<treeCookV1->posC3X<<"    "<<treeCookV1->posC3Y<<"    "<<treeCookV1->posC4X<<"    "<<treeCookV1->posC4Y<<endl;
  if(treeCookV1->posC2Y!=-10000 && treeCookV1->posC3Y!=-10000 && treeCookV1->posC4Y!=-10000){
//cout<<treeCookV1->posC2X<<"    "<<treeCookV1->posC2Y<<"    "<<treeCookV1->posC2Z<<"    "<<treeCookV1->posC3X<<"    "<<treeCookV1->posC3Y<<"    "<<treeCookV1->posC3Z<<"    "<<treeCookV1->posC4X<<"    "<<treeCookV1->posC4Y<<"    "<<treeCookV1->posC4Z<<endl;
    a1=(treeCookV1->posC3X-treeCookV1->posC4X)/(treeCookV1->posC3Z-treeCookV1->posC4Z);
    b1=treeCookV1->posC3X-a1*treeCookV1->posC3Z;
    p1X=a1*p1Z+b1;
    a1P=-1/a1;
    b1P=p1X-a1P*p1Z;


    c1=p1X*p1X+p1Z*p1Z-p2X*p2X;
    c2=treeCookV1->posC2X*p2X-p2X*p2X;
    c3=(2*p1X-p2X-treeCookV1->posC2X)*a1P+2*p1Z-treeCookV1->posC2Z;
    c4=(2*p1X-p2X-treeCookV1->posC2X)*b1P-c1+c2;

    a=-c3+(treeCookV1->posC2X-p2X)*a1P+treeCookV1->posC2Z;
    b=c3*c3-c4+(treeCookV1->posC2X-p2X)*(b1P-2*a1P*treeCookV1->posC2Z)-2*treeCookV1->posC2Z*treeCookV1->posC2Z-c2;
    c=2*c3*c4+c3*treeCookV1->posC2Z*treeCookV1->posC2Z+(treeCookV1->posC2X-p2X)*(a1P*treeCookV1->posC2Z*treeCookV1->posC2Z-2*b1P*treeCookV1->posC2Z)+treeCookV1->posC2Z*treeCookV1->posC2Z*treeCookV1->posC2Z+2*c2*treeCookV1->posC2Z;
    d=c4*c4+c4*treeCookV1->posC2Z*treeCookV1->posC2Z+(treeCookV1->posC2X-p2X)*b1P*treeCookV1->posC2Z*treeCookV1->posC2Z-c2*treeCookV1->posC2Z*treeCookV1->posC2Z;

    alpha=b*c/6/a/a-b*b*b/27/a/a/a-d/2/a;
    beta=c/3/a-b*b/9/a/a;
    delta=alpha*alpha+beta*beta*beta;


    nTrack=0;

    if(delta>=0){ 
      if(alpha+sqrt(delta)>=0 && alpha-sqrt(delta)>=0)
        p0Z=-b/3/a+pow(alpha+sqrt(delta),1./3)+pow(alpha-sqrt(delta),1./3);
      else if(alpha+sqrt(delta)>=0 && alpha-sqrt(delta)<0)
        p0Z=-b/3/a+pow(alpha+sqrt(delta),1./3)-pow(-(alpha-sqrt(delta)),1./3);
      else if(alpha+sqrt(delta)<0 && alpha-sqrt(delta)>=0)
        p0Z=-b/3/a-pow(-(alpha+sqrt(delta)),1./3)+pow(alpha-sqrt(delta),1./3);
      else
        p0Z=-b/3/a-pow(-(alpha+sqrt(delta)),1./3)-pow(-(alpha-sqrt(delta)),1./3);

      p0X=a1P*p0Z+b1P;
      p2Z=(c3*p0Z+c4)/(p0Z-treeCookV1->posC2Z);

      r=sqrt((p1X-p0X)*(p1X-p0X)+(p1Z-p0Z)*(p1Z-p0Z));
      a2=(p2X-treeCookV1->posC2X)/(p2Z-treeCookV1->posC2Z);
      b2=treeCookV1->posC2X-a2*treeCookV1->posC2Z;

      aCond=(p1X-p2X)/(p1Z-p2Z);
      bCond=p1X-aCond*p1Z;
      if(p2Z<=380 && p2Z>-380 && p0Z>(p0X-bCond)/aCond){
        par.clear();
        par.push_back(a1);
        par.push_back(b1);
        par.push_back(r);
        par.push_back(p0Z);
        par.push_back(p0X);
        par.push_back(a2);
        par.push_back(b2);
        vect_par.push_back(par);


        point.clear();
        point.push_back(treeCookV1->posC4Z);
        point.push_back(treeCookV1->posC4X);
        point.push_back(treeCookV1->posC3Z);
        point.push_back(treeCookV1->posC3X);
        point.push_back(p1Z);
        point.push_back(p1X);
        point.push_back(p0Z);
        point.push_back(p0X);
        point.push_back(p2Z);
        point.push_back(p2X);        
        point.push_back(treeCookV1->posC2Z);
        point.push_back(treeCookV1->posC2X);
        vect_point.push_back(point);

        nTrack++;

        if(treeCookV1->event<10084) numTrack++;

      }
    }

    else{
      p0Z1=-b/3/a+2*sqrt(-beta)*cos(acos(alpha/pow(-beta,3./2))/3);
      p0Z2=-b/3/a+2*sqrt(-beta)*cos((acos(alpha/pow(-beta,3./2))+2*TMath::Pi())/3);
      p0Z3=-b/3/a+2*sqrt(-beta)*cos((acos(alpha/pow(-beta,3./2))-2*TMath::Pi())/3);

      p0X1=a1P*p0Z1+b1P;
      p2Z1=(c3*p0Z1+c4)/(p0Z1-treeCookV1->posC2Z);
      p0X2=a1P*p0Z2+b1P;
      p2Z2=(c3*p0Z2+c4)/(p0Z2-treeCookV1->posC2Z);
      p0X3=a1P*p0Z3+b1P;
      p2Z3=(c3*p0Z3+c4)/(p0Z3-treeCookV1->posC2Z);

      r1=sqrt((p2X-p0X1)*(p2X-p0X1)+(p2Z1-p0Z1)*(p2Z1-p0Z1));
      a21=(p2X-treeCookV1->posC2X)/(p2Z1-treeCookV1->posC2Z);
      b21=treeCookV1->posC2X-a21*treeCookV1->posC2Z;

      r2=sqrt((p2X-p0X2)*(p2X-p0X2)+(p2Z2-p0Z2)*(p2Z2-p0Z2));
      a22=(p2X-treeCookV1->posC2X)/(p2Z2-treeCookV1->posC2Z);
      b22=treeCookV1->posC2X-a22*treeCookV1->posC2Z;
      
      r3=sqrt((p2X-p0X3)*(p2X-p0X3)+(p2Z3-p0Z3)*(p2Z3-p0Z3));
      a23=(p2X-treeCookV1->posC2X)/(p2Z3-treeCookV1->posC2Z);
      b23=treeCookV1->posC2X-a23*treeCookV1->posC2Z;

      aCond1=(p1X-p2X)/(p1Z-p2Z1);
      bCond1=p1X-aCond1*p1Z;

      aCond2=(p1X-p2X)/(p1Z-p2Z2);
      bCond2=p1X-aCond2*p1Z;

      aCond3=(p1X-p2X)/(p1Z-p2Z3);
      bCond3=p1X-aCond3*p1Z;

      if(p2Z1<=380 && p2Z1>=-380 && p0Z1>(p0X1-bCond1)/aCond1){
        par.clear();
        par.push_back(a1);
        par.push_back(b1);
        par.push_back(r1);
        par.push_back(p0Z1);
        par.push_back(p0X1);
        par.push_back(a21);
        par.push_back(b21);
        vect_par.push_back(par);

        pL=0.299792458*1.5*r1;
       // s1=sqrt((treeCookV1->posC3Z-p1Z)*(treeCookV1->posC3Z-p1Z)+(treeCookV1->posC3X-p1X)*(treeCookV1->posC3X-p1X));
       // s2=sqrt((treeCookV1->posC2Z-p2Z1)*(treeCookV1->posC2Z-p2Z1)+(treeCookV1->posC2X-p2X)*(treeCookV1->posC2X-p2X));
       // s0=acos(((p1Z-p0Z1)*(p2Z1-p0Z1)+(p1X-p0X1)*(p2X-p0X1))/r1/r1)*r1;
       // pT=pL*fabs(treeCookV1->posC3Y-treeCookV1->posC2Y)/(s1+s2+s0);
        pT=pL*fabs(treeCookV1->posC4Y-treeCookV1->posC3Y)/sqrt((treeCookV1->posC4Z-treeCookV1->posC3Z)*(treeCookV1->posC4Z-treeCookV1->posC3Z)+(treeCookV1->posC4X-treeCookV1->posC3X)*(treeCookV1->posC4X-treeCookV1->posC3X));
        hist->Fill(sqrt(pL*pL+pT*pT));
        if(atan(a21)>0)
          histPVsTheta->Fill(atan(a21),sqrt(pL*pL+pT*pT));
        else
          histPVsTheta->Fill(TMath::Pi()+atan(a21),sqrt(pL*pL+pT*pT));
        histPVsZ->Fill(-b21/a21,sqrt(pL*pL+pT*pT));

        point.clear();
        point.push_back(treeCookV1->posC4Z);
        point.push_back(treeCookV1->posC4X);
        point.push_back(treeCookV1->posC3Z);
        point.push_back(treeCookV1->posC3X);
        point.push_back(p1Z);
        point.push_back(p1X);
        point.push_back(p0Z1);
        point.push_back(p0X1);
        point.push_back(p2Z1);
        point.push_back(p2X);
        point.push_back(treeCookV1->posC2Z);
        point.push_back(treeCookV1->posC2X);
        vect_point.push_back(point);

        nTrack++;

        if(treeCookV1->event<10084) numTrack++;
      }

      if(p2Z2<=380 && p2Z2>=-380 && p0Z2>(p0X2-bCond2)/aCond2){
        par.clear();
        par.push_back(a1);
        par.push_back(b1);
        par.push_back(r2);
        par.push_back(p0Z2);
        par.push_back(p0X2);
        par.push_back(a22);
        par.push_back(b22);
        vect_par.push_back(par);

        pL=0.299792458*1.5*r2;
        pT=pL*fabs(treeCookV1->posC4Y-treeCookV1->posC3Y)/sqrt((treeCookV1->posC4Z-treeCookV1->posC3Z)*(treeCookV1->posC4Z-treeCookV1->posC3Z)+(treeCookV1->posC4X-treeCookV1->posC3X)*(treeCookV1->posC4X-treeCookV1->posC3X));
        hist->Fill(sqrt(pL*pL+pT*pT));
        if(atan(a22)>0)
          histPVsTheta->Fill(atan(a22),sqrt(pL*pL+pT*pT));
        else
          histPVsTheta->Fill(TMath::Pi()+atan(a22),sqrt(pL*pL+pT*pT));
        histPVsZ->Fill(-b22/a22,sqrt(pL*pL+pT*pT));

        point.clear();
        point.push_back(treeCookV1->posC4Z);
        point.push_back(treeCookV1->posC4X);
        point.push_back(treeCookV1->posC3Z);
        point.push_back(treeCookV1->posC3X);
        point.push_back(p1Z);
        point.push_back(p1X);
        point.push_back(p0Z2);
        point.push_back(p0X2);
        point.push_back(p2Z2);
        point.push_back(p2X);
        point.push_back(treeCookV1->posC2Z);
        point.push_back(treeCookV1->posC2X);
        vect_point.push_back(point);

        nTrack++;
        if(treeCookV1->event<10084) numTrack++;
      }

      if(p2Z3<=380 && p2Z3>=-380 && p0Z3>(p0X3-bCond3)/aCond3){
        par.clear();
        par.push_back(a1);
        par.push_back(b1);
        par.push_back(r3);
        par.push_back(p0Z3);
        par.push_back(p0X3);
        par.push_back(a23);
        par.push_back(b23);
        vect_par.push_back(par);

        pL=0.299792458*1.5*r1;
        pT=pL*fabs(treeCookV1->posC4Y-treeCookV1->posC3Y)/sqrt((treeCookV1->posC4Z-treeCookV1->posC3Z)*(treeCookV1->posC4Z-treeCookV1->posC3Z)+(treeCookV1->posC4X-treeCookV1->posC3X)*(treeCookV1->posC4X-treeCookV1->posC3X));
        hist->Fill(sqrt(pL*pL+pT*pT));
        if(atan(a23)>0)
          histPVsTheta->Fill(atan(a23),sqrt(pL*pL+pT*pT));
        else
          histPVsTheta->Fill(TMath::Pi()+atan(a23),sqrt(pL*pL+pT*pT));
        histPVsZ->Fill(-b23/a23,sqrt(pL*pL+pT*pT));

        point.clear();
        point.push_back(treeCookV1->posC4Z);
        point.push_back(treeCookV1->posC4X);
        point.push_back(treeCookV1->posC3Z);
        point.push_back(treeCookV1->posC3X);
        point.push_back(p1Z);
        point.push_back(p1X);
        point.push_back(p0Z3);
        point.push_back(p0X3);
        point.push_back(p2Z3);
        point.push_back(p2X);
        point.push_back(treeCookV1->posC2Z);
        point.push_back(treeCookV1->posC2X);
        vect_point.push_back(point);

        nTrack++;
        if(treeCookV1->event<10084) numTrack++;
      }
    }

    vect_nTrack.push_back(nTrack);

    if(nTrack==2) cout<<treeCookV1->event<<endl;

  }


  return 0;
}

Long_t Det_MWPC::done_trackV2(){
  TCanvas *canv=new TCanvas("Momentum distribution","Momentum distribution",1000,600);
  hist->Draw();
  canv->SaveAs("momentum distribution.pdf");


  TCanvas *canvPVsTheta=new TCanvas("Momentum vs theta","Momentum vs theta",1000,1000);
  histPVsTheta->Draw("colz");
  canvPVsTheta->SaveAs("momentum-vs-theta.pdf");

  TCanvas *canvPVsZ=new TCanvas("Momentum vs z","Momentum vs z",1000,1000);
  histPVsZ->Draw("colz");
  canvPVsZ->SaveAs("momentum-vs-z.pdf");

  return 0;
}

