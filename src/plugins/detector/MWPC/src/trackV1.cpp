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

Long_t Det_MWPC::histos_trackV1(){
  ostringstream sstr_name,sstr_title;

  return 0;
}

Long_t Det_MWPC::startup_trackV1(){
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

static int numTrack=0;
Long_t Det_MWPC::process_trackV1(){
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

static double myfunction(double *x, double *par){
  double xx =x[0];
  double f=0; 
  if(xx>=380) f=par[0]*xx+par[1];
  else if(xx>=par[7] && xx<=380) f=sqrt(par[2]*par[2]-(xx-par[3])*(xx-par[3]))+par[4];
  else f=par[5]*xx+par[6];
  return f;
}

static double myfunction1(double *x, double *par){
  double xx =x[0];
  double f=0;    
  if(xx>=380) f=par[0]*xx+par[1];
  else if(xx>=-par[2]+par[3] && xx<=380) f=sqrt(par[2]*par[2]-(xx-par[3])*(xx-par[3]))+par[4];
  return f;
}

static double myfunction2(double *x, double *par){
  double xx =x[0];
  double f=0;
  if(xx>=par[7]) f=par[5]*xx+par[6];
  else if(xx>=-par[2]+par[3] && xx<par[7]) f=-sqrt(par[2]*par[2]-(xx-par[3])*(xx-par[3]))+par[4];
  return f;
} 

Long_t Det_MWPC::done_trackV1(){


cout<<numTrack<<"  aaaaaaaa"<<endl;


  int num[4]={0};
  for(int i=0;i<vect_nTrack.size();i++){
    if(vect_nTrack[i]==0) num[0]++;
    else if(vect_nTrack[i]==1) num[1]++;
    else if(vect_nTrack[i]==2) num[2]++;
    else if(vect_nTrack[i]==3) num[3]++;
  }
  cout<<num[0]<<"  "<<num[1]<<"  "<<num[2]<<"  "<<num[3]<<endl;


  double parameter[7];
  double pointZX[12];

  TCanvas *canv1Track=new TCanvas("example for 1 track","example for 1 track",1000,1000);
  TF1 *ff=new TF1("ff","x+1700-800",-500,800);
  ff->SetTitle("");
  ff->SetLineColor(0);
  ff->Draw("l");

  TLine *lineC4=new TLine(699.5,1240-360,699.5,1240+360);
  lineC4->SetLineColor(kBlack);
  lineC4->Draw("same");
  TLine *lineC3=new TLine(480.5,1240-320,480.5,1240+320);
  lineC3->SetLineColor(kBlack);
  lineC3->Draw("same");
  TLine *lineC2=new TLine(0-280,629.6,0+280,629.6);
  lineC2->SetLineColor(kBlack);
  lineC2->Draw("same");

  TLine *line1=new TLine(-380,725.01,380,725.01);
  line1->SetLineColor(kBlack);
  line1->Draw("same");
  TLine *line2=new TLine(-380,1120*2-725.01,380,1120*2-725.01);
  line2->SetLineColor(kBlack);
  line2->Draw("same");
  TLine *line3=new TLine(-380,725.01,-380,1120*2-725.01);
  line3->SetLineColor(kBlack);
  line3->Draw("same");
  TLine *line4=new TLine(380,725.01,380,1120*2-725.01);
  line4->SetLineColor(kBlack);
  line4->Draw("same");

  for(int j=0;j<7;j++){
    parameter[j]=vect_par[57][j];
  }
  for(int j=0;j<12;j++){
    pointZX[j]=vect_point[57][j];
  }

  TF1 *func=new TF1("func",myfunction,-400,2000,8);
  TF1 *func1=new TF1("func",myfunction1,-400,2000,8);
  TF1 *func2=new TF1("func",myfunction2,-400,2000,8);
  TMarker *markPoint;
  if(parameter[5]>=0){
    func->FixParameter(0,parameter[0]);
    func->FixParameter(1,parameter[1]);
    func->FixParameter(2,parameter[2]);
    func->FixParameter(3,parameter[3]);
    func->FixParameter(4,parameter[4]);
    func->FixParameter(5,parameter[5]);
    func->FixParameter(6,parameter[6]);
    func->FixParameter(7,pointZX[8]);
    func->DrawClone("lsame");
  }

  else{
    func1->FixParameter(0,parameter[0]);
    func1->FixParameter(1,parameter[1]);
    func1->FixParameter(2,parameter[2]);
    func1->FixParameter(3,parameter[3]);
    func1->FixParameter(4,parameter[4]);
    func1->FixParameter(5,parameter[5]);
    func1->FixParameter(6,parameter[6]);
    func1->FixParameter(7,pointZX[8]);

    func2->FixParameter(0,parameter[0]);
    func2->FixParameter(1,parameter[1]);
    func2->FixParameter(2,parameter[2]);
    func2->FixParameter(3,parameter[3]);
    func2->FixParameter(4,parameter[4]);
    func2->FixParameter(5,parameter[5]);
    func2->FixParameter(6,parameter[6]);
    func2->FixParameter(7,pointZX[8]);

    for(int k=0;k<1000;k++){
      markPoint=new TMarker(-parameter[2]+parameter[3]+(800+parameter[2]-parameter[3])/1000*k,func1->Eval(-parameter[2]+parameter[3]+(800+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();
      markPoint=new TMarker(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]+parameter[2]-parameter[3])/1000*k,func2->Eval(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();

      markPoint=new TMarker(-parameter[2]+parameter[3]+(10+parameter[2]-parameter[3])/1000*k,func1->Eval(-parameter[2]+parameter[3]+(10+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();

      markPoint=new TMarker(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]-2+parameter[2]-parameter[3])/1000*k,func2->Eval(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]-2+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();
    }
  }

  TMarker *mark[6];
  for(int i=0;i<6;i++){
    mark[i]=new TMarker(pointZX[2*i],pointZX[2*i+1],20);
    if(i!=3) mark[i]->SetMarkerColor(kBlue);
    else mark[i]->SetMarkerColor(kBlack);
    mark[i]->Draw();
  }
 
  canv1Track->SaveAs("example for 1 track.pdf");


  TCanvas *canv2Track=new TCanvas("example for 2 track","example for 2 track",1000,1000);
  ff->Draw("l");

  lineC4->Draw("same");
  lineC3->Draw("same");
  lineC2->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");

  for(int j=0;j<7;j++){
    parameter[j]=vect_par[5140][j];
  }
  for(int j=0;j<12;j++){
    pointZX[j]=vect_point[5140][j];
  }

  if(parameter[5]>=0){
    func->FixParameter(0,parameter[0]);
    func->FixParameter(1,parameter[1]);
    func->FixParameter(2,parameter[2]);
    func->FixParameter(3,parameter[3]);
    func->FixParameter(4,parameter[4]);
    func->FixParameter(5,parameter[5]);
    func->FixParameter(6,parameter[6]);
    func->FixParameter(7,pointZX[8]);
    func->DrawClone("lsame");
  }

  else{
    func1->FixParameter(0,parameter[0]);
    func1->FixParameter(1,parameter[1]);
    func1->FixParameter(2,parameter[2]);
    func1->FixParameter(3,parameter[3]);
    func1->FixParameter(4,parameter[4]);
    func1->FixParameter(5,parameter[5]);
    func1->FixParameter(6,parameter[6]);
    func1->FixParameter(7,pointZX[8]);

    func2->FixParameter(0,parameter[0]);
    func2->FixParameter(1,parameter[1]);
    func2->FixParameter(2,parameter[2]);
    func2->FixParameter(3,parameter[3]);
    func2->FixParameter(4,parameter[4]);
    func2->FixParameter(5,parameter[5]);
    func2->FixParameter(6,parameter[6]);
    func2->FixParameter(7,pointZX[8]);

    for(int k=0;k<1000;k++){
      markPoint=new TMarker(-parameter[2]+parameter[3]+(800+parameter[2]-parameter[3])/1000*k,func1->Eval(-parameter[2]+parameter[3]+(800+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();
      markPoint=new TMarker(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]+parameter[2]-parameter[3])/1000*k,func2->Eval(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();

      markPoint=new TMarker(-parameter[2]+parameter[3]+(10+parameter[2]-parameter[3])/1000*k,func1->Eval(-parameter[2]+parameter[3]+(10+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();

      markPoint=new TMarker(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]-2+parameter[2]-parameter[3])/1000*k,func2->Eval(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]-2+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();
    }
  }

  for(int j=0;j<7;j++){
    parameter[j]=vect_par[5141][j];
  }
  for(int j=0;j<12;j++){
    pointZX[j]=vect_point[5141][j];
  }

  if(parameter[5]>=0){
    func->FixParameter(0,parameter[0]);
    func->FixParameter(1,parameter[1]);
    func->FixParameter(2,parameter[2]);
    func->FixParameter(3,parameter[3]);
    func->FixParameter(4,parameter[4]);
    func->FixParameter(5,parameter[5]);
    func->FixParameter(6,parameter[6]);
    func->FixParameter(7,pointZX[8]);
    func->DrawClone("lsame");
  }

  else{
    func1->FixParameter(0,parameter[0]);
    func1->FixParameter(1,parameter[1]);
    func1->FixParameter(2,parameter[2]);
    func1->FixParameter(3,parameter[3]);
    func1->FixParameter(4,parameter[4]);
    func1->FixParameter(5,parameter[5]);
    func1->FixParameter(6,parameter[6]);
    func1->FixParameter(7,pointZX[8]);

    func2->FixParameter(0,parameter[0]);
    func2->FixParameter(1,parameter[1]);
    func2->FixParameter(2,parameter[2]);
    func2->FixParameter(3,parameter[3]);
    func2->FixParameter(4,parameter[4]);
    func2->FixParameter(5,parameter[5]);
    func2->FixParameter(6,parameter[6]);
    func2->FixParameter(7,pointZX[8]);

    for(int k=0;k<1000;k++){
      markPoint=new TMarker(-parameter[2]+parameter[3]+(800+parameter[2]-parameter[3])/1000*k,func1->Eval(-parameter[2]+parameter[3]+(800+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();
      markPoint=new TMarker(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]+parameter[2]-parameter[3])/1000*k,func2->Eval(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();

      markPoint=new TMarker(-parameter[2]+parameter[3]+(10+parameter[2]-parameter[3])/1000*k,func1->Eval(-parameter[2]+parameter[3]+(10+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();

      markPoint=new TMarker(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]-2+parameter[2]-parameter[3])/1000*k,func2->Eval(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]-2+parameter[2]-parameter[3])/1000*k),20);
      markPoint->SetMarkerSize(0.2);
      markPoint->SetMarkerColor(kRed);
      markPoint->Draw();
    }
  }

  for(int i=0;i<6;i++){
    mark[i]=new TMarker(pointZX[2*i],pointZX[2*i+1],20);
    if(i!=3) mark[i]->SetMarkerColor(kBlue);
    else mark[i]->SetMarkerColor(kBlack);
    mark[i]->Draw();
  }

  canv2Track->SaveAs("example for 2 track.pdf");


/*
  TCanvas *canv=new TCanvas("canv","canv",1000,1000);
  TF1 *ff=new TF1("ff","x+1700-800",-500,800);
  ff->SetTitle("");
  ff->SetLineColor(0);
  ff->Draw("l");

  TLine *lineC4=new TLine(699.5,1240-360,699.5,1240+360);
  lineC4->SetLineColor(kBlack);
  lineC4->Draw("same");
  TLine *lineC3=new TLine(480.5,1240-320,480.5,1240+320);
  lineC3->SetLineColor(kBlack);
  lineC3->Draw("same");
  TLine *lineC2=new TLine(0-280,629.6,0+280,629.6);
  lineC2->SetLineColor(kBlack);
  lineC2->Draw("same");

  TLine *line1=new TLine(-380,725.01,380,725.01);
  line1->SetLineColor(kBlack);
  line1->Draw("same");
  TLine *line2=new TLine(-380,1120*2-725.01,380,1120*2-725.01);
  line2->SetLineColor(kBlack);
  line2->Draw("same");
  TLine *line3=new TLine(-380,725.01,-380,1120*2-725.01);
  line3->SetLineColor(kBlack);
  line3->Draw("same");
  TLine *line4=new TLine(380,725.01,380,1120*2-725.01);
  line4->SetLineColor(kBlack);
  line4->Draw("same");


  TF1 *func=new TF1("func",myfunction,-400,2000,8);
  TF1 *func1=new TF1("func",myfunction1,-400,2000,8);
  TF1 *func2=new TF1("func",myfunction2,-400,2000,8);
  TMarker *markPoint;
  for(int i=0;i<vect_par.size();i++){
    for(int j=0;j<7;j++){
      parameter[j]=vect_par[i][j];
    }
    for(int j=0;j<12;j++){
      pointZX[j]=vect_point[i][j];
    }

    if(parameter[5]>=0){
      func->FixParameter(0,parameter[0]);
      func->FixParameter(1,parameter[1]);
      func->FixParameter(2,parameter[2]);
      func->FixParameter(3,parameter[3]);
      func->FixParameter(4,parameter[4]);
      func->FixParameter(5,parameter[5]);
      func->FixParameter(6,parameter[6]);
      func->FixParameter(7,pointZX[8]);
      func->DrawClone("lsame");
    }

    else{
      func1->FixParameter(0,parameter[0]);
      func1->FixParameter(1,parameter[1]);
      func1->FixParameter(2,parameter[2]);
      func1->FixParameter(3,parameter[3]);
      func1->FixParameter(4,parameter[4]);
      func1->FixParameter(5,parameter[5]);
      func1->FixParameter(6,parameter[6]);
      func1->FixParameter(7,pointZX[8]);

      func2->FixParameter(0,parameter[0]);
      func2->FixParameter(1,parameter[1]);
      func2->FixParameter(2,parameter[2]);
      func2->FixParameter(3,parameter[3]);
      func2->FixParameter(4,parameter[4]);
      func2->FixParameter(5,parameter[5]);
      func2->FixParameter(6,parameter[6]);
      func2->FixParameter(7,pointZX[8]);

      for(int k=0;k<1000;k++){
        markPoint=new TMarker(-parameter[2]+parameter[3]+(800+parameter[2]-parameter[3])/1000*k,func1->Eval(-parameter[2]+parameter[3]+(800+parameter[2]-parameter[3])/1000*k),20);
        markPoint->SetMarkerSize(0.2);
        markPoint->SetMarkerColor(kRed);
        markPoint->Draw();
        markPoint=new TMarker(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]+parameter[2]-parameter[3])/1000*k,func2->Eval(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]+parameter[2]-parameter[3])/1000*k),20);
        markPoint->SetMarkerSize(0.2);
        markPoint->SetMarkerColor(kRed);
        markPoint->Draw();


        markPoint=new TMarker(-parameter[2]+parameter[3]+(10+parameter[2]-parameter[3])/1000*k,func1->Eval(-parameter[2]+parameter[3]+(10+parameter[2]-parameter[3])/1000*k),20);
        markPoint->SetMarkerSize(0.2);
        markPoint->SetMarkerColor(kRed);
        markPoint->Draw();

        markPoint=new TMarker(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]-2+parameter[2]-parameter[3])/1000*k,func2->Eval(-parameter[2]+parameter[3]+((340-parameter[6])/parameter[5]-2+parameter[2]-parameter[3])/1000*k),20);
        markPoint->SetMarkerSize(0.2);
        markPoint->SetMarkerColor(kRed);
        markPoint->Draw();
      }

    }


    TMarker *mark[6];
    for(int i=0;i<6;i++){
      mark[i]=new TMarker(pointZX[2*i],pointZX[2*i+1],20);
      if(i!=3) mark[i]->SetMarkerColor(kBlue);
      else mark[i]->SetMarkerColor(kBlack);
      mark[i]->Draw();
    }

  }
  canv->SaveAs("track.pdf");
*/
  return 0;
}

