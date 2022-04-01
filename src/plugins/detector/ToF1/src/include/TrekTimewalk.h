#ifndef __TREK_TIMEWALK__
#define __TREK_TIMEWALK__
#include <string>
#include <TH1D.h>
#include <TDirectory.h>
#include <TList.h>
#include <TH2D.h>
#include <TProfile.h> 
class TrekTimewalk{
 private:
  unsigned int mNBinsX;
  unsigned int mNBinsY;
  double mLowBoundX;
  double mUpBoundX;
  double mLowBoundY;
  double mUpBoundY;
  std::string mName;
  TH2D* mH2;
  TProfile* mProf;
  double mConst;
  double mOrder;
  double mErrorConst;
  double mErrorOrder;
 public:
 TrekTimewalk():mNBinsX(100),mNBinsY(100),mLowBoundX(0),mUpBoundX(4000),mLowBoundY(-20),mUpBoundY(20),mName("unNamed"),mH2(0),mProf(0),mConst(0),mOrder(2),mErrorConst(0),mErrorOrder(0){}
 TrekTimewalk(const std::string& name):mNBinsX(100),mNBinsY(100),mLowBoundX(0),mUpBoundX(4000),mLowBoundY(-20),mUpBoundY(20),mName(name),mH2(0),mProf(0),mConst(0),mOrder(2),mErrorConst(0),mErrorOrder(0){}
  void addData(const double& diff,const double& adc);
  void setName(std::string name);
  bool getValue(double& lambda,double& order,bool writeImg=false);
  void getError(double& errLambda,double& errOrder);
  void setUpperX(const double& range){mUpBoundX=range;}
  void setNBinsX(unsigned n){mNBinsX=n;}
  void setNBinsY(unsigned n){mNBinsY=n;}
};

#endif
