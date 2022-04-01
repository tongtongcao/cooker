#ifndef __MWPC_Event__
#define __MWPC_Event__
#include <string>
#include <TH1D.h>
#include <TDirectory.h>
#include <TList.h>
#include <vector>
using namespace std;

class MwpcEvent{
 private:
  int eventNumber;
  std::string eventName;
  TH1D* eventHist;
  int numBin;
  double lowBound;
  double upBound;
  bool stripTransStat; // =1: Dont need transformation; 0: need transformation
 public:
 MwpcEvent(int nBin, bool stTransStat):numBin(nBin),stripTransStat(stTransStat){} // Constructor for data analysis
 MwpcEvent(int eventNum,int nBin, bool stTransStat):eventNumber(eventNum),numBin(nBin),stripTransStat(stTransStat){} // Constructor for data analysis
 MwpcEvent(int eventNum, std::string name,int nBin,double low,double up,bool stTransStat):eventNumber(eventNum),eventName(name),eventHist(0),numBin(nBin),lowBound(low),upBound(up),stripTransStat(stTransStat){} // Constructor for histogram by event
 int getEventNumber(){return eventNumber;}  // Return event number
 ///////////////////////  Histogram of Strip ADC Readout by Event  //////////////////////////
 int getNBin(){return numBin;}  //return number of bins of histogram == number of strips    
 double getLowBound(){return lowBound;} // return lower of histogram == -0.5
 double getUpBound(){return upBound;} // return upper of histogram  == number of strips - 0.5
 void setHist(); // Define histogram
 void addHistData(int* value); // Set bin content of histogram without pedestal subtraction, gain calibtration and threshold cut
 void addHistData_pedSub(int* value, const double* pdd, const double *sigma); // Set bin content of histogram with pedestal subtraction and gain calibtration
 void addHistData_pedSub_sigmaCut(int* value, const double* pdd, const double* sigma, double nSigma); // Set bin content of histogram with pedestal subtraction,gain calibtration and global threshold cut; Need to specify factor of sigma to determine threshold
 TH1D* getHist(){return eventHist;} // Return histogram 
 void drawHist(); // Draw histogram
 void saveHist(); // Plot histogram
 ///////////////////////  Plot Strip ADC Readout by Event  //////////////////////////
  double* getAdcValue(int* value, const double* pdd);  // Return Adc value per strip after pedetal subtraction 
  double* getAdcValue(int* value, const double* pdd, const double* sigma);  // Return Adc value per strip after pedetal subtraction and gain calibtration 

 ///////////////////////  Data Analysis after pedetal subtraction and global threshold cut; Need to specify factor of sigma to determine threshold //////////////////////////
 double* getAdcValue(int* value, const double* pdd, const double* sigma, double nSigma);  // Return Adc value per strip after pedetal subtraction and global threshold cut; Need to specify factor of sigma to determine threshold
 int getNHit(int* value, const double* pdd, const double* sigma, double nSigma); // Return number of hits per event after pedetal subtraction and global threshold cut; Need to specify factor of sigma to determine threshold
 int* multiplet3Type(int* value, const double* pdd, const double* sigma, double nSigma); // Return number of singlet, doublet, triplet (and greater) per event after pedetal subtraction and global threshold cut; Need to specify factor of sigma to determine threshold
 int* multiplet4Type(int* value, const double* pdd, const double* sigma, double nSigma); // Return number of singelt, doublet, triplet, quartet (and greater) per event after pedetal subtraction and global threshold cut; Need to specify factor of sigma to determine threshold
 int* multipletNType(int* value, const double* pdd, const double* sigma, double nSigma, int nType); // Return number of multiplets per event after pedetal subtraction and global threshold cut; Need to specify factor of sigma to determine threshold; Need to specify number of multiplets
 vector<vector<int> > clusterStripNumber(int* value, const double* pdd, const double* sigma,double nSigma); // Return strip numbers of clusters in order of cluster per event after pedetal subtraction and global threshold cut; Need to specify factor of sigma to determine threshold
 int getNCluster(MwpcEvent* event,int* value, const double* pdd, const double* sigma,double nSigma);  // Return number of clusters per event after pedetal subtraction and global threshold cut; Need to specify factor of sigma to determine threshold
 vector<int> getClusterNStrip(MwpcEvent* event,int* value, const double* pdd, const double* sigma,double nSigma); // Return number of cluster's strips in order of cluster per event after pedetal subtraction and global threshold cut; Need to specify factor of sigma to determine threshold
 double* getClusterTotCharge(MwpcEvent* event, int* value, const double* pdd, const double* sigma, double nSigma); // Return total charge of clusters in order of cluster per event after pedetal subtraction and global threshold cut; Need to specify factor of sigma to determine threshold
 ///////////////////////  Data Analysis after pedetal subtraction and global threshold cut; Need to specify factor of sigma to determine threshold //////////////////////////
 
 ///////////////////////  Data Analysis after pedetal subtraction and specified total charge and threshold cut by cluster //////////////////////////
 vector<vector<int> > clusterStripNumber(MwpcEvent* event, int* value, const double* pdd, const double* sigma, bool statusTotalChargeCut, const double* totalChargeCut, bool statusSingletThreshold, double singletThreshold, bool statusDoubletThresholdTotalChargeCut, double doubletThreshold, double doubletTotalChargeCut); 
 int getNCluster(MwpcEvent* event, int* value, const double* pdd, const double* sigma, bool statusTotalChargeCut, const double* totalChargeCut, bool statusSingletThreshold, double singletThreshold, bool statusDoubletThresholdTotalChargeCut, double doubletThreshold, double doubletTotalChargeCut); 
 vector<int> getClusterNStrip(MwpcEvent* event, int* value, const double* pdd, const double* sigma, bool statusTotalChargeCut, const double* totalChargeCut, bool statusSingletThreshold, double singletThreshold, bool statusDoubletThresholdTotalChargeCut, double doubletThreshold, double doubletTotalChargeCut);
 double* getClusterTotCharge(MwpcEvent* event, int* value, const double* pdd, const double* sigma, bool statusTotalChargeCut, const double* totalChargeCut, bool statusSingletThreshold, double singletThreshold, bool statusDoubletThresholdTotalChargeCut, double doubletThreshold, double doubletTotalChargeCut);
 vector<vector<double> > clusterStripCharge(MwpcEvent* event, int* value, const double* pdd, const double* sigma, bool statusTotalChargeCut, const double* totalChargeCut, bool statusSingletThreshold, double singletThreshold, bool statusDoubletThresholdTotalChargeCut, double doubletThreshold, double doubletTotalChargeCut);

 ///////////////////////  Data Analysis for each cluster //////////////////////////
 vector<double> getEachClusterStCharge(MwpcEvent* event,int* value, const double* pdd,vector<int> stNum);
 double getEachClusterPositionCentroid(vector<int> stNum,vector<double> stCharge);
 double getEachClusterPositionCentroid(MwpcEvent* event,int* value, const double* pdd,vector<int> stNum); 
 double getEachClusterPositionErrorCentroid(vector<double> stCharge);
 double getEachClusterPositionErrorCentroid(MwpcEvent* event,int* value, const double* pdd,vector<int> stNum);  

  virtual ~MwpcEvent();
  ClassDef(MwpcEvent,1);
};

#endif
