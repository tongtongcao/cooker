/*
  The length of all of array variables is 12, which corresponding C2XL, C2XR, C2YL, C2YR, C3XL, C3XR, C3YL, C3YR, C4XL, C4XR, C4YL and C4YR.
  Variables for pairs are introduced after pairing X and Y clusters. In the variable arrays for pairs, there are 6 pairs, such as 0&2, 1&3, 4&6, 5&7, 8&10 and 9&11.
  In the calibration code, all cuts are divided into two layers. The first layer includes threshod, total-charge cut and first-layer flag cut, and the second layer includes pair cut and second-layer flag cut.
  For convinence of futher analysis, the same variables are set for no cut, after threshold and total-charge cuts, after pair cut, and after all of threshold, total-charge, and pair cuts. That is, four sets of variables are specified for four cut contidtions, respectively. 
  For threshold and total-charge cuts, cable swapping is handled, so extra three subset of variables specified for cable swapping are involved.
  The first- and second-layer flag cuts are flagged for each set of variables for Left and Right separately. In the variable arrays, Left corresponds 0, 2, 4, 6, 8, and 10, and Right corresponds 1, 3, 5, 7, 9 and 11. Data type of flag variables is bool, 1 denotes accpetion, and 0 denotes rejection. 
*/


#ifndef __MWPCTREE_H_
#define __MWPCTREE_H_
#include "cookerrawtree.h" // for CRTBase
#include <vector>
using namespace std;
class CRTCaliMwpc:public CRTBase{
  public:
    //// Common variables
    UInt_t run;
    UInt_t event;
    vector<int> adcValueOriginal[12];  // Orignal adc value
    vector<double> adcValuePedSub[12]; // Adc value afer applying pedestal subtraction and C2X strip transformation; Negative value is set as 0 after subtraction; The following variables are assigned based on these adc values
   
    //// Specified variables for no cut
    UInt_t noCutNumCluster[12]; //Number of clusters 
    vector<double> noCutClusterTC[12] ;//Total charge for each cluster; All clusters are saved as a vector
    vector<UInt_t> noCutClusterNumSt[12]; // Number of strips for each cluster; All clusters are saved as a vector
    vector<vector<int> > noCutClusterStNum[12]; // Strip numbers saved as a vector for each clustr; All clusters are saved as a nested vector
    vector<vector<double> > noCutClusterStCharge[12]; // Strip charges saved as a vector for each clustr; All clusters are saved as a nested vector
    vector<double> noCutClusterPosCentroid[12]; //Position determined by centroid method for each cluster; All clusters are saved as a vector
    vector<double> noCutClusterPosErrCentroid[12]; //Position error determined by centroid method for each cluster; All clusters are saved as a vector
    bool noCutLeftFlagFirst; // If there are tracks kept at the left of detector after applying the first layer cuts; Corresponds indices 0, 2, 4, 6, 8 and 10 for array variables
    bool noCutRightFlagFirst; // If there are tracks kept at the right of detector after applying the first layer cuts; Corresponds indices 1, 3, 5, 7, 9 and 11 for array variables
    UInt_t noCutNumPair[12]; // Number of pairs; Number of pairs are the same for relevant indices, such as 0&2, 1&3, 4&6, 5&7, 8&10 and 9&11. 
    vector<UInt_t> noCutPairClusterNum[12]; // X and Y cluster numbers for each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<double> noCutPairClusterTC[12]; // Total charge of X and Y clusters for each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<UInt_t> noCutPairClusterNumSt[12]; // Number of strips for X and Y cluster of each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<vector<int> > noCutPairClusterStNum[12]; // Strip numbers for X and Y clusters of each pair; All pairs are saved as a nested vector; X and Y are saved separately in the array 
    vector<vector<double> > noCutPairClusterStCharge[12]; // Strip charges for X and Y clusters of each pair; All pairs are saved as a nested vector; X and Y are saved separately in the array 
    vector<double> noCutPairPosCentroid[12]; //Position determined by centroid method for all pairs
    vector<double> noCutPairPosErrCentroid[12]; //Position error determined by centroid method for all pairs
    bool noCutLeftFlagSecond; // If there are tracks kept at the left of detector after applying the second layer cuts; Corresponds indices 0, 2, 4, 6, 8 and 10 for arrays
    bool noCutRightFlagSecond; // If there are tracks kept at the right of detector after applying the second layer cuts; Corresponds indices 1, 3, 5, 7, 9 and 11 for arrays
    
    //// Specified variables for threshold and total-charge cuts
    //Cable swapping
    UInt_t thTCCutNumClusterC3XLGap7Swap; //Number of clusters after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 
    vector<double> thTCCutClusterTCC3XLGap7Swap;//Total charge for each cluster; All clusters are saved as a vector after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 
    vector<UInt_t> thTCCutClusterNumStC3XLGap7Swap; // Number of strips for each cluster; All clusters are saved as a vector after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 
    vector<vector<int> > thTCCutClusterStNumC3XLGap7Swap; // Strip numbers saved as a vector for each clustr; All clusters are saved as a nested vector after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 
    vector<vector<double> > thTCCutClusterStChargeC3XLGap7Swap; // Strip charges saved as a vector for each clustr; All clusters are saved as a nested vector after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 
    vector<double> thTCCutClusterPosCentroidC3XLGap7Swap; //Position determined by centroid method for each cluster; All clusters are saved as a vector after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 
    vector<double> thTCCutClusterPosErrCentroidC3XLGap7Swap; //Position error determined by centroid method for each cluster; All clusters are saved as a vector after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 

    UInt_t thTCCutNumClusterC4XLGap7Swap; //Number of clusters after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 
    vector<double> thTCCutClusterTCC4XLGap7Swap;//Total charge for each cluster; All clusters are saved as a vector after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 
    vector<UInt_t> thTCCutClusterNumStC4XLGap7Swap; // Number of strips for each cluster; All clusters are saved as a vector after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 
    vector<vector<int> > thTCCutClusterStNumC4XLGap7Swap; // Strip numbers saved as a vector for each clustr; All clusters are saved as a nested vector after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 
    vector<vector<double> > thTCCutClusterStChargeC4XLGap7Swap; // Strip charges saved as a vector for each clustr; All clusters are saved as a nested vector after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 
    vector<double> thTCCutClusterPosCentroidC4XLGap7Swap; //Position determined by centroid method for each cluster; All clusters are saved as a vector after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 
    vector<double> thTCCutClusterPosErrCentroidC4XLGap7Swap; //Position error determined by centroid method for each cluster; All clusters are saved as a vector after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 

    UInt_t thTCCutNumClusterC4XLGap11Swap; //Number of clusters after switching C4X3 and C4X4 
    vector<double> thTCCutClusterTCC4XLGap11Swap;//Total charge for each cluster; All clusters are saved as a vector after switching C4X3 and C4X4 
    vector<UInt_t> thTCCutClusterNumStC4XLGap11Swap; // Number of strips for each cluster; All clusters are saved as a vector after switching C4X3 and C4X4 
    vector<vector<int> > thTCCutClusterStNumC4XLGap11Swap; // Strip numbers saved as a vector for each clustr; All clusters are saved as a nested vector after switching C4X3 and C4X4
    vector<vector<double> > thTCCutClusterStChargeC4XLGap11Swap; // Strip charges saved as a vector for each clustr; All clusters are saved as a nested vector after switching C4X3 and C4X4 
    vector<double> thTCCutClusterPosCentroidC4XLGap11Swap; //Position determined by centroid method for each cluster; All clusters are saved as a vector after switching C4X3 and C4X4 
    vector<double> thTCCutClusterPosErrCentroidC4XLGap11Swap; //Position error determined by centroid method for each cluster; All clusters are saved as a vector after switching C4X3 and C4X4 
    //Cable swapping


    UInt_t thTCCutNumCluster[12]; //Number of clusters 
    vector<double> thTCCutClusterTC[12];//Total charge for each cluster; All clusters are saved as a vector
    vector<UInt_t> thTCCutClusterNumSt[12]; // Number of strips for each cluster; All clusters are saved as a vector
    vector<vector<int> > thTCCutClusterStNum[12]; // Strip numbers saved as a vector for each clustr; All clusters are saved as a nested vector
    vector<vector<double> > thTCCutClusterStCharge[12]; // Strip charges saved as a vector for each clustr; All clusters are saved as a nested vector
    vector<double> thTCCutClusterPosCentroid[12]; //Position determined by centroid method for each cluster; All clusters are saved as a vector
    vector<double> thTCCutClusterPosErrCentroid[12]; //Position error determined by centroid method for each cluster; All clusters are saved as a vector
    bool thTCCutLeftFlagFirst; // If there are tracks kept at the left of detector after applying the first layer cuts; Corresponds indices 0, 2, 4, 6, 8 and 10 for array variables
    bool thTCCutRightFlagFirst; // If there are tracks kept at the right of detector after applying the first layer cuts; Corresponds indices 1, 3, 5, 7, 9 and 11 for array variables
    UInt_t thTCCutNumPair[12]; // Number of pairs; Number of pairs are the same for relevant indices, such as 0&2, 1&3, 4&6, 5&7, 8&10 and 9&11. 
    vector<UInt_t> thTCCutPairClusterNum[12]; // X and Y cluster numbers for each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<double> thTCCutPairClusterTC[12]; // Total charge of X and Y clusters for each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<UInt_t> thTCCutPairClusterNumSt[12]; // Number of strips for X and Y cluster of each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<vector<int> > thTCCutPairClusterStNum[12]; // Strip numbers for X and Y clusters of each pair; All pairs are saved as a nested vector; X and Y are saved separately in the array 
    vector<vector<double> > thTCCutPairClusterStCharge[12]; // Strip charges for X and Y clusters of each pair; All pairs are saved as a nested vector; X and Y are saved separately in the array 
    vector<double> thTCCutPairPosCentroid[12]; //Position determined by centroid method for all pairs
    vector<double> thTCCutPairPosErrCentroid[12]; //Position error determined by centroid method for all pairs
    bool thTCCutLeftFlagSecond; // If there are tracks kept at the left of detector after applying the second layer cuts; Corresponds indices 0, 2, 4, 6, 8 and 10 for arrays
    bool thTCCutRightFlagSecond; // If there are tracks kept at the right of detector after applying the second layer cuts; Corresponds indices 1, 3, 5, 7, 9 and 11 for arrays

    //// Specified variables for pair cut
    UInt_t pairCutNumCluster[12]; //Number of clusters 
    vector<double> pairCutClusterTC[12] ;//Total charge for each cluster; All clusters are saved as a vector
    vector<UInt_t> pairCutClusterNumSt[12]; // Number of strips for each cluster; All clusters are saved as a vector
    vector<vector<int> > pairCutClusterStNum[12]; // Strip numbers saved as a vector for each clustr; All clusters are saved as a nested vector
    vector<vector<double> > pairCutClusterStCharge[12]; // Strip charges saved as a vector for each clustr; All clusters are saved as a nested vector
    vector<double> pairCutClusterPosCentroid[12]; //Position determined by centroid method for each cluster; All clusters are saved as a vector
    vector<double> pairCutClusterPosErrCentroid[12]; //Position error determined by centroid method for each cluster; All clusters are saved as a vector
    bool pairCutLeftFlagFirst; // If there are tracks kept at the left of detector after applying the first layer cuts; Corresponds indices 0, 2, 4, 6, 8 and 10 for array variables
    bool pairCutRightFlagFirst; // If there are tracks kept at the right of detector after applying the first layer cuts; Corresponds indices 1, 3, 5, 7, 9 and 11 for array variables
    UInt_t pairCutNumPair[12]; // Number of pairs; Number of pairs are the same for relevant indices, such as 0&2, 1&3, 4&6, 5&7, 8&10 and 9&11. 
    vector<UInt_t> pairCutPairClusterNum[12]; // X and Y cluster numbers for each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<double> pairCutPairClusterTC[12]; // Total charge of X and Y clusters for each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<UInt_t> pairCutPairClusterNumSt[12]; // Number of strips for X and Y cluster of each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<vector<int> > pairCutPairClusterStNum[12]; // Strip numbers for X and Y clusters of each pair; All pairs are saved as a nested vector; X and Y are saved separately in the array 
    vector<vector<double> > pairCutPairClusterStCharge[12]; // Strip charges for X and Y clusters of each pair; All pairs are saved as a nested vector; X and Y are saved separately in the array 
    vector<double> pairCutPairPosCentroid[12]; //Position determined by centroid method for all pairs
    vector<double> pairCutPairPosErrCentroid[12]; //Position error determined by centroid method for all pairs
    bool pairCutLeftFlagSecond; // If there are tracks kept at the left of detector after applying the second layer cuts; Corresponds indices 0, 2, 4, 6, 8 and 10 for arrays
    bool pairCutRightFlagSecond; // If there are tracks kept at the right of detector after applying the second layer cuts; Corresponds indices 1, 3, 5, 7, 9 and 11 for arrays

    //// Specified variables for threshold, total-charge and pair cuts
    UInt_t thTCPairCutNumCluster[12]; //Number of clusters 
    vector<double> thTCPairCutClusterTC[12] ;//Total charge for each cluster; All clusters are saved as a vector
    vector<UInt_t> thTCPairCutClusterNumSt[12]; // Number of strips for each cluster; All clusters are saved as a vector
    vector<vector<int> > thTCPairCutClusterStNum[12]; // Strip numbers saved as a vector for each clustr; All clusters are saved as a nested vector
    vector<vector<double> > thTCPairCutClusterStCharge[12]; // Strip charges saved as a vector for each clustr; All clusters are saved as a nested vector
    vector<double> thTCPairCutClusterPosCentroid[12]; //Position determined by centroid method for each cluster; All clusters are saved as a vector
    vector<double> thTCPairCutClusterPosErrCentroid[12]; //Position error determined by centroid method for each cluster; All clusters are saved as a vector
    bool thTCPairCutLeftFlagFirst; // If there are tracks kept at the left of detector after applying the first layer cuts; Corresponds indices 0, 2, 4, 6, 8 and 10 for array variables
    bool thTCPairCutRightFlagFirst; // If there are tracks kept at the right of detector after applying the first layer cuts; Corresponds indices 1, 3, 5, 7, 9 and 11 for array variables
    UInt_t thTCPairCutNumPair[12]; // Number of pairs; Number of pairs are the same for relevant indices, such as 0&2, 1&3, 4&6, 5&7, 8&10 and 9&11. 
    vector<UInt_t> thTCPairCutPairClusterNum[12]; // X and Y cluster numbers for each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<double> thTCPairCutPairClusterTC[12]; // Total charge of X and Y clusters for each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<UInt_t> thTCPairCutPairClusterNumSt[12]; // Number of strips for X and Y cluster of each pair; All pairs are saved as a vector; X and Y are saved separately in the array
    vector<vector<int> > thTCPairCutPairClusterStNum[12]; // Strip numbers for X and Y clusters of each pair; All pairs are saved as a nested vector; X and Y are saved separately in the array 
    vector<vector<double> > thTCPairCutPairClusterStCharge[12]; // Strip charges for X and Y clusters of each pair; All pairs are saved as a nested vector; X and Y are saved separately in the array 
    vector<double> thTCPairCutPairPosCentroid[12]; //Position determined by centroid method for all pairs
    vector<double> thTCPairCutPairPosErrCentroid[12]; //Position error determined by centroid method for all pairs
    bool thTCPairCutLeftFlagSecond; // If there are tracks kept at the left of detector after applying the second layer cuts; Corresponds indices 0, 2, 4, 6, 8 and 10 for arrays
    bool thTCPairCutRightFlagSecond; // If there are tracks kept at the right of detector after applying the second layer cuts; Corresponds indices 1, 3, 5, 7, 9 and 11 for arrays

    void init();   
    CRTCaliMwpc();
    virtual ~CRTCaliMwpc();
    ClassDef(CRTCaliMwpc,1);
};

//// To build the class CRTCookMwpc, a programming is processed in order.
// 1. C2X strip stransformation, pedestal subtraction
// 2. cable swapping for gap 7 of C3XL and gaps 7 & 11 of C4XR
// 3. total-charge and threshold cuts
// 4. gap-dependent gain calibration
// 5. gap-depedent event separation due to dead strip
// 6. Position determination using the centroid method 

////V1: Gap number of MWPCs is determined by Sebastien's table from TRIUMF analysis. In the table, one good event corresponds only one track. 
class CRTCookMwpcV1:public CRTBase{
  public:
    UInt_t run;
    UInt_t event;
    vector<int> adcValueOriginal[12];  // Orignal adc value
    vector<double> adcValuePedSub[12]; // Adc value afer applying pedestal subtraction and C2X strip transformation; Negative value is set as 0 after subtraction; The following variables are assigned based on these adc values

    //Cable Swapping
    UInt_t numClusterC3XLGap7Swap; //Number of clusters after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 
    vector<double> clusterTCC3XLGap7Swap;//Total charge for each cluster; All clusters are saved as a vector after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 
    vector<UInt_t> clusterNumStC3XLGap7Swap; // Number of strips for each cluster; All clusters are saved as a vector after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 
    vector<vector<int> > clusterStNumC3XLGap7Swap; // Strip numbers saved as a vector for each clustr; All clusters are saved as a nested vector after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 
    vector<vector<double> > clusterStChargeC3XLGap7Swap; // Strip charges saved as a vector for each clustr; All clusters are saved as a nested vector after replacing C3X3 and C3X4 by C4X3 and C4X4, respectively 

    UInt_t numClusterC4XLGap7Swap; //Number of clusters after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 
    vector<double> clusterTCC4XLGap7Swap;//Total charge for each cluster; All clusters are saved as a vector after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 
    vector<UInt_t> clusterNumStC4XLGap7Swap; // Number of strips for each cluster; All clusters are saved as a vector after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 
    vector<vector<int> > clusterStNumC4XLGap7Swap; // Strip numbers saved as a vector for each clustr; All clusters are saved as a nested vector after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 
    vector<vector<double> > clusterStChargeC4XLGap7Swap; // Strip charges saved as a vector for each clustr; All clusters are saved as a nested vector after replacing C4X3 and C4X4 by C3X3 and C3X4, respectively 

    UInt_t numClusterC4XLGap11Swap; //Number of clusters after switching C4X3 and C4X4 
    vector<double> clusterTCC4XLGap11Swap;//Total charge for each cluster; All clusters are saved as a vector after switching C4X3 and C4X4 
    vector<UInt_t> clusterNumStC4XLGap11Swap; // Number of strips for each cluster; All clusters are saved as a vector after switching C4X3 and C4X4 
    vector<vector<int> > clusterStNumC4XLGap11Swap; // Strip numbers saved as a vector for each clustr; All clusters are saved as a nested vector after switching C4X3 and C4X4
    vector<vector<double> > clusterStChargeC4XLGap11Swap; // Strip charges saved as a vector for each clustr; All clusters are saved as a nested vector after switching C4X3 and C4X4 
    //Cable Swapping

    UInt_t numCluster[12]; //Number of clusters 
    vector<int>  clusterGapNum[12]; // Gap Number for each cluster; All clusters are saved as a vector 
    vector<double> clusterTC[12];//Total charge for each cluster; All clusters are saved as a vector
    vector<UInt_t> clusterNumSt[12]; // Number of strips for each cluster; All clusters are saved as a vector
    vector<vector<int> > clusterStNum[12]; // Strip numbers saved as a vector for each clustr; All clusters are saved as a nested vector
    vector<vector<double> > clusterStCharge[12]; // Strip charges saved as a vector for each clustr; All clusters are saved as a nested vector

    //After gap-dependent cable swapping and gain calibration, positions are determined
    int  gapNum; // Gap Number of events from TOF2 information 
    double totalCharge[12];
    UInt_t numSt[12];
    vector<int> stNum[12];
    vector<double> stCharge[12];
    double posCentroid[12]; 
    double posErrCentroid[12]; 
    double posChargeRatio[12];
    double posChargeRatioLocal[12];
    double posC2X;
    double posC2Y;
    double posC2Z;
    double posC3X;
    double posC3Y;
    double posC3Z;
    double posC4X;
    double posC4Y;
    double posC4Z;
    int flagLeftRight; // 0:default,event is bad; -1:left; 1:right
    int flagDeadStripCluster; // 0: default; 1: cluster around a dead strip of C4XR
    void init();
    CRTCookMwpcV1();
    virtual ~CRTCookMwpcV1();
    ClassDef(CRTCookMwpcV1,1);
};

class CRTCookMwpcV2:public CRTBase{
  public:
    UInt_t run;
    UInt_t event;
    vector<double> posC2X[12];
    vector<double> posC2Y[12];
    vector<double> posC2Z[12];
    vector<double> posC3X[12];
    vector<double> posC3Y[12];
    vector<double> posC3Z[12];
    vector<double> posC4X[12];
    vector<double> posC4Y[12];
    vector<double> posC4Z[12];
    vector<bool> flagDeadStripClusterC4XR[6];
    void init();
    CRTCookMwpcV2();
    virtual ~CRTCookMwpcV2();
    ClassDef(CRTCookMwpcV2,1);
};

#endif
