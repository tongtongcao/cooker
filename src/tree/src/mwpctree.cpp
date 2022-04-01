#include "mwpctree.h"
void CRTCaliMwpc::init(){
  for(int i=0;i<12;i++){
    adcValueOriginal[i].clear();
    adcValuePedSub[i].clear();
  }

  noCutLeftFlagFirst=0;
  noCutRightFlagFirst=0;
  noCutLeftFlagSecond=0;
  noCutRightFlagSecond=0;
  for(int i=0;i<12;i++){
    noCutClusterTC[i].clear();
    noCutClusterNumSt[i].clear();
    noCutClusterStNum[i].clear();
    noCutClusterStCharge[i].clear();
    noCutClusterPosCentroid[i].clear();
    noCutClusterPosErrCentroid[i].clear();
    noCutPairClusterNum[i].clear();
    noCutPairClusterTC[i].clear();
    noCutPairClusterNumSt[i].clear();
    noCutPairClusterStNum[i].clear();
    noCutPairClusterStCharge[i].clear();
    noCutPairPosCentroid[i].clear();
    noCutPairPosErrCentroid[i].clear();
    noCutNumCluster[i]=0;
    noCutNumPair[i]=0;
  }

  thTCCutLeftFlagFirst=0;
  thTCCutRightFlagFirst=0;
  thTCCutLeftFlagSecond=0;
  thTCCutRightFlagSecond=0;
  for(int i=0;i<12;i++){
    thTCCutClusterTC[i].clear();
    thTCCutClusterNumSt[i].clear();
    thTCCutClusterStNum[i].clear();
    thTCCutClusterStCharge[i].clear();
    thTCCutClusterPosCentroid[i].clear();
    thTCCutClusterPosErrCentroid[i].clear();
    thTCCutPairClusterNum[i].clear();
    thTCCutPairClusterTC[i].clear();
    thTCCutPairClusterNumSt[i].clear();
    thTCCutPairClusterStNum[i].clear();
    thTCCutPairClusterStCharge[i].clear();
    thTCCutPairPosCentroid[i].clear();
    thTCCutPairPosErrCentroid[i].clear();
    thTCCutNumCluster[i]=0;
    thTCCutNumPair[i]=0;
  }

  pairCutLeftFlagFirst=0;
  pairCutRightFlagFirst=0;
  pairCutLeftFlagSecond=0;
  pairCutRightFlagSecond=0;
  for(int i=0;i<12;i++){
    pairCutClusterTC[i].clear();
    pairCutClusterNumSt[i].clear();
    pairCutClusterStNum[i].clear();
    pairCutClusterStCharge[i].clear();
    pairCutClusterPosCentroid[i].clear();
    pairCutClusterPosErrCentroid[i].clear();
    pairCutPairClusterNum[i].clear();
    pairCutPairClusterTC[i].clear();
    pairCutPairClusterNumSt[i].clear();
    pairCutPairClusterStNum[i].clear();
    pairCutPairClusterStCharge[i].clear();
    pairCutPairPosCentroid[i].clear();
    pairCutPairPosErrCentroid[i].clear();
    pairCutNumCluster[i]=0;
    pairCutNumPair[i]=0;
  }

  thTCPairCutLeftFlagFirst=0;
  thTCPairCutRightFlagFirst=0;
  thTCPairCutLeftFlagSecond=0;
  thTCPairCutRightFlagSecond=0;
  for(int i=0;i<12;i++){
    thTCPairCutClusterTC[i].clear();
    thTCPairCutClusterNumSt[i].clear();
    thTCPairCutClusterStNum[i].clear();
    thTCPairCutClusterStCharge[i].clear();
    thTCPairCutClusterPosCentroid[i].clear();
    thTCPairCutClusterPosErrCentroid[i].clear();
    thTCPairCutPairClusterNum[i].clear();
    thTCPairCutPairClusterTC[i].clear();
    thTCPairCutPairClusterNumSt[i].clear();
    thTCPairCutPairClusterStNum[i].clear();
    thTCPairCutPairClusterStCharge[i].clear();
    thTCPairCutPairPosCentroid[i].clear();
    thTCPairCutPairPosErrCentroid[i].clear();
    thTCPairCutNumCluster[i]=0;
    thTCPairCutNumPair[i]=0;
  }

  //Cable swapping
  thTCCutNumClusterC3XLGap7Swap=0;
  thTCCutClusterTCC3XLGap7Swap.clear();
  thTCCutClusterNumStC3XLGap7Swap.clear();
  thTCCutClusterStNumC3XLGap7Swap.clear();
  thTCCutClusterStChargeC3XLGap7Swap.clear();
  thTCCutClusterPosCentroidC3XLGap7Swap.clear();
  thTCCutClusterPosErrCentroidC3XLGap7Swap.clear();

  thTCCutNumClusterC4XLGap7Swap=0;
  thTCCutClusterTCC4XLGap7Swap.clear();
  thTCCutClusterNumStC4XLGap7Swap.clear();
  thTCCutClusterStNumC4XLGap7Swap.clear();
  thTCCutClusterStChargeC4XLGap7Swap.clear();
  thTCCutClusterPosCentroidC4XLGap7Swap.clear();
  thTCCutClusterPosErrCentroidC4XLGap7Swap.clear();

  thTCCutNumClusterC4XLGap11Swap=0;
  thTCCutClusterTCC4XLGap11Swap.clear();
  thTCCutClusterNumStC4XLGap11Swap.clear();
  thTCCutClusterStNumC4XLGap11Swap.clear();
  thTCCutClusterStChargeC4XLGap11Swap.clear();
  thTCCutClusterPosCentroidC4XLGap11Swap.clear();
  thTCCutClusterPosErrCentroidC4XLGap11Swap.clear();
  //Cable swapping
}
CRTCaliMwpc::CRTCaliMwpc(){}
CRTCaliMwpc::~CRTCaliMwpc(){};
ClassImp(CRTCaliMwpc);

void CRTCookMwpcV1::init(){
  for(int i=0;i<12;i++){
    adcValueOriginal[i].clear();
    adcValuePedSub[i].clear();
  }

  for(int i=0;i<12;i++){
    clusterGapNum[i].clear();
    clusterTC[i].clear();
    clusterNumSt[i].clear();
    clusterStNum[i].clear();
    clusterStCharge[i].clear();
    numCluster[i]=0;
  }
  //Cable Swapping
  numClusterC3XLGap7Swap=0;
  clusterTCC3XLGap7Swap.clear();
  clusterNumStC3XLGap7Swap.clear();
  clusterStNumC3XLGap7Swap.clear();
  clusterStChargeC3XLGap7Swap.clear();

  numClusterC4XLGap7Swap=0;
  clusterTCC4XLGap7Swap.clear();
  clusterNumStC4XLGap7Swap.clear();
  clusterStNumC4XLGap7Swap.clear();
  clusterStChargeC4XLGap7Swap.clear();

  numClusterC4XLGap11Swap=0;
  clusterTCC4XLGap11Swap.clear();
  clusterNumStC4XLGap11Swap.clear();
  clusterStNumC4XLGap11Swap.clear();
  clusterStChargeC4XLGap11Swap.clear();
  //Cable Swapping

  //After gap-dependent cable swapping and gain calibration, positions are determined
  gapNum=-1;
  flagDeadStripCluster=0;
  flagLeftRight=0;
  for(int i=0;i<12;i++){
    totalCharge[i]=-10000;
    numSt[i]=-10000;
    stNum[i].clear();
    stCharge[i].clear();
    posCentroid[i]=-10000;
    posErrCentroid[i]=-10000;       
    posChargeRatio[i]=-10000;
    posChargeRatioLocal[i]=-10000;
  }
  posC2X=-10000;
  posC2Y=-10000;
  posC2Z=-10000;
  posC3X=-10000;
  posC3Y=-10000;
  posC3Z=-10000;
  posC4X=-10000;
  posC4Y=-10000;
  posC4Z=-10000;

}
CRTCookMwpcV1::CRTCookMwpcV1(){}
CRTCookMwpcV1::~CRTCookMwpcV1(){};
ClassImp(CRTCookMwpcV1);

void CRTCookMwpcV2::init(){
  for(int i=0;i<12;i++){
    posC2X[i].clear();
    posC2Y[i].clear();
    posC2Z[i].clear();
    posC3X[i].clear();
    posC3Y[i].clear();
    posC3Z[i].clear();
    posC4X[i].clear();
    posC4Y[i].clear();
    posC4Z[i].clear();
  }
  for(int i=0;i<6;i++){
    flagDeadStripClusterC4XR[i].clear();
  }
}
CRTCookMwpcV2::CRTCookMwpcV2(){}
CRTCookMwpcV2::~CRTCookMwpcV2(){};
ClassImp(CRTCookMwpcV2);




