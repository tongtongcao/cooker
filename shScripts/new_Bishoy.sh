#!/bin/bash
pwdd=`pwd`
cd $pwdd

./build/bin/cooker recipes/ToF1/conversion.xml /data/trek/E36/gautam/crnRoot/data2Root/Merged/Run$1MS.root ../cooker/New_Rootfiles/TOF_MWPC/ToF1Con_$1.root 
./build/bin/cooker recipes/ToF2/conversion.xml /data/trek/E36/gautam/crnRoot/data2Root/Merged/Run$1MS.root ../cooker/New_Rootfiles/TOF_MWPC/ToF2Con_$1.root
./build/bin/cooker recipes/TOF/tofCao1.xml ../cooker/New_Rootfiles/TOF_MWPC/ToF1Con_$1.root:../cooker/New_Rootfiles/TOF_MWPC/ToF2Con_$1.root ../cooker/New_Rootfiles/TOF_MWPC/tof$1.root  

#mwpc
./build/bin/cooker recipes/MWPC/mwpcCookedDataV2.xml /data/trek/E36/gautam/crnRoot/data2Root/Merged/Run$1MS.root ../cooker/New_Rootfiles/TOF_MWPC/mwpc1_$1.root 
./build/bin/cooker recipes/MWPC/mwpcCookedDataV21.xml ../cooker/New_Rootfiles/TOF_MWPC/mwpc1_$1.root ../cooker/New_Rootfiles/TOF_MWPC/mwpc2_$1.root

#target
cd /data/trek/E36/gautam/crnRoot/Tools/Target_ToF/targetDataConvert/
python2 run_python.py $1
rm *.so *.d
root -l -q "convertTargetData.C+($1)"
cd $pwdd

cd /data/trek/E36/gautam/crnRoot/tracking_trek/
root -l -q "mergeTargetTofMwpc.C+($1)"
cd $pwdd
./build/bin/cooker recipes/tracking/tracking.xml /data/trek/E36/gautam/crnRoot/tracking_trek/New_Rootfiles/mergedTargetTofMwpcRun$1.root ../cooker/New_Rootfiles/traking_$1.root
./build/bin/cooker $tracking_recipes/basicDataType.xml ../cooker/New_Rootfiles/traking_$1.root:/data/trek/E36/gautam/crnRoot/data2Root/Merged/Run$1MS.root /data/trek/E36/gautam/crnRoot/cooker/New_Rootfiles/basicDataType/basicDataType_$1.root
