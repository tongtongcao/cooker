#!/bin/bash
pwdd=`pwd`
source ../myroot/root_5/root/build/bin/thisroot.sh
cd build/
make clean; rm -r CMake* cmake_install.cmake ; clear
cmake ../ && make -j8 && make install
cd $pwdd
#./tof_m2.sh
#./build/bin/cooker recipes/tracking/tracking.xml /data/trek/E36/gautam/crnRoot/g4pseudo/cookG4/mergedTargetTofMwpc_101.root New_Rootfiles/tracking_101.root
#./build/bin/cooker $tracking_recipes/basicDataType.xml New_Rootfiles/traking_thir_$1.root:/data/trek/E36/gautam/crnRoot/data2Root/Merged/Run$1MS.root New_Rootfiles/basicDataType/basicDataType_thir_$1.root
#./build/bin/cooker $tracking_recipes/basicDataType.xml New_Rootfiles/traking_thir_$1.root:/data/trek/E36/gautam/crnRoot/data2Root/Merged/Run$1MS.root New_Rootfiles/basicDataType/basicDataType_thir_$1.root
#./build/bin/cooker $tracking_recipes/basicDataType.xml New_Rootfiles/traking_thir_$1.root:/data/trek/E36/gautam/crnRoot/data2Root/Merged/Run$1MS.root New_Rootfiles/basicDataType/basicDataType_thir_$1.root
