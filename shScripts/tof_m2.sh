#!/bin/bash
pwdd=`pwd`
#cd build/
#make clean; rm -r CMake* cmake_install.cmake ; clear
#cmake ../ && make -j8 && make install
cd $pwdd
cor_cal="$1"
src_dir="/data/trek/E36/gautam/crnRoot/cooker/src/plugins/analysis/Tof_M2/src"
cd $src_dir
#./run_3p.sh $cor_cal
cd $pwdd

raw_dat="/data/trek/E36/gautam/crnRoot/data2Root/Merged/merged"
raw_ac="/data/trek/E36/gautam/crnRoot/cooker/output_root"
raw_pgc="/data/trek/E36/gautam/crnRoot/cooker/output_root/PgC"
track_root="../../E36_Ana/forCooker/Rootfiles"
cooked_tof="/data/trek/E36/gautam/E36_Ana/cookedtof/Rootfiles"
output_root="New_Rootfiles"
echo $raw_dat
echo $track_root
echo $cooked_tof
echo $output_root
./make_cook.sh
./build/bin/cooker recipes/analysis/Tof_M2/tofm2.xml $raw_pgc/kmu2_PgC.root:$track_root/track_kmu2.root:$cooked_tof/cookedkmu2.root New_Rootfiles/tofM2_kmu2_${cor_cal}.root
root -l New_Rootfiles/tofM2_kmu2_${cor_cal}.root
#./build/bin/cooker recipes/analysis/Tof_M2/tofm2.xml $raw_pgc/kmu2_PgC.root:$track_root/track_kmu2.root:$cooked_tof/cookedkmu2.root New_Rootfiles/tofM2_epkmu2c_${cor_cal}.root
#root -l New_Rootfiles/tofM2_epkmu2c_${cor_cal}.root

#./build/bin/cooker recipes/analysis/Tof_M2/tofm2.xml $raw_pgc/kmu2_PgC.root:$track_root/track_epPR.root:$cooked_tof/cookedepPR.root New_Rootfiles/tofM2_CepPRepCor_${cor_cal}.root
#root -l New_Rootfiles/tofM2_CepPRepCor_${cor_cal}.root
#./build/bin/cooker recipes/analysis/Tof_M2/tofm2.xml $raw_pgc/ke3_PgC.root:$track_root/track_ke3.root:$cooked_tof/cookedke3.root New_Rootfiles/tofM2_ke3_${cor_cal}.root
#root -l New_Rootfiles/tofM2_ke3_${cor_cal}.root
#./build/bin/cooker recipes/analysis/Tof_M2/tofm2.xml /data/trek/E36/gautam/crnRoot/cooker/output_root/PgC/ke3_PgC.root:../../E36_Ana/forCooker/Rootfiles/track_ke3.root:/data/trek/E36/gautam/E36_Ana/cookedtof/Rootfiles/cookedke3.root New_Rootfiles/tofM2_ke3.root
#./build/bin/cooker recipes/analysis/Tof_M2/tofm2.xml /data/trek/E36/gautam/crnRoot/data2Root/Merged/merged/Run_ke3.root:../../E36_Ana/forCooker/Rootfiles/track_ke3.root:/data/trek/E36/gautam/E36_Ana/cookedtof/Rootfiles/cookedke3.root  New_Rootfiles/tofM2_ke3.root
