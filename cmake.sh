#!/bin/tcsh
if ! $?ROOTSYS then
    setenv ROOTSYS /apps/root/5.34
endif
if ! $?Geant4_DIR then
    setenv Geant4_DIR ${HOME}/local/geant4.10.00.p02/lib64/Geant4-10.0.2
endif

set ld_exist=0
if $?LD_LIBRARY_PATH then
    set ld_exist=1
    set old_ld=${LD_LIBRARY_PATH}
    setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
else
    setenv LD_LIBRARY_PATH ${ROOTSYS}/lib
endif

#cmake ~/src/cooker -DCMAKE_INSTALL_PREFIX=~/local/cooker -DQT_QMAKE_EXECUTABLE=${HOME}/local/qt-everywhere-opensource-src-4.8.6/bin/qmake
cmake `dirname $0` #-DQT_QMAKE_EXECUTABLE=${HOME}/local/qt-everywhere-opensource-src-4.8.6/bin/qmake
make
make install
if $ld_exist then
    setenv LD_LIBRARY_PATH ${old_ld}
endif
