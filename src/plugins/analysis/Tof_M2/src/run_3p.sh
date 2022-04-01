#!/bin/bash
dir_par="kmu2_plugins"
#dir_par="epPR_plugins"
#dir_par="tried_Eppeaks"
var1=""
var2=""
if [ "$1" == "Ual" ]; then
var1="Ual"
echo $var1
cp to_run/$dir_par/Tof_M2_$var1.cpp Tof_M2.cpp

else
var2="Al"
echo $var2
echo "$dir_par"
cp to_run/$dir_par/Tof_M2_$var2.cpp Tof_M2.cpp
fi

