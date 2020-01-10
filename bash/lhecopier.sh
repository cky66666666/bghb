#!/bin/bash
lhepath=/mnt/d/MG5_aMC_v2_5_5
destpath=/mnt/d/work/bghb/data
for ((i=1;i<=9;i++))
do
    j=$[$i-1]
    cd $lhepath/bkg_bghb$i/Events/run_01
    gzip -d unweighted_events.lhe.gz
    mv unweighted_events.lhe bkg$j.lhe
    cp bkg$j.lhe $destpath
    rm bkg$j.lhe
done
