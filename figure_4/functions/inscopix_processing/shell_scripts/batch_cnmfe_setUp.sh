#! /bin/env bash


aids="T14
"
recs="recording_20170220_115346
"

#T1
#recording_20161204_142407
#T3
#recording_20161130_173418
#T12
#recording_20170310_125113
#T60
#recording_20190725_145401
#T63
#recording_20190627_122511
#102
#recording_20200721_140825
#103
#recording_20200711_135354
#101
#recording_20200722_152930
#T6
#recording_20161207_154909
#T19
#recording_20170910_135133
#T57
#recording_20190731_123853
#T55
recording_20190627_162644
#T95
#recording_20200722_141305
#T94
#recording_20200717_160625


for aidname in $aids
do
thisaid="$aidname"
for recname in $recs
do
thisrec="$recname"
echo $thisaid
echo $thisrec
    sbatch cnmfe_setUp.sh $thisaid $thisrec 
    sleep 1
done
done








