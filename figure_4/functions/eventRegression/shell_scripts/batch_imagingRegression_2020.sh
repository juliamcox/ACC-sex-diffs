#! /bin/env bash

versions="stay
"

rawFlag="1"


recs="T6/recording_20161207_154909
T19/recording_20170910_135133
T57/recording_20190731_123853
T55/recording_20190627_162644
T95/recording_20200722_141305
T94/recording_20200717_160625
101/recording_20200722_152930
T1/recording_20161204_142407
T3/recording_20161130_173418
T12/recording_20170310_125113
T14/recording_20170220_115346
T60/recording_20190625_132507
T63/recording_20190627_122511
102/recording_20200721_140825
103/recording_20200711_135354
"

method="f
"

numShuff="500
" 

for recname in $recs
do
thisrec="$recname"
for vername in $versions
do
thisver="$vername"
for methodname in $method
do
thismethod="$methodname"
echo $thisrec
echo $vername
    sbatch imagingRegression.sh $thisrec $thisver $rawFlag $thismethod $numShuff
  # bash imagingTest.sh $thisrec $thisver $rawFlag $thismethod $numShuff
    sleep 1
done
done
done