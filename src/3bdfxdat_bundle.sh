#!/bin/bash

nstart=$1
nend=$2
name=$3
fhist=3bdf_bins.dat_$name

drun=run_dir$name
mkdir $drun
cp 3bdf_bins.dat $drun/$fhist

for i in `seq $nstart $nend`;
do
    mv xdat_cut$i $drun;
done

sleep 3;

cd $drun
ls
3bdfxdat $fhist xdat_cut* 
sleep 3;
cp $fhist ../
cd ../

