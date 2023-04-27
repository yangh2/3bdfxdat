#!/bin/bash

pdfin="3bdf_pdf.in"
symbols=`awk 'NR==6{print $0}' POSCAR`;
ntype=`awk 'NR==6{print $0}' POSCAR | wc -w`;

if [ -e $pdfin ];
then
    rm $pdfin
fi

touch $pdfin
ftmp1=`mktemp`;
ftmp2=`mktemp`;
for i in $symbols;
do
    for j in $symbols;
    do
	fname=pdf.$i$j
	awk '{print $2}' $fname > $ftmp1;
	paste  $pdfin $ftmp1 > $ftmp2
	cp $ftmp2 $pdfin
    done
done
awk '{print $1}' $fname > $ftmp1;
paste $ftmp1 $pdfin > $ftmp2
line=`wc -l $ftmp2 | awk '{print $1}'`;
rmin=`head -n 1 $ftmp1`;
rmax=`tail -n 1 $ftmp1`;
echo $ntype > $pdfin;
echo $line $rmin $rmax >> $pdfin;
cat $ftmp2 >> $pdfin

rmin=$1
rmax=$2
dr=$3
init_bins.awk  -v rmin=$rmin -v rc=$rmax -v dr=$dr POSCAR

cat POSCAR | split_xdat.awk;
mv xdat_cut1 3bdf_xdat.in
