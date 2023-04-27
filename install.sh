#!/bin/bash

cwd=`pwd`;

cd src
make cleanall
make

cp  init_bins.awk init_bins.sh split_xdat.awk 3bdfxdat 3bdfxdat_bundle.sh 3bdfxdat_sum 3bdf2s3 3bdf_view_theta 3bdf_print_sup 3bdfMsup pdf2supp hist2supp 3bdf2s3_eps 3bdf2s3_interp.py ../bin/

cd $cwd
