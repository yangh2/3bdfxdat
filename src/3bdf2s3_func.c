#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "3bdf.h"

double VolumeElement();

int Normalization(hist_p_t hist){
  int i, ir1, ir2, idx, nr, ntype;
  nr=hist->nr;
  ntype=hist->ntype;
  long ***bins=hist->bins;
  double rmin, rmax, r1, r2, dr;
  int nbin = hist->nbin;
  rmin=hist->rmin; rmax=hist->rmax;
  dr = hist->dr;
  double *tmins, *tmaxs;
  tmins=hist->tmins; tmaxs = hist->tmaxs;
  
  for (i=0;i<ntype*ntype*ntype;i++){
    bins[i] = (long **) malloc(sizeof(long*)*nbin);
    
    for (ir1=0;ir1<nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
	idx=MapIndex(ir1, ir2);
      	for (it=0;it<nt[idx];it++)
	  bins[i]
      }
  }

  return 0;
}
