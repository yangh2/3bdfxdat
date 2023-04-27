#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hist.h"
#include "3bdf.h"
#include "integrate.h"
#include "random_16807.h"

#include "supp.h"

int PDF2Prob(hist_p_t supp, pdf_p_t pdf, prob_p_t *prob_p){
  int i, j, thresh;
  *prob_p = (prob_p_t) malloc(sizeof(prob_t));
  prob_p_t prob=*prob_p;

  prob->ntype=pdf->ntype;
  prob->dr = pdf->dr;
  prob->nbin = ceil((supp->rmax-supp->rmin)/pdf->dr)+1;
  thresh=floor(supp->rmin/pdf->dr);
  prob->rmin = thresh*pdf->dr;
  prob->rmax = prob->rmin + prob->dr * prob->nbin;
    
  prob->bins = (double **) malloc(sizeof(double*)*prob->ntype*prob->ntype);
  for (i=0;i<prob->ntype*prob->ntype;i++)
    prob->bins[i] = (double *) malloc(sizeof(double)*prob->nbin);
  
  double val, r;
  for (j=0;j<prob->ntype*prob->ntype*prob->ntype;j++){
    prob->bins[j][0] = pdf->bins[j][thresh];
    for (int i=1;i<prob->nbin;i++){
      r = (i+thresh)*prob->dr;
      val = pdf->bins[j][i+thresh]*r;
      prob->bins[j][i] = prob->bins[j][i-1]+val;
    }
    val = prob->bins[j][prob->nbin-1];
    for (int i=0;i<prob->nbin;i++){
      prob->bins[j][i] /= val;
    }
  }
  return 0;
  
}

double SampleRadius(prob_p_t prob, int type_id){
  int left, mid ,right;
  double seed1=(double)random_16807()/myrandom_max;
  double seed2=(double)random_16807()/myrandom_max;
  left = 0; right=prob->nbin-1;
  while (right-left>1){
    mid = (left+right)/2;
    if (prob->bins[type_id][mid] > seed1)
      right=mid;
    else
      left=mid;
  }
  //printf("%lf\n", prob->rmin+prob->dr*(left+seed2));
  return (prob->rmin+prob->dr*(left+seed2));
}

int PDF2Hist(prob_p_t prob, hist_p_t hist){
  int ntype;
  int idx, it, ir1, ir2, ir12, ir23, ir31;
  double r1, r2, r3, t12, t23, t31;
  int t1, t2, t3;

  ntype = prob->ntype;
  
  for (t1=0;t1<prob->ntype;t1++)
    for (t2=0;t2<prob->ntype;t2++)
      for (t3=0;t3<prob->ntype;t3++){
	  r1 = SampleRadius(prob, t1*prob->ntype+t2);
	  r2 = SampleRadius(prob, t1*prob->ntype+t3);
	  r3 = SampleRadius(prob, t3*prob->ntype+t2);

	  if ((r1+r2<r3) || (r1+r3<r2) || (r2+r3<r1)) continue;
	  if ((r1>hist->rmax+hist->dr)||(r1<hist->rmin-hist->dr)) continue;
	  if ((r2>hist->rmax+hist->dr)||(r2<hist->rmin-hist->dr)) continue;
	  if ((r3>hist->rmax+hist->dr)||(r3<hist->rmin-hist->dr)) continue;

	  t12 = (r1*r1+r2*r2-r3*r3)/2/r1/r2;
	  t31 = (r1*r1-r2*r2+r3*r3)/2/r1/r3;
	  t23 = (-r1*r1+r2*r2+r3*r3)/2/r3/r2;
	  
	  ir12 = floor((r1-hist->rmin)/hist->dr);
	  ir23 = floor((r2-hist->rmin)/hist->dr);
	  ir31 = floor((r3-hist->rmin)/hist->dr);
	  
	  ir1=ir12>ir23?ir12:ir23;
	  ir2=ir12<ir23?ir12:ir23;
	  idx=MapIndex(ir1, ir2);
	  if (idx >= hist->nbin) continue;
	  it = floor((t12-hist->tmins[idx])/hist->dt[idx]);
	  if ((it >= hist->nt[idx]) || (it <0)) continue;
	  hist->bins[t1*ntype*ntype+t2*ntype+t3][idx][it]++;
	
	  ir1=ir31>ir23?ir31:ir23;
	  ir2=ir31<ir23?ir31:ir23;
	  idx=MapIndex(ir1, ir2);
	  if (idx >= hist->nbin) continue;;
	  it = floor((t23-hist->tmins[idx])/hist->dt[idx]);
	  if ((it >= hist->nt[idx]) || (it <0)) continue;
	  hist->bins[t1*ntype*ntype+t2*ntype+t3][idx][it]++;

	  ir1=ir12>ir31?ir12:ir31;
	  ir2=ir12<ir31?ir12:ir31;
	  idx=MapIndex(ir1, ir2);
	  if (idx >= hist->nbin) continue;;
	  it = floor((t31-hist->tmins[idx])/hist->dt[idx]);
	  if ((it >= hist->nt[idx]) || (it <0)) continue;
	  hist->bins[t1*ntype*ntype+t2*ntype+t3][idx][it]++;
	  hist->nsample++;
      }
  return 0; 
}

int Hist2Supp(hist_p_t hist, tbf_p_t *supp_p){
  *supp_p = (supp_p_t) malloc(sizeof(supp_t));
  supp_p_t supp=*supp_p;

  /* copy parameters from hist to supp */
  supp->nsample = hist->nsample;
  supp->ntype = hist->ntype;
  supp->rmin = hist->rmin;
  supp->rmax = hist->rmax;
  supp->dr = hist->dr;
  supp->nr = hist->nr;
  supp->nbin = hist->nbin;

  /* allocate memory for 3bdf */
  int nbin=supp->nbin,i;
  supp->tmins = (double *) malloc(sizeof(double)*nbin);
  supp->tmaxs = (double *) malloc(sizeof(double)*nbin);
  supp->dt = (double *) malloc(sizeof(double)*nbin);
  supp->nt = (long *) malloc(sizeof(long)*nbin);

  for (i=0;i<nbin;i++){
    supp->tmins[i] = hist->tmins[i];
    supp->tmaxs[i] = hist->tmaxs[i];
    supp->dt[i] = hist->dt[i];
    supp->nt[i] = hist->nt[i];
  }

  int ntype=supp->ntype, ir1, ir2, idx,nr=supp->nr, it;
  double ***f3s, *vols, r1, r2 ,dr, v0=0;;
  dr = supp->dr;
  f3s = (double ***) malloc(sizeof(double **)*ntype*ntype*ntype);
  vols = (double *) malloc(sizeof(double)*nbin);

  double ntot=0;
  for (i=0;i<ntype*ntype*ntype;i++){
    f3s[i] = (double **) malloc(sizeof(double*)*nbin);
    for (ir1=0;ir1<nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
	idx=MapIndex(ir1, ir2);
	r1=supp->rmin+dr*ir1; r2=supp->rmin+dr*ir2;

	f3s[i][idx] = (double *) malloc(sizeof(double)*supp->nt[idx]);
	vols[idx] = (double)8.0*pi*pi/9
	  *((r1+dr)*(r1+dr)*(r1+dr)-r1*r1*r1)
	  *((r2+dr)*(r2+dr)*(r2+dr)-r2*r2*r2)
	  *supp->dt[idx];
	/* caluclate volume elements and normalize supp */
	if (ir1==ir2)
	  for (it=0;it<supp->nt[idx];it++){
	    f3s[i][idx][it] = (double)hist->bins[i][idx][it]/vols[idx]*2;
	    v0 += vols[idx];
	    ntot += hist->bins[i][idx][it]*2;
	  }
	else 
	  for (it=0;it<supp->nt[idx];it++){
	    f3s[i][idx][it] = (double)hist->bins[i][idx][it]/vols[idx];
	    v0 += vols[idx];
	    ntot += hist->bins[i][idx][it];
	  }
      }
  }

  for (i=0;i<ntype*ntype*ntype;i++){
    for (ir1=0;ir1<nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
	idx=MapIndex(ir1, ir2);
	r1=supp->rmin+dr*ir1; r2=supp->rmin+dr*ir2;

	/* caluclate volume elements and normalize supp */
	for (it=0;it<supp->nt[idx];it++)
	  f3s[i][idx][it] = f3s[i][idx][it]/ntot*v0;
      }
  }

  supp->f3s=f3s;
  supp->vols=vols;
  return 0;
}

