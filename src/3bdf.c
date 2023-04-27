#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hist.h"
#include "3bdf.h"

double pi=3.14159265359;
double kb=8.617333262e-5;	/* eV/K */
double hc=1.23984193e4;		/* eV/angstrom */
double massunit=931.49410242e6; /* ev/c^2 */

int InitTBF(xdat_p_t xdat, hist_p_t hist, tbf_p_t *tbf_p){
  *tbf_p = (tbf_p_t) malloc(sizeof(tbf_t));
  tbf_p_t tbf=*tbf_p;

  /* copy parameters from hist to tbf */
  tbf->nsample = hist->nsample;
  tbf->ntype = hist->ntype;
  tbf->rmin = hist->rmin;
  tbf->rmax = hist->rmax;
  tbf->dr = hist->dr;
  tbf->nr = hist->nr;
  tbf->nbin = hist->nbin;

  /* allocate memory for 3bdf */
  int nbin=tbf->nbin,i;
  tbf->tmins = (double *) malloc(sizeof(double)*nbin);
  tbf->tmaxs = (double *) malloc(sizeof(double)*nbin);
  tbf->dt = (double *) malloc(sizeof(double)*nbin);
  tbf->nt = (long *) malloc(sizeof(long)*nbin);

  for (i=0;i<nbin;i++){
    tbf->tmins[i] = hist->tmins[i];
    tbf->tmaxs[i] = hist->tmaxs[i];
    tbf->dt[i] = hist->dt[i];
    tbf->nt[i] = hist->nt[i];
  }

  int ntype=tbf->ntype, ir1, ir2, idx,nr=tbf->nr, it;
  double ***f3s, *vols, r1, r2 ,dr;
  dr = tbf->dr;
  f3s = (double ***) malloc(sizeof(double **)*ntype*ntype*ntype);
  vols = (double *) malloc(sizeof(double)*nbin);

  double v0 = (4.0/3*tbf->rmax*tbf->rmax*tbf->rmax*pi); v0=v0*v0;
  for (i=0;i<ntype*ntype*ntype;i++){
    f3s[i] = (double **) malloc(sizeof(double*)*nbin);
    for (ir1=0;ir1<nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
	idx=MapIndex(ir1, ir2);
	r1=tbf->rmin+dr*ir1; r2=tbf->rmin+dr*ir2;

	f3s[i][idx] = (double *) malloc(sizeof(double)*tbf->nt[idx]);
	vols[idx] = (double)8.0*pi*pi/9
	  *((r1+dr)*(r1+dr)*(r1+dr)-r1*r1*r1)
	  *((r2+dr)*(r2+dr)*(r2+dr)-r2*r2*r2)
	  *tbf->dt[idx];
	/* caluclate volume elements and normalize tbf */
	if (ir1==ir2)
          for (it=0;it<tbf->nt[idx];it++){
            f3s[i][idx][it] = (double)hist->bins[i][idx][it]/vols[idx]/tbf->nsample
              *xdat->vol*xdat->vol/xdat->natom/(xdat->natom-1)/(xdat->natom-2)*2;
          }
        else 
          for (it=0;it<tbf->nt[idx];it++){
            f3s[i][idx][it] = (double)hist->bins[i][idx][it]/vols[idx]/tbf->nsample
              *xdat->vol*xdat->vol/xdat->natom/(xdat->natom-1)/(xdat->natom-2);
          }

      }
  }

  tbf->f3s=f3s;
  tbf->vols=vols;
  return 0;
}

int FreeTBF(tbf_p_t tbf){
  int i,ntype=tbf->ntype, ir1, ir2, idx,nr=tbf->nr;
  for (i=0;i<ntype*ntype*ntype;i++){
    for (ir1=0;ir1<nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
	idx=MapIndex(ir1, ir2);
	free(tbf->f3s[i][idx]);
      }
    free(tbf->f3s[i]);
  }
  free(tbf->f3s);
  if ( tbf->vols != NULL ) free(tbf->vols);
  free(tbf->tmins);
  free(tbf->tmaxs);
  free(tbf->dt);
  free(tbf->nt);
  free(tbf);
  return 0;
}

int ReadTBF(char *fname, tbf_p_t *tbf_p){
  *tbf_p = (tbf_p_t) malloc(sizeof(tbf_t));
  tbf_p_t tbf=*tbf_p;
  
  double rmin, rmax, dr;		/* range [rmin, rmax], exclude 0 within rmin */
  int ntype, nsample;
  long nbin, *nt, idx;
  double *tmins, *tmaxs, *dt, ***f3s;
  int nr,i;
  int ir1, ir2,it;
  double r1, r2, r3max, r3min;

  FILE *fid = fopen(fname, "r");
  fscanf(fid, "%d", &nsample);
  fscanf(fid, "%d", &ntype);
  fscanf(fid, "%lf %lf %lf", &rmin, &rmax, &dr);

  tbf->nsample = nsample;
  tbf->ntype = ntype;
  tbf->rmin = rmin;
  tbf->rmax = rmax;
  tbf->dr = dr;
  
  nr=ceil((rmax-rmin)/dr);
  f3s = (double ***) malloc(sizeof(double **)*ntype*ntype*ntype);
  /* save memory using symmetry of r1 and r2,
     and only consider r1>=r2 */
  nbin = nr*(nr+1)/2;
  tmins = (double *) malloc(sizeof(double)*nbin);
  tmaxs = (double *) malloc(sizeof(double)*nbin);
  dt = (double *) malloc(sizeof(double)*nbin);
  nt = (long *) malloc(sizeof(long)*nbin);
  for (i=0;i<ntype*ntype*ntype;i++){
    f3s[i] = (double **) malloc(sizeof(double*)*nbin);
    
    for (ir1=0;ir1<nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
	idx=MapIndex(ir1, ir2);
	r1=rmin+dr*ir1; r2=rmin+dr*ir2;
	tmins[idx] = CosThetaMin(rmax, r1, r2);
	tmaxs[idx] = CosThetaMax(rmin, r1, r2);
	r3min=rmin>r1-r2-dr?rmin:r1-r2-dr;
	r3max=rmax<r1+r2+2*dr?rmax:r1+r2+2*dr;
	nt[idx] = 2*ceil((r3max-r3min)/dr);
	dt[idx] = (tmaxs[idx]-tmins[idx]+mepsilon)/nt[idx];

	f3s[i][idx]=(double *) malloc(sizeof(double)*nt[idx]);
      	for (it=0;it<nt[idx];it++)
	  fscanf(fid, "%lf", &(f3s[i][idx][it]));
      }
  }
  fclose(fid);
  tbf->nr = nr;
  tbf->f3s = f3s;
  tbf->vols = NULL;
  tbf->nbin = nbin;
  tbf->tmins = tmins;
  tbf->tmaxs = tmaxs;
  tbf->nt = nt;
  tbf->dt = dt;
  return 0;
}

int PrintTBF(char *fname, tbf_p_t tbf){
  int i, ir1, ir2, idx, it;
  
  FILE *fid = fopen(fname, "w");
  fprintf(fid, "%d\n", tbf->nsample);
  fprintf(fid, "%d\n", tbf->ntype);
  fprintf(fid, "%lf %lf %lf\n", tbf->rmin, tbf->rmax, tbf->dr);
  for (i=0;i<tbf->ntype*tbf->ntype*tbf->ntype;i++)
    for (ir1=0;ir1<tbf->nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
        idx=MapIndex(ir1, ir2);
        //r1=rmin+dr*ir1; r2=rmin+dr*ir2;
        //printf("%lf %lf %lf\n", tbf->tmins[idx], tbf->tmaxs[idx], tbf->dt[idx]);
        for (it=0;it<tbf->nt[idx];it++)
          fprintf(fid, "%.12lf ", tbf->f3s[i][idx][it]);
        fprintf(fid, "\n");
      }
  fclose(fid);
  return 0;
}

int PrintTBFatR(tbf_p_t tbf, double r){
  int i, ir1, ir2, idx, it;
  for (ir1=0;ir1<tbf->nr;ir1++)
    if (tbf->rmin+tbf->dr*ir1>=r)
      break;

  idx=MapIndex(ir1, ir1);

  for (it=0;it<tbf->nt[idx];it++){
    printf("%.12lf ", acos(tbf->tmins[idx]+tbf->dt[idx]*it));
    printf("%.12lf ", sqrt(2*r*r*(1-tbf->tmins[idx]-tbf->dt[idx]*it-tbf->dt[idx])));
    for (i=0;i<tbf->ntype*tbf->ntype*tbf->ntype;i++){
      printf("%.12lf ", tbf->f3s[i][idx][it]);
    }
    printf("\n");
  }

  return 0;
}

int PrintTBFatR1R2(tbf_p_t tbf, double r1, double r2){
  int i, ir1, ir2, idx, it;
  for (ir1=0;ir1<tbf->nr;ir1++)
    if (tbf->rmin+tbf->dr*ir1>=r1)
      break;
  for (ir2=0;ir2<tbf->nr;ir2++)
    if (tbf->rmin+tbf->dr*ir2>=r2)
      break;

  it=ir1>ir2?ir1:ir2;
  ir2=ir1<ir2?ir1:ir2;
  ir1=it;
  idx=MapIndex(ir1, ir2);

  for (it=0;it<tbf->nt[idx];it++){
    printf("%.12lf ", acos(tbf->tmins[idx]+tbf->dt[idx]*it));
    printf("%.12lf ", sqrt(r1*r1+r2*r2-2*r2*r1*(tbf->tmins[idx]+tbf->dt[idx]*it+tbf->dt[idx])));
    for (i=0;i<tbf->ntype*tbf->ntype*tbf->ntype;i++){
      printf("%.12lf ", tbf->f3s[i][idx][it]);
    }
    printf("\n");
  }

  return 0;
}

