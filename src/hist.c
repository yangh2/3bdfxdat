#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hist.h"

double mepsilon=1e-15;

double CosThetaMin(double rc, double r1, double r2){
  double val;
  val=(r1*r1+r2*r2-rc*rc)/2/r1/r2;
  val=val>-1?val:-1;
  return val;
}
double CosThetaMax(double rmin, double r1, double r2){
  double val;
  val=(r1*r1+r2*r2-rmin*rmin)/2/r1/r2;
  val=val<1?val:1;
  return val;
}

long MapIndex(long ir1, long ir2){
  long t;
  if (ir1<ir2){
    t=ir1;
    ir1=ir2;
    ir2=t;
  }
  return (ir1+1)*ir1/2+ir2;  
}

int InitXdat(char *fxdat, xdat_p_t *xdat_p){
  *xdat_p = (xdat_p_t) malloc(sizeof(xdat_t));
  xdat_p_t xdat=*xdat_p;
  int i,j, nlist[100], ntype, natom=0;
  FILE *fid = fopen(fxdat, "r");
  fscanf(fid, "%d", &ntype);
  nlist[0]=0;
  for (i=1;i<=ntype;i++){
    fscanf(fid, "%d", &(nlist[i]));
    natom += nlist[i];
    nlist[i]+=nlist[i-1];
  }
  xdat->ntype = ntype;
  xdat->natom = natom;
  xdat->typemap = (int *) malloc(sizeof(int) *natom);
  for (i=0;i<natom;i++) xdat->typemap[i]=0;
  for (i=1;i<=ntype;i++)
    for (j=nlist[i-1];j<nlist[i];j++)
      xdat->typemap[j]=i-1;
  fscanf(fid, "%lf %lf %lf", &(xdat->a1), &(xdat->a2), &(xdat->a3));
  fscanf(fid, "%lf %lf %lf", &(xdat->b1), &(xdat->b2), &(xdat->b3));
  fscanf(fid, "%lf %lf %lf", &(xdat->c1), &(xdat->c2), &(xdat->c3));

  double a1,a2,a3;
  double b1,b2,b3;
  double c1,c2,c3;
  a1=xdat->a1; b1=xdat->b1; c1=xdat->c1;
  a2=xdat->a2; b2=xdat->b2; c2=xdat->c2;
  a3=xdat->a3; b3=xdat->b3; c3=xdat->c3;
  xdat->vol=(double)(a1*(b2*c3-b3*c2)-a2*(b1*c3-b3*c1)+a3*(b1*c2-b2*c1));

  fclose(fid);
  return 0;
}

int FreeXdat(xdat_p_t xdat){
  free(xdat->typemap);
  free(xdat);
  return 0;
}

int InitHist(char *fname, hist_p_t *hist_p){
  *hist_p = (hist_p_t) malloc(sizeof(hist_t));
  hist_p_t hist=*hist_p;
  
  double rmin, rmax, dr;		/* range [rmin, rmax], exclude 0 within rmin */
  int ntype, nsample;
  long nbin, ***bins, *nt, idx;
  double *tmins, *tmaxs, *dt;
  int nr,i;
  int ir1, ir2,it;
  double r1, r2, r3max, r3min;

  FILE *fid = fopen(fname, "r");
  fscanf(fid, "%d", &nsample);
  fscanf(fid, "%d", &ntype);
  fscanf(fid, "%lf %lf %lf", &rmin, &rmax, &dr);

  hist->nsample = nsample;
  hist->ntype = ntype;
  hist->rmin = rmin;
  hist->rmax = rmax;
  hist->dr = dr;
  
  nr=ceil((rmax-rmin)/dr);
  bins = (long ***) malloc(sizeof(long **)*ntype*ntype*ntype);
  /* save memory using symmetry of r1 and r2,
     and only consider r1>=r2 */
  nbin = nr*(nr+1)/2;
  tmins = (double *) malloc(sizeof(double)*nbin);
  tmaxs = (double *) malloc(sizeof(double)*nbin);
  dt = (double *) malloc(sizeof(double)*nbin);
  nt = (long *) malloc(sizeof(long)*nbin);
  for (i=0;i<ntype*ntype*ntype;i++){
    bins[i] = (long **) malloc(sizeof(long*)*nbin);
    
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

	bins[i][idx]=(long *) malloc(sizeof(long)*nt[idx]);
      	for (it=0;it<nt[idx];it++)
	  bins[i][idx][it]=0;
      }
  }
  fclose(fid);
  hist->nr = nr;
  hist->bins = bins;
  hist->nbin = nbin;
  hist->tmins = tmins;
  hist->tmaxs = tmaxs;
  hist->nt = nt;
  hist->dt = dt;
  return 0;
}

int FreeHist(hist_p_t hist){
  int i, ntype, ir1, ir2, nr, idx;
  nr = hist->nr;
  ntype=hist->ntype;
  free(hist->tmins);
  free(hist->tmaxs);
  free(hist->dt);
  free(hist->nt);
  for (i=0;i<ntype*ntype*ntype;i++){
    for (ir1=0;ir1<nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
	idx=MapIndex(ir1, ir2);
	if (hist->bins[i][idx] != NULL)
	  free(hist->bins[i][idx]);
      }
    free(hist->bins[i]);
  }
  free(hist->bins);
  free(hist);
  return 0;
}

int ReadHist(char *fname, hist_p_t hist){
  double rmin, rmax, dr;		/* range [rmin, rmax], exclude 0 within rmin */
  int ntype, nsample;
  long nbin, ***bins, *nt, idx;
  double *tmins, *tmaxs, *dt;
  int nr,i;
  int ir1, ir2, it;
  double r1, r2, r3max, r3min;

  FILE *fid = fopen(fname, "r");
  fscanf(fid, "%d", &nsample);
  fscanf(fid, "%d", &ntype);
  fscanf(fid, "%lf %lf %lf", &rmin, &rmax, &dr);

  hist->nsample = nsample;
  hist->ntype = ntype;
  hist->rmin = rmin;
  hist->rmax = rmax;
  hist->dr = dr;
  
  nr=ceil((rmax-rmin)/dr);
  /* save memory using symmetry of r1 and r2,
     and only consider r1>=r2 */
  nbin = nr*(nr+1)/2;
  hist->nr = nr;
  bins=hist->bins;
  tmins=hist->tmins;
  tmaxs=hist->tmaxs;
  nt=hist->nt;
  dt=hist->dt;
  for (i=0;i<ntype*ntype*ntype;i++){
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

	for (it=0;it<nt[idx];it++)
	  fscanf(fid, "%ld", &(bins[i][idx][it]));
      }
  }
  fclose(fid);
  return 0;
}

int AddHist(hist_p_t hist1, hist_p_t hist2){
  /* add hist2 to hist1 */
  int ntype;
  long *nt, idx;
  int nr,i;
  int ir1, ir2,it;

  hist1->nsample += hist2->nsample;

  nr = hist1->nr;
  nt = hist1->nt;
  ntype = hist1->ntype;
  for (i=0;i<ntype*ntype*ntype;i++){
    for (ir1=0;ir1<nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
	idx=MapIndex(ir1, ir2);
	for (it=0;it<nt[idx];it++)
	  hist1->bins[i][idx][it]+=hist2->bins[i][idx][it];
      }
  }
  return 0;
}

int PrintHist(char *fname, hist_p_t hist){
  int i, ir1, ir2, idx, it;
  
  FILE *fid = fopen(fname, "w");
  fprintf(fid, "%d\n", hist->nsample);
  fprintf(fid, "%d\n", hist->ntype);
  fprintf(fid, "%lf %lf %lf\n", hist->rmin, hist->rmax, hist->dr);
  for (i=0;i<hist->ntype*hist->ntype*hist->ntype;i++)
    for (ir1=0;ir1<hist->nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
	idx=MapIndex(ir1, ir2);
	//r1=rmin+dr*ir1; r2=rmin+dr*ir2;
	//fprintf(fid, "%lf %lf %lf %lf %lf\n", r1, r2, tmins[idx], tmaxs[idx], dt[idx]);
	for (it=0;it<hist->nt[idx];it++)
	  fprintf(fid, "%ld ", hist->bins[i][idx][it]);
	fprintf(fid, "\n");
      }
  fclose(fid);
  return 0;
}

int Xdat2Hist(char *fxyz, xdat_p_t xdat, hist_p_t hist){
  int i,j,k,ntype,val;
  FILE *fid = fopen(fxyz, "r");
  double a1,a2,a3;
  double b1,b2,b3;
  double c1,c2,c3;
  vec_t *slist, vec1, vec2, vec3, pos1, pos2, pos3;
  int idx, it, ir1, ir2, ir12, ir23, ir31, *typemap=xdat->typemap;
  double r1, r2, r3, t12, t23, t31;
  
  fscanf(fid, "%d", &ntype);
  for (i=0;i<ntype;i++){
    fscanf(fid, "%d", &val);
  }
  fscanf(fid, "%lf %lf %lf", &a1, &a2, &a3);
  fscanf(fid, "%lf %lf %lf", &b1, &b2, &b3);
  fscanf(fid, "%lf %lf %lf", &c1, &c2, &c3);
  slist = (vec_t *) malloc(sizeof(vec_t)*xdat->natom);  
  for (i=0;i<xdat->natom;i++){
    fscanf(fid, "%lf %lf %lf", &(slist[i][0]), &(slist[i][1]), &(slist[i][2]));
  }
  fclose(fid);
  
  for (i=0;i<xdat->natom;i++)
    for (j=0;j<i;j++)
      for (k=0;k<j;k++){
	/*
	  triplets (i,j,k). vec1: i->j; vec2: j->k; vec3: k->i;
	  r1, r2, r3: distances of three pairs in Cartesian Coordinates.
	  t1, t2, t3: cos(theta) at vertex i, j and k.
	  (Be careful of the "-" in cos(theta) if you calculate from vecs.)
	*/
	vec1[0] = slist[j][0]-slist[i][0];
	vec1[1] = slist[j][1]-slist[i][1];
	vec1[2] = slist[j][2]-slist[i][2];
	if (vec1[0] > 0.5) vec1[0]--; if (vec1[0]< -0.5) vec1[0]++;
	if (vec1[1] > 0.5) vec1[1]--; if (vec1[1]< -0.5) vec1[1]++;
	if (vec1[2] > 0.5) vec1[2]--; if (vec1[2]< -0.5) vec1[2]++;
	pos1[0] = a1*vec1[0]+b1*vec1[1]+ c1*vec1[2];
	pos1[1] = a2*vec1[0]+b2*vec1[1]+ c2*vec1[2];
	pos1[2] = a3*vec1[0]+b3*vec1[1]+ c3*vec1[2];
	r1 = sqrt(pos1[0]*pos1[0]+pos1[1]*pos1[1]+pos1[2]*pos1[2]);
	if ((r1>hist->rmax+hist->dr)||(r1<hist->rmin-hist->dr)) continue;
	
	vec2[0] = slist[k][0]-slist[j][0];
	vec2[1] = slist[k][1]-slist[j][1];
	vec2[2] = slist[k][2]-slist[j][2];
	if (vec2[0] > 0.5) vec2[0]--; if (vec2[0]< -0.5) vec2[0]++;
	if (vec2[1] > 0.5) vec2[1]--; if (vec2[1]< -0.5) vec2[1]++;
	if (vec2[2] > 0.5) vec2[2]--; if (vec2[2]< -0.5) vec2[2]++;
	
	pos2[0] = a1*vec2[0]+b1*vec2[1]+ c1*vec2[2];
	pos2[1] = a2*vec2[0]+b2*vec2[1]+ c2*vec2[2];
	pos2[2] = a3*vec2[0]+b3*vec2[1]+ c3*vec2[2];
	r2 = sqrt(pos2[0]*pos2[0]+pos2[1]*pos2[1]+pos2[2]*pos2[2]);
	if ((r2>hist->rmax+hist->dr)||(r2<hist->rmin-hist->dr)) continue;
	
	vec3[0] = slist[i][0]-slist[k][0];
	vec3[1] = slist[i][1]-slist[k][1];
	vec3[2] = slist[i][2]-slist[k][2];
	if (vec3[0] > 0.5) vec3[0]--; if (vec3[0]< -0.5) vec3[0]++;
	if (vec3[1] > 0.5) vec3[1]--; if (vec3[1]< -0.5) vec3[1]++;
	if (vec3[2] > 0.5) vec3[2]--; if (vec3[2]< -0.5) vec3[2]++;
	pos3[0] = a1*vec3[0]+b1*vec3[1]+ c1*vec3[2];
	pos3[1] = a2*vec3[0]+b2*vec3[1]+ c2*vec3[2];
	pos3[2] = a3*vec3[0]+b3*vec3[1]+ c3*vec3[2];
	r3 = sqrt(pos3[0]*pos3[0]+pos3[1]*pos3[1]+pos3[2]*pos3[2]);
	if ((r3>hist->rmax+hist->dr)||(r3<hist->rmin-hist->dr)) continue;
	
	//Be careful of the "-" in cos(theta) if you calculate from vecs.
	t12 = -(pos1[0]*pos2[0]+pos1[1]*pos2[1]+pos1[2]*pos2[2])/r1/r2;
	t23 = -(pos2[0]*pos3[0]+pos2[1]*pos3[1]+pos2[2]*pos3[2])/r2/r3;
	t31 = -(pos3[0]*pos1[0]+pos3[1]*pos1[1]+pos3[2]*pos1[2])/r3/r1;
	
	ir12 = floor((r1-hist->rmin)/hist->dr);
	ir23 = floor((r2-hist->rmin)/hist->dr);
	ir31 = floor((r3-hist->rmin)/hist->dr);
	
	ir1=ir12>ir23?ir12:ir23;
	ir2=ir12<ir23?ir12:ir23;
	idx=MapIndex(ir1, ir2);
	if (idx < hist->nbin){
	  it = floor((t12-hist->tmins[idx])/hist->dt[idx]);
	  //it = (int)((t12-hist->tmins[idx])/hist->dt[idx]);
	  if ((it < hist->nt[idx]) && (it >=0))
	    hist->bins[typemap[i]*ntype*ntype+typemap[j]*ntype+typemap[k]][idx][it]++;
	}
	
	ir1=ir31>ir23?ir31:ir23;
	ir2=ir31<ir23?ir31:ir23;
	idx=MapIndex(ir1, ir2);
	if (idx < hist->nbin){ 
	  it = floor((t23-hist->tmins[idx])/hist->dt[idx]);
	  //it = (int)((t23-hist->tmins[idx])/hist->dt[idx]);
	  if ((it < hist->nt[idx]) && (it >=0))
	    hist->bins[typemap[i]*ntype*ntype+typemap[j]*ntype+typemap[k]][idx][it]++;
	}
	
	ir1=ir12>ir31?ir12:ir31;
	ir2=ir12<ir31?ir12:ir31;
	idx=MapIndex(ir1, ir2);
	if (idx < hist->nbin){
	  it = floor((t31-hist->tmins[idx])/hist->dt[idx]);
	  //it = (int)((t31-hist->tmins[idx])/hist->dt[idx]);	  
	  if ((it < hist->nt[idx]) && (it >=0))
	    hist->bins[typemap[i]*ntype*ntype+typemap[j]*ntype+typemap[k]][idx][it]++;
	}
      }
  hist->nsample++;
  return 0; 
}

int Xdat2Hist_NeighList(char *fxyz, xdat_p_t xdat, hist_p_t hist){
  int i,j,k,ntype,val;
  FILE *fid = fopen(fxyz, "r");
  double a1,a2,a3;
  double b1,b2,b3;
  double c1,c2,c3;
  vec_t *slist, vec1, vec2, vec3, pos1, pos2, pos3;
  int idx, it, ir1, ir2, ir12, ir23, ir31, *typemap=xdat->typemap;
  double r1, r2, r3, t12, t23, t31;
  
  fscanf(fid, "%d", &ntype);
  for (i=0;i<ntype;i++){
    fscanf(fid, "%d", &val);
  }
  fscanf(fid, "%lf %lf %lf", &a1, &a2, &a3);
  fscanf(fid, "%lf %lf %lf", &b1, &b2, &b3);
  fscanf(fid, "%lf %lf %lf", &c1, &c2, &c3);
  slist = (vec_t *) malloc(sizeof(vec_t)*xdat->natom);  
  for (i=0;i<xdat->natom;i++){
    fscanf(fid, "%lf %lf %lf", &(slist[i][0]), &(slist[i][1]), &(slist[i][2]));
  }
  fclose(fid);
  
  /* Create Neighbor List */
  neighl_p_t *neighbors, neigh, neighi, neighj;
  neighbors=(neighl_p_t*) malloc(sizeof(neighl_p_t)*xdat->natom);
  for (i=0;i<xdat->natom;i++){
    neighbors[i]=(neighl_p_t) malloc(sizeof(neighl_t));
    neighbors[i]->next=NULL;
    neighbors[i]->atom_id=i;
    for (j=0;j<i;j++){
	vec1[0] = slist[j][0]-slist[i][0];
	vec1[1] = slist[j][1]-slist[i][1];
	vec1[2] = slist[j][2]-slist[i][2];
	if (vec1[0] > 0.5) vec1[0]--; if (vec1[0]< -0.5) vec1[0]++;
	if (vec1[1] > 0.5) vec1[1]--; if (vec1[1]< -0.5) vec1[1]++;
	if (vec1[2] > 0.5) vec1[2]--; if (vec1[2]< -0.5) vec1[2]++;
	pos1[0] = a1*vec1[0]+b1*vec1[1]+ c1*vec1[2];
	pos1[1] = a2*vec1[0]+b2*vec1[1]+ c2*vec1[2];
	pos1[2] = a3*vec1[0]+b3*vec1[1]+ c3*vec1[2];
	r1 = sqrt(pos1[0]*pos1[0]+pos1[1]*pos1[1]+pos1[2]*pos1[2]);
	if ((r1>hist->rmax+hist->dr)||(r1<hist->rmin-hist->dr)) continue;
	neigh=(neighl_p_t) malloc(sizeof(neighl_t));
	neigh->atom_id = j;
	neigh->next=neighbors[i]->next;
	neighbors[i]->next=neigh;
    }
  }
  
  for (i=0;i<xdat->natom;i++){
    neighi=neighbors[i];
    while (neighi->next != NULL){
      neighi=neighi->next;
      j=neighi->atom_id;
      neighj = neighbors[j];
      while (neighj->next != NULL){
	neighj=neighj->next;
	k=neighj->atom_id;
	/*
	  triplets (i,j,k). vec1: i->j; vec2: j->k; vec3: k->i;
	  r1, r2, r3: distances of three pairs in Cartesian Coordinates.
	  t1, t2, t3: cos(theta) at vertex i, j and k.
	  (Be careful of the "-" in cos(theta) if you calculate from vecs.)
	*/
	vec1[0] = slist[j][0]-slist[i][0];
	vec1[1] = slist[j][1]-slist[i][1];
	vec1[2] = slist[j][2]-slist[i][2];
	if (vec1[0] > 0.5) vec1[0]--; if (vec1[0]< -0.5) vec1[0]++;
	if (vec1[1] > 0.5) vec1[1]--; if (vec1[1]< -0.5) vec1[1]++;
	if (vec1[2] > 0.5) vec1[2]--; if (vec1[2]< -0.5) vec1[2]++;
	pos1[0] = a1*vec1[0]+b1*vec1[1]+ c1*vec1[2];
	pos1[1] = a2*vec1[0]+b2*vec1[1]+ c2*vec1[2];
	pos1[2] = a3*vec1[0]+b3*vec1[1]+ c3*vec1[2];
	r1 = sqrt(pos1[0]*pos1[0]+pos1[1]*pos1[1]+pos1[2]*pos1[2]);
	if ((r1>hist->rmax+hist->dr)||(r1<hist->rmin-hist->dr)) continue;
	
	vec2[0] = slist[k][0]-slist[j][0];
	vec2[1] = slist[k][1]-slist[j][1];
	vec2[2] = slist[k][2]-slist[j][2];
	if (vec2[0] > 0.5) vec2[0]--; if (vec2[0]< -0.5) vec2[0]++;
	if (vec2[1] > 0.5) vec2[1]--; if (vec2[1]< -0.5) vec2[1]++;
	if (vec2[2] > 0.5) vec2[2]--; if (vec2[2]< -0.5) vec2[2]++;
	
	pos2[0] = a1*vec2[0]+b1*vec2[1]+ c1*vec2[2];
	pos2[1] = a2*vec2[0]+b2*vec2[1]+ c2*vec2[2];
	pos2[2] = a3*vec2[0]+b3*vec2[1]+ c3*vec2[2];
	r2 = sqrt(pos2[0]*pos2[0]+pos2[1]*pos2[1]+pos2[2]*pos2[2]);
	if ((r2>hist->rmax+hist->dr)||(r2<hist->rmin-hist->dr)) continue;
	
	vec3[0] = slist[i][0]-slist[k][0];
	vec3[1] = slist[i][1]-slist[k][1];
	vec3[2] = slist[i][2]-slist[k][2];
	if (vec3[0] > 0.5) vec3[0]--; if (vec3[0]< -0.5) vec3[0]++;
	if (vec3[1] > 0.5) vec3[1]--; if (vec3[1]< -0.5) vec3[1]++;
	if (vec3[2] > 0.5) vec3[2]--; if (vec3[2]< -0.5) vec3[2]++;
	pos3[0] = a1*vec3[0]+b1*vec3[1]+ c1*vec3[2];
	pos3[1] = a2*vec3[0]+b2*vec3[1]+ c2*vec3[2];
	pos3[2] = a3*vec3[0]+b3*vec3[1]+ c3*vec3[2];
	r3 = sqrt(pos3[0]*pos3[0]+pos3[1]*pos3[1]+pos3[2]*pos3[2]);
	if ((r3>hist->rmax+hist->dr)||(r3<hist->rmin-hist->dr)) continue;
	
	//Be careful of the "-" in cos(theta) if you calculate from vecs.
	t12 = -(pos1[0]*pos2[0]+pos1[1]*pos2[1]+pos1[2]*pos2[2])/r1/r2;
	t23 = -(pos2[0]*pos3[0]+pos2[1]*pos3[1]+pos2[2]*pos3[2])/r2/r3;
	t31 = -(pos3[0]*pos1[0]+pos3[1]*pos1[1]+pos3[2]*pos1[2])/r3/r1;
	
	ir12 = floor((r1-hist->rmin)/hist->dr);
	ir23 = floor((r2-hist->rmin)/hist->dr);
	ir31 = floor((r3-hist->rmin)/hist->dr);
	
	ir1=ir12>ir23?ir12:ir23;
	ir2=ir12<ir23?ir12:ir23;
	idx=MapIndex(ir1, ir2);
	if (idx < hist->nbin){
	  it = floor((t12-hist->tmins[idx])/hist->dt[idx]);
	  //it = (int)((t12-hist->tmins[idx])/hist->dt[idx]);
	  if ((it < hist->nt[idx]) && (it >=0))
	    hist->bins[typemap[i]*ntype*ntype+typemap[j]*ntype+typemap[k]][idx][it]++;
	}
	
	ir1=ir31>ir23?ir31:ir23;
	ir2=ir31<ir23?ir31:ir23;
	idx=MapIndex(ir1, ir2);
	if (idx < hist->nbin){ 
	  it = floor((t23-hist->tmins[idx])/hist->dt[idx]);
	  //it = (int)((t23-hist->tmins[idx])/hist->dt[idx]);
	  if ((it < hist->nt[idx]) && (it >=0))
	    hist->bins[typemap[i]*ntype*ntype+typemap[j]*ntype+typemap[k]][idx][it]++;
	}
	
	ir1=ir12>ir31?ir12:ir31;
	ir2=ir12<ir31?ir12:ir31;
	idx=MapIndex(ir1, ir2);
	if (idx < hist->nbin){
	  it = floor((t31-hist->tmins[idx])/hist->dt[idx]);
	  //it = (int)((t31-hist->tmins[idx])/hist->dt[idx]);	  
	  if ((it < hist->nt[idx]) && (it >=0))
	    hist->bins[typemap[i]*ntype*ntype+typemap[j]*ntype+typemap[k]][idx][it]++;
	}
      }
    }
  }
  hist->nsample++;

  for (i=0;i<xdat->natom;i++){
    neigh=neighbors[i];
    while (neigh->next != NULL){
      neighi=neigh->next;
      free(neigh);
      neigh=neighi;
    }
  }
  return 0; 
}
