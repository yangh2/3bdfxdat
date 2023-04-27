#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hist.h"
#include "3bdf.h"
#include "integrate.h"

int InitPDF(char *fname, pdf_p_t *pdf_p){
  int i, j;
  *pdf_p = (pdf_p_t) malloc(sizeof(pdf_t));
  pdf_p_t pdf=*pdf_p;

  FILE *fid = fopen(fname, "r");
  fscanf(fid, "%d", &(pdf->ntype));
  fscanf(fid, "%d %lf %lf", &(pdf->nbin), &(pdf->rmin), &(pdf->rmax));
  pdf->dr = (pdf->rmax-pdf->rmin)/(pdf->nbin-1);
  pdf->bins = (double **) malloc(sizeof(double*)*pdf->ntype*pdf->ntype);
  for (i=0;i<pdf->ntype*pdf->ntype;i++)
    pdf->bins[i] = (double *) malloc(sizeof(double)*pdf->nbin);
  double val;
  for (int i=0;i<pdf->nbin;i++){
    fscanf(fid, "%lf", &val);
      for (j=0;j<pdf->ntype*pdf->ntype*pdf->ntype;j++)
	fscanf(fid, "%lf", &(pdf->bins[j][i]));
  }
  fclose(fid);
  return 0;
}

int PrintPDF(char *fname, pdf_p_t pdf){
  FILE *fid = fopen(fname, "w");
  fprintf(fid, "%d\n", pdf->ntype);
  fprintf(fid, "%d %lf %lf", pdf->nbin, pdf->rmin, pdf->rmax);
  for (int i=0;i<pdf->nbin;i++){
    fprintf(fid, "%lf", i*pdf->dr+pdf->rmin);
    for (int j=0;j<pdf->ntype*pdf->ntype*pdf->ntype;j++)
      fprintf(fid, " %lf", pdf->bins[j][i]);
    fprintf(fid, "\n");
  }
  fclose(fid);
  return 0;
}

int FreePDF(pdf_p_t pdf){
  int i;
  for (i=0;i<pdf->ntype*pdf->ntype;i++)
    free(pdf->bins[i]);
  free(pdf->bins);
  free(pdf);
  return 0;
}

int Trilinear(double *xs, double *gs, double *as){
  double x0, y0, z0, x1, y1, z1;
  double c000, c100, c010, c110, c001, c101, c011, c111;
  double dx,dy,dz, dxdydz;
  
  x0 = xs[0]; y0 = xs[1]; z0 = xs[2];
  x1 = xs[3]; y1 = xs[4]; z1 = xs[5];
  dx = x0-x1; dy = y0-y1; dz = z0-z1;
  dxdydz=dx*dy*dz;
  
  c000 = gs[0]; c100 = gs[1]; c010 = gs[2]; c110 = gs[3];
  c001 = gs[4]; c101 = gs[5]; c011 = gs[6]; c111 = gs[7];

  as[0] = (-c000*x1*y1*z1+c001*x1*y1*z0+c010*x1*y0*z1-c011*x1*y0*z0)/dxdydz+(c100*x0*y1*z1-c101*x0*y1*z0-c110*x0*y0*z1+c111*x0*y0*z0)/dxdydz;
  as[1] = (c000*y1*z1-c001*y1*z0-c010*y0*z1+c011*y0*z0)/dxdydz+(-c100*y1*z1+c101*y1*z0+c110*y0*z1-c111*y0*z0)/dxdydz;
  as[2] = (c000*x1*z1-c001*x1*z0-c010*x1*z1+c011*x1*z0)/dxdydz+(-c100*x0*z1+c101*x0*z0+c110*x0*z1-c111*x0*z0)/dxdydz;
  as[3] = (c000*x1*y1-c001*x1*y1-c010*x1*y0+c011*x1*y0)/dxdydz+(-c100*x0*y1+c101*x0*y1+c110*x0*y0-c111*x0*y0)/dxdydz;
  as[4] = (-c000*z1+c001*z0+c010*z1-c011*z0+c100*z1-c101*z0-c110*z1+c111*z0)/dxdydz;
  as[5] = (-c000*y1+c001*y1+c010*y0-c011*y0+c100*y1-c101*y1-c110*y0+c111*y0)/dxdydz;
  as[6] = (-c000*x1+c001*x1+c010*x1-c011*x1+c100*x0-c101*x0-c110*x0+c111*x0)/dxdydz;
  as[7] = (c000-c001-c010+c011-c100+c101+c110-c111)/dxdydz;

  return 0;
}

double LinInterp(pdf_p_t pdf, int type_id, double r){
  int l;
  if (r<pdf->dr/2) return 0;
  l = floor((r-pdf->dr/2)/pdf->dr);
  double val;
  val = pdf->bins[type_id][l]+(pdf->bins[type_id][l+1]-pdf->bins[type_id][l])/pdf->dr*(r-pdf->dr/2-pdf->dr*l);
  return val;
}

double TriLinearIntegral(pdf_p_t pdf, tbf_p_t tbf,
			 int t1, int t2, int t3, /* type */
			 int ir1, int ir2, int it){
  double r1, r2, r3, cost, val, g2_1, g2_2, g2_3;
  int idx=MapIndex(ir1, ir2);
  r1=tbf->rmin+tbf->dr*ir1; r2=tbf->rmin+tbf->dr*ir2;
  cost=tbf->dt[idx]*it+tbf->tmins[idx];
  double xs[6], gs[8], as[8];
  xs[0] = r1;
  xs[1] = r2;
  xs[2] = cost;
  xs[3] = r1+tbf->dr;
  xs[4] = r2+tbf->dr;
  xs[5] = cost+tbf->dt[idx];

  int ntype=pdf->ntype;

  r1=tbf->rmin+tbf->dr*ir1;
  r2=tbf->rmin+tbf->dr*ir2;
  cost=tbf->dt[idx]*it+tbf->tmins[idx];
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
  g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
  g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
  gs[0] = g2_1*g2_2*g2_3;

  r1=tbf->rmin+tbf->dr*(ir1+1);
  r2=tbf->rmin+tbf->dr*ir2;
  cost=tbf->dt[idx]*it+tbf->tmins[idx];
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
  g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
  g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
  gs[1] = g2_1*g2_2*g2_3;

  r1=tbf->rmin+tbf->dr*ir1;
  r2=tbf->rmin+tbf->dr*(ir2+1);
  cost=tbf->dt[idx]*it+tbf->tmins[idx];
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
  g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
  g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
  gs[2] = g2_1*g2_2*g2_3;

  r1=tbf->rmin+tbf->dr*ir1;
  r2=tbf->rmin+tbf->dr*ir2;
  cost=tbf->dt[idx]*(it+1)+tbf->tmins[idx];
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
  g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
  g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
  gs[3] = g2_1*g2_2*g2_3;

  r1=tbf->rmin+tbf->dr*(ir1+1);
  r2=tbf->rmin+tbf->dr*(ir2+1);
  cost=tbf->dt[idx]*it+tbf->tmins[idx];
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
  g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
  g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
  gs[4] = g2_1*g2_2*g2_3;

  r1=tbf->rmin+tbf->dr*(ir1+1);
  r2=tbf->rmin+tbf->dr*ir2;
  cost=tbf->dt[idx]*(it+1)+tbf->tmins[idx];
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
  g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
  g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
  gs[5] = g2_1*g2_2*g2_3;

  r1=tbf->rmin+tbf->dr*ir1;
  r2=tbf->rmin+tbf->dr*(ir2+1);
  cost=tbf->dt[idx]*(it+1)+tbf->tmins[idx];
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
  g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
  g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
  gs[6] = g2_1*g2_2*g2_3;

  r1=tbf->rmin+tbf->dr*(ir1+1);
  r2=tbf->rmin+tbf->dr*(ir2+1);
  cost=tbf->dt[idx]*(it+1)+tbf->tmins[idx];
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
  g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
  g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
  gs[7] = g2_1*g2_2*g2_3;

  double ix0, iy0, iz0, ix1, iy1, iz1;
  //ix0 = 4*pi/3*(pow(xs[3],3)-pow(xs[0],3));
  ix0 = 4*pi/3*tbf->dr*(xs[3]*xs[3]+xs[3]*xs[0]+xs[0]*xs[0]);
  //iy0 = 2*pi/3*(pow(xs[4],3)-pow(xs[1],3));
  iy0 = 2*pi/3*tbf->dr*(xs[4]*xs[4]+xs[4]*xs[1]+xs[1]*xs[1]);
  iz0 = tbf->dt[idx];
  ix1 = pi*tbf->dr*(xs[3]*xs[3]+xs[0]*xs[0])*(xs[3]+xs[0]);
  iy1 = pi/2*tbf->dr*(xs[4]*xs[4]+xs[1]*xs[1])*(xs[4]+xs[1]);
  iz1 = tbf->dt[idx]*(xs[5]+xs[2])/2;

  Trilinear(xs, gs, as);
  val = (as[0]*ix0*iy0*iz0 +as[1]*ix1*iy0*iz0 +as[2]*ix0*iy1*iz0 +as[3]*ix0*iy0*iz1 +as[4]*ix1*iy1*iz0 +as[5]*ix1*iy0*iz1 +as[6]*ix0*iy1*iz1 +as[7]*ix1*iy1*iz1);

  return val;
}

double TriLinearIntegralN(pdf_p_t pdf, tbf_p_t tbf,
			 int t1, int t2, int t3, /* type */
			  int ir1, int ir2, int it, int dim){
  double r1, r2, r3, cost, val, g2_1, g2_2, g2_3;
  int dim_x, dim_y, dim_z;
  
  int idx=MapIndex(ir1, ir2);
  int ntype=pdf->ntype;
  double xs[6], gs[8], as[8];

  double ix0, iy0, iz0, ix1, iy1, iz1;
  double dr, dt;

  val = 0;
  dr = tbf->dr/dim;
  dt = tbf->dt[idx]/dim;
  for (dim_x=0;dim_x<dim;dim_x++)
    for (dim_y=0;dim_y<dim;dim_y++)
      for (dim_z=0;dim_z<dim;dim_z++){
	r1=tbf->rmin+tbf->dr*ir1+dr*dim_x; r2=tbf->rmin+tbf->dr*ir2+dr*dim_y;
	cost=tbf->dt[idx]*it+tbf->tmins[idx] + dt*dim_z;
	
	xs[0] = r1;
	xs[1] = r2;
	xs[2] = cost;
	xs[3] = r1+dr;
	xs[4] = r2+dr;
	xs[5] = cost+dt;

        r1=xs[0];
	r2=xs[1];
	cost=xs[2];
	r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
	g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
	g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
	gs[0] = g2_1*g2_2*g2_3;

	r1=xs[3];
	r2=xs[1];
	cost=xs[2];
	r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
	g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
	g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
	gs[1] = g2_1*g2_2*g2_3;

	r1=xs[0];
	r2=xs[4];
	cost=xs[2];
	r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
	g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
	g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
	gs[2] = g2_1*g2_2*g2_3;
	
	r1=xs[0];
	r2=xs[1];
	cost=xs[5];
	r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
	g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
	g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
	gs[3] = g2_1*g2_2*g2_3;
	
	r1=xs[3];
	r2=xs[4];
	cost=xs[2];
	r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
	g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
	g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
	gs[4] = g2_1*g2_2*g2_3;
	
	r1=xs[3];
	r2=xs[1];
	cost=xs[5];
	r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
	g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
	g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
	gs[5] = g2_1*g2_2*g2_3;

	r1=xs[0];
	r2=xs[4];
	cost=xs[5];
	r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
	g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
	g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
	gs[6] = g2_1*g2_2*g2_3;
	
	r1=xs[3];
	r2=xs[4];
	cost=xs[5];
	r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
	g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
	g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
	gs[7] = g2_1*g2_2*g2_3;
	
	ix0 = 4*pi/3*dr*(xs[3]*xs[3]+xs[3]*xs[0]+xs[0]*xs[0]);
	iy0 = 2*pi/3*dr*(xs[4]*xs[4]+xs[4]*xs[1]+xs[1]*xs[1]);
	iz0 = dt;
	ix1 = pi*dr*(xs[3]*xs[3]+xs[0]*xs[0])*(xs[3]+xs[0]);
	iy1 = pi/2*dr*(xs[4]*xs[4]+xs[1]*xs[1])*(xs[4]+xs[1]);
	iz1 = dt*(xs[5]+xs[2])/2;

	Trilinear(xs, gs, as);
	val += (as[0]*ix0*iy0*iz0 +as[1]*ix1*iy0*iz0 +as[2]*ix0*iy1*iz0 +as[3]*ix0*iy0*iz1 +as[4]*ix1*iy1*iz0 +as[5]*ix1*iy0*iz1 +as[6]*ix0*iy1*iz1 +as[7]*ix1*iy1*iz1);
}

  return val;
}

int FreeIntegral(integral_p_t dest){
  int i;
  for (i=0;i<dest->ntype*dest->ntype*dest->ntype;i++)
    free(dest->bins[i]);
  free(dest->bins);
  free(dest);
  return 0;
}

int PrintIntegral(char *fname, integral_p_t dest){
  int i,ir;
  double r;
  FILE *fid=fopen(fname, "w");
  for (ir=0;ir<dest->nr;ir++){
    r=dest->rmin+dest->dr*ir;
    fprintf(fid, "%lf ", r);
    for (i=0;i<dest->ntype*dest->ntype*dest->ntype;i++)
      fprintf(fid, "%lf ", dest->bins[i][ir]);
    fprintf(fid, "\n");
  }
  fclose(fid);
  return 0;
}

int IntegrateB(tbf_p_t tbf, pdf_p_t pdf, xdat_p_t xdat, integral_p_t *dest_p){
  *dest_p = (integral_p_t)malloc(sizeof(integral_t));
  integral_p_t dest= *dest_p;

  dest->ntype=tbf->ntype;
  dest->nr = ceil(tbf->rmax/tbf->dr);
  dest->rmin = 0;
  dest->rmax = tbf->rmax;
  dest->dr = tbf->dr;
  
  int i,j;
  dest->bins = (double **) malloc(sizeof(double*)*dest->ntype*dest->ntype*dest->ntype);
  for (i=0;i<dest->ntype*dest->ntype*dest->ntype;i++){
    dest->bins[i] = (double *) malloc(sizeof(double)*dest->nr);
    for (j=0;j<dest->nr;j++)
      dest->bins[i][j] = 0;
  }

  int  nr, ntype;
  nr = dest->nr;
  ntype = dest->ntype;

  double rho=(double)xdat->natom/xdat->vol;
  double val, r1, r2, r3, dr=dest->dr, cost, dt, vol;
  int ir1, ir2, it, nt;;
  int t1, t2, t3;
  int p1, p2, p3, pmax;
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++){

	i=t1*ntype*ntype+t2*ntype+t3;
	for (ir1=0;ir1<nr;ir1++)
	  for (ir2=0;ir2<=ir1;ir2++){
	    r1=dest->rmin+dest->dr*ir1; r2=dest->rmin+dest->dr*ir2;
	    nt=ceil(2*r2/dr);
	    dt=(double) 2.0/nt;
	    /* caluclate volume elements and normalize tbf */

	    for (it=0;it<nt;it++){
	      vol = (double)8.0*pi*pi/9 /* volumn element */
		*((r1+dr)*(r1+dr)*(r1+dr)-r1*r1*r1)
		*((r2+dr)*(r2+dr)*(r2+dr)-r2*r2*r2)
		*dt;
	      cost = -1+it*dt;
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	      p1 = floor((r1-pdf->rmin)/pdf->dr);
	      p2 = floor((r2-pdf->rmin)/pdf->dr);
	      p3 = floor((r3-pdf->rmin)/pdf->dr);
	      
	      val = (pdf->bins[t1*ntype+t2][p1]-1)*(pdf->bins[t1*ntype+t3][p2]-1)*(pdf->bins[t3*ntype+t2][p3]-1);

	      
	      p1 = ir1;
	      p2 = ir2;
	      p3 = floor((r3-dest->rmin)/dest->dr);

	      if ( p3 >= nr ) continue;
	      
	      pmax=p1>p2?p1:p2;
	      pmax=pmax>p3?pmax:p3;
	      
	      dest->bins[t1*ntype*ntype+t2*ntype+t3][pmax]+=val*vol;
	    }
	  }
	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[t1*ntype*ntype+t2*ntype+t3][ir1] += dest->bins[t1*ntype*ntype+t2*ntype+t3][ir1-1];
	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[i][ir1] *= (double)1.0/6*rho*rho;
      }
  return 0;
}

int IntegrateB_ext(pdf_p_t pdf, xdat_p_t xdat, integral_p_t *dest_p){
  *dest_p = (integral_p_t)malloc(sizeof(integral_t));
  integral_p_t dest= *dest_p;

  dest->ntype=pdf->ntype;
  dest->nr = pdf->nbin;
  dest->rmin = 0;
  dest->rmax = pdf->rmax;
  dest->dr = pdf->dr;
  
  int i,j;
  dest->bins = (double **) malloc(sizeof(double*)*dest->ntype*dest->ntype*dest->ntype);
  for (i=0;i<dest->ntype*dest->ntype*dest->ntype;i++){
    dest->bins[i] = (double *) malloc(sizeof(double)*dest->nr);
    for (j=0;j<dest->nr;j++)
      dest->bins[i][j] = 0;
  }

  int  nr, ntype;
  nr = dest->nr;
  ntype = dest->ntype;

  double rho=(double)xdat->natom/xdat->vol;
  double val, r1, r2, r3, dr=dest->dr, cost, dt, vol;
  int ir1, ir2, it, nt;;
  int t1, t2, t3;
  int p1, p2, p3, pmax;
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++){

	i=t1*ntype*ntype+t2*ntype+t3;
	for (ir1=0;ir1<nr;ir1++)
	  for (ir2=0;ir2<=ir1;ir2++){
	    r1=dest->rmin+dest->dr*ir1; r2=dest->rmin+dest->dr*ir2;
	    nt=ceil(2*r2/dr);
	    dt=(double) 2.0/nt;
	    /* caluclate volume elements and normalize tbf */

	    for (it=0;it<nt;it++){
	      vol = (double)8.0*pi*pi/9 /* volumn element */
		*((r1+dr)*(r1+dr)*(r1+dr)-r1*r1*r1)
		*((r2+dr)*(r2+dr)*(r2+dr)-r2*r2*r2)
		*dt;
	      cost = -1+it*dt;
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	      p1 = floor((r1-pdf->rmin)/pdf->dr);
	      p2 = floor((r2-pdf->rmin)/pdf->dr);
	      p3 = floor((r3-pdf->rmin)/pdf->dr);
	      
	      val = (pdf->bins[t1*ntype+t2][p1]-1)*(pdf->bins[t1*ntype+t3][p2]-1)*(pdf->bins[t3*ntype+t2][p3]-1);

	      
	      p1 = ir1;
	      p2 = ir2;
	      p3 = floor((r3-dest->rmin)/dest->dr);

	      if ( p3 >= nr ) continue;
	      
	      pmax=p1>p2?p1:p2;
	      pmax=pmax>p3?pmax:p3;
	      
	      dest->bins[t1*ntype*ntype+t2*ntype+t3][pmax]+=val*vol;
	    }
	  }
	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[t1*ntype*ntype+t2*ntype+t3][ir1] += dest->bins[t1*ntype*ntype+t2*ntype+t3][ir1-1];
	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[i][ir1] *= (double)1.0/6*rho*rho;
      }
  return 0;
}

int R3MinMax(tbf_p_t tbf, int ir1, int ir2, int it, double *r3min, double *r3max){
  int idx;
  double r1, r2, cost, dt, r3;
  double dr = tbf->dr;
  *r3min = 10000;
  *r3max = -1;
  
  idx=MapIndex(ir1, ir2);
  dt = tbf->dt[idx];
  r1=tbf->rmin+dr*ir1;
  r2=tbf->rmin+dr*ir2;
  cost=tbf->tmins[idx]+dt*it;
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  *r3min = *r3min<r3?*r3min:r3;
  *r3max = *r3max>r3?*r3max:r3;

  r1=tbf->rmin+dr*(ir1+1);
  r2=tbf->rmin+dr*ir2;
  cost=tbf->tmins[idx]+dt*it;
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  *r3min = *r3min<r3?*r3min:r3;
  *r3max = *r3max>r3?*r3max:r3;

  r1=tbf->rmin+dr*ir1;
  r2=tbf->rmin+dr*(ir2+1);
  cost=tbf->tmins[idx]+dt*it;
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  *r3min = *r3min<r3?*r3min:r3;
  *r3max = *r3max>r3?*r3max:r3;

  r1=tbf->rmin+dr*ir1;
  r2=tbf->rmin+dr*ir2;
  cost=tbf->tmins[idx]+dt*(it+1);
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  *r3min = *r3min<r3?*r3min:r3;
  *r3max = *r3max>r3?*r3max:r3;

  r1=tbf->rmin+dr*(ir1+1);
  r2=tbf->rmin+dr*(ir2+1);
  cost=tbf->tmins[idx]+dt*it;
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  *r3min = *r3min<r3?*r3min:r3;
  *r3max = *r3max>r3?*r3max:r3;

  r1=tbf->rmin+dr*(ir1+1);
  r2=tbf->rmin+dr*ir2;
  cost=tbf->tmins[idx]+dt*(it+1);
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  *r3min = *r3min<r3?*r3min:r3;
  *r3max = *r3max>r3?*r3max:r3;

  r1=tbf->rmin+dr*ir1;
  r2=tbf->rmin+dr*(ir2+1);
  cost=tbf->tmins[idx]+dt*(it+1);
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  *r3min = *r3min<r3?*r3min:r3;
  *r3max = *r3max>r3?*r3max:r3;

  r1=tbf->rmin+dr*(ir1+1);
  r2=tbf->rmin+dr*(ir2+1);
  cost=tbf->tmins[idx]+dt*(it+1);
  r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
  *r3min = *r3min<r3?*r3min:r3;
  *r3max = *r3max>r3?*r3max:r3;

  return 0;
}

int InitEplison(tbf_p_t tbfin, pdf_p_t pdf, tbf_p_t *eps_p){
  *eps_p = (tbf_p_t) malloc(sizeof(tbf_t));
  tbf_p_t eps=*eps_p;

  /* copy parameters from hist to eps */
  eps->nsample = tbfin->nsample;
  eps->ntype = tbfin->ntype;
  eps->rmin = tbfin->rmin;
  eps->rmax = tbfin->rmax;
  eps->dr = tbfin->dr;
  eps->nr = tbfin->nr;
  eps->nbin = tbfin->nbin;

  /* allocate memory for 3bdf */
  int nbin=eps->nbin,i;
  eps->tmins = (double *) malloc(sizeof(double)*nbin);
  eps->tmaxs = (double *) malloc(sizeof(double)*nbin);
  eps->dt = (double *) malloc(sizeof(double)*nbin);
  eps->nt = (long *) malloc(sizeof(long)*nbin);

  for (i=0;i<nbin;i++){
    eps->tmins[i] = tbfin->tmins[i];
    eps->tmaxs[i] = tbfin->tmaxs[i];
    eps->dt[i] = tbfin->dt[i];
    eps->nt[i] = tbfin->nt[i];
  }

  int ntype=eps->ntype, ir1, ir2, idx,nr=eps->nr, it;
  double ***f3s, *vols, r1, r2, r3, dr, cost, val, r3min, r3max;
  double vleft, vright, vlow, vhigh, val_low, val_high;
  int pleft, pright;
  int t1, t2, t3;
  dr = eps->dr;
  f3s = (double ***) malloc(sizeof(double **)*ntype*ntype*ntype);
  vols = (double *) malloc(sizeof(double)*nbin);
  
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++){
	i=t1*ntype*ntype+t2*ntype+t3;
	f3s[i] = (double **) malloc(sizeof(double*)*nbin);
	for (ir1=0;ir1<nr;ir1++)
	  for (ir2=0;ir2<=ir1;ir2++){
	    idx=MapIndex(ir1, ir2);
	    r1=eps->rmin+dr*ir1; r2=eps->rmin+dr*ir2;
	    
	    f3s[i][idx] = (double *) malloc(sizeof(double)*eps->nt[idx]);
	    vols[idx] = (double)8.0*pi*pi/9
	      *((r1+dr)*(r1+dr)*(r1+dr)-r1*r1*r1)
	      *((r2+dr)*(r2+dr)*(r2+dr)-r2*r2*r2)
	      *eps->dt[idx];
	    /* caluclate volume elements and normalize eps */
	    for (it=0;it<eps->nt[idx];it++){
	      vlow = 1; vhigh = 1;
	      
	      pleft = floor((r1-pdf->rmin)/pdf->dr);
	      vleft = pdf->bins[t1*ntype+t2][pleft];
	      pright = floor((r1+tbfin->dr-pdf->rmin)/pdf->dr);
	      vright = pdf->bins[t1*ntype+t2][pright];
	      vlow *= vleft<vright?vleft:vright;
	      vhigh *= vleft>vright?vleft:vright;
	      
	      pleft = floor((r2-pdf->rmin)/pdf->dr);
	      vleft = pdf->bins[t1*ntype+t3][pleft];
	      pright = floor((r2+tbfin->dr-pdf->rmin)/pdf->dr);
	      vright = pdf->bins[t1*ntype+t3][pright];
	      vlow *= vleft<vright?vleft:vright;
	      vhigh *= vleft>vright?vleft:vright;

	      R3MinMax(tbfin, ir1, ir2, it, &r3min, &r3max);
	      pleft = floor((r3min-pdf->rmin)/pdf->dr);
	      vleft = pdf->bins[t1*ntype+t3][pleft];
	      pright = floor((r3max-pdf->rmin)/pdf->dr);
	      vright = pdf->bins[t1*ntype+t3][pright];
	      vlow *= vleft<vright?vleft:vright;
	      vhigh *= vleft>vright?vleft:vright;

	      if ( tbfin->f3s[i][idx][it] < vlow )
		f3s[i][idx][it] = tbfin->f3s[i][idx][it] - vlow;
	      else {
		if ( tbfin->f3s[i][idx][it] > vhigh )
		  f3s[i][idx][it] = tbfin->f3s[i][idx][it] - vhigh;
		else
		  f3s[i][idx][it] = 0;
	      }
	}
      }
  }
  eps->vols=vols;
  eps->f3s=f3s;
  return 0;
}

int InitSuperpos(tbf_p_t tbfin, pdf_p_t pdf, tbf_p_t *tbf_p){
  *tbf_p = (tbf_p_t) malloc(sizeof(tbf_t));
  tbf_p_t tbf=*tbf_p;

  /* copy parameters from hist to tbf */
  tbf->nsample = tbfin->nsample;
  tbf->ntype = tbfin->ntype;
  tbf->rmin = tbfin->rmin;
  tbf->rmax = tbfin->rmax;
  tbf->dr = tbfin->dr;
  tbf->nr = tbfin->nr;
  tbf->nbin = tbfin->nbin;

  /* allocate memory for 3bdf */
  int nbin=tbf->nbin,i;
  tbf->tmins = (double *) malloc(sizeof(double)*nbin);
  tbf->tmaxs = (double *) malloc(sizeof(double)*nbin);
  tbf->dt = (double *) malloc(sizeof(double)*nbin);
  tbf->nt = (long *) malloc(sizeof(long)*nbin);

  for (i=0;i<nbin;i++){
    tbf->tmins[i] = tbfin->tmins[i];
    tbf->tmaxs[i] = tbfin->tmaxs[i];
    tbf->dt[i] = tbfin->dt[i];
    tbf->nt[i] = tbfin->nt[i];
  }

  int ntype=tbf->ntype, ir1, ir2, idx,nr=tbf->nr, it;
  double ***f3s, *vols, r1, r2, r3, dr, cost, val;
  int p1, p2, p3;
  int t1, t2, t3;
  dr = tbf->dr;
  f3s = (double ***) malloc(sizeof(double **)*ntype*ntype*ntype);
  vols = (double *) malloc(sizeof(double)*nbin);
  
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++){
	i=t1*ntype*ntype+t2*ntype+t3;
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
	    for (it=0;it<tbf->nt[idx];it++){
	      cost = tbf->tmins[idx]+it*tbf->dt[idx];
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	      p1 = floor((r1-pdf->rmin)/pdf->dr);
	      p2 = floor((r2-pdf->rmin)/pdf->dr);
	      p3 = floor((r3-pdf->rmin)/pdf->dr);

	      val = pdf->bins[t1*ntype+t2][p1]*pdf->bins[t1*ntype+t3][p2]*pdf->bins[t3*ntype+t2][p3];
	      f3s[i][idx][it] = val;
	}
      }
  }
  tbf->vols=vols;
  tbf->f3s=f3s;
  return 0;
}

int InitSuperpos_Intgrl(tbf_p_t tbfin, pdf_p_t pdf, tbf_p_t *tbf_p){
  *tbf_p = (tbf_p_t) malloc(sizeof(tbf_t));
  tbf_p_t tbf=*tbf_p;

  /* copy parameters from hist to tbf */
  tbf->nsample = tbfin->nsample;
  tbf->ntype = tbfin->ntype;
  tbf->rmin = tbfin->rmin;
  tbf->rmax = tbfin->rmax;
  tbf->dr = tbfin->dr;
  tbf->nr = tbfin->nr;
  tbf->nbin = tbfin->nbin;

  /* allocate memory for 3bdf */
  int nbin=tbf->nbin,i;
  tbf->tmins = (double *) malloc(sizeof(double)*nbin);
  tbf->tmaxs = (double *) malloc(sizeof(double)*nbin);
  tbf->dt = (double *) malloc(sizeof(double)*nbin);
  tbf->nt = (long *) malloc(sizeof(long)*nbin);

  for (i=0;i<nbin;i++){
    tbf->tmins[i] = tbfin->tmins[i];
    tbf->tmaxs[i] = tbfin->tmaxs[i];
    tbf->dt[i] = tbfin->dt[i];
    tbf->nt[i] = tbfin->nt[i];
  }

  int ntype=tbf->ntype, ir1, ir2, idx,nr=tbf->nr, it;
  double ***f3s, *vols, r1, r2, r3, dr, cost, val;
  int p1, p2, p3;
  int t1, t2, t3;
  dr = tbf->dr;
  f3s = (double ***) malloc(sizeof(double **)*ntype*ntype*ntype);
  vols = (double *) malloc(sizeof(double)*nbin);
  
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++){
	i=t1*ntype*ntype+t2*ntype+t3;
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
	    for (it=0;it<tbf->nt[idx];it++){
	      cost = tbf->tmins[idx]+it*tbf->dt[idx];
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	      p1 = floor((r1-pdf->rmin)/pdf->dr);
	      p2 = floor((r2-pdf->rmin)/pdf->dr);
	      p3 = floor((r3-pdf->rmin)/pdf->dr);

	      val = TriLinearIntegral(pdf, tbf, t1, t2, t3, ir1, ir2, it)/vols[idx];
	      f3s[i][idx][it] = val;
	}
      }
  }
  tbf->vols=vols;
  tbf->f3s=f3s;
  return 0;
}

int InitSuperpos_IntgrlN(tbf_p_t tbfin, pdf_p_t pdf, tbf_p_t *tbf_p, int dim){
  *tbf_p = (tbf_p_t) malloc(sizeof(tbf_t));
  tbf_p_t tbf=*tbf_p;

  /* copy parameters from hist to tbf */
  tbf->nsample = tbfin->nsample;
  tbf->ntype = tbfin->ntype;
  tbf->rmin = tbfin->rmin;
  tbf->rmax = tbfin->rmax;
  tbf->dr = tbfin->dr;
  tbf->nr = tbfin->nr;
  tbf->nbin = tbfin->nbin;

  /* allocate memory for 3bdf */
  int nbin=tbf->nbin,i;
  tbf->tmins = (double *) malloc(sizeof(double)*nbin);
  tbf->tmaxs = (double *) malloc(sizeof(double)*nbin);
  tbf->dt = (double *) malloc(sizeof(double)*nbin);
  tbf->nt = (long *) malloc(sizeof(long)*nbin);

  for (i=0;i<nbin;i++){
    tbf->tmins[i] = tbfin->tmins[i];
    tbf->tmaxs[i] = tbfin->tmaxs[i];
    tbf->dt[i] = tbfin->dt[i];
    tbf->nt[i] = tbfin->nt[i];
  }

  int ntype=tbf->ntype, ir1, ir2, idx,nr=tbf->nr, it;
  double ***f3s, *vols, r1, r2, r3, dr, cost, val;
  int p1, p2, p3;
  int t1, t2, t3;
  dr = tbf->dr;
  f3s = (double ***) malloc(sizeof(double **)*ntype*ntype*ntype);
  vols = (double *) malloc(sizeof(double)*nbin);
  
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++){
	i=t1*ntype*ntype+t2*ntype+t3;
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
	    for (it=0;it<tbf->nt[idx];it++){
	      cost = tbf->tmins[idx]+it*tbf->dt[idx];
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	      p1 = floor((r1-pdf->rmin)/pdf->dr);
	      p2 = floor((r2-pdf->rmin)/pdf->dr);
	      p3 = floor((r3-pdf->rmin)/pdf->dr);

	      val = TriLinearIntegralN(pdf, tbf, t1, t2, t3, ir1, ir2, it, dim)/vols[idx];
	      f3s[i][idx][it] = val;
	}
      }
  }
  tbf->vols=vols;
  tbf->f3s=f3s;
  return 0;
}

int InitSuperposScale(tbf_p_t tbfin, pdf_p_t pdf, double scale, tbf_p_t *tbf_p){
  *tbf_p = (tbf_p_t) malloc(sizeof(tbf_t));
  tbf_p_t tbf=*tbf_p;

  /* copy parameters from hist to tbf */
  tbf->nsample = tbfin->nsample;
  tbf->ntype = tbfin->ntype;
  tbf->rmin = tbfin->rmin;
  tbf->rmax = tbfin->rmax;
  tbf->dr = tbfin->dr;
  tbf->nr = tbfin->nr;
  tbf->nbin = tbfin->nbin;

  /* allocate memory for 3bdf */
  int nbin=tbf->nbin,i;
  tbf->tmins = (double *) malloc(sizeof(double)*nbin);
  tbf->tmaxs = (double *) malloc(sizeof(double)*nbin);
  tbf->dt = (double *) malloc(sizeof(double)*nbin);
  tbf->nt = (long *) malloc(sizeof(long)*nbin);

  for (i=0;i<nbin;i++){
    tbf->tmins[i] = tbfin->tmins[i];
    tbf->tmaxs[i] = tbfin->tmaxs[i];
    tbf->dt[i] = tbfin->dt[i];
    tbf->nt[i] = tbfin->nt[i];
  }

  int ntype=tbf->ntype, ir1, ir2, idx,nr=tbf->nr, it;
  double ***f3s, *vols, r1, r2, r3, dr, cost;
  double g2_1, g2_2, g2_3;
  int t1, t2, t3;
  dr = tbf->dr;
  f3s = (double ***) malloc(sizeof(double **)*ntype*ntype*ntype);
  vols = (double *) malloc(sizeof(double)*nbin);
  
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++){
	i=t1*ntype*ntype+t2*ntype+t3;
	f3s[i] = (double **) malloc(sizeof(double*)*nbin);
	for (ir1=0;ir1<nr;ir1++)
	  for (ir2=0;ir2<=ir1;ir2++){
	    idx=MapIndex(ir1, ir2);
	    /*
	      the width of a bin is dr and for each bin its range is (rmin+dr*ir1, rmin+dr*ir1+dr).
	      Factor "scale" determine the acctual position of the bin
	      which is useful for g2g2g2
	     */
	    r1=tbf->rmin+dr*ir1;
	    r2=tbf->rmin+dr*ir2;
	    
	    f3s[i][idx] = (double *) malloc(sizeof(double)*tbf->nt[idx]);
	    vols[idx] = (double)8.0*pi*pi/9
	      *((r1+dr)*(r1+dr)*(r1+dr)-r1*r1*r1)
	      *((r2+dr)*(r2+dr)*(r2+dr)-r2*r2*r2)
	      *tbf->dt[idx];
	    /* The acctual postion of bins 0: left side; 0.5: the middle; 1: the right side; */
	    r1 += dr*scale;
	    r2 += dr*scale;
	    /* caluclate volume elements and normalize tbf */
	    for (it=0;it<tbf->nt[idx];it++){
	      cost = tbf->tmins[idx]+it*tbf->dt[idx] + tbf->dt[idx]*scale;
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);

	      g2_1 = LinInterp(pdf, t1*ntype+t2, r1);
	      g2_2 = LinInterp(pdf, t1*ntype+t3, r2);
	      g2_3 = LinInterp(pdf, t3*ntype+t2, r3);
	      f3s[i][idx][it] = g2_1*g2_2*g2_3;
	}
      }
  }
  tbf->vols=vols;
  tbf->f3s=f3s;
  return 0;
}

int IntegrateFluct(tbf_p_t tbf, tbf_p_t sup, xdat_p_t xdat, integral_p_t *dest_p){
  *dest_p = (integral_p_t) malloc(sizeof(integral_t));
  integral_p_t dest= *dest_p;

  dest->ntype=tbf->ntype;
  dest->nr = tbf->nr;
  dest->rmin = tbf->rmin;
  dest->rmax = tbf->rmax;
  dest->dr = tbf->dr;
  
  int i,j;
  dest->bins = (double **) malloc(sizeof(double*)*dest->ntype*dest->ntype*dest->ntype);
  for (i=0;i<dest->ntype*dest->ntype*dest->ntype;i++){
    dest->bins[i] = (double *) malloc(sizeof(double)*dest->nr);
    for (j=0;j<dest->nr;j++)
      dest->bins[i][j] = 0;
  }

  int  nr, ntype;
  nr = dest->nr;
  ntype = dest->ntype;

  double rho=(double)xdat->natom/xdat->vol; 
  double val, r1, r2, r3, cost;
  int ir1, ir2, it, idx;
  int t1, t2, t3;
  int p1, p2, p3, pmax;
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++){
	i=t1*ntype*ntype+t2*ntype+t3;
	
	for (ir1=0;ir1<nr;ir1++)
	  for (ir2=0;ir2<=ir1;ir2++){
	    idx=MapIndex(ir1, ir2);
	    r1=dest->rmin+dest->dr*ir1; r2=dest->rmin+dest->dr*ir2;
	    /* caluclate volume elements and normalize tbf */

	    for (it=0;it<tbf->nt[idx];it++){
	      val = tbf->f3s[i][idx][it]-sup->f3s[i][idx][it];

	      cost = tbf->tmins[idx]+it*tbf->dt[idx] + tbf->dt[idx];
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	      p1 = ir1;
	      p2 = ir2;
	      p3 = floor((r3-dest->rmin)/dest->dr);

	      if ( p3 >= nr ) continue;

	      pmax=p1>p2?p1:p2;
	      pmax=pmax>p3?pmax:p3;
	      
	      dest->bins[i][pmax]+=val*tbf->vols[idx];
	    }
	  }

	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[i][ir1] += dest->bins[i][ir1-1];
	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[i][ir1] *= (double)1.0/6*rho*rho;
      }
  return 0;
}

int IntegrateInfo(tbf_p_t tbf, tbf_p_t sup, xdat_p_t xdat, integral_p_t *dest_p){
  *dest_p = (integral_p_t) malloc(sizeof(integral_t));
  integral_p_t dest= *dest_p;

  dest->ntype=tbf->ntype;
  dest->nr = tbf->nr;
  dest->rmin = tbf->rmin;
  dest->rmax = tbf->rmax;
  dest->dr = tbf->dr;
  
  int i,j;
  dest->bins = (double **) malloc(sizeof(double*)*dest->ntype*dest->ntype*dest->ntype);
  for (i=0;i<dest->ntype*dest->ntype*dest->ntype;i++){
    dest->bins[i] = (double *) malloc(sizeof(double)*dest->nr);
    for (j=0;j<dest->nr;j++)
      dest->bins[i][j] = 0;
  }

  int  nr, ntype;
  nr = dest->nr;
  ntype = dest->ntype;

  double rho=(double)xdat->natom/xdat->vol;
  double val, r1, r2, r3, cost;
  int ir1, ir2, it, idx;
  int t1, t2, t3;
  int p1, p2, p3, pmax;
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++){
	i=t1*ntype*ntype+t2*ntype+t3;
	
	for (ir1=0;ir1<nr;ir1++)
	  for (ir2=0;ir2<=ir1;ir2++){
	    idx=MapIndex(ir1, ir2);
	    r1=dest->rmin+dest->dr*ir1; r2=dest->rmin+dest->dr*ir2;
	    /* caluclate volume elements and normalize tbf */

	    for (it=0;it<tbf->nt[idx];it++){
	      /* info */
	      if (tbf->f3s[i][idx][it] < mepsilon)
		val=0;
	      else{
		if (sup->f3s[i][idx][it] < mepsilon)
		  val = 0;
		else
		  val = -tbf->f3s[i][idx][it]*log(tbf->f3s[i][idx][it]/sup->f3s[i][idx][it]);
	      }

	      cost = tbf->tmins[idx]+it*tbf->dt[idx] + tbf->dt[idx];
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	      p1 = ir1;
	      p2 = ir2;
	      p3 = floor((r3-dest->rmin)/dest->dr);

	      if ( p3 >= nr ) continue;

	      pmax=p1>p2?p1:p2;
	      pmax=pmax>p3?pmax:p3;
	      
	      /* sinfo x dV */
	      dest->bins[i][pmax]+= val*tbf->vols[idx];
	    }
	  }

	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[i][ir1] += dest->bins[i][ir1-1];
	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[i][ir1] *= (double)1.0/6*rho*rho;
      }
  return 0;
}

int IntegrateFluctFromEpsilon(tbf_p_t eps, xdat_p_t xdat, integral_p_t *dest_p){
  *dest_p = (integral_p_t) malloc(sizeof(integral_t));
  integral_p_t dest= *dest_p;

  dest->ntype=eps->ntype;
  dest->nr = eps->nr;
  dest->rmin = eps->rmin;
  dest->rmax = eps->rmax;
  dest->dr = eps->dr;
  
  int i,j;
  dest->bins = (double **) malloc(sizeof(double*)*dest->ntype*dest->ntype*dest->ntype);
  for (i=0;i<dest->ntype*dest->ntype*dest->ntype;i++){
    dest->bins[i] = (double *) malloc(sizeof(double)*dest->nr);
    for (j=0;j<dest->nr;j++)
      dest->bins[i][j] = 0;
  }

  int  nr, ntype;
  nr = dest->nr;
  ntype = dest->ntype;

  double rho=(double)xdat->natom/xdat->vol; 
  double val, r1, r2, r3, cost;
  int ir1, ir2, it, idx;
  int t1, t2, t3;
  int p1, p2, p3, pmax;
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++){
	i=t1*ntype*ntype+t2*ntype+t3;
	
	for (ir1=0;ir1<nr;ir1++)
	  for (ir2=0;ir2<=ir1;ir2++){
	    idx=MapIndex(ir1, ir2);
	    r1=dest->rmin+dest->dr*ir1; r2=dest->rmin+dest->dr*ir2;
	    /* caluclate volume elements and normalize tbf */

	    for (it=0;it<eps->nt[idx];it++){
	      val = eps->f3s[i][idx][it];

	      cost = eps->tmins[idx]+it*eps->dt[idx] + eps->dt[idx];
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	      p1 = ir1;
	      p2 = ir2;
	      p3 = floor((r3-dest->rmin)/dest->dr);

	      if ( p3 >= nr ) continue;

	      pmax=p1>p2?p1:p2;
	      pmax=pmax>p3?pmax:p3;
	      
	      dest->bins[i][pmax]+=val*eps->vols[idx];
	    }
	  }

	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[i][ir1] += dest->bins[i][ir1-1];
	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[i][ir1] *= (double)1.0/6*rho*rho;
      }
  return 0;
}

int IntegrateInfoFromEpsilon(tbf_p_t tbf, tbf_p_t eps, xdat_p_t xdat, integral_p_t *dest_p){
  *dest_p = (integral_p_t) malloc(sizeof(integral_t));
  integral_p_t dest= *dest_p;

  dest->ntype=eps->ntype;
  dest->nr = eps->nr;
  dest->rmin = eps->rmin;
  dest->rmax = eps->rmax;
  dest->dr = eps->dr;
  
  int i,j;
  dest->bins = (double **) malloc(sizeof(double*)*dest->ntype*dest->ntype*dest->ntype);
  for (i=0;i<dest->ntype*dest->ntype*dest->ntype;i++){
    dest->bins[i] = (double *) malloc(sizeof(double)*dest->nr);
    for (j=0;j<dest->nr;j++)
      dest->bins[i][j] = 0;
  }

  int  nr, ntype;
  nr = dest->nr;
  ntype = dest->ntype;

  double rho=(double)xdat->natom/xdat->vol;
  double val, r1, r2, r3, cost;
  int ir1, ir2, it, idx;
  int t1, t2, t3;
  int p1, p2, p3, pmax;
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++){
	i=t1*ntype*ntype+t2*ntype+t3;
	
	for (ir1=0;ir1<nr;ir1++)
	  for (ir2=0;ir2<=ir1;ir2++){
	    idx=MapIndex(ir1, ir2);
	    r1=dest->rmin+dest->dr*ir1; r2=dest->rmin+dest->dr*ir2;
	    /* caluclate volume elements and normalize tbf */

	    for (it=0;it<eps->nt[idx];it++){
	      /* info */
	      if (tbf->f3s[i][idx][it] < mepsilon)
		val=0;
	      else{
		if (tbf->f3s[i][idx][it]-eps->f3s[i][idx][it] < mepsilon)
		  val = 0;
		else
		  val = -tbf->f3s[i][idx][it]*log(tbf->f3s[i][idx][it]/(tbf->f3s[i][idx][it]-eps->f3s[i][idx][it]));
	      }

	      cost = eps->tmins[idx]+it*eps->dt[idx] + eps->dt[idx];
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	      p1 = ir1;
	      p2 = ir2;
	      p3 = floor((r3-dest->rmin)/dest->dr);

	      if ( p3 >= nr ) continue;

	      pmax=p1>p2?p1:p2;
	      pmax=pmax>p3?pmax:p3;
	      
	      /* sinfo x dV */
	      dest->bins[i][pmax]+= val*eps->vols[idx];
	    }
	  }

	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[i][ir1] += dest->bins[i][ir1-1];
	for (ir1=1;ir1<nr;ir1++)
	  dest->bins[i][ir1] *= (double)1.0/6*rho*rho;
      }
  return 0;
}

int CheckSuperposMinus(char *fname, tbf_p_t tbf, tbf_p_t sup){
  int i, ir1, ir2, idx, it;

  FILE *fid = fopen(fname, "w");
  fprintf(fid, "%d\n", tbf->nsample);
  fprintf(fid, "%d\n", tbf->ntype);
  fprintf(fid, "%lf %lf %lf\n", tbf->rmin, tbf->rmax, tbf->dr);
  for (i=0;i<tbf->ntype*tbf->ntype*tbf->ntype;i++)
    for (ir1=0;ir1<tbf->nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
        idx=MapIndex(ir1, ir2);
        for (it=0;it<tbf->nt[idx];it++)
          fprintf(fid, "%.12lf ", tbf->f3s[i][idx][it]-sup->f3s[i][idx][it]);
        fprintf(fid, "\n");
      }
  fclose(fid);
  return 0;
}

int CheckSuperposDivide(char *fname, tbf_p_t tbf, tbf_p_t sup){
  int i, ir1, ir2, idx, it;

  FILE *fid = fopen(fname, "w");
  fprintf(fid, "%d\n", tbf->nsample);
  fprintf(fid, "%d\n", tbf->ntype);
  fprintf(fid, "%lf %lf %lf\n", tbf->rmin, tbf->rmax, tbf->dr);
  for (i=0;i<tbf->ntype*tbf->ntype*tbf->ntype;i++)
    for (ir1=0;ir1<tbf->nr;ir1++)
      for (ir2=0;ir2<=ir1;ir2++){
        idx=MapIndex(ir1, ir2);
        for (it=0;it<tbf->nt[idx];it++)
	  if (sup->f3s[i][idx][it] > mepsilon)
	    fprintf(fid, "%.12lf ", tbf->f3s[i][idx][it]/sup->f3s[i][idx][it]);
	  else
	    fprintf(fid, "%.12lf ", (double)0);
        fprintf(fid, "\n");
      }
  fclose(fid);
  return 0;
}
