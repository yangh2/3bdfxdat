#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef double vec_t[3];

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

int main(int argc, char *argv[]){
  const double pi=3.14159265359;
  const double kb=8.617333262e-5;	/* eV/K */
  const double hc=1.23984193e4;		/* eV/angstrom */
  const double massunit=931.49410242e6; /* ev/c^2 */
  int i,j,k,t;
  char fhist[100];

  double rmin, rc;
  int nr1, nr2, ncost;
  double **nbins, **g3_info, **g3_fluct, val,sum=0;
  double **s3_info_r, **s3_fluct_r;
  double **pdf_bins, *pdf_r, pdf_rmin, pdf_rmax;
  int pdf_nr;
  int ntype, nstep;

  strcpy(fhist, argv[1]);
  FILE *fid = fopen(fhist, "r"); /* 3 body distribution function */
  fscanf(fid, "%d", &nstep);
  fscanf(fid, "%d", &ntype);
  fscanf(fid, "%lf %lf", &rmin, &rc);
  fscanf(fid, "%d %d %d", &nr1, &nr2, &ncost);
  nbins = (double **) malloc(sizeof(double *)*nr1*nr2*ncost);
  for (i=0;i<nr1*nr2*ncost;i++)
    nbins[i] = (double *) malloc(sizeof(double)*ntype*ntype*ntype);
  g3_info = (double **) malloc(sizeof(double *)*nr1*nr2*ncost);
  for (i=0;i<nr1*nr2*ncost;i++)
    g3_info[i] = (double *) malloc(sizeof(double)*ntype*ntype*ntype);
  g3_fluct = (double **) malloc(sizeof(double *)*nr1*nr2*ncost);
  for (i=0;i<nr1*nr2*ncost;i++)
    g3_fluct[i] = (double *) malloc(sizeof(double)*ntype*ntype*ntype);
  
  s3_info_r = (double **) malloc(sizeof(double *)*nr1);
  for (i=0;i<nr1;i++)
    s3_info_r[i] = (double *) malloc(sizeof(double)*ntype*ntype*ntype);
  s3_fluct_r = (double **) malloc(sizeof(double *)*nr1);
  for (i=0;i<nr1;i++)
    s3_fluct_r[i] = (double *) malloc(sizeof(double)*ntype*ntype*ntype);
  
  for (i=0;i<nr1*nr2*ncost;i++)
    for (j=0;j<ntype*ntype*ntype;j++){
      fscanf(fid, "%lf", &(nbins[i][j]));
      g3_info[i][j]=0;
      g3_fluct[i][j]=0;
    }
  fclose(fid);

  for (i=0;i<nr1;i++)
    for (j=0;j<ntype*ntype*ntype;j++){
      s3_info_r[i][j]=0;
      s3_fluct_r[i][j]=0;
    }
  
  for (j=0;j<ntype*ntype*ntype;j++){
    sum=0;
    for (i=0;i<nr1*nr2*ncost;i++)
      sum += nbins[i][j];
    /* for (i=0;i<nr1*nr2*ncost;i++) */
    /*   nbins[i][j] /= sum; */
  }
  
  fid = fopen("3bdf_pdf.in", "r"); /* pair distribution function */
  fscanf(fid, "%d %lf %lf", &pdf_nr, &pdf_rmin, &pdf_rmax);
  pdf_r = (double *) malloc(sizeof(double)*pdf_nr);
  pdf_bins = (double **) malloc(sizeof(double *)*pdf_nr);
  for (i=0;i<pdf_nr;i++)
    pdf_bins[i] = (double *) malloc(sizeof(double)*ntype*ntype);
  
  for (i=0;i<pdf_nr;i++){
    fscanf(fid, "%lf", &(pdf_r[i]));
    for (j=0;j<ntype*ntype;j++)
      fscanf(fid, "%lf", &(pdf_bins[i][j]));
  }
  fclose(fid);

  double rho, temperature;	/* rho: 1/angstrom^3; temperature: K */
  int natom;
  int *typelist,tmp;
  double *masslist;
  fid = fopen("3bdf_bins.ini", "r");
  fscanf(fid, "%lf %d", &rho, &natom);
  typelist = (int *) malloc(sizeof(int)*ntype);
  for (i=0;i<ntype;i++) fscanf(fid, "%d", &(typelist[i]));
  fclose(fid);
  fid = fopen("Trun", "r");
  fscanf(fid, "%lf", &temperature);
  fclose(fid);
  fid = fopen("Mass", "r");
  fscanf(fid, "%d", &tmp);
  for (i=0;i<ntype;i++) fscanf(fid, "%d", &tmp);
  masslist = (double *) malloc(sizeof(double)*ntype);
  for (i=0;i<ntype;i++) fscanf(fid, "%lf", &(masslist[i]));
  fclose(fid);
  
  /* See: https://pubs.acs.org/doi/full/10.1021/acs.jpcb.7b10723 */
  double sideal, s1, *s2, *s3, *s3info, *s3fluct;
  double **s2map, **s3map_fluct, **s3map_info;
  double r1,r2,cost,v1,v2,dv1,dv2;
  double dr1=(rc-rmin)/nr1, dr2=(rc-rmin)/nr2, dcost=(double)2.0/ncost;
  s2 = (double *) malloc(sizeof(double)*ntype*ntype);
  for (i=0;i<ntype*ntype;i++) s2[i]=1.0/2*(1-rho*4.0/3*pi*pow(rmin,3)); /* const 0.5 in Eq.(10) */
  s2map = (double **) malloc(sizeof(double *)*nr1);
  for (i=0;i<nr1;i++)
    s2map[i] = (double *) malloc(sizeof(double)*ntype*ntype);
  s3 = (double *) malloc(sizeof(double)*ntype*ntype*ntype);
  s3fluct = (double *) malloc(sizeof(double)*ntype*ntype*ntype);
  s3info = (double *) malloc(sizeof(double)*ntype*ntype*ntype);
  //for (i=0;i<ntype*ntype*ntype;i++) s3[i]=1.0/6*(1-rho*rho*4.0/3*pi*pow(rmin,3)*4.0/3*pi*pow(rmin,3)); /* const 1.0/6 in Eq.(11) */
  for (i=0;i<ntype*ntype*ntype;i++) s3[i]=1.0/6;
  for (i=0;i<ntype*ntype*ntype;i++) s3fluct[i]=0;
  for (i=0;i<ntype*ntype*ntype;i++) s3info[i]=0;
  s3map_fluct = (double **) malloc(sizeof(double *)*nr1*nr2*ncost);
  for (i=0;i<nr1*nr2*ncost;i++)
    s3map_fluct[i] = (double *) malloc(sizeof(double)*ntype*ntype*ntype);
  s3map_info = (double **) malloc(sizeof(double *)*nr1*nr2*ncost);
  for (i=0;i<nr1*nr2*ncost;i++)
    s3map_info[i] = (double *) malloc(sizeof(double)*ntype*ntype*ntype);

  double lambda=sqrt(hc*hc/(2*pi*kb*massunit*temperature));
  sideal=5.0/2; s1=3.0/2;	/* Eq.(3) and Eq.(5) */
  for (i=0;i<ntype;i++){
    sideal -= log(rho*typelist[i]/natom*pow(lambda/sqrt(masslist[i]),3));
    s1 -= log(rho*typelist[i]/natom*pow(lambda/sqrt(masslist[i]),3));
  }
  
  for (t=0;t<ntype*ntype*ntype;t++){
    for (i=0;i<nr1;i++)
      for (j=0;j<nr2;j++)
	for (k=0;k<ncost;k++){
	  r1 = i*dr1+rmin;
	  r2 = (i+1)*dr1+rmin;
	  dv1 = 4*pi*(pow(r2,3)-pow(r1,3))/3;
	  r1 = j*dr2+rmin;
	  r2 = (j+1)*dr2+rmin;
	  dv2 = 2*pi*(pow(r2,3)-pow(r1,3))/3;
	  if ( nbins[i*nr2*ncost+j*ncost+k][t]<=0) continue;
	  nbins[i*nr2*ncost+j*ncost+k][t] = nbins[i*nr2*ncost+j*ncost+k][t]/nstep/(dv1*rho)/(dv2*dcost*rho)/natom;
	}
  }

  double **nbins_r1;
  nbins_r1 = (double **) malloc(sizeof(double *)*nr1*nr2*ncost);
  for (i=0;i<nr1*nr2*ncost;i++){
    nbins_r1[i] = (double *) malloc(sizeof(double)*ntype*ntype);
    for (j=0;j<ntype*ntype;j++)
      nbins_r1[i][j]=0;
  }
  int t1,t2,t3;
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++)
	for (i=0;i<nr1;i++)
	  for (j=0;j<nr2;j++)
	    for (k=0;k<ncost;k++){
	      r1 = (j)*dr2+rmin;
	      r2 = (j+1)*dr2+rmin;
	      dv1 = 2*pi*(pow(r2,3)-pow(r1,3))/3;
	      nbins_r1[i][t1*ntype+t2] += nbins[i*nr2*ncost+j*ncost+k][t1*ntype*ntype+t2*ntype+t3]*(dv1*dcost)/(4*pi/3*pow(rc,3)-4*pi/3*pow(rmin,3));
	    }

  /* See Eq.(10) */
  double g2, g2h, g2l, valh, vall, r1h,r1l;
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (i=1;i<nr1;i++){
	r1h = i*dr1+rmin;
	g2h = nbins_r1[i][t1*ntype+t2];
	r1l = (i-1)*dr1+rmin;
	g2l = nbins_r1[i-1][t1*ntype+t2];
	g2 = (g2l+g2h)/2;
	if ( g2l <=0)
	  vall= 0.5*(g2l-1)*rho;
	else
	  vall = 0.5*(g2l-1 -g2l*log(g2l))*rho;
	if ( g2h <=0)
	  valh= 0.5*(g2h-1)*rho;
	else
	  valh = 0.5*(g2h-1 -g2h*log(g2h))*rho;
	s2map[i][t1*ntype+t2]=(vall+valh)/2;
	s2[t1*ntype+t2] += (vall*4*pi*r1l*r1l*dr1+valh*4*pi*r1h*r1h*dr1)/2;
      }
  
  //================================================================================
  // S3 calculation
  // https://en.wikipedia.org/wiki/Trilinear_interpolation
  //================================================================================
  /* See Eq.(11) */
  double g3, g2bond1, g2bond2, g2bond3;
  double r3, g2r1, g2r2, weight1, weight2;
  int ifloor, iceil, imax;
  double pdf_dr=(pdf_rmax-pdf_rmin)/pdf_nr;
  int pnr = ceil((rc-pdf_rmin)/pdf_dr);
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++)
	for (i=0;i<pnr;i++)
	  for (j=0;j<pnr;j++)
	    for (k=0;k<ncost;k++){
	      r1 = (i+0.5)*pdf_dr+pdf_rmin;
	      r2 = (j+0.5)*pdf_dr+pdf_rmin;
	      cost = (k)*dcost-1;

	      if ( r1>=rc-dr1 ) continue;
	      if ( r2>=rc-dr2 ) continue;
	      if ((r1>rmin-dr1) && (r2>rmin-dr2)) continue;
	      
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	      if ( r3<=pdf_rmin + pdf_dr) continue;
	      
	      ifloor=floor(pdf_nr*(r1-pdf_rmin)/(pdf_rmax-pdf_rmin)); /* linear interpolation for g2 */
	      iceil=ceil(pdf_nr*(r1-pdf_rmin)/(pdf_rmax-pdf_rmin));
	      if ((ifloor < 0 ) || (ifloor >= pdf_nr)) continue;
	      if ((iceil < 0 ) || (iceil >= pdf_nr)) continue;
	      g2r1 = (pdf_rmax-pdf_rmin)/pdf_nr*ifloor+pdf_rmin;
	      g2r2 = (pdf_rmax-pdf_rmin)/pdf_nr*iceil+pdf_rmin;
	      weight1 = (g2r2-r1)/(g2r2-g2r1);
	      weight2 = (r1-g2r1)/(g2r2-g2r1);
	      g2bond1 = weight1*pdf_bins[ifloor][t1*ntype+t2]+weight2*pdf_bins[iceil][t1*ntype+t2];

	      ifloor=floor(pdf_nr*(r2-pdf_rmin)/(pdf_rmax-pdf_rmin));
	      iceil=ceil(pdf_nr*(r2-pdf_rmin)/(pdf_rmax-pdf_rmin));
	      if ((ifloor < 0 ) || (ifloor >= pdf_nr)) continue;
	      if ((iceil < 0 ) || (iceil >= pdf_nr)) continue;
	      g2r1 = (pdf_rmax-pdf_rmin)/pdf_nr*ifloor+pdf_rmin;
	      g2r2 = (pdf_rmax-pdf_rmin)/pdf_nr*iceil+pdf_rmin;
	      weight1 = (g2r2-r2)/(g2r2-g2r1);
	      weight2 = (r2-g2r1)/(g2r2-g2r1);
	      //printf("%lf %lf\n", g2r1, g2r2);
	      g2bond2 = weight1*pdf_bins[ifloor][t1*ntype+t3]+weight2*pdf_bins[iceil][t1*ntype+t3];

	      ifloor=floor(pdf_nr*(r3-pdf_rmin)/(pdf_rmax-pdf_rmin));
	      iceil=ceil(pdf_nr*(r3-pdf_rmin)/(pdf_rmax-pdf_rmin));
	      if ((ifloor < 0 ) || (ifloor >= pdf_nr)) continue;
	      if ((iceil < 0 ) || (iceil >= pdf_nr)) continue;
	      g2r1 = (pdf_rmax-pdf_rmin)/pdf_nr*ifloor+pdf_rmin;
	      g2r2 = (pdf_rmax-pdf_rmin)/pdf_nr*iceil+pdf_rmin;
	      weight1 = (g2r2-r3)/(g2r2-g2r1);
	      weight2 = (r3-g2r1)/(g2r2-g2r1);
	      g2bond3 = weight1*pdf_bins[ifloor][t3*ntype+t2]+weight2*pdf_bins[iceil][t3*ntype+t2];

	      val =(double)1.0/6*(-(g2bond1*g2bond2+g2bond2*g2bond3+g2bond3*g2bond1)+(g2bond1+g2bond2+g2bond3)-1);
	      s3fluct[t1*ntype*ntype+t2*ntype+t3] += val*rho*rho*(4*pi*pow(r1,2)*pdf_dr*2*pi*pow(r2,2)*pdf_dr*dcost);

	      imax=i>j?i:j;
	      r3 = (imax+0.5)*pdf_dr+pdf_rmin;
	      imax=floor((r3-rmin)/dr1);
	      if (imax>=nr1)continue;
	      if (imax<0)
		s3_fluct_r[0][t1*ntype*ntype+t2*ntype+t3] += val*rho*rho*(4*pi*pow(r1,2)*pdf_dr*2*pi*pow(r2,2)*pdf_dr*dcost);
	      else
		s3_fluct_r[imax][t1*ntype*ntype+t2*ntype+t3] += val*rho*rho*(4*pi*pow(r1,2)*pdf_dr*2*pi*pow(r2,2)*pdf_dr*dcost);
	    }
  
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++)
	for (i=0;i<nr1;i++)
	  for (j=0;j<nr2;j++)
	    for (k=0;k<ncost;k++){
	      r1 = (i+0.5)*dr1+rmin;
	      r2 = (j+0.5)*dr2+rmin;
	      cost = (k+0.5)*dcost-1;
	      r3 = sqrt(r1*r1+r2*r2-2*r1*r2*cost);
	      
	      g3 = nbins[i*nr2*ncost+j*ncost+k][t1*ntype*ntype+t2*ntype+t3];

	      ifloor=floor(pdf_nr*(r1-pdf_rmin)/(pdf_rmax-pdf_rmin)); /* linear interpolation for g2 */
	      iceil=ceil(pdf_nr*(r1-pdf_rmin)/(pdf_rmax-pdf_rmin));
	      g2r1 = (pdf_rmax-pdf_rmin)/pdf_nr*ifloor+pdf_rmin;
	      g2r2 = (pdf_rmax-pdf_rmin)/pdf_nr*iceil+pdf_rmin;
	      weight1 = (g2r2-r1)/(g2r2-g2r1);
	      weight2 = (r1-g2r1)/(g2r2-g2r1);
	      if ( ifloor==iceil)
		g2bond1 = pdf_bins[ifloor][t1*ntype+t2];
	      else
		g2bond1 = weight1*pdf_bins[ifloor][t1*ntype+t2]+weight2*pdf_bins[iceil][t1*ntype+t2];

	      ifloor=floor(pdf_nr*(r2-pdf_rmin)/(pdf_rmax-pdf_rmin));
	      iceil=ceil(pdf_nr*(r2-pdf_rmin)/(pdf_rmax-pdf_rmin));
	      g2r1 = (pdf_rmax-pdf_rmin)/pdf_nr*ifloor+pdf_rmin;
	      g2r2 = (pdf_rmax-pdf_rmin)/pdf_nr*iceil+pdf_rmin;
	      weight1 = (g2r2-r2)/(g2r2-g2r1);
	      weight2 = (r2-g2r1)/(g2r2-g2r1);
	      //printf("%lf %lf\n", g2r1, g2r2);
	      if ( ifloor==iceil)
		g2bond2 = pdf_bins[ifloor][t1*ntype+t3];
	      else
		g2bond2 = weight1*pdf_bins[ifloor][t1*ntype+t3]+weight2*pdf_bins[iceil][t1*ntype+t3];

	      ifloor=floor(pdf_nr*(r3-pdf_rmin)/(pdf_rmax-pdf_rmin));
	      iceil=ceil(pdf_nr*(r3-pdf_rmin)/(pdf_rmax-pdf_rmin));
	      g2r1 = (pdf_rmax-pdf_rmin)/pdf_nr*ifloor+pdf_rmin;
	      g2r2 = (pdf_rmax-pdf_rmin)/pdf_nr*iceil+pdf_rmin;
	      weight1 = (g2r2-r3)/(g2r2-g2r1);
	      weight2 = (r3-g2r1)/(g2r2-g2r1);
	      if ( ifloor==iceil)
		g2bond3 = pdf_bins[ifloor][t3*ntype+t2];
	      else
		g2bond3 = weight1*pdf_bins[ifloor][t3*ntype+t2]+weight2*pdf_bins[iceil][t3*ntype+t2];

	      g3_fluct[i*nr2*ncost+j*ncost+k][t1*ntype*ntype+t2*ntype+t3]=(double)1.0/6*(g3-(g2bond1*g2bond2+g2bond2*g2bond3+g2bond3*g2bond1)+(g2bond1+g2bond2+g2bond3)-1);
	      //g3_fluct[i*nr2*ncost+j*ncost+k][t1*ntype*ntype+t2*ntype+t3]*=4*pi*r1*r1*2*pi*r2*r2*rho*rho;
	      
	      
	      if (( g3 <= 0.0001) || ( g2bond1 <= 0.0001)|| ( g2bond2 <= 0.0001) || ( g2bond3 <= 0.0001))
		g3_info[i*nr2*ncost+j*ncost+k][t1*ntype*ntype+t2*ntype+t3] = 0;
	      else
		g3_info[i*nr2*ncost+j*ncost+k][t1*ntype*ntype+t2*ntype+t3] = -1.0/6*g3*log(g3/(g2bond1*g2bond2*g2bond3));
	      //g3_info[i*nr2*ncost+j*ncost+k][t1*ntype*ntype+t2*ntype+t3] = -1.0/6*g3*log(g3);
	      //g3_info[i*nr2*ncost+j*ncost+k][t1*ntype*ntype+t2*ntype+t3]*=4*pi*r1*r1*2*pi*r2*r2*rho*rho;
	      //s3map[i*nr2*ncost+j*ncost+k][t1*ntype*ntype+t2*ntype+t3] = val*rho*rho;
	      //s3[t1*ntype*ntype+t2*ntype+t3] += val*rho*rho *4*pi*r1*r1*dr1 * 2*pi*r2*r2*dr2*dcost;
	      //s3map_r1[i][t1*ntype*ntype+t2*ntype+t3] += val*rho*rho *4*pi*r1*r1*dr1 * 2*pi*r2*r2*dr2*dcost;
	    }


  double xs[6], as[8], gs[8];
  double ix0,iy0,iz0, ix1,iy1,iz1;
  for (t1=0;t1<ntype;t1++)
    for (t2=0;t2<ntype;t2++)
      for (t3=0;t3<ntype;t3++)
	for (i=1;i<nr1;i++)
	  for (j=1;j<nr2;j++)
	    for (k=1;k<ncost;k++){
	      xs[0] = (i-1)*dr1+rmin; /* xs[]: x0, y0, z0, x1, y1, z1 */
	      xs[1] = (j-1)*dr2+rmin;
	      xs[2] = (k-1)*dcost-1;
	      xs[3] = (i+0)*dr1+rmin;
	      xs[4] = (j+0)*dr2+rmin;
	      xs[5] = (k+0)*dcost-1;

	      //ix0 = 4*pi/3*(pow(xs[3],3)-pow(xs[0],3));
	      ix0 = 4*pi/3*dr1*(xs[3]*xs[3]+xs[3]*xs[0]+xs[0]*xs[0]);
	      //iy0 = 2*pi/3*(pow(xs[4],3)-pow(xs[1],3));
	      iy0 = 2*pi/3*dr2*(xs[4]*xs[4]+xs[4]*xs[1]+xs[1]*xs[1]);
	      iz0 = dcost;
	      ix1 = pi*dr1*(xs[3]*xs[3]+xs[0]*xs[0])*(xs[3]+xs[0]);
	      iy1 = pi/2*dr2*(xs[4]*xs[4]+xs[1]*xs[1])*(xs[4]+xs[1]);
	      iz1 = dcost*(xs[5]+xs[2])/2;
	      
	      gs[0] = g3_fluct[(i-1)*nr2*ncost+(j-1)*ncost+(k-1)][t1*ntype*ntype+t2*ntype+t3];
	      gs[1] = g3_fluct[(i-0)*nr2*ncost+(j-1)*ncost+(k-1)][t1*ntype*ntype+t2*ntype+t3];
	      gs[2] = g3_fluct[(i-1)*nr2*ncost+(j-0)*ncost+(k-1)][t1*ntype*ntype+t2*ntype+t3];
	      gs[3] = g3_fluct[(i-0)*nr2*ncost+(j-0)*ncost+(k-1)][t1*ntype*ntype+t2*ntype+t3];
	      gs[4] = g3_fluct[(i-1)*nr2*ncost+(j-1)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];
	      gs[5] = g3_fluct[(i-0)*nr2*ncost+(j-1)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];
	      gs[6] = g3_fluct[(i-1)*nr2*ncost+(j-0)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];
	      gs[7] = g3_fluct[(i-0)*nr2*ncost+(j-0)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];

	      Trilinear(xs, gs, as);
	      s3map_fluct[(i-0)*nr2*ncost+(j-0)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3] = rho*rho*(as[0]*ix0*iy0*iz0 +as[1]*ix1*iy0*iz0 +as[2]*ix0*iy1*iz0 +as[3]*ix0*iy0*iz1 +as[4]*ix1*iy1*iz0 +as[5]*ix1*iy0*iz1 +as[6]*ix0*iy1*iz1 +as[7]*ix1*iy1*iz1);
	      //for(tmpi=0;tmpi<8;tmpi++) printf("%e ", as[tmpi]);printf("\n");
	      //for(tmpi=0;tmpi<6;tmpi++) printf("%e ", xs[tmpi]);
	      //printf("%e %e %e %e %e\n",  ix0, iy0, iz0, as[0], rho*rho*as[0]*ix0*iy0*iz0);
	      gs[0] = g3_info[(i-1)*nr2*ncost+(j-1)*ncost+(k-1)][t1*ntype*ntype+t2*ntype+t3];
	      gs[1] = g3_info[(i-0)*nr2*ncost+(j-1)*ncost+(k-1)][t1*ntype*ntype+t2*ntype+t3];
	      gs[2] = g3_info[(i-1)*nr2*ncost+(j-0)*ncost+(k-1)][t1*ntype*ntype+t2*ntype+t3];
	      gs[3] = g3_info[(i-0)*nr2*ncost+(j-0)*ncost+(k-1)][t1*ntype*ntype+t2*ntype+t3];
	      gs[4] = g3_info[(i-1)*nr2*ncost+(j-1)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];
	      gs[5] = g3_info[(i-0)*nr2*ncost+(j-1)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];
	      gs[6] = g3_info[(i-1)*nr2*ncost+(j-0)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];
	      gs[7] = g3_info[(i-0)*nr2*ncost+(j-0)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];

	      Trilinear(xs, gs, as);
	      s3map_info[(i-0)*nr2*ncost+(j-0)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3] = rho*rho*(as[0]*ix0*iy0*iz0 +as[1]*ix1*iy0*iz0 +as[2]*ix0*iy1*iz0 +as[3]*ix0*iy0*iz1 +as[4]*ix1*iy1*iz0 +as[5]*ix1*iy0*iz1 +as[6]*ix0*iy1*iz1 +as[7]*ix1*iy1*iz1);

	      s3fluct[t1*ntype*ntype+t2*ntype+t3] += s3map_fluct[(i-0)*nr2*ncost+(j-0)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];
	      s3info[t1*ntype*ntype+t2*ntype+t3]  += s3map_info[(i-0)*nr2*ncost+(j-0)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];

	      imax=i>j?i:j;
	      s3_info_r[imax][t1*ntype*ntype+t2*ntype+t3]  += s3map_info[(i-0)*nr2*ncost+(j-0)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];
	      s3_fluct_r[imax][t1*ntype*ntype+t2*ntype+t3] += s3map_fluct[(i-0)*nr2*ncost+(j-0)*ncost+(k-0)][t1*ntype*ntype+t2*ntype+t3];
	    }
  //================================================================================
  //================================================================================
  
  fid =fopen("3bdf_pdf.dat", "w");
  for (i=0;i<nr1;i++){
    r1 = (i+0.5)*dr1+rmin;
    fprintf(fid, "%.6e ", r1);
    for (t=0;t<ntype*ntype;t++)
      fprintf(fid, "%.6e ", nbins_r1[i][t]);
    fprintf(fid, "\n");
  }
  fclose(fid);
  fid =fopen("3bdf_pdf2s2_part.dat", "w");
  for (i=0;i<nr1;i++){
    r1 = (i+0.5)*dr1+rmin;
    fprintf(fid, "%.6e ", r1);
    for (t=0;t<ntype*ntype;t++)
      fprintf(fid, "%.6e ", s2map[i][t]);
    fprintf(fid, "\n");
  }
  fclose(fid);
  
  fid =fopen("3bdf_normalized.dat", "w");
  for (i=0;i<nr1*nr2*ncost;i++){
    for (t=0;t<ntype*ntype*ntype;t++)
      fprintf(fid, "%.6e ", nbins[i][t]);
    fprintf(fid, "\n");
  }
  fclose(fid);
  fid =fopen("3bdf_g3_fluct.dat", "w");
  for (i=0;i<nr1*nr2*ncost;i++){
    for (t=0;t<ntype*ntype*ntype;t++)
      fprintf(fid, "%.6e ", g3_fluct[i][t]);
    fprintf(fid, "\n");
  }
  fclose(fid);
  fid =fopen("3bdf_g3_info.dat", "w");
  for (i=0;i<nr1*nr2*ncost;i++){
    for (t=0;t<ntype*ntype*ntype;t++)
      fprintf(fid, "%.6e ", g3_info[i][t]);
    fprintf(fid, "\n");
  }
  fclose(fid);
  fid =fopen("3bdf_s3map_fluct.dat", "w");
  for (i=0;i<nr1*nr2*ncost;i++){
    for (t=0;t<ntype*ntype*ntype;t++)
      fprintf(fid, "%.6e ", s3map_fluct[i][t]);
    fprintf(fid, "\n");
  }
  fclose(fid);
  fid =fopen("3bdf_s3map_info.dat", "w");
  for (i=0;i<nr1;i++){
    for (t=0;t<ntype*ntype*ntype;t++)
      fprintf(fid, "%.6e ", s3map_info[i][t]);
    fprintf(fid, "\n");
  }
  fclose(fid);

  printf("Sideal = %.6e [kB/at]\n", sideal);
  printf("S1 = %.6e [kB/at]\n", s1);
  printf("S2 = \n");
  for (i=0;i<ntype*ntype;i++)
    printf("%.6e ", s2[i]);
  printf("\n");
  printf("Weighted S2 = \n");
  for (i=0;i<ntype;i++)
    for (j=0;j<ntype;j++)
      printf("%.6e ", s2[i*ntype+j]*typelist[i]*typelist[j]/natom/natom);
  printf("\n");
  printf("S3 = \n");
  for (i=0;i<ntype*ntype*ntype;i++){
    s3[i] += s3info[i] + s3fluct[i];
    printf("%.6e ", s3[i]);
  }
  printf("\n");
  printf("S3info = \n");
  for (i=0;i<ntype*ntype*ntype;i++){
    printf("%.6e ", s3info[i]);
  }
  printf("\n");
  printf("S3fluct = \n");
  for (i=0;i<ntype*ntype*ntype;i++){
    printf("%.6e ", s3fluct[i]);
  }
  printf("\n");
  printf("Weighted S3 = \n");
  for (i=0;i<ntype;i++)
    for (j=0;j<ntype;j++)
      for (k=0;k<ntype;k++)
	printf("%.6e ", s3[i*ntype*ntype+j*ntype+k]*typelist[i]*typelist[j]*typelist[k]/natom/natom/natom);
  printf("\n");
  
  for (j=0;j<ntype*ntype*ntype;j++){
      s3info[j]=0;
      s3fluct[j]=0;
    }
  
  fid =fopen("3bdf_s3_fluct.dat", "w");
  for (i=0;i<nr1;i++){
    for (t=0;t<ntype*ntype*ntype;t++)
      s3fluct[t]+=s3_fluct_r[i][t];
    for (t=0;t<ntype*ntype*ntype;t++)
      fprintf(fid, "%.6e ", s3fluct[t]);
    fprintf(fid, "\n");
  }
  fclose(fid);
  fid =fopen("3bdf_s3_info.dat", "w");
  for (i=0;i<nr1;i++){
    for (t=0;t<ntype*ntype*ntype;t++)
      s3info[t]+=s3_info_r[i][t];
    for (t=0;t<ntype*ntype*ntype;t++)
      fprintf(fid, "%.6e ", s3info[t]);
    fprintf(fid, "\n");
  }
  fclose(fid);

  for (i=0;i<nr1*nr2*ncost;i++)
    free(nbins_r1[i]);
  free(nbins_r1);

  for (i=0;i<nr1*nr2*ncost;i++)
    free(nbins[i]);
  free(nbins);
}
