#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "hist.h"
#include "3bdf.h"

typedef struct pair_dist_func_type{
  int ntype, nbin;
  double rmin, rmax;
  double dr;

  double **bins;
}pdf_t;

typedef pdf_t *pdf_p_t;
typedef pdf_t prob_t;
typedef prob_t *prob_p_t;

typedef struct integral_type{
  int ntype, nr;
  double rmin, rmax;
  double dr;
  
  double **bins;
}integral_t;

typedef integral_t *integral_p_t;

int InitPDF(char *fname, pdf_p_t *pdf_p);

int PrintPDF(char *fname, pdf_p_t pdf);

int FreePDF(pdf_p_t pdf);

int FreeIntegral(integral_p_t dest);

int PrintIntegral(char *fname, integral_p_t dest);

/* Integrate Function B: B=(g2-1)*(g2-1)*(g2-1); same distance with 3bdf */
int IntegrateB(tbf_p_t tbf, pdf_p_t pdf, xdat_p_t xdat, integral_p_t *dest_p);
/* Integrate Function B: B=(g2-1)*(g2-1)*(g2-1); extended to half of the box*/
int IntegrateB_ext(pdf_p_t pdf, xdat_p_t xdat, integral_p_t *dest_p);

double TriLinearIntegral(pdf_p_t pdf, tbf_p_t tbf,
			 int t1, int t2, int t3, /* type */
			 int ir1, int ir2, int it);
double TriLinearIntegralN(pdf_p_t pdf, tbf_p_t tbf,
			 int t1, int t2, int t3, /* type */
			  int ir1, int ir2, int it, int dim);

/* Initialize Superposition Function g2*g2*g2 */
int InitSuperpos(tbf_p_t tbfin, pdf_p_t pdf, tbf_p_t *tbf_p);
int InitSuperposScale(tbf_p_t tbfin, pdf_p_t pdf, double scale, tbf_p_t *tbf_p);
int InitSuperpos_Intgrl(tbf_p_t tbfin, pdf_p_t pdf, tbf_p_t *tbf_p);
int InitSuperpos_IntgrlN(tbf_p_t tbfin, pdf_p_t pdf, tbf_p_t *tbf_p, int dim);

int R3MinMax(tbf_p_t tbf, int ir1, int ir2, int it, double *r3min, double *r3max);
int InitEplison(tbf_p_t tbfin, pdf_p_t pdf, tbf_p_t *eps_p);

int IntegrateFluct(tbf_p_t tbf, tbf_p_t sup, xdat_p_t xdat, integral_p_t *dest_p);
int IntegrateFluctFromEpsilon(tbf_p_t eps, xdat_p_t xdat, integral_p_t *dest_p);

int IntegrateInfo(tbf_p_t tbf, tbf_p_t sup, xdat_p_t xdat, integral_p_t *dest_p);
int IntegrateInfoFromEpsilon(tbf_p_t tbf, tbf_p_t eps, xdat_p_t xdat, integral_p_t *dest_p);

/* print g3-g2g2g2 */
int CheckSuperposMinus(char *fname, tbf_p_t tbf, tbf_p_t sup);
/* print g3/g2g2g2 */
int CheckSuperposDivide(char *fname, tbf_p_t tbf, tbf_p_t sup);

#endif
