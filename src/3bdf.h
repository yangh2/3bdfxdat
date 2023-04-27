#ifndef BDFH
#define BDFH

#include "hist.h"

typedef struct threebodyfunc_t{
  int nsample, ntype;
  double rmin, rmax, dr;
  int nr;
  /*
    "nbin" is the number of bins for r1,r2 leaving out cos(theta).
    "nt" is the number of bins for cos(theta).
    "f3s" is the normalized 3-body distribution function. Its first index denote types
    of species, the second denote r1,r2 and the third one is cos(theta).
    "vols" is an array of volume elements.
   */
  double ***f3s;
  double *vols;
  long nbin, *nt;
  /*
    "tmins" store minimum of cos(theta) for given r1 and r2
    "tmaxs" store maximum of cos(theta) for given r1 and r2
    "dt" store their interval
  */
  double *tmins, *tmaxs, *dt;
} tbf_t;

typedef tbf_t *tbf_p_t;

typedef tbf_t supp_t;
typedef tbf_t *supp_p_t;

extern double pi;
extern double kb;
extern double hc;
extern double massunit;

/*
  Initialize a three body distribution function
  and normalize it with volume elements
*/
int InitTBF(xdat_p_t xdat, hist_p_t hist, tbf_p_t *tbf_p);

/* free a three body distribution function */
int FreeTBF(tbf_p_t tbf);

int ReadTBF(char *fname, tbf_p_t *tbf);

int PrintTBF(char *fname, tbf_p_t tbf);

int PrintTBFatR(tbf_p_t tbf, double r);

int PrintTBFatR1R2(tbf_p_t tbf, double r1, double r2);
#endif
