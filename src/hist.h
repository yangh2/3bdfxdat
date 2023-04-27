#ifndef HIST
#define HIST

/* define a data type for xdatcar */
typedef struct xdatcar_t{
  int natom, ntype;
  int *typemap;
  double vol;
  double a1,a2,a3;
  double b1,b2,b3;
  double c1,c2,c3;
}xdat_t;

typedef xdat_t *xdat_p_t;

typedef struct histogram_t{
  int nsample, ntype;
  double rmin, rmax, dr;
  int nr;
  /*
    "nbin" is the number of bins for r1,r2 leaving out cos(theta).
    "nt" is the number of bins for cos(theta).
    "bins" store counts of triplet. Its first index denote types
    of species, the second denote r1,r2 and the third one is
    cos(theta).
   */
  long nbin, ***bins, *nt;
  /*
    "tmins" store minimum of cos(theta) for given r1 and r2
    "tmaxs" store maximum of cos(theta) for given r1 and r2
    "dt" store their interval
  */
  double *tmins, *tmaxs, *dt;
} hist_t;

typedef hist_t *hist_p_t;

typedef struct neighborlist_type{
  int atom_id;
  struct neighborlist_type *next;
}neighl_t;
typedef neighl_t *neighl_p_t;

typedef double vec_t[3];
extern double mepsilon;

/* Determine the minimum of Cos(theta) for given r1 and r2 */
double CosThetaMin(double rc, double r1, double r2);

/* Determine the minimum of Cos(theta) for given r1 and r2 */
double CosThetaMax(double rmin, double r1, double r2);

/* map 2D index (ir1, ir2) into  1D index*/
long MapIndex(long ir1, long ir2);

/* Initilize a xyz or xdat file */
int InitXdat(char *fxdat, xdat_p_t *xdat_p);

/* free xdat */
int FreeXdat(xdat_p_t xdat);
  
/*
  Initilize a histogram: read data from input file "fname"
   and allocate memory for arrays.
 */
int InitHist(char *fname, hist_p_t *hist_p);

/* free up hist structure data */
int FreeHist(hist_p_t hist);

/* read hist from "fname". Dont allocate memory */
int ReadHist(char *fname, hist_p_t hist);

/* add hist2 to hist1 */
int AddHist(hist_p_t hist1, hist_p_t hist2);

/* print a hist to file "fname" */
int PrintHist(char *fname, hist_p_t hist);

/* read data from a xdatcar and count number of triplets*/
int Xdat2Hist(char *fxyz, xdat_p_t xdat, hist_p_t hist);

/*
  Read data from a xdatcar and count number of triplets
  Using Neighborlist to reduce unwanted pairs;
 */
int Xdat2Hist_NeighList(char *fxyz, xdat_p_t xdat, hist_p_t hist);

#endif
