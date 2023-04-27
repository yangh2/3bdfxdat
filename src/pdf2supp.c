#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hist.h"
#include "3bdf.h"
#include "integrate.h"
#include "random_16807.h"

#include "supp.h"

int main(int argc, char *argv[]){

  int ntot=atoi(argv[2]);
  
  pdf_p_t pdf;
  InitPDF("3bdf_pdf.in", &pdf);
  
  hist_p_t hist;
  InitHist(argv[1], &hist);

  prob_p_t prob;
  PDF2Prob(hist, pdf, &prob);

  getseed();

  int nloop=1000;

  while (hist->nsample < nloop*ntot)
    PDF2Hist(prob, hist);
  hist->nsample/=nloop;
  
  //PrintPDF("tmp", prob);
  PrintHist(argv[1], hist);
  
  FreePDF(pdf);
  FreeHist(hist);
}
