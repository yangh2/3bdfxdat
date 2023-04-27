#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hist.h"
#include "3bdf.h"
#include "integrate.h"
#include "random_16807.h"

#include "supp.h"

int main(int argc, char *argv[]){
  hist_p_t hist;
  InitHist(argv[1], &hist);
  ReadHist(argv[1], hist);

  tbf_p_t supp;
  Hist2Supp(hist, &supp);
  
  PrintTBF("supp.dat", supp);
  
  FreeHist(hist);
  FreeTBF(supp);
}
