#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hist.h"

int main(int argc, char *argv[]){
  int fi;
  hist_p_t hist;

  InitHist(argv[1], &hist);
  
  xdat_p_t xdat;
  InitXdat(argv[2], &xdat);
  
  for (fi=2;fi<argc;fi++)
    //Xdat2Hist(argv[fi], xdat, hist);
    Xdat2Hist_NeighList(argv[fi], xdat, hist);

  PrintHist(argv[1], hist);

  FreeHist(hist);
  FreeXdat(xdat);
  printf("DONE\n");
}
