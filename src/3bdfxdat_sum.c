#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hist.h"

int usage(char *argv){
  printf("Usage:\n");
  printf("%s fout fins\n", argv);
  printf("Sum histograms in fins and write to fout.\n");
  return 0;
}

int main(int argc, char *argv[]){
  if ( argc <= 2 ){
    usage(argv[0]);
    exit(0);
  }
  
  char fout[100];
  strcpy(fout, argv[1]);

  int k;
  char fhist[100];

  hist_p_t hist, hist_tmp;
  strcpy(fhist, argv[2]);
  InitHist(fhist, &hist);
  ReadHist(fhist, hist);
  InitHist(fhist, &hist_tmp);
  
  for (k=3;k<argc;k++){
    strcpy(fhist, argv[k]);
    ReadHist(fhist, hist_tmp);
    AddHist(hist, hist_tmp);
  }

  PrintHist(fout, hist);
  FreeHist(hist);
  FreeHist(hist_tmp);
}
