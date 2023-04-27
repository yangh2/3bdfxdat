#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hist.h"
#include "3bdf.h"

int usage(char *argv){
  printf("Usage: print angular distribution at R1, R2\n");
  printf("%s fname R1 R2\n", argv);
  printf("output\n");
  printf("theta[0,PI] R3[A] f\n");
  return 0;
}

int main(int argc, char *argv[]){
  tbf_p_t tbf;

  if ( argc != 4){
    usage(argv[0]);
    exit(0);
  }
  
  ReadTBF(argv[1], &tbf);
  PrintTBFatR1R2(tbf, atof(argv[2]), atof(argv[3]));
  
  FreeTBF(tbf);
}
