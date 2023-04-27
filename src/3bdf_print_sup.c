#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hist.h"
#include "3bdf.h"
#include "integrate.h"

int usage(char *argv){
  printf("Usage:\n");
  printf("%s dim fin fout\n", argv);
  printf("val is an integer indicate integration grid dimxdimxdim. dim=1, 2, 3, ...\n");
  printf("Read g3 from fin, g2 from '3bdf_pdf.in' and write g2g2g2 to fout\n");
  return 0;
}

int main(int argc, char *argv[]){
  if (argc != 4){
    usage(argv[0]);
    exit(0);
  }
  
  tbf_p_t tbf;
  ReadTBF(argv[2], &tbf);

  pdf_p_t pdf;
  InitPDF("3bdf_pdf.in", &pdf);
  tbf_p_t sup;
  //InitSuperposScale(tbf, pdf, atof(argv[1]), &sup);
  InitSuperpos_IntgrlN(tbf, pdf, &sup, atoi(argv[1]));
  PrintTBF(argv[3], sup);
  
  FreeTBF(tbf);
  FreePDF(pdf);
  FreeTBF(sup);
}
