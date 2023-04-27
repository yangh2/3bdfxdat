#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hist.h"
#include "3bdf.h"
#include "integrate.h"

int usage(char *argv){
  printf("Usage:\n");
  printf("%s val fin fout\n", argv);
  printf("val is in [0,1]. 0: use the left side of 3bdf bins for pdf calcualtion. 0.5: use the middle point of 3bdf bins for pdf calcualtion. 1: use the right side of 3bdf bins for pdf calculation.\n");
  printf("Read g3 from fin, g2 from '3bdf_pdf.in' and write g3-g2g2g2 to fout\n");
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
  InitSuperpos_Intgrl(tbf, pdf, &sup);
  
  CheckSuperposMinus(argv[3], tbf, sup);
  
  FreeTBF(tbf);
  FreePDF(pdf);
  FreeTBF(sup);
}
