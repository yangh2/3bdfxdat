#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hist.h"
#include "3bdf.h"
#include "integrate.h"

int usage(char *argv){
  printf("Usage:\n");
  printf("%s fin\n", argv);
  printf("Normalize a histogram from fin and write 3bdf to '3bdf.dat'\n");
  printf("Write entropy integrations\n");
  return 0;
}

int main(int argc, char *argv[]){
  if ( argc != 2 ){
    usage(argv[0]);
    exit(0);
  }
  
  xdat_p_t xdat;
  InitXdat("3bdf_xdat.in", &xdat);
  
  hist_p_t hist;
  InitHist(argv[1], &hist);
  ReadHist(argv[1], hist);
  
  tbf_p_t tbf;
  InitTBF(xdat, hist, &tbf);
  PrintTBF("3bdf.dat", tbf);
  FreeHist(hist);

  pdf_p_t pdf;
  InitPDF("3bdf_pdf.in", &pdf);
  
  integral_p_t B;
  IntegrateB(tbf, pdf, xdat, &B);
  PrintIntegral("S3B.dat", B);
  FreeIntegral(B);

  tbf_p_t eps;
  InitEplison(tbf, pdf, &eps);
  PrintTBF("epsl.dat", eps);
  
  integral_p_t Sfluct;
  IntegrateFluctFromEpsilon(eps, xdat, &Sfluct);
  PrintIntegral("S3Fluct.dat", Sfluct);
  FreeIntegral(Sfluct);
  integral_p_t Sinfo;
  IntegrateInfoFromEpsilon(tbf, eps, xdat, &Sinfo);
  PrintIntegral("S3Info.dat", Sinfo);
  FreeIntegral(Sinfo);
  
  FreeXdat(xdat);
  FreeTBF(tbf);
  FreePDF(pdf);
  FreeTBF(eps);
}
