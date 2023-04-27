#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hist.h"
#include "3bdf.h"
#include "integrate.h"

int usage(char *argv){
  printf("Usage:\n");
  printf("%s fin dim\n", argv);
  printf("Normalize a histogram from fin and write 3bdf to '3bdf.dat'\n");
  printf("dim is an integer indicate integration grid dimxdimxdim. dim=1, 2, 3, ...\n");
  printf("Write entropy integrations\n");
  return 0;
}

int main(int argc, char *argv[]){
  if ( argc != 3 ){
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
  tbf_p_t sup;
  //InitSuperposScale(tbf, pdf, atof(argv[2]), &sup);
  InitSuperpos_IntgrlN(tbf, pdf, &sup, atoi(argv[2]));
  
  integral_p_t B;
  IntegrateB(tbf, pdf, xdat, &B);
  PrintIntegral("S3B.dat", B);
  FreeIntegral(B);
  IntegrateB_ext(pdf, xdat, &B);
  PrintIntegral("S3B-ext.dat", B);
  FreeIntegral(B);

  integral_p_t Sfluct;
  IntegrateFluct(tbf, sup, xdat, &Sfluct);
  PrintIntegral("S3Fluct.dat", Sfluct);
  FreeIntegral(Sfluct);
  integral_p_t Sinfo;
  IntegrateInfo(tbf, sup, xdat, &Sinfo);
  PrintIntegral("S3Info.dat", Sinfo);
  FreeIntegral(Sinfo);
  
  FreeXdat(xdat);
  FreeTBF(tbf);
  FreePDF(pdf);
  FreeTBF(sup);
}
