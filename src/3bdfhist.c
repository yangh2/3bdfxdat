#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[]){
  char fsamples[100], fhist[100];
  strcpy(fsamples, argv[1]);
  strcpy(fhist, argv[2]);

  int i;
  double rc;
  double r1, r2, cost;
  double dr1, dr2, dcost;
  int nr1, nr2, ncost;
  int ir1, ir2, icost;
  int *nbins;

  FILE *fid = fopen(fhist, "r");
  fscanf(fid, "%lf", &rc);
  fscanf(fid, "%d", &nr1);
  fscanf(fid, "%d", &nr2);
  fscanf(fid, "%d", &ncost);
  nbins = (int *) malloc(sizeof(int)*nr1*nr2*ncost);

  for (i=0;i<nr1*nr2*ncost;i++)
    fscanf(fid, "%d", &(nbins[i]));
  fclose(fid);

  dr1=rc/nr1; dr2=rc/nr2; dcost=2.0/ncost;
  
  fid = fopen(fsamples, "r");
  while (!feof(fid)){
    fscanf(fid, "%lf %lf %lf", &r1, &r2, &cost);
    ir1 = (int)r1/dr1;
    ir2 = (int)r2/dr2;
    icost = (int)cost/dcost;
    nbins[ir1*ir2*icost+ir2*icost+icost]++;
  }
  fclose(fid);

  fid = fopen(fhist, "w");
  fprintf(fid, "%lf\n", rc);
  fprintf(fid, "%d %d %d\n", nr1, nr2, ncost);
  for (i=0;i<nr1*nr2*ncost;i++)
    fprintf(fid, "%d\n", nbins[i]);
  fclose(fid);
}
