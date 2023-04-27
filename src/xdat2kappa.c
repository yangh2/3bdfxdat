#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <unistd.h>
#include <string.h>

typedef double mat33_t[3][3]; 	/* type of 3x3 matrix */
typedef double vec3_t[3];

const double epsilon=1e-6;

int Usage(char* cmd){
  printf("Usage %s: -f xyzfile XDATCARs [-s symmetryfile] [-r rcut] [-m sample_max] [-i nskip]\n", cmd);
  /* printf("-f xyzfile: xyz like file\n" ); */
  /* printf("-s symmetryfile: symmetry operation. Use 'aflow-sym' to generate."); */
  /* printf("-m sample_max: sample maximum."); */
  /* printf("-i nskip: skip n sample runs."); */
  return 0;
}

void reverse(char str[], int length){
  char tmp;
  int start = 0;
  int end = length -1;
  while (start < end){
    tmp=*(str+start);
    *(str+start)=*(str+end);
    *(str+end)=tmp;
    start++;
    end--;
  }
}

char* itoa(int num, char* str, int base)
{
  int i = 0;
  int isNegative = 0;
  
  /* Handle 0 explicitely, otherwise empty string is printed for 0 */
  if (num == 0)
    {
      str[i++] = '0';
      str[i] = '\0';
      return str;
    }
  
  // In standard itoa(), negative numbers are handled only with 
  // base 10. Otherwise numbers are considered unsigned.
  if (num < 0 && base == 10)
    {
      isNegative = 1;
      num = -num;
    }

  // Process individual digits
  while (num != 0)
    {
      int rem = num % base;
      str[i++] = (rem > 9)? (rem-10) + 'a' : rem + '0';
      num = num/base;
    }
  
  // If number is negative, append '-'
  if (isNegative)
    str[i++] = '-';
  
  str[i] = '\0'; // Append string terminator
  
  // Reverse the string
  reverse(str, i);
  
  return str;
}

int main(int argc, char * argv[]){
  const double kb=8.614e-5;	    /* eV/T */
  
  int i,j,k;
  FILE *fid, *fidw;
  const int linelen=1000;
  char line[linelen], fname[100], strtmp[100];

  int natom, count;
  int nmax=2;
  int n2mean=0, nmean2=0, ncount=0;
  int n2mean1=0, nmean21=0, ncount1=0;
  int n2means=0, nmean2s=0, ncounts=0;
  int ***npercell;
  
  int max_sample = 10000000, interv=1;
  int c;
  while ((c = getopt (argc, argv, "n:")) != -1)
    switch (c)
      {
      case 'n':
	nmax=atoi(optarg);
	break;
      default:
	Usage(argv[0]);
        abort ();
      }
  npercell = (int ***) malloc(sizeof(int**)*nmax);
  for (i=0;i<nmax;i++){
    npercell[i] = (int **) malloc(sizeof(int*)*nmax);
    for (j=0;j<nmax;j++)
      npercell[i][j] = (int*)malloc(sizeof(int)*nmax);
  }
  
  double temperature, beta;
  double basis[3][3];
  fid = fopen(argv[optind], "r");
  //for (i=0;i<2;i++) fgets(line, linelen, fid);
  for (i=0;i<2;i++) fgets(line, linelen, fid);
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      fscanf(fid, "%lf", &(basis[i][j]));
  natom=0;
  double vol=fabs(basis[0][0]*basis[1][1]*basis[2][2]-basis[0][0]*basis[1][2]*basis[2][1]-basis[0][1]*basis[1][0]*basis[2][2]+
		  basis[0][1]*basis[1][2]*basis[2][0]+basis[0][2]*basis[1][0]*basis[2][1]-basis[0][2]*basis[1][1]*basis[2][0]);
  for (i=0;i<2;i++) fgets(line, linelen, fid);
  fscanf(fid, "%d", &natom);
  fclose(fid);

  fid = fopen("Trun", "r");
  fscanf(fid, "%lf", &temperature);
  fclose(fid);
  beta = 1/(temperature * kb);

  //printf("\nCalculating density matrix elements ...\n");
  int iargc;
  double dx,dy,dz;
  int nx, ny, nz;
  for (iargc=optind;iargc<argc;iargc++){
    //printf("%s ", argv[iargc]);
    fid = fopen(argv[iargc], "r");
    for (i=0;i<7;i++) fgets(line, linelen, fid);
    while (! feof(fid)){
      fgets(line, linelen, fid);
      if (feof(fid)) break;
      for (i=0;i<nmax;i++)
	for(j=0;j<nmax;j++)
	  for(k=0;k<nmax;k++)
	    npercell[i][j][k]=0;
      for (i=0;i<natom;i++){
	fscanf(fid, "%lf %lf %lf", &dx, &dy, &dz);
	nx = (int)floor(dx*nmax);
	ny = (int)floor(dy*nmax);
	nz = (int)floor(dz*nmax);
	if ( nx  == nmax ) nx --;
	if ( ny  == nmax ) ny --;
	if ( nz  == nmax ) nz --;
	npercell[nx][ny][nz]++;
	fgets(line, linelen, fid);
      }
      n2mean1 += npercell[0][0][0]*npercell[0][0][0];
      nmean21 += npercell[0][0][0];
      ncount1 ++;
      for (i=0;i<nmax;i++)
	for(j=0;j<nmax;j++)
	  for(k=0;k<nmax;k++){
	    //printf("%d\n", npercell[i][j][k]);
	    n2mean += npercell[i][j][k]*npercell[i][j][k];
	    nmean2 += npercell[i][j][k];
	    ncount ++;
	    if ((i==nmax-1) && (j==nmax-1) && (k==nmax-1)) break;
	    n2means += npercell[i][j][k]*npercell[i][j][k];
	    nmean2s += npercell[i][j][k];
	    ncounts ++;
	  }

    }
    fclose(fid);
  }

  double aev2bar=1.6021e4;
  double n2=((double)n2mean)/ncount;
  double nn=((double)nmean2)/ncount;
  double nd=n2-nn*nn;
  double kappa=nd*vol/natom/(natom)/kb/temperature*nmax*nmax*nmax;
  printf("number density: %lf %lf %lf\n", n2, nn, nd);
  printf("compressibility kappa [A3/eV]: %lf\n", kappa);
  printf("compressibility kappa [1/kbar]: %e\n", kappa/aev2bar);

  n2=((double)n2means)/ncounts;
  nn=((double)nmean2s)/ncounts;
  nd=n2-nn*nn;
  kappa=nd*vol/natom/(natom)/kb/temperature*nmax*nmax*nmax;
  printf("number density: %lf\n", nd);
  printf("compressibility kappa [A3/eV]: %lf\n", kappa);
  printf("compressibility kappa [1/kbar]: %e\n", kappa/aev2bar);

  n2=((double)n2mean1)/ncount1;
  nn=((double)nmean21)/ncount1;
  nd=n2-nn*nn;
  kappa=nd*vol/natom/(natom)/kb/temperature*nmax*nmax*nmax;
  printf("number density: %lf\n", nd);
  printf("compressibility kappa [A3/eV]: %lf\n", kappa);
  printf("compressibility kappa [1/kbar]: %e\n", kappa/aev2bar);

  return 0;
}
