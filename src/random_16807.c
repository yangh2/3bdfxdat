
#include <stdlib.h>
#include <time.h>

#include "random_16807.h"

long SEED;

long getseed(){
  
  SEED =time(NULL);	
  if (SEED) 
    return 0;
  else
    return 1;
}

long random_16807(){
  long FRUIT;
  long q;
  long r;
  
  q = myrandom_max/a; 
  r = myrandom_max % a;	
  FRUIT = a*(SEED % q) - r*(SEED/q);	
  if (FRUIT < 0)	
    FRUIT = FRUIT + myrandom_max;
  SEED = FRUIT;	
  return FRUIT;
}

