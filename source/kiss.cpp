/* Kiss.c
 * The function to generate random number
 */

//doh Added stdlib.h
#include "stdlib.h"
#include "limits.h"
//doh Added kiss.h and moved seed structure there
#include "kiss.h"
#include <iostream>

using namespace std;

seed_type seed;

/* set the seeds of RNG */
void kisset(unsigned long ii,
	    unsigned long jj,
	    unsigned long kk,
		bool exclaim)
{
  seed.i = ii;     /* Initiate the RNG using */
  seed.j = jj;     /* user-input seeds */
  seed.k = kk;
  if(exclaim)
	cout << "\n kisset (" << ii << ", " << jj << ", " << kk << ")" << endl;
}

void kisseed(unsigned long &ii, unsigned long &jj, unsigned long &kk){
	ii = seed.i;
	jj = seed.j;
	kk = seed.k;
}

void kissprint(){
	cout << "\n kissprint (" << seed.i << ", " << seed.j << ", " << seed.k << ")" << endl;
}

/* the real part of RNG */
long kiss()
{
  seed.j = seed.j ^ (seed.j<<17);
  seed.k = (seed.k ^(seed.k<<18)) & 0x7fffffff;
  return (seed.i = 69069*seed.i + 23606797) +
    (seed.j ^= (seed.j>>15)) + (seed.k ^= (seed.k>>13));
}

/* generates random number in [0,1) */
double rkiss()
{
 double r;
 r=((double) kiss())/(2*(LONG_MAX+1.0)) + 0.5;
 return(r);
}

/* generates random number in (0,1) */
double r2kiss()
{
 double r;
 r=((double) kiss())/(2*(LONG_MAX+1.0)) + 0.5;
 while(r<=0.0)
       {
       r=((double) kiss())/(2*(LONG_MAX+1.0)) + 0.5;
       }
 return(r);
}



