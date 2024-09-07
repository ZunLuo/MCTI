#include <stdio.h>
#define IA3 16807
#define IM3 2147483647
#define AM3 (1.0/IM3)
#define IQ3 127773
#define IR3 2836
#define NTAB3 32
#define NDIV3 (1+(IM3-1)/NTAB3)
#define EPS3 1.2e-7
#define RNMX3 (1.0-EPS3)

double ran3(int *idum)
{
	int j;
	int k;
	//static
	static int iy=0;
	static int iv[NTAB3];
	
	double temp;
    //printf("ran3_seed=%i\n",*idum);
	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB3+7;j>=0;j--) {
			k=(*idum)/IQ3;
			*idum=IA3*(*idum-k*IQ3)-IR3*k;
			if (*idum < 0) *idum += IM3;
			if (j < NTAB3) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ3;
	*idum=IA3*(*idum-k*IQ3)-IR3*k;
	if (*idum < 0) *idum += IM3;
	j=iy/NDIV3;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM3*iy) > RNMX3){
        //printf("ran3=%f\n",RNMX3);
        return RNMX3;
    }
	else{
        //printf("ran3=%f\n",temp);
        return temp;
    }
}

#undef IA3
#undef IM3
#undef AM3
#undef IQ3
#undef IR3
#undef NTAB3
#undef NDIV3
#undef EPS3
#undef RNMX3
/* (C) Copr. 1986-92 Numerical Recipes Software )1!. */
