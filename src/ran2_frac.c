#define IA2 16807
#define IM2 2147483647
#define AM2 (1.0/IM2)
#define IQ2 127773
#define IR2 2836
#define NTAB2 32
#define NDIV2 (1+(IM2-1)/NTAB2)
#define EPS2 1.2e-7
#define RNMX2 (1.0-EPS2)

double ran2(int *idum)
{
	int j;
	int k;
	//static
	static int iy=0;
	static int iv[NTAB2];
	
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB2+7;j>=0;j--) {
			k=(*idum)/IQ2;
			*idum=IA2*(*idum-k*IQ2)-IR2*k;
			if (*idum < 0) *idum += IM2;
			if (j < NTAB2) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ2;
	*idum=IA2*(*idum-k*IQ2)-IR2*k;
	if (*idum < 0) *idum += IM2;
	j=iy/NDIV2;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM2*iy) > RNMX2) return RNMX2;
	else return temp;
}

#undef IA2
#undef IM2
#undef AM2
#undef IQ2
#undef IR2
#undef NTAB2
#undef NDIV2
#undef EPS2
#undef RNMX2
/* (C) Copr. 1986-92 Numerical Recipes Software )1!. */
