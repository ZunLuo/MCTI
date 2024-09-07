#define IA4 16807
#define IM4 2147483647
#define AM4 (1.0/IM4)
#define IQ4 127773
#define IR4 2836
#define NTAB4 32
#define NDIV4 (1+(IM4-1)/NTAB4)
#define EPS4 1.2e-7
#define RNMX4 (1.0-EPS4)

double ran4(int *idum)
{
	int j;
	int k;
	//static
	static int iy=0;
	static int iv[NTAB4];
	
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB4+7;j>=0;j--) {
			k=(*idum)/IQ4;
			*idum=IA4*(*idum-k*IQ4)-IR4*k;
			if (*idum < 0) *idum += IM4;
			if (j < NTAB4) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ4;
	*idum=IA4*(*idum-k*IQ4)-IR4*k;
	if (*idum < 0) *idum += IM4;
	j=iy/NDIV4;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM4*iy) > RNMX4) return RNMX4;
	else return temp;
}

#undef IA4
#undef IM4
#undef AM4
#undef IQ4
#undef IR4
#undef NTAB4
#undef NDIV4
#undef EPS4
#undef RNMX4
/* (C) Copr. 1986-92 Numerical Recipes Software )1!. */
