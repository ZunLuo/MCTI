#include <math.h>
#include <stdio.h>
#define PI 3.141592654

double poidev(double xm, int *idum)
{
	double gammln(double xx);
	double ran1(int *idum);
	static double sq,alxm,g,oldm=(-1.0);
	double em,t,y;

    //printf("xm,idum=%f,%i\n",xm,*idum);
	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran1(idum);
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran1(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran1(idum) > t);
	}
	return em;
}

#undef PI
/* (C) Copr. 1986-92 Numerical Recipes Software )1!. */
