
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nrutil.h"
void creattraps(int *seed,int *num_trap,int ny,int nsp, int nmax,\
        double c,double beta,double ***ntrap,int *ny_cal)
{
  //creat ntrap traps along one column
  //sp[i][0] the number of trap in the ith pixel;
  //sp[i][j] (c,1+c) the height of the jth trap in the ith pixel;
  double ran3(int *idum);
  void sort(unsigned long n, double arr[]);
  float tmp,betare,height;
  int i,nyt,ntmp,j,l,k,*num_trap_all;
  double *tmpv;
  tmpv = dvector(1,10000);
  num_trap_all = ivector(0,ny);

  //if(nmax>50){printf(" the upper limit of trap in each pixel is too large, nmax=%d\n",nmax); getchar();}
  
  
  betare=1./beta;
  for(l=0;l<nsp;l++){
    for(i=0;i<ny;i++)ntrap[l][i][0]=0;
  }
  for(i=1;i<=100;i++)tmpv[i] = 0;
  for(i=0;i<ny;i++)num_trap_all[i] = 0;
  for(l=0;l<nsp;l++){
    for(i=0;i<num_trap[l];i++){
      nyt=floor(ran3(&seed[l])*ny);
      ntrap[l][nyt][0]++;
      num_trap_all[nyt]++;
    }
  }
  //printf("num_trap_all:%i,%i,%i,%i,%i,%i,%i,%i\n",num_trap_all[0],num_trap_all[1],num_trap_all[2],num_trap_all[3],num_trap_all[4],num_trap_all[5],num_trap_all[6],num_trap_all[7]);
  //printf("trap allocated\n");
  k=0;
  for(i=0;i<ny;i++)
  {
    if(num_trap_all[i]==0) continue;
    else ny_cal[k]=i;
    k++;
  }
  if(k!=ny)ny_cal[k] = 999999;
  //printf("marker array generated\n");
  i = 0;
  while(ny_cal[i]!=999999&&i<ny){
    nyt = ny_cal[i];
//    printf("%i\n",nyt);
    for(l=0;l<nsp;l++){
      ntmp=ntrap[l][nyt][0];
//      printf("ntmp=%i\n",ntmp);
      if(ntmp==1){
        height=ran3(&seed[l])-c;
        if(height<=0){ntrap[l][nyt][1]=0;}
        else{
          ntrap[l][nyt][1]=pow(height,betare);
        }
        ntrap[l][nyt][ntmp+1]=999999;
//        printf("row %i nsp %i traph=%f,%f\n",nyt,l,ntrap[l][nyt][1],ntrap[l][nyt][2]);
      }
      if(ntmp>1){
        if(ntmp>nmax)ntrap[l][nyt][0]=ntmp=nmax; // upper limite of trap in each pixel is nmax
        for(j=1;j<=ntmp;j++)tmpv[j]=ran3(&seed[l])-c;
        sort(ntmp, tmpv);
        //printf("row %i nsp %i traph=",nyt,l);
        for(j=1;j<=ntmp;j++){
	        height=tmpv[j];
            if(height<=0){ntrap[l][nyt][j]=0;}
	        else{ntrap[l][nyt][j]=pow(height,betare);}
//	        printf("%f,",ntrap[l][nyt][j]);
        }
//        printf("\n");
      }
      ntrap[l][nyt][ntmp+1]=999999;
    }
    i++;
  }
  //printf("ny_cal:%i,%i,%i,%i,%i,%i,%i,%i\n",ny_cal[0],ny_cal[1],ny_cal[2],ny_cal[3],ny_cal[4],ny_cal[5],ny_cal[6],ny_cal[7]);
  //for(i=4000;i<=4010;i++)printf("%i ",ny_cal[i]);
  //printf("\n");
  free_dvector(tmpv,1,100);
  free_ivector(num_trap_all,0,ny);
}
