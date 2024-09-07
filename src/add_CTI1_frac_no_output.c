#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nrutil.h"
#include <time.h>

double poidev(double x, int *idum);
double poidev1(double x, int *idum);
void creattraps(int *seed,int *num_trap,int ny,int nsp,int nmax,double c,\
        double beta,double ***ntrap,int *ny_cal);
void addCTI(double *a0,int ny,int noverscan,int nsp, int nmax, double w,\
        double *t, int oversample, double roversample, double ***ntrap,\
        int *ny_cal, double *acti, int *release_seed, double *alpha,\
        int *cap_seed, int inj_flag, double *injn, int *inj_seed);
//int read_fits_2D(const char *argv,float *galaxy,int imagesize);
//int write_fits_2D(const char *argv,float **stars,int *dim);

void CTI_simul(double **image, int nx, int ny, int noverscan, int nsp,\
        double *rho_trap, double *t, double beta, double w, double c,\
        int nmax, int oversample, int *trap_seeds, int release_seed,\
        double *cap_prob, int cap_seed,double **imagecti, int inj_flag,\
        double **injc, int dark_flag, double Tdark, double Tpix,\
        int inj_seed){
  int ntotal,m;
  printf("image parameters: nx=%i ny=%i noverscan=%i\n",nx,ny,noverscan);
  printf("rho_trap ");
  for(m=0;m<nsp;m++)printf("rho%i=%f,",m,rho_trap[m]);
  printf("\nt ");
  for(m=0;m<nsp;m++)printf("t%i=%f,",m,t[m]);
  printf("\nvolume parameter beta=%f,w=%f,c=%f\n",beta,w,c);
  printf("nsp=%i,nmax=%i,oversample=%i\n",nsp,nmax,oversample);
  printf("trap_seeds = ");
  for(m=0;m<nsp;m++)printf("%i,",trap_seeds[m]);
  printf("\nrelease_seed = %i\n",release_seed);
  printf("cap_seed = %i\n",cap_seed);
  printf("cap_prob = ");
  for(m=0;m<nsp;m++)printf("%f,",cap_prob[m]);
  printf("\n");
  if(inj_flag==1){
    printf("injection parameters: dark_flag=%i Tdark=%f Tpix=%f inj_seed=%i\n",\
            dark_flag,Tdark,Tpix,inj_seed);
  }
  printf("begin CTI simulation\n");
  double ***ntrap;
  double **sp,**injn,*alpha;
  double tmp, prob_sum;
  double *a0,*acti,roversample=1./oversample;
  printf("roversample=%f\n",roversample);
  int ntrap_tmp,*num_trap;
  int dim[2],*ny_cal;
  int i,j,k,l;
  ntotal=ny+noverscan;
  ntrap = d3tensor(0,nsp,0,ny,0,(nmax+1)*oversample);  
  sp = dmatrix(0,ny,0,(nmax+1)*oversample);
  injn = dmatrix(0,nx,0,ny);
  a0=dvector(0,ny);
  acti=dvector(0,ntotal);
  ny_cal = ivector(0,ny+1);
  num_trap = ivector(0,nsp);
  alpha = dvector(0,nsp);
  prob_sum = 0;
  
  for(i=0;i<=ny;i++)ny_cal[i] = 999999;
  for(i=0;i<nsp;i++)prob_sum = prob_sum+cap_prob[i];
  for(i=0;i<nsp;i++)alpha[i] = cap_prob[i]/prob_sum;
  if(inj_flag==1){
    for(i=0;i<nx;i++)for(j=0;j<ny;j++)injn[i][j]=injc[i][j]*Tpix;
  }

  for(k=0;k<nx;k++){
    printf("column k=%i/%i\r",k,nx);
    for(i=0;i<ny;i++)ny_cal[i]=0;
    for(l=0;l<nsp;l++){
        tmp = rho_trap[l]*ny*oversample;
        num_trap[l] = poidev(tmp,&trap_seeds[l]);
        //printf("tmp=%f,num_trap[%i]=%i,trap_seeds[l]=%i\n",tmp,l,num_trap[l],trap_seeds[l]);
    }
    creattraps(trap_seeds,num_trap,ny,nsp,nmax*oversample,c,beta,ntrap,ny_cal);
    for(i=0;i<ntotal;i++)acti[i]=0;
    if(dark_flag==0)for(i=0;i<ny;i++)a0[i]=image[k][i];
    else for(i=0;i<ny;i++)a0[i]=image[k][i]+poidev1(injc[k][i]*Tdark*oversample,&inj_seed)*roversample;
    //if(k==nx-1){
      //  printf("a0[1]=%f\n",a0[1]);
      //  printf("ntrap0=%f,%f,%f,%f,%f\n",ntrap[0][0][0],ntrap[0][0][1],ntrap[0][0][2],ntrap[0][0][3],ntrap[0][0][4]);
      //  printf("ntrap1=%f,%f,%f,%f,%f\n",ntrap[1][0][0],ntrap[1][0][1],ntrap[1][0][2],ntrap[1][0][3],ntrap[1][0][4]);
      //  printf("ntrap2=%f,%f,%f,%f,%f\n",ntrap[2][0][0],ntrap[2][0][1],ntrap[2][0][2],ntrap[2][0][3],ntrap[2][0][4]);
    //}
    addCTI(a0,ny,noverscan,nsp,nmax,w,t,oversample,roversample,ntrap,\
            ny_cal,acti,&release_seed,alpha,&cap_seed,inj_flag,injn[k],\
            &inj_seed); 
    for(i=0;i<ntotal;i++)imagecti[k][i]=acti[i];
  }
  printf("\nCTI simulation finished!\n");
  free_dvector(a0,0,ny);
  free_dvector(acti,0,ntotal);
  free_d3tensor(ntrap,0,nsp,0,ny+1,0,(nmax+1)*oversample);
  free_dmatrix(sp,0,ny,0,(nmax+1)*oversample);
  free_dmatrix(injn,0,nx,0,ny);
  free_ivector(ny_cal,0,ny);
  free_ivector(num_trap,0,nsp);
  free_dvector(alpha,0,nsp);
}

void addCTI(double *a0,int ny,int noverscan, int nsp, int nmax, double w,\
        double *t,int oversample,double roversample,double ***ntrap,\
        int *ny_cal, double *acti,int *release_seed,double *alpha,\
        int *cap_seed, int inj_flag, double *injn,int *inj_seed){
    double ran2(int *idum);
    double height,wre,*f,*h,*ntrap_sp,rnd,flow,flowo;
    double htrap,*alpha_arr,rnd_sum,rnd_cum,rel_rnd;
    int **watermark,*ws,wsp;
    int i,j,k,l,lmax,m,sp,ntraptot,ntotal,mark,nyc,trap_num,*ls,*sps;
    double inject_num,*inject_nums,inj_num,inject_sum=0;

    f = dvector(0,nsp);
    for(sp=0;i<nsp;i++)f[i] = 1-exp(-1/t[i]);
    watermark = imatrix(0,ny,0,nsp);
    inject_nums = dvector(0,ny);
    h = dvector(0,nsp);
    ls = ivector(0,nsp*nmax*oversample);
    sps = ivector(0,nsp*nmax*oversample);
    alpha_arr = dvector(0,nmax*nsp*oversample);
    for(i=0;i<ny;i++)for(sp=0;sp<nsp;sp++)watermark[i][sp] = 0;
    wre = 1/w;
    ntotal=ny+noverscan;
    mark = 0;
    if(inj_flag == 0)for(j=0;j<=ny;j++)inject_nums[j] = 0;
    for(i=0;i<ny+noverscan;i++){
        if(i<ny){
            flow = a0[i];
            if(mark<0)mark=0;
            while(ny_cal[mark]<=i&&ny_cal[mark]!=999999)mark++;
            mark = mark-1;
        }
        else flow=0;
        //ny_cal record the rows with traps, rows with no traps will be skipped
        flowo = flow*oversample;
        if(inj_flag==1){
            if(mark==-1){
                inject_num = 0;
                for(l=0;l<=i;l++)inject_num = inject_num + poidev1(injn[l]*oversample,inj_seed);
                acti[i] = (flowo+inject_num)*roversample;
                inject_sum = inject_sum+inject_num;
                continue;
            }
            for(j=mark;j>=0;j--){
                nyc = ny_cal[j];
                inject_num = 0;
                if(j==0){
                    if(ny_cal[j+1]<=i)lmax= ny_cal[j+1];
                    else lmax = i+1;
                    for(l=0;l<lmax;l++){
                        inj_num = poidev1(injn[l]*oversample,inj_seed);
                        inject_num = inject_num + inj_num;
                    }
                }
                else if(j!=mark){
                    for(l=nyc;l<ny_cal[j+1];l++){
                        inj_num  = poidev1(injn[l]*oversample,inj_seed);
                        inject_num = inject_num + inj_num;
                    }
                }
                else if(i<ny){
                    for(l=nyc;l<=i;l++){
                        inj_num = poidev1(injn[l]*oversample,inj_seed);
                        inject_num = inject_num + inj_num;
                    }
                }
                else{ 
                    for(l=nyc;l<ny;l++){
                        inj_num = poidev1(injn[l]*oversample,inj_seed);
                        inject_num = inject_num + inj_num;
                    }
                }
                inject_nums[j] = inject_num;
                inject_sum = inject_sum+inject_num;
            }
        }
        for(j=mark;j>=0;j--){
            flowo = flowo+inject_nums[j];
            height = flowo*wre*roversample;
            nyc = ny_cal[j];  //the jth row with trap
            ws = watermark[nyc];
            ntraptot = 0;
            for(sp=0;sp<nsp;sp++){
                h[sp] = ntrap[sp][nyc][ws[sp]];
                if(ws[sp]==0)h[sp] = 0;
                else if(h[sp]>1) h[sp] -= 1;
                ntraptot = ntraptot+ntrap[sp][nyc][0];
            }   
            if(flowo<ntraptot&&flowo!=0){
                trap_num = 0;
                rnd_sum = 0;
                m = 0;
                for(sp=0;sp<nsp;sp++){//release
                    l = 0;
                    htrap = 0;
                    ntrap_sp = ntrap[sp][nyc];
                    while(l<ntrap[sp][nyc][0]){
                        htrap = ntrap_sp[l+1];
                        if(htrap>1)htrap--;
                        if(height<htrap)break;
                        if(ntrap_sp[l+1]<=1){
                            sps[trap_num] = sp;
                            ls[trap_num] = l+1;
                            rnd_sum = rnd_sum+alpha[sp];
                            alpha_arr[trap_num] = alpha[sp];
                            trap_num++;
                        }
                        l++;
                    }
                    if(height>h[sp])ws[sp]=l;
                    for(k=l;k<=ntrap_sp[0]-1;k++){
                        if(ntrap_sp[k+1]>1){
                            rel_rnd = ran2(release_seed);
                            if(f[sp]>rel_rnd){
                                flowo++;
                                ntrap_sp[k+1]--;
                            }
                        }
                    }
                }
                while(trap_num>0&&flowo>0){
                    rnd = ran2(cap_seed)*rnd_sum;
                    m = 0;
                    rnd_cum=alpha_arr[0];
                    while(rnd>rnd_cum){
                        rnd_cum = rnd_cum+alpha_arr[m+1];
                        m++;
                    }
                    sp = sps[m];
                    l = ls[m];
                    ntrap[sp][nyc][l] = ntrap[sp][nyc][l]+1;
                    flowo--;
                    rnd_sum = rnd_sum-alpha[sp];
                    alpha_arr[m] = 0;
                    trap_num--;
                }
                continue;
            }
            for(sp=0;sp<nsp;sp++){
                ntrap_sp = ntrap[sp][nyc];
                if(ntrap_sp[0]!=0){
                    if(height>h[sp]){//trap only
                        for(k=1;k<=ntrap_sp[0];k++){
                            if(ntrap_sp[k]>1)continue;
                            if(ntrap_sp[k]>height)break;
                            flowo--;
                            ntrap_sp[k] += 1;
                        }
                        ws[sp] = k-1;
                    }
                    else{
                        for(k=1;k<=ntrap_sp[0];k++){
                            if(ntrap_sp[k]>1){
                                rnd = ran2(release_seed);
                                if(f[sp]>rnd){//release
                                    flowo++;
                                    ntrap_sp[k] -= 1;
                                }
                            }
                            if(ntrap_sp[k]<height&&ntrap_sp[k]<1){//trap
                                flowo--;
                                ntrap_sp[k] += 1;
                            }
                        }
                    }
                }
            }
        }
        acti[i] = flowo*roversample;
    }
    free_imatrix(watermark,0,ny,0,nsp);
    free_dvector(inject_nums,0,ny);
    free_dvector(f,0,nsp);
    free_dvector(alpha_arr,0,nsp*nmax);
    free_ivector(sps,0,nsp*nmax);
    free_ivector(ls,0,nsp*nmax);
}

