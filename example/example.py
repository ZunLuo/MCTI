import numpy as np
from astropy.io import fits
from CTI_modeling_frac import CTI_sim

if __name__ =='__main__':
    #Set parameters
    nx,ny,noverscan,nsp,nmax = 200,200,84,3,10
    beta,w,c = 0.478,84700,0
    oversample = 1
    t = np.array([0.74,7.7,37],dtype=np.float64)
    rho_trap = np.array([0.6,1.6,1.4],dtype=np.float64)
    trap_seeds = np.array([0,100,1000],dtype=np.int32)
    cap_prob = np.array([1,1,1],dtype=np.float64)
    inj_flag = 1
    release_seed,cap_seed,inj_seed = 10000,100000,13
    dark_flag = 1
    Tpix,Tdark = 0.03,100
    #Read input image
    image = fits.getdata("input/demo.fits").astype(np.float64)
    #Read charge injection map
    injc = fits.getdata("input/demo_darkc.fits").astype(np.float64)
    #Run MCTI
    image_cti = CTI_sim(image,nx,ny,noverscan,nsp,rho_trap,t,beta,\
            w,c,nmax,oversample,trap_seeds,release_seed,cap_prob,cap_seed,\
            inj_flag,inj_seed,injc,Tpix,dark_flag,Tdark)  
    #Write output image to FITS file
    fits.writeto("output/demo_CTI_dark_inj.fits",data=image_cti,overwrite=True)
