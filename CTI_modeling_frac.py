from ctypes import CDLL, POINTER, c_int, c_double,c_float,c_long,c_char_p
from numpy.ctypeslib import ndpointer
import numpy.ctypeslib as clb
import numpy as np
from astropy.io import fits
from scipy.stats import randint
from glob import glob
from datetime import datetime
import os

lib_path = os.path.dirname(os.path.realpath(__file__))

lib_path += "/add_CTI1_frac_no_output.so"
print(lib_path)
lib = CDLL(lib_path)
CTI_simul = lib.__getattr__('CTI_simul')
CTI_simul.argtypes = [POINTER(POINTER(c_double)),c_int,c_int,c_int,\
                    c_int,POINTER(c_double),POINTER(c_double),\
                    c_double,c_double,c_double,c_int,c_int,\
                    POINTER(c_int),c_int,POINTER(c_double),c_int,\
                    POINTER(POINTER(c_double)),c_int,\
                    POINTER(POINTER(c_double)),c_int,c_double,\
                    c_double,c_int]
        
def numpy_matrix_to_int_pointer(arr):
    int_pointer_array = (POINTER(c_int)*arr.shape[0])()
    for i in range(arr.shape[0]):
        arr1 = np.array(arr[i].copy().tolist(),dtype=np.int32)
        int_pointer_array[i] = np.ctypeslib.as_ctypes(arr1)
    return int_pointer_array
def numpy_matrix_to_float_pointer(arr):
    float_pointer_array = (POINTER(c_float)*arr.shape[0])()
    for i in range(arr.shape[0]):
        arr1 = np.array(arr[i].copy().tolist(),dtype=np.float32)
        float_pointer_array[i] = np.ctypeslib.as_ctypes(arr1)
    return float_pointer_array
def numpy_matrix_to_double_pointer(arr):
    double_pointer_array = (POINTER(c_double)*arr.shape[0])()
    for i in range(arr.shape[0]):
        arr1 = np.array(arr[i].copy().tolist(),dtype=np.float64)
        double_pointer_array[i] = np.ctypeslib.as_ctypes(arr1)
    return double_pointer_array
def pointer_to_numpy_matrix(arr_pointer,row,col):
    arr = np.zeros((row,col))
    for i in range(row):
        for j in range(col):
            arr[i,j] = arr_pointer[i][j]
    return arr
def CTI_sim(im,nx,ny,noverscan,nsp,rho_trap,t,beta,w,c,nmax,\
        oversample,trap_seeds,release_seed,cap_prob=None,\
        cap_seed=0,inj_flag=0,inj_seed=0,injc=None,Tpix=None,\
        dark_flag=0,Tdark=None):
    if cap_prob is None:
        cap_prob = np.ones(len(rho_trap))
    if np.max(rho_trap)*oversample>65535:
        raise Exception("rho_trap*ny must be less than 65535")
    if len(rho_trap)!=len(t) or len(rho_trap)!=len(trap_seeds) or len(rho_trap)!=len(cap_prob):
        raise Exception("rho_trap,t and trap_seeds must share the same length")
    image = im.T.astype(np.float64)
    nx_c,ny_c,noverscan_c,nsp_c,nmax_c,oversample_c = c_int(nx),\
        c_int(ny),c_int(noverscan),c_int(nsp),c_int(nmax),\
        c_int(oversample)
    ntotal = ny+noverscan
    beta_c,w_c,c_c = c_double(beta),c_double(w),c_double(c)
    t_p = np.ctypeslib.as_ctypes(t)
    rho_trap_p = np.ctypeslib.as_ctypes(rho_trap)
    cap_prob_p = np.ctypeslib.as_ctypes(cap_prob)
    cap_seed_c = int(cap_seed)
    image_p = numpy_matrix_to_double_pointer(image)
    trap_seeds1 = trap_seeds.astype(np.int32)
    trap_seeds_p = np.ctypeslib.as_ctypes(trap_seeds1)
    release_seed_c = c_int(release_seed)
    image_cti = np.zeros((nx,ntotal))
    image_cti = image_cti.astype(np.float64)
    image_cti_p = numpy_matrix_to_double_pointer(image_cti)
    if inj_flag==0:
        injc1 = np.array([[0]])
        Tpix = 0
        inj_seed = 0
    else:
        if injc is None:
            raise Exception("inj must be given if inj_flag==1")
        injc1 = injc.T
        if injc1.shape[0]!=nx or injc1.shape[1]!=ny:
            raise Exception("injc must have the same shape as input image")
        if Tpix is None:
            raise Exception("Tpix must be given when inj_flag==1")
    if dark_flag==1 and Tdark is None:
        raise Exception("Tdark must be given when dark_flag==1")
    inj_flag_c = c_int(inj_flag)
    injc1_p = numpy_matrix_to_double_pointer(injc1)
    inj_seed_p = c_int(inj_seed)
    Tpix_c = c_double(Tpix)
    if dark_flag==0:
        Tdark = 0
    Tdark_c = c_double(Tdark)
    dark_flag_c = int(dark_flag)
    print(datetime.now())
    CTI_simul(image_p,nx_c,ny_c,noverscan_c,nsp_c,rho_trap_p,\
            t_p,beta_c,w_c,c_c,nmax_c,oversample_c,trap_seeds_p,\
            release_seed_c,cap_prob_p,cap_seed,image_cti_p,inj_flag_c,\
            injc1_p,dark_flag_c,Tdark_c,Tpix_c,inj_seed_p)
    print(datetime.now())
    image_cti_result = np.zeros((nx,ntotal))
    for i in range(nx):
        for j in range(ntotal):
            image_cti_result[i,j] = image_cti_p[i][j]
    return image_cti_result.T

if __name__ =='__main__':
    nx,ny,noverscan,nsp,nmax = 4608,4616,84,3,10
    beta,w,c = 0.478,84700,0
    oversample = 100
    t = np.array([0.74,7.7,37],dtype=np.float64)
    rho_trap = np.array([0.6,1.6,1.4],dtype=np.float64)
    #rho_trap = np.array([0.,0.,0.],dtype=np.float64)
    trap_seeds = np.array([0,100,1000],dtype=np.int32)
    cap_prob = np.array([1,1,1],dtype=np.float64)
    inj_flag = 0
    release_seed,cap_seed,inj_seed = 10000,100000,13
    dark_flag = 0
    Tpix,Tdark = 0.03,100
    image = fits.getdata("input/im_test.fits").astype(np.float64)
    injc = fits.getdata("input/darkc.fits").astype(np.float64)
    image_cti = CTI_sim(image,nx,ny,noverscan,nsp,rho_trap,t,beta,\
            w,c,nmax,oversample,trap_seeds,release_seed,cap_prob,cap_seed,\
            inj_flag,inj_seed,injc,Tpix,dark_flag,Tdark)  
    print(image_cti)
    fits.writeto("output/im_test_CTI_100_py.fits",data=image_cti,overwrite=True)
