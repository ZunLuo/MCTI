# MCTI
Monte-Carlo Simulation for Charge Transfer Inefficiency
# I. installation
1. requirement: gcc compiler, astropy, numpy, scipy
2. change $install_path in install.sh
3. run: bash install.sh
# II. usage
add "from CTI_modeling_frac import CTI_sim" to your python manuscript
CTI_sim(image,nx,ny,noverscan,nsp,rho_trap,t,beta,w,c,nmax,oversample,\
        trap_seeds,release_seed,cap_prob,cap_seed,inj_flag,inj_seed,injc,\
        Tpix,dark_flag_Tdark):
1. image parameters
image: 2D float64 numpy array, input image with size (ny,nx)
nx: int32, column number of the image
ny: int32, row number of the image
noverscan: int32, number of parallel overscan pixels
2. trap parameters
nsp: int32, number of trap species
rho_trap: float64 numpy array, trap densities
t: float64 numpy array, trap release timescales(with unit of pix), with the same order as rho_trap
3. CCD parameters
beta: float64, power exponent of electron cloud volume 
w: float64, fullwell capacity
c: float64, notch channel depth
The volume of the electron cloud V is determined by V/Vmax = max{(N/w-c)^beta,0}\
where Vmax is the maximum volume allowed in a pixel, N is the electron count in the pixel
4. calculation parameters
nmax: int32, maximum number of traps allowed to be generated in a single pixel
oversample: int32, number of trap will be multiplied by oversample, while one trap will trap/release 1/oversample electrons
5. random seeds
trap_seed: int32 numpy array, with the same order as rho_trap, random seeds to generate traps
release_seed: int32, random seed to determine trap release process
6. capture parameters
cap_prob: float64 numpy array, default [1,1,1], with the same order as rho_trap, when the number of free electrons is not enough for trapping,\
        each free electron will be assigned to a random trap, this array gives the relative capture probability of each type of trap
cap_seed: int32, default 0, random seed to determine the capture process when free electron number is not enough
7. charge injection parameters (only for inj_flag=1)
inj_flag: int32, defalut 0, set to 1 if allow charge injection
inj_seed: int32, default 0, random seed for charge injection
injc: 2D float64 numpy array, default None, with the same size as image, charge injection current for each pixel
Tpix: float64, default None, dwell time, injc*Tpix is the number of charge injected to each pixel in a transfer
dark_flag: int32 default 0, set to 1 if the charge injection is from dark current
Tdark: float64, default None, integration time, injc*Tdark will be added to the input image before CTI is added
# III. See example for a quick start
