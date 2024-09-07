# MCTI
Monte-Carlo Simulation for Charge Transfer Inefficiency
# Installation
1. Requirement: gcc compiler, astropy, numpy, scipy
2. Change $install_path in install.sh
3. Run: bash install.sh
# Usage
Add "from CTI_modeling_frac import CTI_sim" to your python manuscript.\
Use function CTI_sim to add CTI to your image.\
CTI_sim(image,nx,ny,noverscan,nsp,rho_trap,t,beta,w,c,nmax,oversample,\
        trap_seeds,release_seed,cap_prob,cap_seed,inj_flag,inj_seed,injc,\
        Tpix,dark_flag_Tdark):
1. __Image parameters__\
(1) *image: 2D float64 numpy array*\
    Input image with size (ny,nx)\
   \
(2) *nx: int32*\
    Column number of the image\
   \
(3) *ny: int32*\
    Row number of the image\
   \
(4) *noverscan: int32*\
    Number of parallel overscan pixels\
   
2. __Trap parameters__\
(1) *nsp: int32*\
    Number of trap species\
   \
(2) *rho_trap: float64 numpy array*\
    Trap densities\
   \
(3) *t: float64 numpy array*\
    Trap release timescales (in unit of pix), with the same order as rho_trap\
   
3. __CCD parameters__\
The volume of the electron cloud V is determined by: $V/V_{\mathrm{max}} = \mathrm{max}(((N/w-c)^{\beta},0)$, where $V_{\mathrm{max}}$ is the maximum volume allowed in a pixel, $N$ is the electron count in the pixel.\
(1) *beta: float64*\
    Power exponent of electron cloud volume \
   \
(2) *w: float64*\
    Fullwell capacity\
   \
(3) *c: float64*\
    Notch channel depth\
   
4. __Calculation parameters__\
(1) *nmax: int32*\
Maximum of trap number allowed to be generated in a single pixel for one trap species\
\
(2) *oversample: int32*\
Number of trap will be multiplied by oversample, while each trap will trap or release 1/oversample electrons at a time.\

5. __Random seeds__\
(1) *trap_seed: int32 numpy array*\
Random seeds to generate traps\
\
(2) *release_seed: int32*\
Random seed to determine trap release process\

6. __Capture parameters__\
(1) *cap_prob: float64 numpy array, default None*\
With the same order as rho_trap, when the number of free electrons is not enough for trapping, each free electron will be assigned to a random trap, this array gives the relative capture probability of each trap species. If not given all trap species will have the same trap probabilities.\
\
(2) *cap_seed: int32, default 0*\
Random seed to determine the capture process when free electron number is not enough.\

7. __Charge injection parameters__ (only for inj_flag=1)\
(1) *inj_flag: int32, defalut 0*\
Set to 1 if allow charge injection\
\
(2) *inj_seed: int32, default 0*\
Random seed for charge injection\
\
(3) *injc: 2D float64 numpy array, default None*\
Charge injection current for each pixel\
\
(4) *Tpix: float64, default None*\
Dwell time, injc times Tpix is the number of charge injected to each pixel in a transfer\
\
(5) *dark_flag: int32 default 0*\
Set to 1 if the charge injection is from dark current\
\
(6) *Tdark: float64, default None*\
Integration time, injc times Tdark will be added to the input image before CTI is added\

# See example for a quick start
