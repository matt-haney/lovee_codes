
The LOVEE package:

A set of 8 Matlab scripts and 7 functions to forward model and invert 
a collection of Love-wave phase and/or group velocities of any order modes

Below are descriptions for the scripts, the functions, and the examples.

SCRIPTS

make_synthetic_ex2_rayleigh:     writes "velocity_values.txt", "velocity_values_errs.txt",
                         	  "frequency_values.txt", "mode_values.txt", and "vtype_values.txt"
                        	 calls "raylee_lysmer.m"

make_synthetic_ex2_love:     	 writes "velocity_values.txt", "velocity_values_errs.txt",
                         	  "frequency_values.txt", "mode_values.txt", and "vtype_values.txt"
                        	 calls "lovee_lysmer.m"

make_synthetic_ex1:     	 writes "modx_phase_vels.ascii" and "modx_freqs.ascii"
                        	 calls "lovee_lysmer.m"

make_initial_model_dix_ex2: 	 writes "vp_init.txt", "vs_init.txt", "rho_init.txt",
                         	  "vpf.txt", "rhof.txt",
                         	  "grid_values_solid.txt", "grid_values_fluid.txt", 
                         	  and "input_params.txt"
                        	 no calls

make_initial_model_dix_ex1: 	 writes "velocity_values.txt", "velocity_values_errs.txt",
                         	  "frequency_values.txt", "mode_values.txt", "vtype_values.txt", 
                         	  "vp_init.txt", "vs_init.txt", "rho_init.txt",
                         	  "vpf.txt", "rhof.txt",
                         	  "grid_values_solid.txt", "grid_values_fluid.txt", 
                         	  and "input_params.txt"
                        	 no calls

lovee_invert:          		 reads "velocity_values.txt", "velocity_values_errs.txt",
                         	  "frequency_values.txt", "mode_values.txt", "vtype_values.txt", 
                         	  "vp_init.txt", "vs_init.txt", "rho_init.txt",
                         	  "vpf.txt", "rhof.txt",
                         	  "grid_values_solid.txt", "grid_values_fluid.txt", 
                         	  and "input_params.txt"
                        	 no write
                        	 calls "lovee_sensitivity.m", "lovee_lysmer.m", 
                         	  "linvers.m", and "check_nans.m"

plot_results_dix_ex1:   	 no read/write/calls
                        	 to be run immediately after lovee_invert

plot_results_dix_ex2:      	 no read/write/calls
                        	 to be run immediately after lovee_invert


FUNCTIONS

lovee_lysmer:       INPUT
                        Nn      number of nodes in solid part of model
                        hv      vector of grid spacings for solid (meters)
                        f       frequency (Hz)
                        modn    which mode (1=fundamental, 2=first overtone, etc) 
                        vsv     S-wave velocity model in solid, a vector (m/s)
                        rhov    density model in solid, a vector (kg/m^3)

                    OUTPUT
                        kk      wavenumber for the Love wave at this 
                                frequency
                        vpk     phase velocity for the Love wave at 
                                this frequency
                        vgk     group velocity for the Love wave at 
                                this frequency
                        ev      horizontal displacement 
                                eigenfunction (mode shapes)

                    DEPENDENCIES
                        none

raylee_lysmer:      INPUT
                        Nn      number of nodes in solid part of model
                        Nnf     number of nodes in fluid part of model
                        hv      vector of grid spacings for solid (meters)
                        hvf     vector of grid spacings for fluid (meters)
                        f       frequency (Hz)
                        modn    which mode (1=fundamental, 2=first overtone, etc) 
                        vsv     S-wave velocity model in solid, a vector (m/s)
                        vpv     P-wave velocity model in solid, a vector (m/s)
                        rhov    density model in solid, a vector (kg/m^3)
                        vpv     P-wave velocity model in fluid, a vector (m/s)
                        rhov    density model in fluid, a vector (kg/m^3)

                    OUTPUT
                        kk      wavenumber for the Rayleigh wave at this 
                                frequency
                        vpk     phase velocity for the Rayleigh wave at 
                                this frequency
                        vgk     group velocity for the Rayleigh wave at 
                                this frequency
                        ev      vertical and horizontal displacement 
                                eigenfunctions (mode shapes), note these 
                                are scrambled
                    DEPENDENCIES
                        "stoneley_vel.m" - a version is included here although 
                        it isn't used in the examples since it addresses the 
                        case of a water layer

lovee_sensitivity:  INPUT
                        Nn           number of nodes in solid part of model
                        hv           vector of grid spacings for solid (meters)
                        f            frequency (Hz)
                        modn         vector of mode numbers (1=fundamental, 2=first overtone, etc) 
                        vsv          S-wave velocity model in solid, a vector (m/s)
                        rhov         density model in solid, a vector (kg/m^3)
                        vflg         vector of phase or group flag (0=phase, 1=group)

                    OUTPUT
                        U            modeled velocities (group or phase 
                                     depending on vflg) over the entire 
                                     frequency range
                        snsmf_vstotf group or phase velocity sensitivity 
                                     kernel for Vs (again, depending on vflg)
                        snsmf_htotf  group or phase velocity sensitivity 
                                     kernel (again, depending on vflg)
                                     for an interface in the layering 
                                     changing its depth

                    DEPENDENCIES
                        "lovee_lymser.m"

linvers:            INPUT
                        U_data      velocity data to be inverted
                        U           modeled velocity data
                        snsmf_vstot the jacobian or kernel matrix
                        mcmisr      the inverse square root of the model 
                                    covariance matrix
                        dcmisr      the inverse square root of the data 
                                    covariance matrix
                        Nn          number of elements 
                        vsv         the current S-wave wave velocity model
                        vsg         the initial guess for S-wave velocity

                    OUTPUT
                        dvs         the velocity update

                    DEPENDENCIES
                        none

check_nans          INPUT
                        U             modeled velocity data
                        U_data        velocity data to be inverted
                        fks           vector of frequencies
                        modn          vector of mode numbers
                        vflg          vector flags indicating group or phase
                        snsmf_vstotf  Vs sensitivity kernel

                    OUTPUT
                        Ur            modeled velocity data with NaNs removed
                        U_datar       velocity data to be inverted with NaNs removed
                        fksr          vector of frequencies with NaNs removed
                        fksri         vector of original frequency indicies 
                        modnr         vector of mode numbers with NaNs removed
                        vflgr         vector flags indicating group or phase with NaNs removed
                        snsmf_vstotfr Vs sensitivity kernel with NaNs removed

                    DEPENDENCIES
                        none

convsm_1d is an additional utility function which zero-phase convolutional 
smoothing to an input vector


The codes write and read the following text files:

velocity_values.txt:      The measured Love velocities
velocity_values_errs.txt: Error bars on the measurements
frequency_values.txt:     Frequencies at which the measurements are made
mode_values.txt:          Mode number (Fundamental=1, First overtone=2, etc)
vtype_values.txt:         Vector of velocity type, either phase (0) or group (1)
vp_init.txt               Initial P-wave velocity model in solid (for Love-waves, this is not used)
vs_init.txt               Initial S-wave velocity model in solid
rho_init.txt              Initial density model in solid
vpf.txt                   P-wave velocity in fluid layer (for Love-waves, this is not used)
rhof.txt                  Density in fluid layer (for Love-waves, this is not used)
grid_values_solid.txt     Finite element grid in solid layer
grid_values_fluid.txt     Finite element grid in fluid layer (for Love-waves, this is not used)
input_params.txt          See description below

The file "input_params.txt" is automatically generated by the make_initial_model scripts. 
This file contains these quantities:

% flag for fixed poisson's ratio (0=no,1=yes)
% smoothness scale (m)
% a priori model standard deviation factor
% maximum number of updates (iterations)
% number of measurements
% number of elements in solid part of model
% number of elements in fluid part of model
% lower chi squared window
% higher chi squared window

Note that the first and seventh quantities, the ones with poisson's ratio and elements in 
the fluid part, are not used for Love wave inversion. They are kept in the input file so that 
files can be exchanged with the codes for Rayleigh-Scholte wave inversion

The lower and higher chi-squared values define a window of chi-squared when the code can terminate
the iterative perturbational inversion. 

The smoothness scale is a regularization parameter that is typically less than the maximum sensitivity depth, 
usually one-tenth of that depth.

The a prior model standard deviation factor is usually set to 2 and this factor multiplies the average data 
standard deviation to give the model standard deviation.


EXAMPLES

Place all 8 scripts and 7 functions into the same directory. 

For the first inversion example, execute the 
following scripts in order

>> make_synthetic_ex1
>> make_initial_model_dix_ex1
>> lovee_invert
>> plot_results_ex1

This generates Figures 1, 2, and 3 in the paper.


For the second inversion example, execute the 
following scripts in order

>> make_synthetic_ex2_rayleigh
>> make_initial_model_dix_ex2
>> make_synthetic_ex2_love
>> lovee_invert
>> plot_results_ex2

This generates Figures 4, 5, and 6 in the paper.

Matt Haney and Victor Tsai
1 October 2019




