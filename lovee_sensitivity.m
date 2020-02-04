%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% lovee_sensitivity.m
%
% PROGRAMMERS:
% Matt Haney and Victor Tsai
%
% Last revision date:
% 1 October 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is distributed as part of the source-code package 
%                   lovee_inversion_codes 
% that accompanies Haney and Tsai (2019). The package can be downloaded 
% from the Geophysics source-code archive at 
%                   http://software.seg.org/2019/00XX/index.html
% Use of this code is subject to acceptance of the terms and conditions
% that can be found at http://software.seg.org/disclaimer.txt 
% Copyright 2019 by The Society of Exploration Geophysicists (SEG)
% Reference:
% Haney, M. M., Tsai, V. C. (2019) Perturbational and nonperturbational 
% inversion of Love-wave velocities, Geophysics, in press
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program lovee_sensitivity is a Matlab function to compute phase and 
% group velocity sensitivity kernels of Love waves over a 
% range of frequencies.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
% Nn            number of elements in solid part of model
% hv            vector of grid spacings in solid (meters)
% fks           vector of frequencies (Hz)
% modnv         vector of mode numbers (1=fundamental)
% vsv           shear velocity model, a vector (m/s)
% rhov          density model in solid, a vector (kg/m^3)
% vflg          vector of phase or group flag (=0 for phase, =1 for group)
%
% Output:
% U             modeled velocities (group or phase depending on vflg) over 
%               the entire frequency range
% snsmf_vstotf  group or phase velocity sensitivity kernel (again, 
%               depending on vflg)
% snsmf_htotf   sensitivity kernel for phase  or group velocity due to an interface
%               within the layering changing its depth
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uf, snsmf_vstotf, snsmf_htotf] = ...
    lovee_sensitivity(Nn,vsv,rhov,fks,hv,modnv,vflg)

countr = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% augment the frequency vector if group kernels are needed, this triples
% the size of the frequency vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fks_orig = fks;
modnv_orig = modnv;
if (vflg(1) == 1)
    fks = [ fks_orig(1)*.999 fks_orig(1) fks_orig(1)*1.001];
    modnv = [ modnv_orig(1) modnv_orig(1) modnv_orig(1)];
else
    fks = fks_orig(1);
    modnv = modnv_orig(1);
end
for ii=2:length(fks_orig)
    
    if (vflg(ii) == 1)
        fks = [fks fks_orig(ii)*.999 fks_orig(ii) fks_orig(ii)*1.001];
        modnv = [ modnv modnv_orig(ii) modnv_orig(ii) modnv_orig(ii)];        
    else
        fks = [fks fks_orig(ii)];
        modnv = [ modnv modnv_orig(ii)];        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate phase velocities and eigenfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vp = zeros(1,length(fks));
U = zeros(1,length(fks));
xm = zeros(Nn,length(fks));

for f=fks

countr = countr + 1;
[kk, vpk, vgk, x] = lovee_lysmer(Nn,vsv,rhov,f,hv,modnv(countr));

vp(countr) = vpk;
U(countr) = vgk;
xm(:,countr) = x;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct phase sensitivity kernels
%
% these are derivative matrices and are extremely sparse, so the 
% necessary matrix-vector multiplications are hardcoded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snsmf = zeros(Nn,length(fks));
snsmfrho = zeros(Nn,length(fks));

% for all depths (except the bottom element) and frequencies
for ii=1:(Nn-1)
    
    h = hv(ii);

    countr = 0;
for f=fks
    
    countr = countr + 1;

    % density sensitivity
    snsmfrho(ii,countr) = -(vp(countr)/(2*U(countr)))*...
                           (xm(1*ii+0,countr)*(h/2)*xm(1*ii+0,countr)+ ...
                            xm(1*ii+1,countr)*(h/2)*xm(1*ii+1,countr))*...
                            rhov(ii);                                                                 
    
    % mu sensitivity
    % k^2 term
    snsmf(ii,countr) = (1/(2*vp(countr)*U(countr)))*...
                       (xm(1*ii+0,countr)*(h/3)*xm(1*ii+0,countr)+ ...
                        xm(1*ii+1,countr)*(h/3)*xm(1*ii+1,countr)+ ...
                        xm(1*ii+0,countr)*(h/6)*xm(1*ii+1,countr)+ ...
                        xm(1*ii+1,countr)*(h/6)*xm(1*ii+1,countr))*...
                        (rhov(ii)*(vsv(ii)^2));
    % k^0 term                                                                        
    snsmf(ii,countr) = snsmf(ii,countr) + ...
                       (vp(countr)/(2*((2*pi*f)^2)*U(countr)))*...
                       (xm(1*ii+0,countr)*(1/h)*xm(1*ii+0,countr)+ ...
                        xm(1*ii+1,countr)*(1/h)*xm(1*ii+1,countr)+ ...
                        xm(1*ii+0,countr)*(-1/h)*xm(1*ii+1,countr)+ ...
                        xm(1*ii+1,countr)*(-1/h)*xm(1*ii+0,countr))*...
                        (rhov(ii)*(vsv(ii)^2)); 
                    
                    

    % thickness sensitivity
    % omega^2 term                    
    snsmfh(ii,countr) = -(vp(countr)/(2*U(countr)))*...
                           (xm(1*ii+0,countr)*(rhov(ii)/2)*xm(1*ii+0,countr)+ ...
                            xm(1*ii+1,countr)*(rhov(ii)/2)*xm(1*ii+1,countr))*...
                            h;
                       
    % k^2 term                        
    smod = (rhov(ii)*(vsv(ii)^2));                     
    snsmfh(ii,countr) = snsmfh(ii,countr) + ...
                       (1/(2*vp(countr)*U(countr)))*...
                       (xm(1*ii+0,countr)*(smod/3)*xm(1*ii+0,countr)+ ...
                        xm(1*ii+1,countr)*(smod/3)*xm(1*ii+1,countr)+ ...
                        xm(1*ii+0,countr)*(smod/6)*xm(1*ii+1,countr)+ ...
                        xm(1*ii+1,countr)*(smod/6)*xm(1*ii+1,countr))*...
                        h;
                    
    % k^0 term                 
    snsmfh(ii,countr) = snsmfh(ii,countr) + ...
                       (vp(countr)/(2*((2*pi*f)^2)*U(countr)))*...
                       (xm(1*ii+0,countr)*(smod)*(-1/(h^2))*xm(1*ii+0,countr)+ ...
                        xm(1*ii+1,countr)*(smod)*(-1/(h^2))*xm(1*ii+1,countr)+ ...
                        xm(1*ii+0,countr)*(-smod)*(-1/(h^2))*xm(1*ii+1,countr)+ ...
                        xm(1*ii+1,countr)*(-smod)*(-1/(h^2))*xm(1*ii+0,countr))*...
                        h;                                     
        
end

end


% special case for the bottom element
ii=Nn;
h = hv(ii);
countr = 0;
for f=fks
    
    countr = countr + 1;

    % density sensitivity
    snsmfrho(ii,countr) = -(vp(countr)/(2*U(countr)))*...
                           (xm(1*ii-1,countr)*(h/2)*xm(1*ii-1,countr))*...
                            rhov(ii);
        
    % mu sensitivity
    % k^2 term    
    snsmf(ii,countr) = (1/(2*vp(countr)*U(countr)))*...
                       (xm(1*ii+0,countr)*(h/3)*xm(1*ii+0,countr))*...
                        (rhov(ii)*(vsv(ii)^2));    
    % k^0 term                
    snsmf(ii,countr) = snsmf(ii,countr) + ...
                       (vp(countr)/(2*((2*pi*f)^2)*U(countr)))*...
                       (xm(1*ii+0,countr)*(1/h)*xm(1*ii+0,countr))*...
                        (rhov(ii)*(vsv(ii)^2));
                    
    % thickness sensitivity
    % omega^2 term 
    snsmfh(ii,countr) = -(vp(countr)/(2*U(countr)))*...
                           (xm(1*ii-1,countr)*(rhov(ii)/2)*xm(1*ii-1,countr))*...
                            h;
                        
    % k^2 term                          
    smod = (rhov(ii)*(vsv(ii)^2));                     
    snsmfh(ii,countr) = snsmfh(ii,countr) + ...
                       (1/(2*vp(countr)*U(countr)))*...
                       (xm(1*ii+0,countr)*(smod/3)*xm(1*ii+0,countr))*...
                        h;  
                    
    % k^0 term                      
    snsmfh(ii,countr) = snsmfh(ii,countr) + ...
                       (vp(countr)/(2*((2*pi*f)^2)*U(countr)))*...
                       (xm(1*ii+0,countr)*(smod)*(-1/(h^2))*xm(1*ii+0,countr))*...
                        h;                     
                                         
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the shear velocity phase sensitivity kernel for frequencies of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snsmf_vs = 2*snsmf;

% make a vector of the frequencies of interest
if (vflg(1) == 0)
    vfi(1) = 1;
else
    vfi(1) = 2; 
end

for ii=2:length(vflg)
    
    if (vflg(ii) == 1 & vflg(ii-1) == 0)
        
        vfi(ii) = vfi(ii-1) + 2;
        
    elseif (vflg(ii) == 1 & vflg(ii-1) == 1)
        
        vfi(ii) = vfi(ii-1) + 3;
        
    elseif (vflg(ii) == 0 & vflg(ii-1) == 0)
        
        vfi(ii) = vfi(ii-1) + 1;
        
    else
        
        vfi(ii) = vfi(ii-1) + 2;
        
    end
    
end

% compute group kernels and change relative perturbations to absolute
countr = 0;
for f=fks_orig
    
    countr = countr + 1;
    
    if (vflg(countr) == 0)
        
    Uf(countr) = vp(vfi(countr));    
    % change phase kernel for absolute perturbation in the model and data
    % instead of relative perturbation
    snsmf_vstotf(:,countr) = transpose(vp(vfi(countr))*transpose(snsmf_vs(:,vfi(countr)))*diag(vsv.^-1));

    % change phase kernel with respect to thickness to absolute perturbation
    snsmf_htotf(:,countr) = transpose(vp(vfi(countr))*transpose(snsmfh(:,vfi(countr)))*diag(hv.^-1));        
        
        
    else
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the shear velocity group sensitivity kernel,
% obtained using the method of Rodi et al. BSSA (1975)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snsmf_vstot(:,countr) = snsmf_vs(:,vfi(countr)) + ...
                        ((U(vfi(countr))/vp(vfi(countr)))*(2*pi*fks(vfi(countr)))*...
                        (snsmf_vs(:,vfi(countr)+1)-snsmf_vs(:,vfi(countr)-1))/...
                        (1*(fks(vfi(countr)+1)-fks(vfi(countr)-1))*2*pi));

% for absolute perturbations                    
snsmf_vstotf(:,countr) = transpose(vp(vfi(countr))*transpose(snsmf_vstot(:,countr))*diag(vsv.^-1));


% group velocity sensitivity kernel for changes in element thickness
snsmf_htot(:,countr) = snsmfh(:,vfi(countr)) + ...
                        ((U(vfi(countr))/vp(vfi(countr)))*(2*pi*fks(vfi(countr)))*...
                        (snsmfh(:,vfi(countr)+1)-snsmfh(:,vfi(countr)-1))/...
                        (1*(fks(vfi(countr)+1)-fks(vfi(countr)-1))*2*pi));

% for absolute perturbations                     
snsmf_htotf(:,countr) = transpose(U(vfi(countr))*transpose(snsmf_htot(:,countr))*diag(hv.^-1));

% decimate the group velocity for passback
if (isnan(U(vfi(countr)+1)) | isnan(U(vfi(countr)-1)))
    Uf(countr) = NaN;
else
    Uf(countr) = U(vfi(countr));
end
        
        
    end
    
end




