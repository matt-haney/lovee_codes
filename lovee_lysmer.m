%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% lovee_lysmer.m
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
% Program lovee_lysmer is a Matlab function that computes Love 
% wave phase and group velocities and mode shapes for a given model based 
% on the finite element method of Lysmer (BSSA, 1970). It performs the 
% forward problem at a single frequency.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
% Nn            number of elements in solid part of model
% hv            vector of grid spacings for solid (meters)
% f             frequency (Hz)
% modn          which mode (1=fundamental, 2=first overtone, etc)
% vsv           shear velocity model in solid, a vector (m/s)
% rhov          density model in solid, a vector (kg/m^3)
%
% Output:
% kk            wavenumber for the Love wave at this frequency
% vpk           phase velocity for the Love wave at this frequency
% vgk           group velocity for the Love wave at this frequency
% ev            horizontal displacement eigenfunction (mode shape)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kk, vpk, vgk, ev] = ... 
    lovee_lysmer(Nn,vsv,rhov,f,hv,modn)

% make mu
muv = rhov.*vsv.*vsv;

% make angular frequency
omga = 2*pi*f;

% initialize some local matrices
L1 = sparse(2,2);
L3 = sparse(2,2);
M1 = sparse(2,2);

% initialize the global matrix
Ka1 = sparse(Nn,Nn);
Ka3 = sparse(Nn,Nn);
M = sparse(Nn,Nn);

% for all elements
for ii=1:Nn
   
    % grab grid interval of current element
    h = hv(ii);
    
    % grab material properties of current element
    mu = muv(ii);
    
    % make elemental mass matrix
    M1 = sparse(2,2);
    M1(1,1) = h*rhov(ii)/2;
    M1(2,2) = h*rhov(ii)/2;
    
    % make elemental stiffness matrices
    L1 = sparse(2,2);
    L3 = sparse(2,2);
    
    % some alternate variables from Lysmer
    alph = mu/6;
    bet = mu/6;
    
    % the 4 entries of the 2x2 elemental stiffness matrices of Lysmer
    
    L1(1,1) = 2*alph*h; 
    L3(1,1) = (6*bet/h);
    
    L1(1,2) = alph*h;
    L3(1,2) = -(6*bet/h);
    
    L1(2,1) = L1(1,2);
    L3(2,1) = L3(1,2);
    
    L1(2,2) = L1(1,1);
    L3(2,2) = L3(1,1);
    
    % assemble mass and stiffness matrices from elemental matrices
    if (ii == Nn)
    M((1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)) = ...
        M((1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)) + M1(1:1,1:1);
    Ka1((1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)) = ...
        Ka1((1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)) + L1(1:1,1:1);
    Ka3((1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)) = ...
        Ka3((1*(ii-1)+1):(1*ii),(1*(ii-1)+1):(1*ii)) + L3(1:1,1:1);
    else
    M((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) = ...
        M((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) + M1;
    Ka1((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) = ...
        Ka1((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) + L1;
    Ka3((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) = ...
        Ka3((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) + L3;
    end
    
end

% a lower bound on the Love wave speed
lspd = min(vsv);

% find the eigenvalue closest to the upper-bound eigenvalue
mn = modn;
opts.disp = 0;
[xpp,dpp]=eigs((omga*omga*M)-Ka3,Ka1,mn,(omga/lspd)^2,opts); 

% test for parasitic or artificial numerical modes and skip them
for mm=1:mn
    % test for whether parasitic or not
    if (max(diff(xpp(:,mm))/max(abs(xpp(:,mm)))) > 0.5)    
        pmi(mm) = 1;
    else 
        pmi(mm) = 0;
    end
end

    
% while not enough modes have been found, repeat
wcountr = 0;
while ((sum(pmi == 0) ~= modn) & wcountr < 10*modn)
    
    % solve eigenproblem again for more modes
    mn = modn+sum(pmi);
    opts.disp = 0;
    [xpp,dpp]=eigs((omga*omga*M)-Ka3,Ka1,mn,(omga/lspd)^2,opts); 
        
    % test new modes for whether they are parasitic or not
    for mm=1:(modn+sum(pmi))    
        if (max(diff(xpp(:,mm))/max(abs(xpp(:,mm)))) > 0.5)     
            pmi(mm) = 1;
        else 
            pmi(mm) = 0;
        end
    end  
    wcountr = wcountr + 1;
end

if (wcountr == 10*modn)
    error('Insufficiently dense grid: Too many parasitic modes encountered. Densify grid and re-run.');
else
end

% organize the non parasitic modes
countrm = 0;
dp = zeros(modn,modn);
xp = zeros(Nn,modn);
for mm=1:length(pmi)
   if(pmi(mm) == 0)
       countrm = countrm + 1;
       dp(countrm,countrm) = dpp(mm,mm);
       xp(:,countrm) = xpp(:,mm);
   else
   end
end

% pick the mode of interest    
x = xp(:,modn);
d = dp(modn,modn);

% normalize the eigenfunction
fctr = (1/(transpose(x(1:1:(1*Nn)))*M*x(1:1:(1*Nn))));
ev = x(1:1:(1*Nn))*sqrt(fctr);

% the wavenumber
kk = sqrt(d(1,1));

% the phase velocity 
vpk = omga/kk;
  
% the group velocity 
vgk = ((transpose(x(1:1:(Nn)))*...
      (2*sqrt(d(1,1))*Ka1)*x(1:1:(1*Nn)))/(2*omga))*...
      (1/(transpose(x(1:1:(1*Nn)))*M*x(1:1:(1*Nn))));

% test if it's a guided mode
a = abs(ev);
% depths of solid model
hs = [0 cumsum(hv)];
% index of half depth of model 
[srtv srti] = sort(abs([0 cumsum(hv)]-(sum(hv)/2)));
s3 = min(srti(1:2));
% integral of absolute value of mode over bottom half of model
dum2 = sum(a(s3:Nn)'.*hv(s3:Nn));
% integral of absolute value of mode over top half of model
dum3 = sum(a(1:s3)'.*hv(1:s3));
% index of sensitivity depth
[srtvv srtii] = sort(abs([0 cumsum(hv)]-(modn*.25*(vpk/f))));

% test if it is a guided mode
if (dum3/dum2 < 3)    
    
    % not a guided mode
    vpk = NaN;
    vgk = NaN;
    kk = NaN;
    ev = NaN(Nn,1);
    
elseif ((2*modn*.25*(vpk/f))/(sum(hv)) > 1)    
    
    error('Insufficiently deep grid: Base of model less than twice the sensitivity depth from the top. Extend grid in depth and re-run. ');
        
elseif (sum((vpk/f)./hv(1:min(srtii(1:2))) < 5) > 0)
    
    error('Insufficiently dense grid: Number of elements per wavelength less than 5 above the depth of sensitivity. Densify grid and re-run.');
    
else    
    
    % the mode is acceptable, a guided mode
    
end





