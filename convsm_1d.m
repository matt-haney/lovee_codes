%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% convsm_1d.m
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
% Program convsm_1d is a Matlab function that smooths a vector using 
% convolution with a comb filter.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input
%
% 1. vector A: either row or column
%
% 2. scalar nx: number of data points on one side of the smoothing window
%    
% Output
%
% 1. vector Asm: a smoothed version of the vector A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Notes
%
% The window of the convolutional smoother is symmetric; for example, for 
% nx=2, with an input array f and an output array F, the convolutional 
% smoother gives       
%                                                           
% F_(i) = (1/5)*[f_(i-2)+f_(i-1)+f_(i)+f_(i+1)+f_(i+2)]    
%                                                           
% nx=2 means that the sum ranges from (i-2) to (i+2).       
%                                                           
% Since the smoothing is done using Fourier transforms, edge effects are 
% treated by padding at the edges before smoothing and then removing the 
% padding. The padding is done by mirroring the data closest to the edges.                                                
%                                                           
% For a nice demo, do the following - smooth a small vector using a 
% smoothing window which is 2*2+1 = 5 samples long: 
%                                                           
% smthd = convsm_1d(1:6,2)                                     
%                                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Asm=convsm_1d(A,nx)

% see if the input vector is a row or a column
% if it is a row, change it to a column, do everything, 
% then change it back to a row at the end for output
ln = size(A);
isflip = 0; % a flag which says if vector was flipped
if (ln(1) == 1)
    A = transpose(A);
    isflip = 1;
else
end

% get the length of the data series
ln = size(A);
Nx = ln(1);

% initialize a padded data vector
Ap = zeros(1,Nx+2*nx);

% pad the data the vector
Ap((nx+1):(nx+Nx)) = A;

% impose reflective boundary conditions for the padded data
Ap(1:nx) = Ap((2*nx+1):-1:(nx+2));
Ap((nx+Nx+1):(Nx+2*nx)) = Ap((nx+Nx-1):-1:Nx);

% initialize a padded filter 
Nxp = Nx+2*nx;
Ktx = zeros(1,Nxp); 

% define a vector with length equal to the padded vector
kx = [1:Nxp];

% build smoothing filter in the frequency domain
if (ceil(Nx/2) == floor(Nx/2))
% even number of samples 
    Ktx = (1/(2*nx+1))*...
          sin(.5*(2*nx+1)*(kx-((Nxp+2)/2)+eps)*((2*pi)/Nxp))./...
          sin(.5*(kx-((Nxp+2)/2)+eps)*((2*pi)/Nxp));
else
% odd number of samples 
    Ktx = (1/(2*nx+1))*...
          sin(.5*(2*nx+1)*(kx-((Nxp+1)/2)+eps)*((2*pi)/Nxp))./...
          sin(.5*(kx-((Nxp+1)/2)+eps)*((2*pi)/Nxp));
end

% apply the smoothing filter 
Ap = real(ifftn(ifftshift(fftshift(fftn(Ap)).*Ktx)));

% unpad the smoothed vector
Asm = zeros(1,Nx);
Asm = Ap((nx+1):(nx+Nx)); 

% if originally handed a row vector, flip to give row back
if (isflip)
else
    Asm = transpose(Asm);
end


