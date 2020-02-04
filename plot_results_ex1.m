%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% plot_results_ex1.m
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
% Program plot_results_ex1 is a Matlab script to make plots from the 
% output of program lovee_invert when processing the MODX synthetic 
% example. It is to be run immediately after running lovee_invert.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% true model
vsv_true = [ 194 194 270 270 367 367 485 485 603 603 740 740];
hsst = [0 1.9 2.1 4.2 4.4 6.7 6.9 9.5 9.7 12.7 12.9 100];

% plot data comparisons
figure
fsize = 16;
plot(fks,U_data,'bo','LineWidth',2,'MarkerSize',6); axis([min(fks) max(fks) 100 700]); hold on
plot(fks,U_guess,'ro','LineWidth',2,'MarkerSize',6);
plot(fks,U,'ko','LineWidth',2,'MarkerSize',6);
set(gca,'Fontsize',fsize,'FontWeight','bold');
ylabel(' Velocity (m/s) '); xlabel(' Frequency (Hz) ');
title(' MODX data (blue), initial (red), and final update (black) ');

orient landscape
print(gcf,'-dpsc','-r300','-opengl','Data_space_modx.ps');
    
% plot model comparisons
figure
fsize = 16;
maxv = 1000; 
plot(vsv_true,hsst,'b-','LineWidth',4); axis([0 maxv 0 max(hss/8)]); 
axis ij; hold on
plot(vsv_guess,hss,'r--','LineWidth',4); axis([0 maxv 0 max(hss/8)]); 
axis ij; hold on
plot(vsv_update((nupdat),:),hss,'k--','LineWidth',4);
set(gca,'Fontsize',fsize,'FontWeight','bold');
ylabel(' Depth (m) '); xlabel(' Shear velocity (m/s) ');
title(' MODX model (blue), initial (red), and final update (black) ');

orient landscape
print(gcf,'-dpsc','-r300','-opengl','Model_space_modx.ps');

% sensitivity kernel of final update

% interpolate onto regular grid
snsmf_vstoti = zeros(length([0:min(h):sum(h)]),Nf);
for ii=1:Nf
snsmf_vstoti(:,ii) = interp1(hs,snsmf_vstot(:,ii),[0:min(h):sum(h)],'linear');
end

figure
fsize = 16;
[qq zz] = size(snsmf_vstoti);
smax = round(max(max(snsmf_vstoti/min(h)))*10)/10;
smin = round(min(min(snsmf_vstoti/min(h)))*10)/10;
imagesc(fks,hs,snsmf_vstoti/min(h)); colormap('jet'); caxis([smin smax]); hold on
axis([min(fks) max(fks) 0 100])
hh = colorbar('EastOutside','FontSize',fsize,...
    'FontWeight','bold','Ytick',[smin:((smax-smin)/4):smax],'Ylim',[smin smax]);
label = sprintf(' Sensitivity (m^{-1}) ');
set(get(hh,'YLabel'),'String',label,'FontSize',fsize,'FontWeight','bold')
set(gca,'Fontsize',fsize,'FontWeight','bold');
xlabel(' Frequency (Hz) '); ylabel(' Depth (m) '); 
title(' V_{S} kernel for fundamental mode ')

orient landscape
print(gcf,'-dpsc','-r300','-opengl','Kernel_modx.ps');

