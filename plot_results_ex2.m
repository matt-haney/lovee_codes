%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% plot_results_ex2.m
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
% Program plot_results_dix_ex2 is a Matlab script to make plots from the 
% output of program lovee_invert when processing the crustal synthetic 
% example. It is to be run immediately after running lovee_invert.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
vpvsr = 1.7321;
Nn_orig = 240;
% make a three layered model
layrth1 = 5; % thickness in elements, layrth1*h = thickness in meters
layrth2 = 10; % thickness in elements, layrth2*h = thickness in meters
layrth3 = 50;
% the true model
vplay1 = 4000; vslay1 = vplay1/vpvsr; 
vplay2 = 3696; vslay2 = vplay2/vpvsr; 
vplay3 = 4500; vslay3 = vplay3/vpvsr; 
vplay4 = 6000; vslay4 = vplay4/vpvsr; 
vsv_true = [vslay1*ones(1,layrth1) vslay2*ones(1,layrth2) ...
            vslay3*ones(1,layrth3) ...
            vslay4*ones(1,(Nn_orig-(layrth1+layrth2+layrth3)))];
hsst = cumsum(ones(1,Nn_orig)*250)-125;

% plot data comparisons
figure
fsize = 16;
plot(fks,U_data,'bo','LineWidth',2,'MarkerSize',6); axis([min(fks) max(fks) 1500 3400]); hold on
plot(fksr_guess,U_guess,'ro','LineWidth',2,'MarkerSize',6);
plot(fks,U,'ko','LineWidth',2,'MarkerSize',6);
set(gca,'Fontsize',fsize,'FontWeight','bold');
ylabel(' Velocity (m/s) '); xlabel(' Frequency (Hz) ');
title(' Data (blue), initial (red), and final update (black) ');

orient landscape
print(gcf,'-dpsc','-r300','-opengl','Data_space_dix_love.ps');
    
% plot model comparisons
figure
fsize = 16;
maxv = 1000; 
plot(vsv_true,0.001*hsst,'b-','LineWidth',4); axis([1500 4000 0 20]); 
axis ij; hold on
plot(vsv_guess,0.001*hss,'r--','LineWidth',4); axis([1500 4000 0 20]); 
axis ij; hold on
plot(vsv_update((nupdat),:),0.001*hss,'k--','LineWidth',4);
set(gca,'Fontsize',fsize,'FontWeight','bold');
ylabel(' Depth (km) '); xlabel(' Shear velocity (m/s) ');
title(' Model (blue), initial (red), and final update (black) ');

orient landscape
print(gcf,'-dpsc','-r300','-opengl','Model_space_dix_love.ps');

% sensitivity kernel of final update

% interpolate onto regular grid
snsmf_vstoti = zeros(length([0:min(h):sum(h)]),Nf);
for ii=1:Nf
snsmf_vstoti(:,ii) = interp1(hs,snsmf_vstot(:,ii),[0:min(h):sum(h)],'linear');
end

figure
fsize = 12;
subplot(1,2,1)
[qq zz] = size(snsmf_vstoti);
smax = round(max(max(snsmf_vstoti/min(h)))*10000)/10;
smin = round(min(min(snsmf_vstoti/min(h)))*10000)/10;
imagesc(fks,hs/1000,1000*snsmf_vstoti(:,1:81)/min(h)); colormap('jet'); caxis([smin smax]); hold on
axis([.1 .9 0 20])
hh = colorbar('EastOutside','FontSize',fsize,...
    'FontWeight','bold','Ytick',[smin:((smax-smin)/5):smax],'Ylim',[smin smax]);
label = sprintf(' Sensitivity (km^{-1}) ');
set(get(hh,'YLabel'),'String',label,'FontSize',fsize,'FontWeight','bold')
set(gca,'Fontsize',fsize,'FontWeight','bold');
xlabel(' Frequency (Hz) '); ylabel(' Depth (km) '); 
title(' V_{S} kernel for fundamental mode ')
%
subplot(1,2,2)
imagesc(fks(104:end),hs/1000,1000*snsmf_vstoti(:,104:end)/min(h)); colormap('jet'); caxis([smin smax]); hold on
axis([fks(104) .9 0 20])
hh = colorbar('EastOutside','FontSize',fsize,...
    'FontWeight','bold','Ytick',[smin:((smax-smin)/5):smax],'Ylim',[smin smax]);
label = sprintf(' Sensitivity (km^{-1}) ');
set(get(hh,'YLabel'),'String',label,'FontSize',fsize,'FontWeight','bold')
set(gca,'Fontsize',fsize,'FontWeight','bold');
xlabel(' Frequency (Hz) '); ylabel(' Depth (km) '); 
title(' V_{S} kernel for first overtone ')

orient landscape
print(gcf,'-dpsc','-r300','-opengl','Kernel_dix_love.ps');

