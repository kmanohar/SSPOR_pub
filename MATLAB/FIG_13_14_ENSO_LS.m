% Code to reproduce Figs 13-14 ocean SST/ENSO sensors and L2
% reconstructions
% Data source: NOAA OI SST V2 data provided by the NOAA/OAR/ESRL
% PSD, Boulder, Colorado, USA, from their Web site at http://www.esrl.noaa.gov/psd/
% Author: Krithika Manohar
% Target: January 1997-May 1998 El Nino event

%ncdisp('sst.wkmean.1990-present.nc')
clear; close all; clc
datpath = '../DATA/';
figpath = '../figures/';
[ Lat, Lon, time, mask, sst ] = read_data_enso( [datpath,'sst.wkmean.1990-present.nc'],...
    [datpath,'lsmask.nc'] );


% each element of time array is a new week, in units of days
t0 = datetime(1800,1,1,0,0,0) + days(time(1));
tfin = datetime(1800,1,1,0,0,0) + days(time(end));

[m,n,p] = size(sst);
N = m*n;
X = zeros(N,p);

M = length(mask(mask==1));
Y = zeros(M,length(time));
for i=1:length(time)
   snapshot = reshape(sst(:,:,i),N,1);   
   Y(:,i) = snapshot(mask==1);
end


% train on first 16 years
Iord = 1:length(time);

Itrain = Iord(1:52*16);
Itest = Iord(~ismember(Iord,Itrain));

Train = Y(:,Itrain);
meansst = mean(Train,2);
Train = bsxfun(@minus,Train,meansst);

[Psi,S,V] = svd(Train,'econ');
[m,n] = size(Train);
sing = diag(S);
thresh = optimal_SVHT_coef(n/m,0)*median(sing);
r_opt = length(sing(sing>=thresh));

% select validation snapshot
x = Y(:,Itest(1))-meansst;
bounds = [min(x+meansst) max(x+meansst)];
display_fig(x+meansst,mask,[],[]);
printFormattedEPS(gca,gcf,[figpath,'FIG_enso_true.eps'],[]);

% seed random number generator for reproducibility
rng(729);

% target ranks
R = [100 200 r_opt 600];

%% Reconstruction with random vs. qr sensors
for r = R(1:end-1)
    %% POD approximation with r eigenssts
    
    xproj = Psi(:,1:r)*(Psi(:,1:r)'*x);
    close all; display_fig(xproj+meansst,mask,[],bounds);    
    printFormattedEPS(gca,gcf,[figpath,'FIG_enso_proj',num2str(r),'.eps'],[]);
    
    
    %% Random reconstruction with r sensors
    
    sensors = randperm(m,r);
    
    close all; display_sensors(x+meansst,mask,sensors);    
    printFormattedEPS(gca,gcf,[figpath,'FIG_enso_randmask',num2str(r),'.eps'],[]);
    
    xls = Psi(:,1:r)*(Psi(sensors,1:r)\x(sensors));
    close all; display_fig(xls+meansst,mask,[],bounds);
    printFormattedEPS(gca,gcf,[figpath,'FIG_enso_rand_',num2str(r),'.eps'],[]);
    
    
    %% QDEIM with r QR sensors
    
    [~,~,pivot] = qr(Psi(:,1:r)','vector');
    sensors = pivot(1:r);
    
    close all; display_sensors(x+meansst,mask,sensors);    
    printFormattedEPS(gca,gcf,[figpath,'FIG_enso_opt_mask',num2str(r),'.eps'],[]);
    
    xls = Psi(:,1:r)*(Psi(sensors,1:r)\x(sensors));
    close all; display_fig(xls+meansst,mask,[],bounds);
    printFormattedEPS(gca,gcf,[figpath,'FIG_enso_opt_',num2str(r),'.eps'],[]);
    
end