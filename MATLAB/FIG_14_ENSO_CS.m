% Code to reproduce Fig 14 ocean SST/ENSO Compressed sensing reconstruction
% Data source: NOAA OI SST V2 data provided by the NOAA/OAR/ESRL
% PSD, Boulder, Colorado, USA, from their Web site at http://www.esrl.noaa.gov/psd/
% Author: Krithika Manohar
% Target: January 1997-May 1998 El Nino event

%ncdisp('sst.wkmean.1990-present.nc')
clear; close all; clc
datpath = '../DATA/';
figpath = '../FIGURES/';
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
% display_fig(x+meansst,mask,[],[]);
% printFormattedEPS(gca,gcf,[figpath,'FIG_enso_true.eps'],[]);

% seed random number generator for reproducibility
rng(729);

% target ranks
R = [100 200 r_opt 600];


%% Compressed sensing with r sensors (same random sensors used for random
   % reconstruction in FIG_13_14_ENSO_LS.m 
tic;
% first obtain sensors
sens1 = randperm(m,R(1));
sens2 = randperm(m,R(2));
sens3 = randperm(m,R(3));

% form large Theta matrices all at once
Theta1 = zeros(R(1),n);
Theta2 = zeros(R(2),n);
Theta3 = zeros(R(3),n);
Isp = speye(m);
for i = 1:m
    psi = idct(Isp(:,i));
    Theta1(:,i) = psi(sens1);
    Theta2(:,i) = psi(sens2);
    Theta3(:,i) = psi(sens3);
end
toc
%% r = 100 sensors

r = R(1);
y = x(sens1)+meansst(sens1);

cvx_begin ;
variable a(m);
minimize(norm(a,1)) ;
subject to 
Theta1*a  == y;
cvx_end;

xcs = idct(a);
close all; display_fig(xcs,mask,[],bounds);
printFormattedEPS(gca,gcf,[figpath,'FIG_enso_cs_',num2str(r),'.eps'],[]);

%% r = 200 sensors

r = R(2);
y = x(sens2)+meansst(sens2);

cvx_begin ;
variable a(m);
minimize(norm(a,1)) ;
subject to 
Theta2*a  == y;
cvx_end;

xcs = idct(a);
close all; display_fig(xcs,mask,[],bounds);
printFormattedEPS(gca,gcf,[figpath,'FIG_enso_cs_',num2str(r),'.eps'],[]);

%% r = 302 sensors

r = R(3);
y = x(sens3)+meansst(sens3);

cvx_begin ;
variable a(m);
minimize(norm(a,1)) ;
subject to 
Theta3*a  == y;
cvx_end;

xcs = idct(a);
close all; display_fig(xcs,mask,[],bounds);
printFormattedEPS(gca,gcf,[figpath,'FIG_enso_cs_',num2str(r),'.eps'],[]);
