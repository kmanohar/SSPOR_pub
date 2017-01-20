% Code to reproduce Fig 12 ocean SST/ENSO singular value spectrum and modes
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
timeavg = mean(Train,2);
Train = bsxfun(@minus,Train,timeavg);

[Psi,S,V] = svd(Train,'econ');
[m,n] = size(Train);
sing = diag(S);
thresh = optimal_SVHT_coef(n/m,0)*median(sing);
r_opt = length(sing(sing>=thresh));

R = [100 200 r_opt 600];

%% print incremental eigenssts to file

for r=R

    display_fig(Psi(:,r),mask,[],[]);
    % Uncomment to save images to file
    %printFormattedEPS(gca,gcf,[figpath,'FIG_eigensst',num2str(r),'.eps'],[]);

end

%% plot singular values
close all;
plot(sing,'.','color',.7*[1 1 1]);
hold on
plot(sing(1:r_opt),'b.');

plot(R, sing(R),'ro');
plot(R, sing(R),'r.');
set(gca,'yscale','log','xlim',[0 length(sing)+1]);
grid on

set(gca,'ylim',[1 10^5]);

% Uncomment to save spectrum to file
%printFormattedEPS(gca,gcf,[figpath,'FIG_ENSO_spectrum_new.eps'],...
%    [100 100 800 350]);
