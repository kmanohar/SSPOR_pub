% Code to reproduce Figure 10 L2 reconstructions of Yale B faces
clear; close all; clc

datpath = '../DATA/';
figpath = '../figures/';

load([datpath,'YaleB_32x32.mat']);

%===========================================
[nSmp,nFea] = size(fea);
for i = 1:nSmp
     fea(i,:) = fea(i,:) ./ max(1e-12,norm(fea(i,:)));
end
%===========================================
%Scale the features (pixel values) to [0,1]%
%===========================================
maxValue = max(max(fea));
fea = fea/maxValue;
%===========================================

X = fea';
meanface = mean(X,2);
X = X-repmat(meanface, 1,size(X,2)); % mean centered data

% 64 images of each person
% seed random number generator for predictable sequence
rng(729); 
trainIdx = [];
for i=1:37
   idx = randperm(64,32);
   trainIdx =  [trainIdx, i*idx];
end

Iord = 1:size(X,2);
testIdx = Iord(~ismember(Iord,trainIdx));
XTrain = X(:,trainIdx);

[Psi,S,V] = svd(XTrain,'econ');
[m,n] = size(XTrain);
sing = diag(S);
sing = sing(sing>1e-13);
thresh = optimal_SVHT_coef(m/n,0)*median(sing);

r_opt = length(sing(sing>=thresh));

% select training image
x = X(:,testIdx(1));
print_face(x+meanface,[figpath,'FIG_10_true']);

for r = [50 100 r_opt 300]
    %% Approximation with r eigenfaces
    
    xproj = Psi(:,1:r)*Psi(:,1:r)'*x;
    print_face(xproj+meanface,[figpath,'FIG_10_proj_',num2str(r)]);
    
    
    %% Random reconstruction with r sensors
    
    sensors = randperm(m,r);
    
    mask = zeros(size(x));
    mask(sensors)  = x(sensors)+meanface(sensors);
    print_face(mask,[figpath,'FIG_10_rand_mask_',num2str(r)]);
    
    xls = Psi(:,1:r)*(Psi(sensors,1:r)\x(sensors));
    print_face(xls+meanface,[figpath,'FIG_10_rand_',num2str(r)]);
    
    
    %% QDEIM with r QR sensors
    
    [~,~,pivot] = qr(Psi(:,1:r)','vector');
    sensors = pivot(1:r);
    
    mask = zeros(size(x));
    mask(sensors)  = x(sensors)+meanface(sensors);
    print_face(mask,[figpath,'FIG_10_qr_mask_',num2str(r)]);
    
    xls = Psi(:,1:r)*(Psi(sensors,1:r)\x(sensors));
    print_face(xls+meanface,[figpath,'FIG_10_qr_',num2str(r)]);
    
end