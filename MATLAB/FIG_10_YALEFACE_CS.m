% Code to reproduce Figure 10 compressed sensing reconstruction of Yale B
% faces
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

r = length(sing(sing>=thresh));

% select training image
x = X(:,testIdx(1))+meanface;

R = [100 r 300 600];

for i=1:4
    sensors = randperm(m,R(i));
    
    mask = zeros(size(x));
    mask(sensors) = x(sensors);
    print_face(mask,[figpath,'FIG_1_cs_mask_',num2str(R(i)),'.eps']);
    
    [~,xcs] = compressedsensingF( x, sensors, m);
    print_face(xcs,[figpath,'FIG_1_cs_',num2str(R(i)),'.eps']);
    
end

