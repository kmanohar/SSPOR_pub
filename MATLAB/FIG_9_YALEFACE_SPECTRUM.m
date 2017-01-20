% Code to reproduce Figure 9, singular values of Yale B dataset
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

% Preliminary threshold
sing = sing(sing>=1e-13);
% Optimal truncation threshold
thresh = optimal_SVHT_coef(m/n,0)*median(sing);

r = length(sing(sing>=thresh))

% select training image
x = X(:,testIdx(1));

R = [50 100 r 300];

for i=1:4
    print_face(Psi(:,R(i)),[figpath,'FIG_eigenface',num2str(R(i)),'.eps']);
end
%% plot singular values

close all
semilogy(sing,'.','color',.7*[1 1 1]);
hold on
plot(sing(1:r),'b.');
plot(R, sing(R),'ro');
plot(R, sing(R),'r.');
set(gca,'yscale','log','xlim',[0 length(sing)+1]);
grid on

% Uncomment to save spectrum
%printFormattedEPS(gca,gcf,[figpath,'FIG_spectrum.eps'],[100 100 800 350]);