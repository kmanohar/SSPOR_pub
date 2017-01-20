% Code to reproduce Figure 11 reconstruction convergence plot for YaleB
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

[Phi,S,V] = svd(XTrain,'econ');
[m,n] = size(XTrain);

idx = randperm(length(testIdx),100);
Y = X(:,testIdx(idx)) ;

% target ranks
R = [5 10:10:100];

Eopt = zeros(100,length(R));
Eqr = zeros(100,length(R));
Eproj = zeros(100,length(R));
Erand = zeros(100*25,length(R));

rng(729);
%% Regression onto random vs. optimal sensors
for k = 1:length(R)
    r = R(k);
    
    % Convex optimization sensor selection
    zhat = sens_sel_approxnt(Phi(:,1:r),2*r);
    [zloc,~] = sens_sel_locr(Phi(:,1:r),2*r,zhat);
    sens = find(zloc>.1);
    assert(length(sens)==2*r);
    
    % QR sensor selection
    [~,~,pivot] = qr(Phi(:,1:r)*Phi(:,1:r)','vector');
    pivot = pivot(1:2*r);
    
    zqr = zeros(size(zhat)); zqr(pivot) = 1;
    [zloc,~] = sens_sel_locr(Phi(:,1:r),2*r,zqr);
    pivot = find(zloc>.1);
    assert(length(pivot)==2*r);
    
    ct = 1;
    for ii =1:100
        x = Y(:,ii);
        
        %Projection onto r eigenfaces
        xproj = Phi(:,1:r)*(Phi(:,1:r)'*x);        
        Eproj(ii,k) = norm(xproj-x)/norm(x);
 
        % Convex optimized sensors      
        xopt = Phi(:,1:r)*(pinv(Phi(sens,1:r))*x(sens));        
        Eopt(ii,k) = norm(xopt-x)/norm(x);    
        
        % QR optimized sensors      
        xqr = Phi(:,1:r)*(pinv(Phi(pivot,1:r))*x(pivot));        
        Eqr(ii,k) = norm(xqr-x)/norm(x);   
        
        % Random sensors gappy with 2x as many sensors as modes
        for jj =1:25
            rsens = randperm(m,2*r);            
            xr = Phi(:,1:r)*(pinv(Phi(rsens,1:r))*x(rsens));
            Erand(ct,k) = norm(xr-x)/norm(x);
            ct = ct+1;
        end
    end
    
end


save([datpath,'yale_conv_r2x_small'],'R','Eqr','Eproj','Eopt','Erand');


%%
close all
load([datpath 'yale_conv_r2x_small'])
aboxplot({Eproj*100,Eopt*100,Eqr*100,Erand*100},'labels',num2str(R'));
set(gca,'ylim',[0 100]); grid on
legend('POD','convex','QR sensors','random','Location','sw');

% Uncomment to save plot to file
%printFormattedEPS(gca,gcf,[figpath,'FIG_YALE_CONV_r2x_small.eps'],[]);

