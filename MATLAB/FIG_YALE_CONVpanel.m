% facial reconstruction template
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

meanface = mean(X,2);
X = X-repmat(meanface,1,size(X,2)); % mean centered data

% proper orthogonal decomposition
[Psi,S,V] = svd(X,'econ');
eigenfaces = Psi(:,1:10);

[m,n] = size(XTrain);
imgtrain = XTrain(:,1:10)+repmat(meanface,1,10);
maxval = max(imgtrain(:));
for i=1:10
    imgtrain(:,i) = imgtrain(:,i)/maxval;
end
imgtrain = reshape(imgtrain,32,10*32);
imagesc(imgtrain), axis image,  colormap(gray), axis off
brighten(0.5)
%printFormattedEPS(gca,gcf,[figpath,'yale_train.eps'],get(gcf,'Position'));

for i=1:10
    eigenfaces(:,i) = eigenfaces(:,i)/max(abs(eigenfaces(:,i)));
end
imagesc(reshape(eigenfaces,32,10*32)), axis image,  colormap(gray), axis off
%brighten(0.5)
%printFormattedEPS(gca,gcf,[figpath,'eigenface.eps'],get(gcf,'Position'));


idx = randperm(length(testIdx),100);
Y = X(:,testIdx(idx)) ;

% target ranks
R = 10:10:300;

Eopt = zeros(100,length(R));
Eqr = zeros(100,length(R));
Eproj = zeros(100,length(R));
Erand = zeros(100*25,length(R));

rng(729);
%% Regression onto random vs. optimal sensors
for k = 1:length(R)
    r = R(k);
    
    % L1 optimized sensor selection
%     zhat = sens_sel_approxnt(Psi(:,1:r),r);
%     [zloc,~] = sens_sel_locr(Psi(:,1:r),r,zhat);
%     sens = find(zloc>.1);
%     assert(length(sens)==r);
            
    % QR sensor selection
    [~,~,pivot] = qr(Psi(:,1:r)*Psi(:,1:r)','vector');
    pivot = pivot(1:r);
    
    ct = 1;
    for ii =1:100
        x = Y(:,ii);
        
        %Projection onto r eigenfaces
        xproj = Psi(:,1:r)*(Psi(:,1:r)'*x);        
        Eproj(ii,k) = norm(xproj-x)/norm(x);
 
        % L1 optimized sensors      
        %xopt = Psi(:,1:r)*(pinv(Psi(sens,1:r))*x(sens));
%         xopt = Psi(:,1:r)*(Psi(sens,1:r)\x(sens));
%         Eopt(ii,k) = norm(xopt-x)/norm(x);    
        
        % QR sensors      
        xqr = Psi(:,1:r)*(pinv(Psi(pivot,1:r))*x(pivot));
        %xqr = Psi(:,1:r)*(Psi(pivot,1:r)\x(pivot));
        Eqr(ii,k) = norm(xqr-x)/norm(x);   
        
        % Random sensors gappy with 2x as many sensors as modes
%         for jj =1:25
%             rsens = randperm(m,r);            
%             %xr = Psi(:,1:r)*(pinv(Psi(rsens,1:r))*x(rsens));
%             xr = Psi(:,1:r)*(Psi(rsens,1:r)\x(rsens));
%             Erand(ct,k) = norm(xr-x)/norm(x);
%             ct = ct+1;
%         end
    end
    
end


save([datpath,'yale_conv'],'R','Eproj','Eopt','Erand','Eqr');


%%

aboxplot({Eproj,Eopt,Eqr,Erand},'labels',num2str(R'));
legend('proj','opt','qr','rand');
set(gca,'yscale','log');

printFormattedEPS(gca,gcf,[figpath,'FIG_YALE_CONV.eps'],[]);