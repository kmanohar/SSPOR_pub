% Code to reproduce Figure 8: Cylinder flow example
clear; close all; clc

datpath = '../DATA/';
figpath = '../figures/';

load([datpath,'ALL.mat']);


X = VORTALL;
meanvort = mean(X,2);
X = X-repmat(meanvort, 1,size(X,2)); % mean centered data

% 64 images of each person
% seed random number generator for predictable sequence
rng(729); 

trainIdx = 1:100;

Iord = 1:size(X,2);
testIdx = Iord(~ismember(Iord,trainIdx));
XTrain = X(:,trainIdx);

[Psi,S,V] = svd(XTrain,'econ');
[m,n] = size(XTrain);


semilogy(diag(S)/sum(diag(S)),'-');
set(gca,'XTickLabelRotation',45,'YTick',logspace(-8,0,9));
grid on

% Uncomment to save spectrum
% export_fig([figpath 'ibpm_spect.pdf'],'-transparent');

% Uncomment to save individual modes
% for mode=[1 2 3 4 41 42 43 44 45 46]
%    close; plotCylinderNoSave(Psi(:,mode));
%    caxis([min(Psi(:,mode)) max(Psi(:,mode))]);
%    axis off; box on;
%    set(gcf,'color','none');
%    export_fig([figpath 'ibpm_mode' num2str(mode) '.eps']);
% %    printFormattedEPS(gca,gcf,[figpath 'ibpm_mode' num2str(mode)],...
% %        get(gcf,'Position'));   
% 
% end

Y = X(:,testIdx);

% target ranks
R = 5:5:100;

%Eopt = zeros(100,length(R));
Eqr = zeros(50,length(R));
Eproj = zeros(50,length(R));
Erand = zeros(50*25,length(R));

rng(729);
%% Regression onto random vs. optimized QR sensors
for k = 1:length(R)
    r = R(k);
    
            
    % QR sensor selection
    [~,~,pivot] = qr(Psi(:,1:r)','vector');
    pivot = pivot(1:r);
    
    ct = 1;
    for ii =1:50
        x = Y(:,ii);
        
        %Projection onto r eigenfaces
        xproj = Psi(:,1:r)*(Psi(:,1:r)'*x);        
        Eproj(ii,k) = norm(xproj-x)/norm(x);
 

        % QR sensors      
        xqr = Psi(:,1:r)*(Psi(pivot,1:r)\x(pivot));
        Eqr(ii,k) = norm(xqr-x)/norm(x);   
        
        % Random sensor reconstruction with as many sensors as modes
        for jj =1:25
            rsens = randperm(m,r+1);            
            xr = Psi(:,1:r)*(Psi(rsens,1:r)\x(rsens));
            Erand(ct,k) = norm(xr-x)/norm(x);
            ct = ct+1;
        end
    end
    
end


save([datpath,'ibpm_conv'],'R','Eproj','Erand','Eqr');

%% display convergence plot

close all;
aboxplot({Eproj,Eqr,Erand},'labels',num2str(R'));
legend('POD','QR sensors','random sensors');
set(gca,'yscale','log'); grid on
set(gca,'xTickLabelRotation',45)

%export_fig([figpath,'FIG_IBPM_CONV.eps'],'-transparent');
%printFormattedEPS(gca,gcf,[figpath,'FIG_IBPM_CONV.eps'],get(gcf,'Position'));
