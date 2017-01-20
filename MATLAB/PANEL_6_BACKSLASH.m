% Code to reproduce Fig 6, backslash reconstruction of mountain and coffee
% images
clear; close all; clc
datpath = '../DATA/';
figpath = '../figures/';
set(0, 'defaultlinemarkersize', 5); 
set(0, 'defaultaxeslinewidth', 0.9);

rng(729); % seed random number generator for reproducibility

data1 = imread([datpath 'mountain.jpg']);
data2 = imread([datpath 'coffeesq.jpg']);

subplot(2,1,1),imshow(data1);
subplot(2,1,2),imshow(data2);

data1 = rgb2gray(imresize(data1,[64 64]));
data2 = rgb2gray(imresize(data2,[64 64]));

n = 64*64;

x1 = reshape(im2double(data1),n,1);
x2 = reshape(im2double(data2),n,1);

% Determine sparsity level in Fourier space
coef1 = dct(x1);
coef2 = dct(x2);

figure(1)
subplot(2,1,1), imshow(reshape(x1,64,64));
subplot(2,1,2), imshow(reshape(x2,64,64));
export_fig([figpath 'panel_backslash_true_a.eps'],'-transparent');

figure(2)
ind1 = find(abs(coef1)>.1);
ind2 = find(abs(coef2)>.1);

subplot(2,1,1), plot(1:n,coef1,'b-',ind1,coef1(ind1),'ro');
set(gca,'xlim',[0 n]);
subplot(2,1,2), plot(1:n,coef2,'b-',ind2,coef2(ind2),'ro');
set(gca,'xlim',[0 n]);

export_fig([figpath 'panel_backslash_true_b.eps'],'-transparent');
%% K-sparse IDCT reconstruction

a1hat = zeros(size(coef1));
a2hat = zeros(size(coef1));

ind1 = find(abs(coef1)>.1);
ind2 = find(abs(coef2)>.1);

a1hat(ind1) = coef1(ind1);
a2hat(ind2) = coef2(ind2);

x1hat = idct(a1hat);
x2hat = idct(a2hat);

figure(1)
subplot(2,1,1), imshow(reshape(x1hat,64,64));
subplot(2,1,2), imshow(reshape(x2hat,64,64));
export_fig([figpath 'panel_backslash_sparse_a.eps'],'-transparent');

subplot(2,1,1), plot(1:n,a1hat,'b-');
set(gca,'xlim',[0 n]);
subplot(2,1,2), plot(1:n,a2hat,'b-');
set(gca,'xlim',[0 n]);

export_fig([figpath 'panel_backslash_sparse_b.eps'],'-transparent');
%% p = K
% Compute universal Fourier basis

p = max(length(ind1),length(ind2));

Isp = speye(n);

% Compute Theta
Theta = zeros(p,n);
ek = zeros(n,1);
sensors = randperm(n,p);
for ii = 1:n
    ek(ii) = 1;
    psi = idct(ek);
    ek(ii) = 0;
    Theta(:,ii) = psi(sensors);
end

a1 = Theta\x1(sensors);
a2 = Theta\x2(sensors);

ind = find(abs(a1)>0);

% Verify same nonzero coefficients for both images
assert(isequal(find(abs(a1)>0),find(abs(a2)>0)))

subplot(2,1,1), plot(1:n,a1,'b-',ind,a1(ind),'ro');
set(gca,'xlim',[1000 1200]);
subplot(2,1,2), plot(1:n,a2,'b-',ind,a2(ind),'ro');
set(gca,'xlim',[1000 1200]);
export_fig([figpath 'panel_backslash_p_K_b.eps'],'-transparent');

x1hat = idct(a1);
x2hat = idct(a2);

subplot(2,1,1), imshow(reshape(x1hat,64,64));
subplot(2,1,2), imshow(reshape(x2hat,64,64));

export_fig([figpath 'panel_backslash_p_K_a.eps'],'-transparent');
%% p = 2000
% Compute universal Fourier basis

p = 2000;

% Compute Theta
Theta = zeros(p,n);
ek = zeros(n,1);
sensors = randperm(n,p);
for ii = 1:n
    ek(ii) = 1;
    psi = idct(ek);
    ek(ii) = 0;
    Theta(:,ii) = psi(sensors);
end

a1 = Theta\x1(sensors);
a2 = Theta\x2(sensors);

ind = find(abs(a1)>0);

% Verify same nonzero coefficients for both images
assert(isequal(find(abs(a1)>0),find(abs(a2)>0)))

subplot(2,1,1), plot(1:n,a1,'b-',ind,a1(ind),'ro');
set(gca,'xlim',[1000 1200]);
subplot(2,1,2), plot(1:n,a2,'b-',ind,a2(ind),'ro');
set(gca,'xlim',[1000 1200]);
export_fig([figpath 'panel_backslash_p_2000_b.eps'],'-transparent');

x1hat = idct(a1);
x2hat = idct(a2);

subplot(2,1,1), imshow(reshape(x1hat,64,64));
subplot(2,1,2), imshow(reshape(x2hat,64,64));

export_fig([figpath 'panel_backslash_p_2000_a.eps'],'-transparent');