clear all, close all, clc

figure

X = imread('cappuccino','jpeg');  
[nx,ny] = size(X);
S = fft2(X); % X = Psi*s
Ssort = sort(abs(S(:)));
thresh = Ssort(ceil(.95*nx*ny));
bigind = abs(S)>thresh; % big coeffs
Strunc = S.*bigind; % keep big coeffs
Xrecon = uint8(ifft2(Strunc));


imshow(X)
axis off
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../figures/FIG_X_compressa');

F = log(abs(fftshift(S))+1);  % put FFT on log-scale
imshow(mat2gray(F),[]);
axis off
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../figures/FIG_X_compressb');

Flow = log(abs(fftshift(Struncate))+1);  % put FFT on log-scale
imshow(mat2gray(Flow),[]);
axis off
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../figures/FIG_X_compressc');

imshow(Xrecon)
axis off
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../figures/FIG_X_compressd');