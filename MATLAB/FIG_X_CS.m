clear all, close all, clc;

N = 4096;
t = linspace(0, 1, N); 
x = cos(2*pi*37*t) + cos(2*pi*420*t) + cos(2*pi*711*t); 

%% Randomly sample signal
p = 256; % num. samples, p=N/16
perm = randperm(N,p);
y = x(perm); % compressed measurement

%% Solve compressed sensing problem
Psi = dct(eye(N, N));  % build Psi
Theta = Psi(perm, :);  % Measure Psi

cvx_begin; % L1-minimization with CVX
    variable s(N); 
    minimize( norm(s,1) ); 
    subject to 
        Theta*s == y';
cvx_end;

xrecon = idct(s); % reconstruct x

% CS using matching pursuit
s = cosamp(Theta,y',10,1.e-10,10); 
xrecon = idct(s);      % reconstruct full signal


%% Plot
time_window = [900 1100]/4096;

figure
subplot(2,2,2)
freq = N/(N)*(0:N);  %create the x-axis of frequencies in Hz
L = 1:floor(N/2);  % only plot the first half of freqs
xt = fft(x); % Fourier transformed signal
PSD = xt.*conj(xt)/N;  % Power spectral density
plot(freq(L),PSD(L),'k', 'LineWidth', 2);
axis([0 1024 0 1200]);

subplot(2,2,1)
plot(t, x, 'k', 'LineWidth', 2);
hold on
plot(perm/N,y,'bo','LineWidth',3);
axis([time_window -3 3]);

subplot(2,2,3)
plot(t, xrecon, 'b', 'LineWidth', 2);
ylim([-2 2])
axis([time_window -3 3]);

subplot(2,2,4)
xtrecon = fft(xrecon,N);  % computes the (fast) discrete fourier transform
PSDrecon = xtrecon.*conj(xtrecon)/N;  % Power spectrum (how much power in each freq)
plot(freq(L),PSDrecon(L),'b', 'LineWidth', 2);
axis([0 1024 0 1200]);

set(gcf,'Position',[100 100 600 400])
set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', '../figures/FIG_X_CS');