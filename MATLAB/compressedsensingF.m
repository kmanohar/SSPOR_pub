function [ a, xhat ] = compressedsensingF( x, sens, n )
%compressedsensingF Standard recovery in fourier domain
%   Detailed explanation goes here

tic;
y = x(sens);

%___THETA___
% NOTE: Avoid calculating Psi (nxn) directly to avoid memory issues.
ek = zeros(n,1);
Theta = zeros(length(sens),n);
for ii = 1:n
    ek(ii) = 1;
    psi = idct(ek);
    Theta(:,ii) = psi(sens);
    ek(ii) = 0;
end

cvx_begin ;
variable a(n);
minimize(norm(a,1)) ;
subject to 
Theta*a  == y;
cvx_end;


toc

xhat = idct(a);
end

