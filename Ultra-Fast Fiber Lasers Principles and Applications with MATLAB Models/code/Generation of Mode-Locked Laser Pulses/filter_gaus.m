function Eout = filter_gaus(Ein,f3dB,n)
% filter_gaus(Ein,f3dB,n)
% filter the input signal with the n order gaussian filter
% written by Lam Quoc Huy

% https://en.wikipedia.org/wiki/Gaussian_filter
% T(f) = exp(-log(sqrt(2))*(2/f3dB/Ts/N)^2n*(k.^2n))

global Ts;
N = size(Ein,1);
% the k element in the Ein corresponds to the frequency of
%	(k-N/2)/Ts/N
k = (1:N) -1;
k((N/2+1):N) = k((N/2+1):N) - N;    % something like fftshift, which means it is in the frequency domain
% k = 0, 1, 2, ... , N/2-1(511), -N/2(-512), -N/2+1(-511), ... -2, -1
k = k';
n = 2*n;
% Eout = Ein.*exp(-log(sqrt(2))*(2/f3dB/Ts/N)^n*(k.^n));
% n order gaussian filter
temp = log(sqrt(2))*(2/f3dB/Ts/N)^n;    % log=ln
Eout = Ein.*exp(-temp*(k.^n));	% n order gaussian filter VPI
end
