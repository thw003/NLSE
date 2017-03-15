function Eout = filter_gaus(Ein,f3dB,n)
% Eq. (4.9)?
% filter_gaus(Ein,f3dB,n)
% filter the input signal with the n order gaussian filter
% written by Lam Quoc Huy

% T(f) = exp(-log(sqrt(2))*(2/f3dB/Ts/N)^2n*(k.^2n))

global Ts;
N = size(Ein,1);
% the k element in the Ein corresponds to the frequency of
%	(k-N/2)/Ts/N
k = (1:N) -1;
k((N/2+1):N) = k((N/2+1):N) - N;
% k = 0, 1, 2, ... , N/2-1(511), -N/2(-512), -N/2+1(-511), ... -2, -1
k = k';
n = 2*n;
% Eout = Ein.*exp(-log(sqrt(2))*(2/f3dB/Ts/N)^n*(k.^n));
% n order gaussian filter
temp = log(sqrt(2))*(2/f3dB/Ts/N)^n;
Eout = Ein.*exp(-temp*(k.^n));	% n order gaussian filter VPI
end
