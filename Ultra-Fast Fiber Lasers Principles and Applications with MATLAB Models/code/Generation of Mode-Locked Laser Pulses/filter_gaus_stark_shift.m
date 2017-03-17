function Eout = filter_gaus_stark_shift(Ein,f3dB,shift)
% filter_gaus(Ein,f3dB,n)
% filter the input signal with the n order gaussian filter
% written by Lam Quoc Huy

% https://en.wikipedia.org/wiki/Gaussian_filter
% T(f) = exp(-log(sqrt(2))*(2/f3dB/Ts/N)^2n*(k.^2n))

% k = (1:N)-N/2
% T = (exp(-log(sqrt(2))*(2/f3dB/Ts/N)^2*((k-shift).^2))+exp(-log(sqrt(2))*(2/f3dB/Ts/N)^2*((k+shift).^2)))/2 *(2/f3dB*sqrt(log(2)/pi));

global Ts;
N = size(Ein,1);
% the k element in the Ein corresponds to the frequency of
%	(k-N/2)/Ts/N
k = (1:N) -1;   % index in the frequency domain
k((N/2+1):N) = k((N/2+1):N) - N;    % something like fftshift, which means it is in the frequency domain
% k = 0, 1, 2, ... , N/2-1(511), -N/2(-512), -N/2+1(-511), ... -2, -1
k = k';
% Eout = Ein.*exp(-log(sqrt(2))*(2/f3dB/Ts/N)^n*(k.^n));
% n order gaussian filter
temp = log(2)*(2/f3dB)^2;
gauss = exp(-temp*(((k-shift)/(Ts*N)).^2))+exp(-temp*(((k+shift)/(Ts*N)).^2));
gauss = gauss/max(gauss);
Eout = Ein.*gauss;
% Ts*N is the range of the time domain, 1/(Ts*N) is the refinement in the frequency domain
% k/(Ts*N) is the frequency x-axis
end
