function [Atout] = GaussLPfilt(t,At,N,f3dB)
% global BitRate

% if nargin == 3
% f3dB = 0.7*BitRate;
% end
S = length(t);
dt = abs(t(2)-t(1));
fin = fftshift(1/dt/1*(-S/2:S/2-1)/(S));
% k = (1:S)-1;
% k(S/2+1:S) = k(S/2+1:S) - S;
% % k = k';
% fin = k/dt/S;

Afin = fft(At);

T = exp(-log(sqrt(2))*(fin/f3dB).^(2*N)); % Ham truyen bo loc
fout = fin;
Afout = T.*Afin;
Atout = ifft(Afout);
% figure(4);
% plot(fin,T);
