function [Ato] = ssfm(Ati,h)
% Symmetrized Split-Step Fourier Method
% used for single channel.
% Input:
% L = Length of fiber
% h = Variable simulation step
% Ati = Input field in the time domain
%
% Output:
%		Ato = Output field in the time domain
%
%		written by N D Nhan - PTIT

global beta2
global beta3
global w
global alpha
global gam

% c = 3e8;
% nin = 1.5;
% vg = c/nin;
% Tr = 50/c/nin;
% Nc =size(Ati,1);
% Tw = Nc*1e-13;

D = -i/2*beta2.*(i*w).^2+1/6*beta3.*(i*w).^3-alpha/2; % linear operator
N1 = i*gam.*(abs(Ati).^2); % nonlinear operator
N2 = N1;
%Propagation in the first half dispersion region, z to z +h/2
At1 = ifft(exp(h/2.*D).*fft(Ati)); %
% At1 = ifft(exp(1/Tr*h/2.*D).*fft(Ati)); %
% =======================================================
% Iteration for the nonlinear phase shift (2 iterations)
% =======================================================
for m = 1:4
	At1temp = ifft(exp(h/2.*D).*fft(Ati));
	At2 = exp(h/2*(N1+N2)).*At1temp;400	
	At3 = ifft(exp(h/2.*D).*fft(At2));
	%At3 = At3.';
	N2 = i*gam.*(abs(At3).^2);
	% At1temp = ifft(exp(1/Tr*h/2.*D).*fft(Ati));
	% At2 = exp(1/Tr*h/2*(N1+N2)).*At1temp;
	% At3 = ifft(1/Tr*exp(h/2.*D).*fft(At2));
	% %At3 = At3.';
	% N2 = i*gam.*(abs(At3).^2);
end
At4 = exp(h/2.*(N1+N2)).*At1;
% Propagation in the second Dispersion region, z +h/2 to z + h
Ato = ifft(exp(h/2.*D).*fft(At4));
% Ato = Ato.*exp(j*Tr/Tw);

