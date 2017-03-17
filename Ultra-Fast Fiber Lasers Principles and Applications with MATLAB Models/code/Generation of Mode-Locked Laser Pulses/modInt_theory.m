function Eout = modInt_theory(Ein,m,fm)
% modInt_theory(Ein,m,fm)
% intensity modulator model
% written by Lam Quoc Huy

% modulator parameters
%	chirp alpha factor: alpha
%	extinction ratio: esilon (dB)
% modulation parameters
%	modulation index: m
%	modulation frequency: fm
global Ts;
N = size(Ein,1);
k = (1:N)'; % no something like fftshift, which means that it is in the time domain.
k = k-N/2;  % move the shape to the center, otherwise only half the pulse will appear. 
% -511
% -510
% ...
% 512
Eout = Ein.*exp(-m/4*(2*pi*fm)^2*(k*Ts).^2);% Eq. (2.12)
                                            % k is the index, k*Ts is time-axit
end
