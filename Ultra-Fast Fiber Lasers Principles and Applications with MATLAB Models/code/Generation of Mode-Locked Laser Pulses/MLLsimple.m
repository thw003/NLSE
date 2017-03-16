% MLLsimple.m
% simple model of mode locked laser
% written by Saiyu Luo
clear;
clc;

global Ts;			% sampling period
global Fcar;		% carrier frequency (optical frequency)
c_const = 3e8;		% speed of light

lamda = 1550e-9;	% m
Fcar = c_const/lamda;

Ts = 0.1e-12;		% 0.1 ps
N = 1024;			% number of samples in a block. Tblk = N * Ts = 102.4 ps

% Amplifier parameters:
  GssdB = 20;		% (dB)
  PoutsatdB = 10;	% (dBm)
  % NF = 8;			% (dB)

% filter bandwidth
  lamda3dB = 1e-9;	% m
  f3dB = lamda3dB*(1e11/0.8e-9);

% modulator parameters
  alpha = -0.07;
  epsilon = 40;		% (dB) extinction ratio
  
% modulation parameters
  m = 0.5;			% modulation index
  fm = 10e9;		% modulation frequency

% Loss
  loss = 10;			% dB
  atten = 1/10^(loss/20);

% generate an initial block of signal Ein
% Generate an N-by-1 matrix of complex white Gaussian noise having power -40 dBW. 
Ein = wgn(N,1,-40,'complex');

Eout = Ein;
Eo = Ein;
% N_pass = 500;
N_pass = 50;
for ii = 1:N_pass
    fprintf('----------------------------\n', ii);
    fprintf('pass %d begin\n', ii);
	[Eo,G] = AmpSimpNonoise(Eo,GssdB,PoutsatdB); % no noise
	Eo = fft(Eo);
	Eo = filter_gaus(Eo,f3dB,1);
	Eo = ifft(Eo);
	% Eo = modInt(Eo(1:N),alpha,epsilon,m,fm,0.5);
	Eo = modInt_theory(Eo(1:N),m,fm);
	Eo = Eo*atten;
	if mod(ii,N_pass/50)==0
		Eout = [Eout, Eo];
	end
    fprintf('pass %d end\n', ii);
    fprintf('----------------------------\n', ii);
end
Eout = Eout/atten;
close all

% -------------- Display the results ---------
% mesh (abs(Eout'),'edgecolor','black','meshstyle','row', 'facecolor','none');
Iout = Eout.*conj(Eout);
mesh (Iout','meshstyle','row','facecolor','none');
axis tight;
% set(gca,'XTick',tt_mark);
% set(gca,'XTickLabel',tt_tick);
% set(gca,'XDir','reverse');
xlabel('T (0.1ps)');
% set(gca,'YTick',yy_mark);
% set(gca,'YTickLabel',yy_tick);
ylabel('Pass number');
zlabel('intensity (W)');

% N1 = size(Eout,2);
N1 = 1;
dPhi = angle(Eout(2:N,N1)) - angle(Eout(1:N-1,N1));
figure (2);
plot(fftshift(dPhi));

% return the Full Width at Half Maximum of the pulse x
Tp = fwhm(Iout(:,N1))*Ts;
pulse_alpha = 2*log(2)/(Tp^2);
pulse_beta = (dPhi(N/2+100) - dPhi(N/2-100))/200/Ts/Ts;
%  chirp = pulse_beta/pulse_alpha

Kmag = 8;
Nplot = 100;
Eoutfreq = fft(Eout(:,N1),N*Kmag);
Ioutfreq = Eoutfreq.*conj(Eoutfreq)/(N*Kmag)^2;

figure(3);
ind = (- Nplot/2 : Nplot/2)';
freq = ind/Ts/N/Kmag;
ind = mod((ind + N*Kmag),N*Kmag)+1;
plot(freq,Ioutfreq(ind));

n=1;
n = 2*n;

Tfil = exp(-log(2)*(2/f3dB*freq).^n);	% n order gaussian filter VPI
Tfil = Tfil *max(Ioutfreq(ind));
hold on
plot(freq,Tfil,'r');

% plot the gaussian fit curve
% gaussFit(Iout(:,N1));

%  pulseBW = fwhm(Ioutfreq(ind))/Ts/N/Kmag
%  Tp = fwhm(Iout(:,N1))*Ts
%  TBP = pulseBW*Tp
