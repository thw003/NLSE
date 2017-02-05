% This code solves the NLS equation with the split-step method
% 	idu/dz - sgn(beta2)/2 dˆ2u/d(tau)ˆ2 + Nˆ2*|u|ˆ2*u = 0
% Written by Govind P. Agrawal in March 2005 for the NLFO book

%---Specify input parameters
clear all; %
distance = input('Enter fiber length (in units of L_D) ='); %
beta2 = input('dispersion: 1 for normal, -1 for anomalous'); %
N = input ('Nonlinear parameter N = '); % Soliton order
mshape = input ('m = 0 for sech, m > 0 for super-Gaussian = ');
chirp0 = 0;	% input pulse chirp (default value)

%---set simulation parameters
nt = 1024; Tmax = 32; 	
% FFT points and window size
step_num = round(20*distance*N^2); 	
% No. of z steps
deltaz = distance/step_num; % step size in z
dtau = (2*Tmax)/nt;	% step size in tau

%---tau and omega arrays
tau = (-nt/2:nt/2-1)*dtau; 	
% temporal grid
omega = (pi/Tmax) * [(0:nt/2-1) (-nt/2:-1)]; % frequency grid

%---Input Field profile
if mshape==0 	 % soliton sech shape
	uu = sech(tau).*exp(-0.5i* chirp0* tau.^2);
else % super-Gaussian
	uu = exp(-0.5*(1+1i*chirp0).*tau.^(2*mshape));
end

%---Plot input pulse shape and spectrum
temp = fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi); % spectrum
figure; 	subplot(2,1,1);
	plot (tau, abs(uu).^2, '--k'); hold on;
	axis([-5 5 0 inf]);
	xlabel('Normalized Time');
	ylabel('Normalized Power');
	title('Input and Output Pulse Shape and Spectrum');
subplot(2,1,2);
	plot (fftshift(omega)/(2*pi), abs(temp).^2, '--k'); hold on;
	axis([-.5 .5 0 inf]);
	xlabel('Normalized Frequency');
	ylabel('Spectral Power');
%---Store dispersive phase shifts to speedup code
dispersion = exp(0.5i*beta2*omega.^2*deltaz);	% phase factor
hhz = 1i*N^2*deltaz; 	 % nonlinear phase factor

% ********* [ Beginning of MAIN Loop] ***********
% scheme: 1/2N -> D -> 1/2N; first half step nonlinear
temp = uu.*exp(abs(uu).^2.*hhz/2);	% note hhz/2
for n=1:step_num
	 f_temp = ifft(temp).*dispersion;
	 uu = fft(f_temp);
	 temp = uu.*exp(abs(uu).^2.*hhz);
end
uu = temp.*exp(-abs(uu).^2.*hhz/2); 	 % Final field
temp = fftshift(ifft(uu)).* (nt*dtau)/sqrt(2*pi); % Final spectrum
% *************** [ End of MAIN Loop ] **************

%----Plot output pulse shape and spectrum
subplot(2,1,1)
	plot (tau, abs(uu).^2, '-k')
subplot(2,1,2)
	plot(fftshift(omega)/(2*pi), abs(temp).^2, '-k')
