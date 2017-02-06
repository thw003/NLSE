% This code solves the NLS equation with the split-step method
% 	idu/dz - sgn(beta2)/2 dˆ2u/d(tau)ˆ2 + Nˆ2*|u|ˆ2*u = 0
% Written by Govind P. Agrawal in March 2005 for the NLFO book

%---Specify input parameters
clear all; %
% distance = input('Enter fiber length (in units of L_D) = '); %
distance = 1e3;	fprintf('fiber length (in units of L_D) = %g\n', distance);
% 一般浮点数用%f或者%e
% %f显示小数形式
% %e显示科学计数法形式
% %g是两者的综合，会根据数据选择适当的显示方式
% beta2 = input('dispersion: 1 for normal, -1 for anomalous = '); %
beta2 = -2;	fprintf('dispersion: 1 for normal, -1 for anomalous = %g\n', beta2); %
% N = input ('Nonlinear parameter N = '); % Soliton order
N = 1;	fprintf('Nonlinear parameter N = %d\n', N); % Soliton order
mshape = 0;	fprintf('m = 0 for sech, m > 0 for super-Gaussian = %g\n', mshape);
chirp0 = 0;	% input pulse chirp (default value)

%---set simulation parameters9
nt = 2^14; Tmax = 2^9; 	
step_num = round(20*distance*N^2); 	
% No. of z steps
deltaz = distance/step_num; % step size in z
dtau = (2*Tmax)/nt;	% step size in tau

%---tau and omega arrays
tau = (-nt/2:nt/2-1)*dtau; 	
% temporal grid
omega = (pi/Tmax) * [(0:nt/2-1) (-nt/2:-1)]; % frequency grid(already shifted, to be coincidence with ifft(uu))

% Note that the definitions of the FFT and IFFT 
% in Matlab are the opposite to ours. Therefore,
% in the code listed below we use IFFT for the FFT and vice versa.

% very time we do a fft() or ifft(), fftshift() should follow,
% very time we use omega to plot, fftshift() should follow.

%---Input Field profile
if mshape==0 	 % soliton sech shape
	uu = sech(tau).*exp(-0.5i* chirp0* tau.^2);
else % super-Gaussian
	uu = exp(-0.5*(1+1i*chirp0).*tau.^(2*mshape));
end

%---Plot input pulse shape and spectrum
% temp = fftshift(fft(uu)).*(nt*dtau)/sqrt(2*pi); % spectrum
temp = fftshift(ifft(uu)).*(2*Tmax)/sqrt(2*pi); % spectrum
% temp = fftshift(ifft(uu)).*(nt*dtau)/sqrt2*pi); % spectrum
% fft that uses rad/s as frequency unit will lead to scaling in amplitude.
% If we want to display the results in the frequency domain, elimination should be done. As shown in the above line.
% If we want to ifft back, just keep it. As shown in the MAIN loop.

figure; 	
subplot(2,1,1);
	plot (tau, abs(uu).^2, '--k'); hold on;
	axis([-5 5 0 inf]);
	xlabel('Normalized Time');
	ylabel('Normalized Power');
	title('Input and Output Pulse Shape and Spectrum');
subplot(2,1,2);
	plot(fftshift(omega)/(2*pi), abs(temp).^2, '--k'); hold on;
	axis([-.5 .5 0 inf]);
	xlabel('Normalized Frequency');
	ylabel('Spectral Power');
%---Store dispersive phase shifts to speedup code
dispersion = exp(0.5i*beta2*omega.^2*deltaz);	% phase factor
% obtained from Eq. (2.4.2) by replacing the operator ∂/∂T by −iω , and ω is the frequency in the Fourier domain. 
hhz = 1i*N^2*deltaz; 	 % nonlinear phase factor

% ********* [ Beginning of MAIN Loop] ***********
% scheme: 1/2N -> D -> 1/2N -> 1/2N -> D -> ...; 
% 1/2N -> D -> N -> D -> N -> ... -> D -> N -> -1/2N

% first half step nonlinear
temp = uu.*exp(abs(uu).^2.*hhz/2);	% note hhz/2
for n=1:step_num
	 f_temp = ifft(temp).*dispersion;
	 uu = fft(f_temp);
	 temp = uu.*exp(abs(uu).^2.*hhz);	% (B.1)
end
uu = temp.*exp(-abs(uu).^2.*hhz/2); 	 % Final field, step back hhz/2
% temp = fftshift(ifft(uu)).* (nt*dtau)/sqrt(2*pi); % Final spectrum
temp = fftshift(ifft(uu)).* (2*Tmax)/sqrt(2*pi); % Final spectrum
% *************** [ End of MAIN Loop ] **************

%----Plot output pulse shape and spectrum
subplot(2,1,1)
plot(tau, abs(uu).^2, '-k')
subplot(2,1,2)
plot(fftshift(omega)/(2*pi), abs(temp).^2, '-k')

