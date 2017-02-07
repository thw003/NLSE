%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S O L I T O N O P T I C A L C O M M U N I C A T I O N S %
%
%
% global h0
% solitonpairsnoise.m
clear;
close all
clc;
global Ts; % sampling period
global Fcar; % carrier frequency (optical frequency)
global Vw;
global w;
global tstep;
% cla reset
h0 = figure(1);
c_const = 3e8 ;
ng = 1.47;
lamda = 1.550e-6 ;
% =======================================================
% Fiber parameters
Length_SMF = 80;
% Length of a span in meter
nz = Length_SMF/100;
alpha_indB = 0.2*1e-3; % attenuation (dB/km) --> dB/m
D_SMF = 1*17e-6; % GVD (ps/nm.km); if anomalous dispersion(for compensation),D is negative
n2_SMF = 1*2.6e-20; % nonlinear index (m^2/W)
Aeff_SMF = 55e-12; % effective area (m^2)
% Slope Dispersion
S_SMF = 0.06e3; % ps/(nm^2.km) in SI unit

% CALCULATED QUANTITIES of SMF

alpha_loss_SMF = 1*log(10)*alpha_indB/10; % alpha (1/m)
beta2_SMF = -D_SMF*lamda^2/(2*pi*c_const); % beta2 (ps^2/km);
% -----------------------------------------------------
% beta 3 can be calculated from the Slope Dispersion (S) as follows:]
beta3_SMF = (S_SMF - (4*pi*c_const/(lamda^3))*beta2_SMF)/ ...
((2*pi*c_const/(lamda^2))^2);
% -------------------------------------------------------370	
gamma_SMF = 2*pi*n2_SMF/(lamda*Aeff_SMF); % nonlinearity coef (km^-1.W^-1)
% =======================================================
% filter bandwidth
lamda3dB = 4.0e-9;%1.2e-9; % 1.2nm
f3dB = lamda3dB*(1e11/0.8e-9); % Hz
% Amplifier parameters:
GssdB = 21; % (dB)
PoutsatdB = 13.5; % (dBm)
NF = 5; % (dB)
% EDF fiber parameters
Length_EDF = 15; % Length of a span in meter
alpha_indB = 0.5*1e-3; % attenuation (dB/km) --> dB/m
D_EDF = 15e-6; % GVD (ps/nm.km); if anomalous dispersion(for compensation),D is negative
n2_EDF = 1*2.6e-20; % nonlinear index (m^2/W)
Aeff_EDF = 45e-12; % effective area (m^2)
% Slope Dispersion
S_EDF = 0.06e3; % ps/(nm^2.km) in SI unit
% CALCULATED QUANTITIES of EDF
alpha_loss_EDF = 1*log(10)*alpha_indB/10; % alpha (1/m)
beta2_EDF = -D_EDF*lamda^2/(2*pi*c_const) % beta2 (ps^2/km);
% -----------------------------------------------------
% beta 3 can be calculated from the Slope Dispersion (S) as follows:]
beta3_EDF = (S_EDF - (4*pi*c_const/(lamda^3))*beta2_EDF)/ ...
((2*pi*c_const/(lamda^2))^2);
% -----------------------------------------------------
gamma_EDF = 2*pi*n2_EDF/(lamda*Aeff_EDF); % nonlinearity coef (km^-1.W^-1)
% =======================================================
% Average parameters
% Loss
loss = 11.5; % dB
atten = 1/10^(loss/20);
% to test dispersion let TotDisp=-2e-6, SegNum=2.5*Zo/h
% run again with U1 removed and TotDisp=2e-6
Length_total = Length_SMF+Length_EDF;
TotDisp = (D_SMF*Length_SMF+D_EDF*Length_EDF)/Length_total;
% TotDisp = 17*1e-6 ;

ro = 3e-6 ;
n2 = 3.2*1e-20 ;
Beta2 = -(TotDisp*lamda^2)/(2*pi*c_const);
Beta3 = 0 ;
Gamma = (gamma_SMF*Length_SMF+gamma_EDF*Length_EDF)/Length_total ;
% Gamma = (2*pi*n2)/(lamda*pi*ro^2) ;
alphadB = 0;
Alpha=log(10)*alphadB/(10*1000); % = field loss per meter
% =======================================================
% modulator parameters
Vpi = 6; % Volt
Vbias = 0; % Volt
Vm = 1.8; % Volt
m = Vm/Vpi; % modulation index
fm = 1e9; % modulation frequency
p0 = 0*pi/2;
% =======================================================
% Soliton pulse parameters
Norder = 1; % Order of soliton
Tfwhm = 6.0e-12;
To = Tfwhm/1.763; % Pulsewidth of soliton
Tb = 58.5e-12; % Pulses distance
teta = 2*pi/4; % Relative phase difference between 2 solitons
r = 1; % Amplitude ratio
Ci = 0;% Chirping factor of pulse
qo = Tb/2/To; % Space between 2 solitons
qr = Tb/2/To/2 % Space between 2 solitons

Ld = To^2/abs(Beta2)
Zo = pi*Ld/2 ;
Po = Norder^2/(Gamma*Ld);
P0 = 130e-3;
% Pav = 3*Po/(2*qr)
% =======================================================
Ns = 2^13;
Ts = 0.05e-12;
tstep = Ts*(-Ns/2:1:Ns/2-1);
w = fftshift((2*pi/((Ns-1)*Ts))*(-Ns/2:1:Ns/2-1));
% Usol = sqrt(Po).*sech(tstep./To).*exp(-j*teta).*exp ...
% (-j*Ci*tstep.^2./2/To.^2); %
% Usol = sqrt(P0).*(sech(tstep./To-qo)+r.*sech(r.*(tstep./ ...
% To+qo)).*exp(j*teta)).*exp(-j*Ci*tstep.^2./2/To.^2) ; %
% Usol = sqrt(P0).*(sech(tstep./ ...
% To-qo).*exp(+j*teta)+r.*sech(r.*(tstep./To+qo)).*exp ...
% (-j*teta)).*exp(-j*Ci*tstep.^2./2/To.^2) ; %
% Usol = sqrt(Po).*(r*sech(tstep./ ...
% To-qo).*exp(+j*teta)+sech(tstep./To-0).*exp(- ...
% j*teta)+r.*sech(r.*(tstep./To+qo)).*exp(+j*teta)).*exp ...
% (-j*Ci*tstep.^2./2/To.^2) ; %
Usol = sqrt(P0).*(r*sech(tstep./To-1*qo).*exp(+j*teta)+ ...
r*sech(tstep./To-qo/3).*exp(-j*teta/1)...
+sech(tstep./To+qo/3).*exp(+j*teta/1)+r.*sech(r.*(tstep./ ...
To+1*qo)).*exp(-j*teta)).*exp(-j*Ci*tstep.^2./2/To.^2) ;
% Usol = sqrt(P0).*(r*sech(tstep./To-1*qo).*exp(+j*teta)+ ...
% r*sech(tstep./To-qo/2).*exp(-j*teta/1)+ r*sech(tstep./ ...
% To-0).*exp(j*teta/1)...
% +sech(tstep./To+qo/2).*exp(-j*teta/1)+r.*sech(r.*(tstep./ ...
% To+1*qo)).*exp(j*teta)).*exp(-j*Ci*tstep.^2./2/To.^2) ; %	
% Eo = Usol;
plot(tstep,abs(Usol).^2)
grid
U2(1,:) = abs(Usol).^2 ;
Up(1,:) = Usol*atten;
pause(0.1)
% =======================================================
N_pass = 200;
for ii = 2:N_pass
	% h1 = waitbar((ii-1)/N_pass,'The program is running...');
	[Usol,G] = AmpSimpNoise(Usol,GssdB,PoutsatdB,NF);
	[Usol,A3d,z]=hconst(Usol,Length_EDF,5,1*beta2_EDF,1*beta3_EDF,alpha_loss_EDF,1*gamma_EDF,0); % Fiber section #1
	Usol = GaussLPfilt(tstep,Usol,1,f3dB);
	[Usol,A3d,z]=hconst(Usol,Length_SMF,5,1*beta2_SMF,1*beta3_SMF,alpha_loss_SMF,1*gamma_SMF,0); % Fiber section #2
	Usol = pha_mod(Usol(1:Ns),Vm,Vbias,Vpi,fm,p0); % Phase modulator
	Usol = Usol*atten;
	% -----------------------------------------------------
	ind = mod(ii,20);
	if ind==0
		in = round(ii/20) + 1;
		U2(in,:)=abs(Usol).^2 ;
		Up(in,:) = Usol;
	end
% % close(h1);
end
Pav = mean(abs(Usol).^2);
PavdB = 10*log10(Pav)
Up = Up./atten;
% =======================================================
NN = size(Up,1);
figure(1);
plot(tstep,abs(Usol).^2)
colormap('default');
mesh (real(abs(Up).^2),'meshstyle','row','facecolor','none');
% waterfall(ARR2);
% waterfall(abs(Up).^2);
view(17.5,42);
axis tight;
% delete solpair.mat
% save soldoub9db06md Up tstep Po TotDisp Gamma Ld To Ci r teta Tfwhm Tb Vm Vpi NF PoutsatdB GssdB
Pav = 3*max(max(U2))/(2*qr);
Pav1 = mean(U2(NN,:));

figure(2);
PP = angle(Up(1,:));
Pha_tinput = unwrap(PP);
ff = -diff(Pha_tinput)/Ts;
chrate = diff(ff)/Ts;
subplot(221);
plot (tstep*1e12,Pha_tinput);grid;
% plot(tstep(1:Ns-1)*1e12,ff);grid;
% plot(tstep(2:Ns-1)*1e12,chrate);grid;
PP = angle(Up(5,:));
Pha_tinput = unwrap(PP);
ff = -diff(Pha_tinput)/Ts;
chrate = diff(ff)/Ts;
subplot(222);
plot (tstep*1e12,Pha_tinput);grid;
% plot(tstep(1:Ns-1)*1e12,ff);grid;
% plot(tstep(2:Ns-1)*1e12,chrate);grid;
in = round(NN/2);
PP = angle(Up(in,:));
Pha_tinput = unwrap(PP);
ff = -diff(Pha_tinput)/Ts;
chrate = diff(ff)/Ts;
subplot(223);
plot (tstep*1e12,Pha_tinput);grid;
% plot (tstep(1:Ns-1)*1e12,ff);grid;
% plot (tstep(2:Ns-1)*1e12,chrate);grid;
PP = angle(Up(NN,:));
Pha_tinput = unwrap(PP);
ff = -diff(Pha_tinput)/Ts;
chrate = diff(ff)/Ts;
subplot(224);
plot (tstep*1e12,Pha_tinput);grid;
% plot (tstep(1:Ns-1)*1e12,ff);grid;
% plot (tstep(2:Ns-1)*1e12,chrate);grid;
figure(3);
title('The phase evolution of soliton pairs');
subplot(221);
plot(tstep*1e12,U2(1,:));grid;
in = round(NN/4);
subplot(222);
plot(tstep*1e12,U2(in,:));grid;
subplot(223);
in = round(NN/2);
plot(tstep*1e12,U2(in,:));grid;
subplot(224);
plot(tstep*1e12,U2(NN,:));grid;
figure(4);
Kmag = 1;
Nplot = 100;
Uf = Up.';
Eoutfreq = fft(Uf,Ns); %(:,N1)
Eoutfreq1 = fft(Uf,Ns*Kmag); %(:,N1)
Ioutfreq = Eoutfreq1.*conj(Eoutfreq1)/(Ns*Kmag)^2;
ind = (- Nplot/2 : Nplot/2)';
freq = ind/Ts/Ns/Kmag;
ind = mod((ind + Ns*Kmag),Ns*Kmag)+1;
title('The spectrum evolution of pulse');
subplot(221)
% plot(freq,Ioutfreq(ind,1));
plot(freq,10*log10(Ioutfreq(ind,1))+30);
xlabel('Freq (Hz)');
ylabel('P (W)');
subplot(222)
% in = round(N1/4);
% plot(freq,Ioutfreq(ind,5));
plot(freq,10*log10(Ioutfreq(ind,5))+30);
xlabel('Freq (Hz)');
ylabel('P (W)');
subplot(223)
in = round(NN/2);
% plot(freq,Ioutfreq(ind,in));
plot(freq,10*log10(Ioutfreq(ind,in))+30);
xlabel('Freq (Hz)');
ylabel('P (W)');
subplot(224)
% plot(freq,Ioutfreq(ind,NN));
plot(freq,10*log10(Ioutfreq(ind,NN))+30);
xlabel('Freq (Hz)');
ylabel('P (W)');
figure(5);grid;
Pha = Vw*pi./Vpi; % Phase change in time
dff = -diff(Pha)/Ts;
subplot(211);plot(tstep*1e12,Pha);xlabel('Time (ps)');ylabel('Phase (rad)');axis tight;
subplot(212);plot(tstep(1:Ns-1)*1e12,dff*1e-12);xlabel('Time (ps)');ylabel('Chirping (THz)');
