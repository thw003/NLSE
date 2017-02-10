function Eout = hpha_mod(Ein,Vm,Vbias,Vpi,fm,st,ph0)
% phase modulator parameters
% m: phase modulation index
% bias: DC phase shift (rad)
% modulation frequency: fm
% ph0: initial phase
%% written by Nguyen Duc Nhan
global tstep;
global Ts;
global Vw;
% N = size(Ein,1);
% k = (1:N)';
% tstep = Ts*(k-N/2);
mrad = Vm/Vpi*pi;
Norder = 1;
% Vw = Vbias+Vm*cos(2*pi*fm*(tstep-st)+ph0)+0.27*Vm*cos
(2*pi*2*fm*(tstep-st)+ph0-0.6*pi/1)...
% +0.001*Vm*cos(2*pi*3*fm*(tstep-st)+ph0-0.6*pi/1)+0.0001*Vm
*cos(2*pi*4*fm*(tstep-st)+ph0-0.6*pi/1);
Vw = Vbias+Norder*Vm*cos(2*pi*Norder*fm*(tstep-st)+ph0+
0.0*pi)+0.48*Norder*Vm*cos(2*pi*2*Norder*fm*
(tstep-st)+2*ph0-1.4*pi/1)...
+0.003*Norder*Vm*cos(2*pi*3*Norder*fm*(tstep-st)+
3*ph0-0.1*pi/1)+0.0001*Norder*Vm*cos(2*pi*4*fm*
(tstep-st)+4*ph0-0.1*pi/1);
Â�
Eout = Ein.*exp(j*Vw);
function Vout = synth_sig(Vbias,Vac,fm,t,initphase,opt)
% Synth_sig is a function to synthesize an arbitrary waveforms
% to generate the signal driving a phase modulator
% Vbias - DC voltage
% Vac - amplitude of the ac component
% fm - modulation frequency
% t - vector of times
Nharm = 38;
period = 1/fm;
y = 0;
if opt == 1
step = 2;
Vdc = Vbias;
for ii = -Nharm:step:Nharm
if ii == 0
y = y + 0;
else
ai = -2/(pi*ii)^2;
y = y + ai*exp(j*2*pi*ii*(t-initphase)/period);
end
end
Va = Vac*2*y;
elseif opt == 2
step = 1;
Vdc = Vbias + 0;%15/20;
Va = Vac/2;
for ii = -Nharm:step:Nharm
if ii == 0
y = y + 0;
else
ai = -2/(pi*ii)^2;
y = y + ai*exp(j*2*pi*ii*(t-initphase)/period);
end
endAppendix B	
401
Va = Va*2*y;
else
% Vdc= Vbias;
% pp = initphase*2*pi/period;
% ff = cos(2*pi*fm*t+pp);
% Va = Vac*exp(ff);
Vdc = Vbias;
pp
= initphase*2*pi/period;
mm
= 0.02;
Tfwhm= mm*period;
ti = t;
% ff = Vac*cos(2*pi*fm*t+pp);
% ti = rem(t,period) ; % fraction of time within BitPeriod
n = fix(t/period)+1 ; % extract input sequence number
Va = Vac*exp(-log(2)*1/2*(2*ti./Tfwhm).^(2*1));
end
Vout = Vdc + Va;
