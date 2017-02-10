function Eout = pha_mod(Ein,Vm,Vbias,Vpi,fm,ph0)
% phase modulator parameters
% m: phase modulation index
% bias: DC phase shift (rad)
% modulation frequency: fm
% ph0: initial phase
% % written by Nguyen Duc Nhan

global tstep;
global Ts;
global Vw;

N = size(Ein,1);
% k = (1:N)';
% tstep = Ts*(k-N/2);

T0 = 1/fm;
modeph = 0;
if modeph == 0
	Vw = Vbias+Vm*cos(2*pi*fm*tstep+ph0);
elseif modeph == 1
	Vw = real(synth_sig(Vbias,Vm,fm,tstep,ph0*T0/(2*pi),2));
else
	Vw = real(synth_sig(Vbias,Vm,fm,tstep,ph0*T0/(2*pi),1));
end
% delta_phi = pi/4*(2- bias*2 - ext*Vm);
Phi = Vw*pi/Vpi;
Eout = Ein.*exp(j*Phi);
% Eout = Ein.*Vm;
end
