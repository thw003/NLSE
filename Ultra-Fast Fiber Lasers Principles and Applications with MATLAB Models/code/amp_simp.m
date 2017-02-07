function [Eout,gain] = amp_simp(Ein,GssdB,PoutsatdB,NF)
% amp_simp(Ein,GssdB,PoutsatdB,NF)
% simple model of optical amplifier. The model includes the gain
% saturation and the ASE noise
% written by Lam Quoc Huy

% Amplifier parameters:
%	small signal gain: GssdB (dB)
%	output saturation power: PoutsatdB (dBm)
%	noise figure: NF (dB)
%
% Simulation parameters:
%	Gain Tolerance: tol, used as the threshold value to exit the gain calculation loop.
%
% The input is a column vector containing block N samples of the optical signal sampling at the
% rate 1/Ts
% The output is calculated using	if step >0
	step = -step/2;
	end
	err = -er
%	Eout = Ein*sqrt(G) + Enoise
% where: G is the saturated gain
%		 G = Gss*exp(-(G−1)Pin/Psat) (eq1)
%	  Enoise is the complex noise with the noise power
%		 Pase = (10ˆ(NF/10)) * (G−1)hf/2 * BW
%		 BW = 1/Ts

global Ts;	 % sampling period
global Fcar;
h_plan = 6.626e-34;
