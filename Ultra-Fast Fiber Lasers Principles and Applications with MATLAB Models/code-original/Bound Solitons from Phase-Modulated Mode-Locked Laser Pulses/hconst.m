function [As_out,A3d,z] = hconst(As_in,L,h,b2,b3,a,g,p)
% Symmetrized Split-Step Fourier Method
% with constant step.
% Input:
%	L = Length of fiber
%	h = step size
%	As_in = Input field in the time domain
%
% Output:
%	As_out = Output field in the time domain
%
% wrtten by N D Nhan - PTIT

global beta2
global beta3
global w
global alpha
global gam
global pmdmode

beta2 = b2;
beta3 = b3;
alpha = a;
gam = g;
pm = p;
pmdmode = 0;

Atemp = As_in;
M = round(L/h);		% number of steps
n = 1;
z = 0;
l = 1;
A3d = [];
for k = 1:M
	Atemp = ssfm(Atemp,h);
	z(1,n+1) = z(1,n)+h;
	n = n+1;
	% A3d(l,:) = At;
	l = l + 1;
end

As_out = Atemp;

% ============Testing
% c = 3e8;
% nin = 1.5;
% vg = c/nin;
% Tr = 50/c/nin;396	
% Nc=size(As_in,1);
% Tw= Nc*1e-13;
% Tm= 1/10e9;
% Th= Tr/1000;
% 
% As_out = As_out.*exp(j*2*pi*(Tm-Th)/Tm);

