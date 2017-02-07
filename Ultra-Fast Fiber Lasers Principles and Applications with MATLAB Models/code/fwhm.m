function width = fwhm(x)
% return the Full Width at Half Maximum of the pulse x
% written by Lam Quoc Huy

nsize = size(x,1);
% half_peak = max(x)/2;
[peak ind] = max(x);
half_peak = peak/2;		% 3 dB
for iii=1:nsize-1
	if (x(iii)<=half_peak) && (x(iii+1)>half_peak)
		break;
	end
end
for jjj=ind:nsize-1
	if (x(jjj)>=half_peak) && (x(jjj+1)<half_peak)
		break;
	end
end
width = jjj-iii;
end
