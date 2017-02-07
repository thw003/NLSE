% solplot.m
% plotting the soliton – soliton pairs
load solpair5
Ts = tstep(2)-tstep(1);
Ns = size(Up,2);
ARR2 = abs(Up).^2;
figure(1);
colormap('default')
mesh (real(abs(Up).^2),'meshstyle','row','facecolor','none');
%waterfall(abs(Up).^2);
view(17.5,42);
axis tight;
NN = size(Up,1);
figure(2);
PP = angle(Up(1,:));
Pha_tinput = unwrap(PP);
ff = -diff(Pha_tinput)/Ts;
chrate = diff(ff)/Ts;
subplot(221);
% plot (tstep*1e12,Pha_tinput);grid;
% plot (tstep(1:Ns-1)*1e12,ff);grid;
plot (tstep(2:Ns-1)*1e12,chrate);grid;
PP = angle(Up(5,:));
Pha_tinput = unwrap(PP);
ff = -diff(Pha_tinput)/Ts;
chrate = diff(ff)/Ts;
subplot(222);
% plot (tstep*1e12,Pha_tinput);grid;
% plot (tstep(1:Ns-1)*1e12,ff);grid;
plot (tstep(2:Ns-1)*1e12,chrate);grid;
in = round(NN/2);
PP = angle(Up(in,:));
Pha_tinput = unwrap(PP);
ff = -diff(Pha_tinput)/Ts;
chrate = diff(ff)/Ts;398	
subplot(223);
% plot (tstep*1e12,Pha_tinput);grid;
% plot (tstep(1:Ns-1)*1e12,ff);grid;
plot (tstep(2:Ns-1)*1e12,chrate);grid;
PP = angle(Up(NN,:));
Pha_tinput = unwrap(PP);
ff = -diff(Pha_tinput)/Ts;
chrate = diff(ff)/Ts;
subplot(224);
% plot (tstep*1e12,Pha_tinput);grid;
% plot (tstep(1:Ns-1)*1e12,ff);grid;
plot (tstep(2:Ns-1)*1e12,chrate);grid;
figure(3);
title('The phase evolution of soliton pairs');
subplot(221);
plot (tstep*1e12,ARR2(1,:));grid;
subplot(222);
plot (tstep*1e12,ARR2(5,:));grid;
subplot(223);
in = round(NN/2);
plot (tstep*1e12,ARR2(in,:));grid;
subplot(224);
plot (tstep*1e12,ARR2(NN,:));grid;
figure(4);
Kmag = 1;
Nplot = 100;
Uf = Up.';
%(:,N1)
Eoutfreq = fft(Uf,Ns);
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
