%Test bed for basic masking signal algorithm development
clear all;close all
%% Set data recording parameters
Fs  = 1e3;%Sample rate
dt = 1/Fs;
Nsamp  = 2^10;%number of samples
Ts  = Nsamp/Fs;
t   = ([1:Nsamp]-1)/Fs;

%% Make data.  Start with basic Gabor wave packet riding on lower frequency sinusoid.

f0 = 50;
f1 = 5;
Tsig = Ts/30;
t0 = Ts/2;
phi0 = 0;
A0 = .1;
A1 = 4;
SigNoise = A0*2;
x0  = A0/Tsig/sqrt(2*pi)*exp( -(t-t0).^2/Tsig^2/2).*...
    cos(2*pi*(f0.*(t-t0) + phi0)) +...
    A1*cos(2*pi*f1*(t-t0));
    
x = x0 + SigNoise*randn(size(t));
figure;plot(t,x,t,x0)
xlim([min(t) max(t)]);


sig1 =  struct('val',x0,'period',dt);
cwtS1 = cwtft(sig1,'scales',1./[1:1:200]);
scales = cwtS1.scales;
MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
freq = 1./(scales.*MorletFourierFactor);
figure;imagesc(t,freq,abs(cwtS1.cfs));set(gca,'YDir','normal');
xlabel('time (sec)'); ylabel('Pseudo-frequency');
title('Morlet spectrogram')
% set(gca,'YTick',[0:10:150])
grid on
ylim([0 100])

%% Decompose with standard EMD
maxIMF = 4;
IMFx0 = emdos(x0,'method','emd');
figure;plot_imf(IMFx0(1:maxIMF,:),t,'CA1')
IMFx = emdos(x,'method','emd');
figure;plot_imf(IMFx(1:maxIMF,:),t,'CA1')

%% Add a masking signal to bias amplitude of higher frequency component
As = 2*A0;
s = cos(2*pi*f0*(t-t0));
y0 = x0 + s;
y = x + s;
figure;plot(t,x0,t,y0)
IMFy0 = emdos(y0,'method','emd');
IMFy = emdos(y,'method','emd');
IMFy0(1,:) = IMFy0(1,:) - s;
IMFy(2,:) = IMFy(2,:) - s;
figure;plot_imf(IMFy0(1:maxIMF,:),t,'CA1')
IMFy = emdos(y,'method','emd');
figure;plot_imf(IMFy(1:maxIMF,:),t,'CA1')
