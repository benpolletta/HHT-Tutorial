%Test bed for basic masking signal algorithm development
clear all;close all
%% Set data recording parameters
Fs  = 1e3;%Sample rate
dt = 1/Fs;
Nsamp  = 2^10;%number of samples
Ts  = Nsamp/Fs;
t   = ([1:Nsamp]-1)/Fs;

%% Make data.  Start with basic Gabor wave packet riding on lower frequency sinusoid.

f0 = 50
f1 = 8;
Tsig = Ts/30;
t0 = Ts/2;
phi0 = 0;
A0 = .1;
A1 = 4;
SigNoise = A0*3;
x01 = A0/Tsig/sqrt(2*pi)*exp( -(t-t0).^2/Tsig^2/2).*...
    cos(2*pi*(f0.*(t-t0) + phi0));
x02 = A1*cos(2*pi*f1*(t-t0));
x0  =  x01 + x02;
    
    
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
figure;plot_imf(IMFx0(1:maxIMF,:),t,'x0')
IMFx = emdos(x,'method','emd');
figure;plot_imf(IMFx(1:maxIMF,:),t,'x0 + noise')
MSE(1) = sum((x01-IMFx(3,:)).^2)*dt;
%% Add a masking signal to bias amplitude of higher frequency component
As = 15*A0;
fm = f0
s = As*cos(2*pi*fm*(t-t0)+pi);
y0Plus = x0 + s;
y0Minus = x0 - s;
figure;plot(t,x0,t,y0Plus, t,y0Minus)
IMFy0Plus = emdos(y0Plus,'method','emd');
IMFy0Minus = emdos(y0Minus,'method','emd');
IMFy01 = zeros(size(IMFy0Plus));
IMFy01(1,:) = IMFy0Plus(1,:) - s;
IMFy01(2:end,:) = IMFy0Plus(2:end,:);
[rP cP] = size(IMFy0Plus);
[rM cM] = size(IMFy0Minus);
rmin = min([rP rM]);
IMFy02 = (IMFy0Plus(1:rmin,:) + IMFy0Minus(1:rmin,:))/2;

figure;plot_imf(IMFy01,t,'Masking1')
figure;plot_imf(IMFy02,t,'Masking2')

%% Masking in presence of noise
sig1 =  struct('val',x,'period',dt);
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


%plot multitaper spectrum
TW = 3;
Hs=spectrum.mtm(TW);
figure;psd(Hs,x,'Fs',Fs)

fNyq = Fs/2;
NFFT = 2^nextpow2(Nsamp); % Next power of 2 from length of y
W = 2*abs(fft(x,NFFT))/Nsamp;
f = fNyq*linspace(0,1,NFFT/2+1);
figure(2);semilogy(f, W(1:NFFT/2+1))
hold on
% Low-pass filter below burst frequency 
fp = (f0 + 10)/fNyq;fst = (f0 + 15)/fNyq;Ap = 1;Ast = 20;
h=fdesign.lowpass('Fp,Fst,Ap,Ast',fp,fst,Ap,Ast);
d=design(h,'equiripple'); %Lowpass FIR filter
y=filtfilt(d.Numerator,1,x); %zero-phase filtering
Y = 2*abs(fft(y,NFFT))/Nsamp;
figure(2);semilogy(f,Y(1:NFFT/2+1),'g');
IMFxFilt = emdos(y,'method','emd');
figure;plot_imf(IMFxFilt,t,'Low-pass Filtered')
MSE(2) = sum((x01-IMFxFilt(1,:)).^2)*dt;

yy = y + s;
IMFyMask = emdos(yy,'method','emd');
IMFyMask(1,:) = IMFyMask(1,:) - s;
figure;plot_imf(IMFyMask,t,'low-pass and Masked')
MSE(3) = sum((x01-IMFyMask(1,:)).^2)*dt;

fst1 = (f0 - 30)/fNyq;
fst2 = (f0 + 20)/fNyq;
fp1 = (f0 - 10)/fNyq;
fp2 = (f0 + 10)/fNyq;
Ast1 = 20;Ast2 = 8;Ap = 1;

h = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',fst1,fp1,fp2,fst2,Ast1,Ap,Ast2);
d=design(h,'FIR'); %bandpass FIR filter
yband=filtfilt(d.Numerator,1,y); %zero-phase filtering
MSE(4) = sum((x01-yband).^2)*dt;

Y = 2*abs(fft(yband,NFFT))/Nsamp;
figure(2);semilogy(f,Y(1:NFFT/2+1),'r');
hold off
legend('raw','low-pass','band-pass')
figure;plot(t,yband)
IMFxFilt = emdos(yband,'method','emd');
figure;plot_imf(IMFxFilt,t,'band-pass Filtered')
MSE(5) = sum((x01-IMFxFilt(1,:)).^2)*dt;

%% Now need to estimate phase

