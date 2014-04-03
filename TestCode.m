clear all;close all
%% Generate sample data
Fs  = 1e3;
Nsamp  = 1e3;%1e4+1;
Ts  = Fs*Nsamp;
t   = ([1:Nsamp]-1)/Fs;

wavespec.form   = 'tri';%waveform type
wavespec.param  = 1/4;%starting phase of cycle
wavespec.freq   = 10*ones(1,Nsamp);%Hz
wavespec.amp    = 1;
wavespec.freqmod= [];
wavespec.ampmod = [1 1];
Tscale = .1;
t0 = .5;
A = sqrt(1/Tscale)*exp(-pi*(t-t0).^2/Tscale);
[wave1, subwaves]= dg_mkwave(wavespec, Fs);

wavespec.form   = 'pink';%waveform type
fcutoff = 10;
wavespec.param  = 1/2;%starting phase of cycle
wavespec.freq   = 20*ones(1,Nsamp);%Hz
wavespec.amp    = .5;
wavespec.freqmod= [];
wavespec.ampmod = [];
[wave2, subwaves]= dg_mkwave(wavespec, Fs);

wave =  A'.*wave1 + wave2;
figure(1);plot(t,wave)

%% Plot spectrogram
NFFT = 2^nextpow2(Nsamp); % Next power of 2 from length of y
W = fft(wave,NFFT)/Nsamp;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure(2);plot(f, 2*abs(W(1:NFFT/2+1)))

%% Plot CWT of data
sig1 =  struct('val',wave,'period',1/Fs);
maxFrq = 100;
cwtS1 = cwtft(sig1,'scales',1./[1:1:maxFrq]);
scales = cwtS1.scales;
MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
freq = 1./(scales.*MorletFourierFactor);

figure(3);
imagesc(t,freq,abs(cwtS1.cfs));set(gca,'YDir','normal');
xlabel('time (sec)'); ylabel('Pseudo-frequency');
title('Morlet spectrogram CA1 LFP \theta')
% ylim([0 150]);set(gca,'YTick',[0:10:150]);
grid on
colormap gray
colorbar

%% Perform standard EMD
IMF = emdos(wave,'method','emd');
[nIMF nt] = size(IMF);
figure(4);plot_imf(IMF,t,'triangle')

figure(5);plot(t, sum(IMF(4:5,:),1))

%% Identify Significant modes 


%% Now high-pass filter and redo EMD
fc = 50/Fs/2;%cutoff frequency
h=fdesign.lowpass('N,Fc',20,fc);
d=design(h,'FIR'); %Lowpass FIR filter
y=filtfilt(d.Numerator,1,wave); %zero-phase filtering
Y = 2*abs(fft(y,NFFT))/Nsamp;
figure;plot(f,Y(1:NFFT/2+1))
IMF = emdos(wave,'method','emd');
[nIMF nt] = size(IMF);
figure;plot(t, sum(IMF(5,:),1))

