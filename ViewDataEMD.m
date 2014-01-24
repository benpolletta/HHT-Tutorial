clear all;close all;
%% Data
load('2sample_CA1_DG_LFP_10s')
Nt = length(sample_CA1_theta);
fs = 1e3;
dt = 1/fs;
t = [0:Nt-1]*dt;

CA1theta = sample_CA1_theta;
DGLFP = sample_DG_complexLFP;
clear sample_DG_complexLFP sample_CA1_theta
figure;
subplot(211);plot(t,CA1theta);
xlabel('time (sec');title('CA1 LFP');axis tight
subplot(212);plot(t,DGLFP);axis tight
xlabel('time (sec');title('DG LFP');axis tight

%% Gabor spectrograms
L       = 10;
NFFT    = 2^L;
tw      = -NFFT/2+1:NFFT/2;
sigma   = .2;%[sec]
sigSamp = sigma*fs;
w       = sqrt(sqrt(2)/sigSamp)*exp(-pi*tw.*tw/sigSamp/sigSamp);
overlap = NFFT-1;
[SpecCA1 T F]=spectrogram(CA1theta,w,overlap,NFFT,fs);%deriv gaussian windowed spectrogram
figure;imagesc(F,T,abs(SpecCA1));set(gca,'YDir','normal');
ylim([0 150])
set(gca,'YTick',[0:10:150])
grid on
xlabel('time (sec)');ylabel('Hz')
title('Gabor spectrogram CA1 LFP \theta')

[SpecDG F T]=spectrogram(DGLFP,w,overlap,NFFT,fs);%deriv gaussian windowed spectrogram
figure;imagesc(T,F,abs(SpecDG));set(gca,'YDir','normal');
ylim([0 150])
set(gca,'YTick',[0:10:150])
grid on
xlabel('time (sec)');ylabel('Hz')
title('Gabor spectrogram DG LFP')

%% Morlet wavelet spectrograms
sig1 =  struct('val',CA1theta,'period',dt);
cwtS1 = cwtft(sig1,'scales',1./[1:1:200]);
scales = cwtS1.scales;
MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
freq = 1./(scales.*MorletFourierFactor);
figure;imagesc(t,freq,abs(cwtS1.cfs));set(gca,'YDir','normal');
xlabel('time (sec)'); ylabel('Pseudo-frequency');
title('Morlet spectrogram CA1 LFP \theta')
set(gca,'YTick',[0:10:150])
grid on
ylim([0 150])

sig1 =  struct('val',DGLFP,'period',dt);
cwtS1 = cwtft(sig1,'scales',1./[1:1:200]);
scales = cwtS1.scales;
MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
freq = 1./(scales.*MorletFourierFactor);
figure;imagesc(t,freq,abs(cwtS1.cfs));set(gca,'YDir','normal');
xlabel('time (sec)'); ylabel('Pseudo-frequency');
title('Morlet spectrogram DG LFP \theta')
set(gca,'YTick',[0:10:150])
grid on
ylim([0 150])

%% HHT IF diagram

%Do EMD using conventional algorithm
IMFCA1 = emdos(CA1theta,'method','emd');
IMFDG  = emdos(DGLFP,'method','emd');

figure;plot_imf(IMFCA1,t,'CA1')
figure;plot_imf(IMFDG,t,'DG')

%plot CWT magnitude of IMFs
nIMF = 7;
for j = 1:nIMF
    %     subplot(nIMF,1,j)
    %     psd(Hs,IMFCA1(j,:),'Fs',fs)
    sig1 =  struct('val',IMFDG(j,:),'period',dt);
    cwtS1 = cwtft(sig1,'scales',1./[1:1:200]);
    scales = cwtS1.scales;
    MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
    freq = 1./(scales.*MorletFourierFactor);
    figure;imagesc(t,freq,abs(cwtS1.cfs));set(gca,'YDir','normal');
    xlabel('time (sec)'); ylabel('Pseudo-frequency');
    title('Morlet spectrogram DG LFP \theta')
    set(gca,'YTick',[0:10:150])
    grid on
    ylim([0 150])
end
for j = 1:nIMF
    %     subplot(nIMF,1,j)
    %     psd(Hs,IMFCA1(j,:),'Fs',fs)
    sig1 =  struct('val',IMFCA1(j,:),'period',dt);
    cwtS1 = cwtft(sig1,'scales',1./[1:1:200]);
    scales = cwtS1.scales;
    MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
    freq = 1./(scales.*MorletFourierFactor);
    figure;imagesc(t,freq,abs(cwtS1.cfs));set(gca,'YDir','normal');
    xlabel('time (sec)'); ylabel('Pseudo-frequency');
    title('Morlet spectrogram DG LFP \theta')
    set(gca,'YTick',[0:10:150])
    grid on
    ylim([0 150])
end

% %plot multitaper spectra of each IMF
% TW = 5;
% Hs=spectrum.mtm(TW);
% figure;
% for j = 1:nIMF
%     subplot(nIMF,1,j)
%     psd(Hs,IMFDG(j,:),'Fs',fs)
% end

% return
% %Plot instantaneous frequency of each IMF
% xIMFCA1   = hilbert(IMFCA1');
% xIMFDG    = hilbert(IMFDG');
% PhiIMFCA1 = angle(xIMFCA1');
% PhiIMFDG  = angle(xIMFDG);
% MagIMFCA1 = abs(xIMFCA1);figure;plot(t,MagIMFCA1)
% MagIMFDG  = abs(xIMFDG);
% IF_IMFCA1 = diff(PhiIMFCA1);
% IF_IMFDG  = diff(PhiIMFDG);
% figure;plot(t(2:end),IF_IMFCA1(:,1:3))
% figure;plot(t(2:end),IF_IMFDG(:,1:3))
