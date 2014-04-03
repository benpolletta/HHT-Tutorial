function [Morlet,freq,t]=IMF_morlet(IMF,sampling_freq,fig_opt)

datalength = size(IMF,2);

t = (1:datalength)/sampling_freq;

nIMF = min(size(IMF,1),20);

for j = 1:nIMF
    
    sig1 =  struct('val',IMF(j,:),'period',1/sampling_freq);
    cwtS1 = cwtft(sig1,'scales',1./[1:1:200]);
    Morlet(:,:,j) = cwtS1.cfs;
    scales = cwtS1.scales;
    MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
    freq = 1./(scales.*MorletFourierFactor);
    
    if fig_opt == 1
        
        figure;imagesc(t,freq,abs(cwtS1.cfs));set(gca,'YDir','normal');
        xlabel('time (sec)'); ylabel('Pseudo-frequency');
        title('Morlet spectrogram DG LFP')
        set(gca,'YTick',[0:10:150])
        grid on
        ylim([0 150])
        
    end
    
end