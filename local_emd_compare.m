function local_emd_compare(sampling_freq,data,smooth_vec,localEMDfunc)
% I presume:
%   smooth_vec is a vector of smoothing widths to try

datalength = length(data);

t=(1:datalength)/sampling_freq;

no_windows = length(smooth_vec);

colors = [linspace(1,0,no_windows)' abs(linspace(1,0,no_windows)'-.5) ...
    linspace(0,1,no_windows)'];

IMF_temp = cell(no_windows+1,1);

size(IMF_temp)

no_IMFs = zeros(no_windows+1,1);

IMF_temp{1} = memd_emd(data);

no_IMFs(1) = size(IMF_temp{1},1);

[max_amp,~] = plot_imf_1axis(IMF_temp{1},t,'Non-Local', [], [], [], data);

for w = 1:no_windows
    
    IMF_temp{w+1} = memd_emd(data,struct('localEMDfunc',localEMDfunc,...
        'localEMDparam',smooth_vec(w)));
    
    plot_imf_1axis(IMF_temp{w+1}, t, ['Smoothing Window = ', ...
        num2str(smooth_vec(w)), ' ms.'], [], [], [], data);
    
    no_IMFs(w+1) = size(IMF_temp{w+1},1);
    
end

max_no_IMFs = max(no_IMFs);

IMF = zeros(max_no_IMFs, datalength, no_windows+1);

for w=1:(no_windows+1)
    
    [r,c] = size(IMF_temp{w});
    
    IMF(1:r,1:c,w) = IMF_temp{w};
    
end

IMF_diff = IMF - repmat(IMF(:,:,1),[1 1 no_windows+1]);
IMF_diff = IMF_diff(:,:,2:end);

handle = figure(); hold on

for w = 1:no_windows
    
    plot_imf_1axis(IMF_diff(:,:,w),t,'',max_amp,handle,colors(w,:));
    
end
title('Different Smoothing Widths Overlaid');

set(handle,'ColorMap',colors)
colorbar('YTick',(1:no_windows)+.5,'YTickLabel',fliplr(smooth_vec))

[IMF1_morlet,f,~] = IMF_morlet(IMF(:,:,1),1000,0);

[r,c] = subplot_size(max_no_IMFs);

for w = 1:no_windows
    
    figure()
   
    [IMFw_morlet,~,~] = IMF_morlet(IMF(:,:,w+1),1000,0);
    
    for i = 1:max_no_IMFs
       
        subplot(r,c,i)
        
        imagesc(t,f,abs(IMF1_morlet(:,:,i))-abs(IMFw_morlet(:,:,i)))
        axis xy
        colorbar
        title(sprintf('\\DeltaIMF %d s=%d', i, smooth_vec(w)))
        
    end
    
end



    

