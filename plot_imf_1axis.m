function [max_amp, handle] = plot_imf_1axis(imf, t, suffixe, max_amp ,hF, ...
    linecolor, rawsignal)

if nargin<2
    t = linspace(0,1,size(imf,2));
    suffixe = '';
elseif nargin<3
    suffixe = '';
end
if nargin < 4
    max_amp = [];
end
if nargin < 5
    hF = [];
end
if nargin < 6 || isempty(linecolor)
    linecolor = [0 0 0];
end
if nargin < 7
    rawsignal = [];
end

M = size(imf,1);

if abs(t-1)<=eps
    t = 1:size(imf,2);
end

if isempty(max_amp)
    max_amp=max(max(imf')-min(imf'));
end

if isempty(hF)
    figure;
else
    hA = findobj(hF, 'Type', 'axes');
    set(hA,'NextPlot', 'add');
end
hold on;
title(suffixe);
xlabel('Time')
ylabel('Number')

for i=1:M
    if ~isempty(rawsignal)
        plot(t,rawsignal+(i-1/2)*max_amp, 'Color', (linecolor + [1 1 1])/2);
    end
    handle = plot(t,imf(M-i+1,:)+(i-1/2)*max_amp, 'Color', linecolor);
    hold on;
end

set(gca,'YTick',([1:M]-1/2)*max_amp,'YTickLabel',M:-1:1)
xlim([min(t) max(t)])

end