function max_amp = plot_imf_1axis(imf,t,suffixe,max_amp,hF)

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


M = size(imf,1);

if abs(t-1)<=eps
    t = 1:size(imf,2);
end

if isempty(max_amp)
    max_amp=max(max(imf')-min(imf'));
end

if isempty(hF)
    figure;
    linecolor = 'b';
else
    hA = findobj(hF, 'Type', 'axes');
    set(hA,'NextPlot', 'add');
    linecolor = 'r';
end
hold on;
title(suffixe);
xlabel('Time')
ylabel('Number')

for i=1:M
    plot(t,imf(M-i+1,:)+(i-1/2)*max_amp, 'Color', linecolor);
    hold on;
end

set(gca,'YTick',([1:M]-1/2)*max_amp,'YTickLabel',M:-1:1)

end