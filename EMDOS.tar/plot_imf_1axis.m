function plot_imf_1axis(imf,t,suffixe)

if nargin<2
    t = linspace(0,1,size(imf,2));
    suffixe = '';
elseif nargin<3
    suffixe = '';
end

M = size(imf,1);

if abs(t-1)<=eps
    t = 1:size(imf,2);
end

max_amp=max(max(imf')-min(imf'));

figure;
hold on;
title(suffixe);
xlabel('Time')
ylabel('Number')

for i=1:M
    plot(t,imf(M-i+1,:)+(i-1/2)*max_amp);
    hold on;
end

set(gca,'YTick',([1:M]-1/2)*max_amp,'YTickLabel',M:-1:1)

end