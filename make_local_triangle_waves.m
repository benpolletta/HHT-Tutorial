function [wave] = make_local_triangle_waves(step_size)

fs = 2e3;
dt = 1/fs;
N = 10001;
t = (0:N-1)*dt;

% Low frequency signal is constant amplitude and constant frequency.
amp_low = 1;
freq_low = 15;

wavespec(1) = struct( ...
    'form', 'tri', ...
    'param', 0.7, ...
    'freq', freq_low*ones(1,N), ...
    'amp', amp_low, ...
    'freqmod', [], ...
    'ampmod', [] );

% High frequency signal is constant amplitude, but steps closer to low
% frequency signal in frequency, on timescale step_size.
amp_hi = .1;

freq_hi = zeros(1,N);
no_steps = floor(N/step_size);
freq_steps = linspace(0.05, 0.95, no_steps);
for s = 1:(no_steps-1)
    freq_hi((s-1)*step_size+(1:step_size)) = freq_low/freq_steps(s);
end
freq_hi(((no_steps-1)*step_size+1):end) = freq_low/0.05;
% freq_hi_reflected = [fliplr(freq_hi(1:step_size)) freq_hi fliplr(freq_hi(1:step_size))];
% freq_hi_smoothed = conv(freq_hi_reflected, ones(1,step_size/10), 'same');
% freq_hi = freq_hi_smoothed((step_size+1):(end-step_size));

figure; plot(t, [freq_hi' freq_low*ones(N,1)])

wavespec(2) = struct( ...
    'form', 'tri', ...
    'param', 0.7, ...
    'freq', freq_hi, ...
    'amp', amp_hi, ...
    'freqmod', [], ...
    'ampmod', [] );

% Noise.
wavespec(3) = struct( ...
    'form', 'pink', ...
    'param', exp(-2*pi*2/fs), ...
    'freq', [], ...
    'amp', 0.001, ...
    'freqmod', [], ...
    'ampmod', []);

% Making waves.
[wave, subwaves] = dg_mkwave(wavespec, fs);
figure; plot_imf(subwaves', t, 'subwaves');
figure; plot(t, wave)



