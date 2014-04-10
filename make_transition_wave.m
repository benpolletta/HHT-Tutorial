function [wave] = make_transition_wave(transition_index)

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

% Noise.
wavespec(2) = struct( ...
    'form', 'pink', ...
    'param', exp(-2*pi*0.1/fs), ...
    'freq', [], ...
    'amp', .01, ...
    'freqmod', [], ...
    'ampmod', []);

% Making waves.
[~, subwaves] = dg_mkwave(wavespec, fs);
subwaves(:,1) = sin(subwaves(:,1));
subwaves(:,2) = subwaves(:,2) - subwaves(transition_index+1,2);
subwaves(1:transition_index,2) = 0;
wave = sum(subwaves,2);

figure; plot_imf(subwaves', t, 'Subwaves');
figure; plot(t, wave)



