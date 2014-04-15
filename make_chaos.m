function [mixedwave, subwaves] = make_chaos(fspacing, freqlim)
% Makes a fake LFP by mixing randomly AMed and FMed sine waves whose center
% frequencies f are logarithmically spaced.  Center amplitudes are 1/f.  AM
% and FM are independent of each other. Modulation depths are tastefully
% adjusted to yield vaguely realistic-looking components.  A pink noise
% component with a low cutoff equal to the highest center frequency
% provides a unifying sauce to round out the dish.  10 s long series of 2/s
% bursts sampled at 1000 Hz.
%INPUTS
% fspacing: scalar ratio between successive center frequencies.
% freqlim: two-element vector where freqlim(1) is the the center frequency
%   of the lowest sinusoid and freqlim(2) is a limit which will not be
%   exceeded by the center frequency of the highest sinusoid.
%OUTPUTS
% mixedwave: you guessed it
% subwaves: in subwaves X samples format

numsamp = 10001;
Fs = 1000;

% convert everytghin to log units to construct the frequency grid.
logfspacing = log(fspacing);
logfreqlim = log(freqlim);
logf = logfreqlim(1) : logfspacing : logfreqlim(2);
f = exp(logf);

for fidx = 1:length(f)
    lowcut = f(fidx) / 2;
    putburst = double(rand(numsamp,1) < 0.15*f(fidx)/Fs);
    burstdur = 2 * round(3*Fs/f(fidx)) + 1;
    wavespec((fidx-1)*2 + 1) = struct( 'form', 'pink', ...
        'param', exp(-2*pi*lowcut/Fs), ...
        'freq', ones(numsamp,1), ...
        'amp', 1, ...
        'freqmod', [], ...
        'ampmod', [] ); %#ok<AGROW>
    wavespec((fidx-1)*2 + 2) = struct( 'form', 'sin', ...
        'param', 0, ...
        'freq', f(fidx), ...
        'amp', 1/f(fidx) * conv(putburst, hanning(burstdur), 'same'), ...
        'freqmod', [(fidx-1)*2 + 1, 0.03 * fspacing * f(fidx)], ...
        'ampmod', [] ); %#ok<AGROW>
end
wavespec(end+1) = struct( 'form', 'pink', ...
    'param', exp(-2*pi*f(fidx)/Fs), ...
    'freq', 1, ...
    'amp', 0.1 * 1/f(fidx), ...
    'freqmod', [], ...
    'ampmod', [] );
[~, subwaves] = dg_mkwave(wavespec, Fs);
subwaves = subwaves(:, [(1:length(f)) * 2, length(wavespec)])';
mixedwave = sum(subwaves, 1);
