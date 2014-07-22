function [wave, subwaves] = dg_mkwave(wavespec, Fs)
% Creates a waveform by summing an arbitrary number of components specified
% by <wavespec>.  The components are created one at a time, and previously
% created components can be used as modulation sources for later
% components.  Does not do any sanity checks on specified values, so
% stupid results, crashes, etc. can potentially result.  The number of
% samples created is equal to the length of the vector(s) in
% <wavespec.freq>.
%INPUTS
% wavespec: a vector of struct with the following fields:
%   form - can be 'sin', 'sqr', 'tri', 'table', 'pink'.  The value 'table'
%       indicates that the waveform is specified by interpolating values
%       from a table.  The value 'pink' specifies pink noise.
%   param - for 'sin' this is the starting phase in cycles; for 'sqr'
%       and 'tri', it is the duty cycle; for 'table' it is a vector of
%       sample values for one full cycle. For 'pink' it is the gain to be
%       applied to the previous value when integrating the signal (1 =
%       perfect integrator, 0 = no integration at all, so noise is white).
%       Looked on as a single-pole highpass filtering of perfectly
%       integrated white noise, the half-power frequency is fcutoff = 
%       -1/(2*pi/(Fs*ln(<param>))); conversely, <param> =
%       exp(-2*pi*fcutoff/Fs).
%   freq - a vector of sample-by-sample instantaneous frequency in Hz, or a
%       scalar for constant frequency.  There must be at least one
%       component where 'freq' is a vector, and its length specifies the
%       number of samples to be produced.  All vector values must be the
%       same length.
%   amp - a vector of sample-by-sample instantaneous amplitude, or a scalar
%       for constant amplitude.
%   freqmod - a two-element vector where the first is an index into
%       <wavespec> specifying a component to use as a modulation source,
%       and the second is a modulation gain by which to multiply that
%       component before adding it to <freq> to determine the final
%       frequency. If the modulation gain is zero, or the field is empty,
%       there is no modulation.  The modulation source must have a lower
%       index than the current channel so that the modulation waveform will
%       have a value when the current waveform is computed.
%   ampmod - same as 'freqmod', but for amplitude.
% Fs - sampling frequency in Hz.
%OUTPUTS
% wave: column vector containing the summed waveform.
% subwaves: the separate components, in samples X components format.
%NOTES
% Does not make any attempt to cope with roundoff problems at extremely
% large numbers of samples.
%   The fractional precision of timing of the last two samples is given by
% (eps / (numsamp/(numsamp-1) - 1)), so with eps = 2.22044604925031e-016,
% you would be down to 1 part in 1000 precision when numsamp =
% -1/(1/(eps*1000 + 1) - 1) = 4.503600e+012, which is 142 years at 1 kHz.
% No worries.

%$Rev: 189 $
%$Date: 2014-01-24 18:48:11 -0500 (Fri, 24 Jan 2014) $
%$Author: dgibson $

% (c) 2014 Daniel J. Gibson, all rights reserved.
% This code may not be used or copied without written permission of the
% author.  Please email dgibson@mit.edu.

% Find number of samples.
nsamp = 1;
for cidx = 1:length(wavespec)
    if length(wavespec(cidx).freq) > 1
        if nsamp == 1
            nsamp = length(wavespec(cidx).freq);
        elseif nsamp ~= length(wavespec(cidx).freq)
            error('dg_mkwave:nsamp', ...
                'Component %d has an inconsistent number of samples.', ...
                cidx);
        end
        wavespec(cidx).freq = reshape(wavespec(cidx).freq, [], 1);
        wavespec(cidx).amp = reshape(wavespec(cidx).amp, [], 1);
    end
end

subwaves = zeros(nsamp, length(wavespec));
Ts = 1/Fs;
for cidx = 1:length(wavespec)
    if isempty(wavespec(cidx).freqmod) || wavespec(cidx).freqmod(2) == 0
        modidx = 1;
        modgain = 0;
    else
        modidx = wavespec(cidx).freqmod(1);
        modgain = wavespec(cidx).freqmod(2);
    end
    % <phase> is in cycles.
    if isempty(wavespec(cidx).freq)
        wavespec(cidx).freq = ones(nsamp, 1);
    elseif isscalar(wavespec(cidx).freq)
        wavespec(cidx).freq = wavespec(cidx).freq * ones(nsamp, 1);
    end
    phase = cumsum((wavespec(cidx).freq + ...
        modgain * subwaves(:, modidx)) * Ts);
    if isempty(wavespec(cidx).ampmod) || wavespec(cidx).ampmod(2) == 0
        modidx = 1;
        modgain = 0;
    else
        modidx = wavespec(cidx).ampmod(1);
        modgain = wavespec(cidx).ampmod(2);
    end
    amp = wavespec(cidx).amp + modgain * subwaves(:, modidx);
    switch wavespec(cidx).form
        case 'sin'
            subwaves(:, cidx) = amp .* ...
                sin(2*pi*(wavespec(cidx).param + phase));
        case 'sqr'
            isupstate = mod(phase,1) < wavespec(cidx).param;
            subwaves(isupstate, cidx) = 1;
            subwaves(~isupstate, cidx) = -1;
            subwaves(:, cidx) = amp .* subwaves(:, cidx);
        case 'tri'
            isupstate = mod(phase,1) < wavespec(cidx).param;
            startsamp = [1
                find(isupstate(2:end) & ~isupstate(1:end-1)) + 1];
            endsamp = find(~isupstate(2:end) & isupstate(1:end-1));
            if length(endsamp) < length(startsamp)
                endsamp(end+1) = nsamp; %#ok<AGROW>
            end
            for k = 1:length(startsamp)
                subwaves(startsamp(k):endsamp(k), cidx) = reshape( ...
                    linspace(-1, 1, endsamp(k) - startsamp(k) + 1), ...
                    [], 1 );
                if k < length(startsamp)
                    subwaves(endsamp(k)+1:startsamp(k+1), cidx) = reshape( ...
                        linspace(1, -1, startsamp(k+1) - endsamp(k)), ...
                        [], 1 );
                else
                    subwaves(endsamp(k)+1:nsamp, cidx) = ...
                        reshape( linspace(1, -1, nsamp - endsamp(k)), ...
                        [], 1 );
                end
            end
            subwaves(:, cidx) = amp .* subwaves(:, cidx);
        case 'table'
            phase = mod(phase,1);
            startsamp = [1; find(phase(2:end) < phase(1:end-1)) + 1];
            for k = 1:length(startsamp)
                if k < length(startsamp)
                    endsamp = startsamp(k+1);
                else
                    endsamp = nsamp;
                end
                subwaves(startsamp(k):endsamp, cidx) = interp1(reshape( ...
                    linspace(0, 1, length(wavespec(cidx).param)), [], ...
                    1 ), reshape(wavespec(cidx).param, [], 1 ), ...
                    phase(startsamp(k):endsamp) );
            end
            subwaves(:, cidx) = amp .* subwaves(:, cidx);
        case 'pink'
            noise = randn(nsamp,1);
            if wavespec(cidx).param == 0
                subwaves(:, cidx) = amp .* noise;
            elseif wavespec(cidx).param == 1
                subwaves(:, cidx) = amp .* cumsum(noise);
            else
                subwaves(1, cidx) = noise(1);
                for k = 2:nsamp
                    subwaves(k, cidx) = ...
                        wavespec(cidx).param * subwaves(k-1, cidx) ...
                        + noise(k);
                end
                subwaves(:, cidx) = amp .* subwaves(:, cidx);
            end
        otherwise
            warning( 'dg_mkwave:form', ...
                'Unknown waveform "%s"; setting component to zero.', ...
                wavespec(cidx).form );
            subwaves(:, cidx) = 0;
    end
end

% And finally...
wave = sum(subwaves, 2);
