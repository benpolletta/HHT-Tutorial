function [A, F] = dg_AMFMdecomp(imf, method, verboseflag)
%[A, F] = dg_AMFMdecomp(imf)  Norden Huang's empirical AM-FM decomposition
% as described in Huang et al
%INPUTS
% imf: a row vector containing an individual IMF, as defined by Liang et al
% method: optional param, see Matlab 'interp1' function.  Default =
%   'spline'.

%$Rev: 97 $
%$Date: 2011-01-14 21:06:08 -0500 (Fri, 14 Jan 2011) $
%$Author: dgibson $

%REFERENCES
%
% Liang H, Bressler SL, Buffalo EA, Desimone R, and Fries P, "Empirical
% mode decomposition of field potentials from macaque V4 in visual spatial
% attention", Biol. Cybern. (2005) 92: 380-392
%
% Huang NE, Wu Z, Long SR, Arnold KC, Chen X, Blank K, "On Instaneous
% Frequency", Advances in Adaptive Data Analysis Vol. 1, No. 2 (2009)
% 177?229

if nargin < 2
    method = 'spline';
end
if nargin < 3
    verboseflag = false;
end

pospeaks = dg_findpks(imf);
negpeaks = dg_findpks(-imf);
peaks = union(pospeaks, negpeaks);
if isempty(peaks)
    error('dg_AMFMdecomp:nopeaks', ...
        'There are no extrema in the IMF');
end
if abs(imf(1)) > abs(imf(peaks(1)))
    startval = abs(imf(1));
else
    startval = abs(imf(peaks(1)));
end
if abs(imf(end)) > abs(imf(peaks(end)))
    endval = abs(imf(end));
else
    endval = abs(imf(peaks(end)));
end
A = interp1([0 peaks length(imf)+1], ...
    [startval abs(imf(peaks)) endval], ...
    1:length(imf), method);
F = imf ./ A;
loopcnt = 1;
newpospeaks = dg_findpks(F);
newnegpeaks = dg_findpks(-F);
if ~isequal(newpospeaks, pospeaks) || ~isequal(newnegpeaks, negpeaks)
    posadded = setdiff(newpospeaks, pospeaks);
    poslost = setdiff(pospeaks, newpospeaks);
    negadded = setdiff(newnegpeaks, negpeaks);
    neglost = setdiff(negpeaks, newnegpeaks);
    pospeaks = newpospeaks;
    negpeaks = newnegpeaks;
    peaks = union(pospeaks, negpeaks);
end
while any(abs(abs(F(peaks)) - 1) > eps)
    % If values at start or end are causing trouble, include them in the
    % new envelope; otherwise, propagate the values from the nearest peak.
    if abs(F(1)) > abs(F(peaks(1)))
        startval = abs(F(1));
    else
        startval = abs(F(peaks(1)));
    end
    if abs(F(end)) > abs(F(peaks(end)))
        endval = abs(F(end));
    else
        endval = abs(F(peaks(end)));
    end
    e = interp1([1 peaks length(F)+1], ...
        [startval abs(F(peaks)) endval], 1:length(F), method);
    F = F ./ e;
    A = A .* e;
    newpospeaks = dg_findpks(F);
    newnegpeaks = dg_findpks(-F);
    if ~isequal(newpospeaks, pospeaks) || ~isequal(newnegpeaks, negpeaks)
        pospeaks = newpospeaks;
        negpeaks = newnegpeaks;
        peaks = union(pospeaks, negpeaks);
    end
    peakvals{loopcnt} = reshape(F(peaks), 1, []);
    loopcnt = loopcnt + 1;
    if loopcnt > 100
        warning('dg_AMFMdecomp:loopcnt', ...
            'Quitting due to loopcnt = %d', loopcnt);
        return
    end
end
if verboseflag
    fprintf('normal termination of dg_AMFMdecomp after %d iterations\n', ...
        loopcnt);
end



