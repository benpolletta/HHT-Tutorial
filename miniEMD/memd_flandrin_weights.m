function w = memd_flandrin_weights(r, envmax, envmin, alpha)
% Weighting function for local EMD which uses the criterion that sifting
% must continue in any section where there are negative frequencies coming
% out of the Hilbert transform.  The local part of the Flandrin criterion
% is also applied to ensure that at least some sifting is done on all
% sections of signal.
%   <width> is the half-width of the Hanning window used to smooth <w>.  I
% think it might be unavoidable to have an a priori scale setting parameter
% like this in any local EMD weighting function.

width = 10;

amp = mean(abs(envmax-envmin))/2;
sx = abs((envmin+envmax)/2)./amp;
isbad = sx > 2*alpha;
if sum(isbad) > length(isbad)/2         % Checking that more than half of indices are bad.
    w = ones(size(r));
else
    w = memd_smoothBinSeries(isbad, width);
end

