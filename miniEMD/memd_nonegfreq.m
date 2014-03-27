function w = memd_nonegfreq(r, envmax, envmin, alpha, width)
% Weighting function for local EMD which uses the criterion that sifting
% must continue in any section where there are negative frequencies coming
% out of the Hilbert transform.
%   <width> is the half-width of the Hanning window used to smooth <w>.  I
% think it might be unavoidable to have an a priori scale setting parameter
% like this in any local EMD weighting function.

hfunc = hilbert(r);
th = unwrap(angle(hfunc));
freq = [NaN diff(th)];
% first segments containing freq < 0 and fill any short gaps between them:
isbad = freq < 0;
if sum(isbad) > length(isbad)/2
    w = ones(size(r));
else
    % O'barf: I have no idea how to create a nicely smoothed waveform from
    % a random rectangular wave without doing a nonlinear iteration of some
    % kind.
    hw = hanning(2*width+1)';
    w = conv(double(freq < 0), hw, 'same');
end

