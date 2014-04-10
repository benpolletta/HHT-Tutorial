function w = benp_emd_local(r, envmax, envmin, alpha, smoothing)
% The guts extraced from memd_emd_local on the benp branch, as modified
% from version of Thu, Mar 27, 2014  2:39:24 PM to work with non-crashing
% version of memd_boundary_conditions.
if isempty(smoothing)
    smoothing = 1000;
end
envmoy = (envmin+envmax)/2;
amp = abs(envmax-envmin)/2;
sx = abs(envmoy)./amp;
w = sx > 2*alpha;

if sum(w) > length(w)/2 % If more than half of indices are bad, forget local.
    w = ones(size(r));
else
    w = memd_smoothBinSeries(w, smoothing);
end

% Flandrin
envmoy = w.*envmoy;
amp = mean(abs(envmax-envmin))/2; % Half of mean difference of max. and min. envelopes is the mean amp. of the signal.
sx = abs(envmoy)./amp; % Divide underlying trend by this mean amplitude at each point.
if ~(mean(sx > alpha) > 0.05 || any(sx > 10*alpha))
    % stop sifting
    w = zeros(size(r));
end
