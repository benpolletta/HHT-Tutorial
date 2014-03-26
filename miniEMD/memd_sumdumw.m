function w = memd_sumdumw(r, envmax, envmin, alpha, crit)
% A dumb test function to determine whether the 'localEMDfunc',
% 'localEMDparam' mechanism works.  <crit> works the same as 'stop' in
% 'memd_emd'.

switch(crit)
    case 'f'
        disp('Flandrin et al. criterion');
        envmoy = (envmin+envmax)/2;
        amp = mean(abs(envmax-envmin))/2;
        sx = abs(envmoy)./amp;
        stop_sift = ~(mean(sx > alpha) > 0.05 | any(sx > 10*alpha));
    case 'h'
        disp('Hello, Ferd.  This is the Huang.');
        envmoy = (envmin+envmax)/2;
        nr = r - envmoy;
        stop_sift = norm(nr-r)/(norm(r)+eps) < alpha;
end

if stop_sift
    w = zeros(size(r));
else
    w = ones(size(r));
end

