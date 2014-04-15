function [fspacevals, err] = chaosSeparation

warning off rate_imfs:unassigned
fspacevals = 10 .^ ((1:50)/100);
for fidx = 1:length(fspacevals)
    [mixedwave, subwaves] = make_chaos(fspacevals(fidx), [10 100]);
    imfs = memd_emd(mixedwave);
    err(fidx) = rate_imfs(imfs, subwaves);
    fprintf('Done %d/%d\n', fidx, length(fspacevals));
end
