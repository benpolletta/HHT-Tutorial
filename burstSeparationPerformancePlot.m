function [a2vals, f2vals, err] = burstSeparationPerformancePlot
f1 = 15;
f2vals = (30:90)/3;
a2vals = logspace(-2, 2, 100);
err = NaN(length(f2vals), length(a2vals));
for f2idx = 1:length(f2vals)
    for a2idx = 1:length(a2vals)
        [mixedwave, subwaves] = make_bursts(f1, f2vals(f2idx), ...
            a2vals(a2idx));
        imfs = memd_emd(mixedwave);
        err(f2idx, a2idx) = rate_imfs(imfs, subwaves);
        fprintf('Done [%d %d]\n', f2idx, a2idx);
    end
end
figure;
imagesc(log(a2vals), f2vals/15, err);
axis xy;
colorbar;
xlabel('log amplitude');
ylabel('freq ratio');
