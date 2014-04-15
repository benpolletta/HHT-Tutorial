function err = rate_imfs(imfs, waves)
%err = rate_imfs(imfs, waves)
% Breaks <imfs> into 1000-sample pieces.  If the number of samples is not
% an integral multiple of 1000, then any extra samples are ignored. For
% each 1000-sample piece, computes the zero-lag covariance of each IMF with
% each channel in <waves>.  The IMF-wave pairs are then sorted by
% decreasing covariance value, and waves are assigned to IMFs according to
% the sorted pairs.  If either member of any pair has already been
% assigned, that pair is ignored.  Any IMFs that have still not been
% assigned a matching wave at the end of this process are ignored. Finally,
% the squared error is computed at each sample between each of the
% remaining IMFs and its matching wave, summed over samples, and divided by
% the sum of the squares of the points in <waves> (again ignoring any extra
% samples not included in the 1000-point pieces), and returned as <err>.
% Note that this formula strongly penalizes the splitting of single
% subwaves into multiple IMFs.
%INPUTS
% imfs: in IMFs X samples form.
% waves: in subwaves X samples form, same size as <imfs>.
%OUTPUTS
% err: scalar.

if ~isequal(size(waves, 2), size(imfs, 2))
    error('rate_imfs:size', ...
        '<imfs> and <waves> must be the same length.');
end
pcsz = 1000;
numpieces = fix(size(imfs,2) / pcsz);
squaredError = 0;
for piecenum = 1:numpieces
    sampidx = pcsz*(piecenum-1) + (1:pcsz);
    xcovs = zeros(size(imfs,1), size(waves,1));
    for imfnum = 1:size(imfs,1)
        for wavenum = 1:size(waves,1)
            xcovs(imfnum, wavenum) = xcov(imfs(imfnum, sampidx), ...
                waves(wavenum, sampidx), 0, 'coeff');
        end
    end
    matchedwave = zeros(size(imfs,1), 1); % index into <waves> for each IMF
    assigned = false(size(waves,1), 1); % true if the wave is already assigned
    [~, xcovidx] = sort(xcovs(:), 'descend');
    for xcovidx2 = 1:length(xcovidx)
        [imfnum, wavenum] = ind2sub(size(xcovs), xcovidx(xcovidx2));
        if matchedwave(imfnum) || assigned(wavenum)
            continue
        end
        matchedwave(imfnum) = wavenum;
        assigned(wavenum) = true;
        if all(matchedwave) || all(assigned)
            if ~all(assigned)
                warning('rate_imfs:unassigned', ...
                    'Waves %s were not assigned to any IMF in piece %d', ...
                    mat2str(find(~assigned)), piecenum);
            end
            break
        end
    end
    for imfnum = reshape(find(matchedwave~=0), 1, [])
        squaredError = squaredError + ...
            sum( ( imfs(imfnum, sampidx) - ...
            waves(matchedwave(imfnum), sampidx) ).^2 );
    end
end
err = squaredError / sum(sum(waves(:, 1:pcsz*numpieces).^2));

