function y = memd_smoothBinSeries(x, width)
% Converts a logical vector (e.g. representing some property of a time
% series) to a series of floats with a smoothly graded transition from
% regions of of ones to regions of zeros.  One half of a Hanning window is
% used for this purpose.  Any strings of <false> that are of length 4*width
% or less are simply replaced with ones.
%INPUTS
% x: a logical row vector.
% width: an integer specifying the half-width of the Hanning window.
%OUTPUTS
% y: a double vector of same size as <x>.

y = double(x);
firstfalseidx = find([~x(1) x(1:end-1) & ~x(2:end)]);
firsttrueidx = find([x(1) ~x(1:end-1) & x(2:end)]);

% Match lengths, using "length(x) + 1" to denote that the run of same
% values goes all the way to the end of the vector.  Note that if either
% length is strictly greater than the other one, then x must both start and
% end with that value.
if length(firstfalseidx) > length(firsttrueidx)
    firsttrueidx(end+1) = length(x) + 1;
end
if length(firsttrueidx) > length(firstfalseidx)
    firstfalseidx(end+1) = length(x) + 1;
end

% Find lengths strings of <false>:
if firstfalseidx(1) == 1
    % We have runs of false marked by firstfalseidx, firsttrueidx.
    falselengths = firsttrueidx - firstfalseidx;
else
    % We have runs of true marked by firsttrueidx, firstfalseidx.  That
    % means the runs of false are in between.
    falselengths = firsttrueidx(2:end) - firstfalseidx(1:end-1);
    % ...and it also means that there might be a trailing run of false at
    % the end that we didn't compute yet:
    if firstfalseidx(end) <= length(x)
        falselengths(end+1) = length(x) - firstfalseidx(end) + 1;
    end
end

% Fill in the short strings of <false> and remove from index lists.
% <isshortidx> is an index into <falselengths>, whose length is either
% length(firstfalseidx) or one shorter.
isshortidx = find(falselengths <= 4 * width);
false2delete = zeros(size(isshortidx));
true2delete = zeros(size(isshortidx));
for k = 1:length(isshortidx)
    if firstfalseidx(1) == 1
        % We have runs of false marked by firstfalseidx, firsttrueidx.
        % Whether the last run of false runs to the end or not, it ends at
        % firsttrueidx(end) - 1, so there is nothing to worry about.
        y(firstfalseidx(isshortidx(k)) : ...
            firsttrueidx(isshortidx(k)) - 1) = 1;
        false2delete(k) = isshortidx(k);
        true2delete(k) = isshortidx(k);
    else
        % We have runs of true marked by firsttrueidx, firstfalseidx.  That
        % means that a run of false that goes all the way to the end will
        % need special handling.  It also means that length(falselengths)
        % is length(firstfalseidx) - 1.
        if isshortidx(k) < length(falselengths)
            % Not the last run of false, nothing to worry about.
            y(firstfalseidx(isshortidx(k)) : ...
                firsttrueidx(isshortidx(k) + 1) - 1) = 1;
            false2delete(k) = isshortidx(k);
            true2delete(k) = isshortidx(k) + 1;
        else
            % isshortidx(k) is the very last entry in falselengths.
            % Since we are dealing with runs of true, then either
            % this run of false goes all the way to the end, in which
            % case firstfalseidx(end) <= length(x); or there is another run
            % of true that goes all the way to the end, in which case
            % firstfalseidx(end) = length(x) + 1.
            if firstfalseidx(end) <= length(x)
                y(firstfalseidx(isshortidx(k)) : end) = 1;
                % instead of deleting the marker for the start of the run
                % of false that we are filling in, we move the marker off
                % the end of the vector:
                firstfalseidx(end) = length(x) + 1;
            else
                % there is another run of true
                y(firstfalseidx(isshortidx(k)) : ...
                    firsttrueidx(isshortidx(k) + 1) - 1) = 1;
                false2delete(k) = isshortidx(k);
                true2delete(k) = isshortidx(k) + 1;
            end
        end
    end
end
false2delete(false2delete==0) = [];
true2delete(true2delete==0) = [];
firstfalseidx(false2delete) = [];
firsttrueidx(true2delete) = [];

% Now smooth the transitions:
hw = hanning(2*width+1)';
for k = 1:length(firstfalseidx)
    if firstfalseidx(k) ~= 1
        % Add a right-half hanning starting at firstfalseidx(k):
        yidx = firstfalseidx(k) : min(length(y), firstfalseidx(k) + width - 1);
        hwidx = width + (2 : length(yidx) + 1);
        y(yidx) = hw(hwidx);
    end
    if firsttrueidx(k) <= length(x)
        % Add a left-half hanning ending at firsttrueidx(k) - 1:
        yidx = max(1, firsttrueidx(k) - width) : firsttrueidx(k) - 1;
        hwidx = width - length(yidx) + 1 : width;
        y(yidx) = hw(hwidx);
    end
end


