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

% match lengths:
if length(firstfalseidx) > length(firsttrueidx)
    firsttrueidx(end+1) = length(x) + 1;
end
if length(firsttrueidx) > length(firstfalseidx)
    firstfalseidx(end+1) = length(x) + 1;
end

% Find lengths strings of <false>:
if firstfalseidx(1) == 1
    falselengths = firsttrueidx - firstfalseidx;
else
    % if firstfalseidx(1) ~= 1, then firsttrueidx(1) == 1.
    falselengths = firsttrueidx(2:end) - firstfalseidx(1:end-1);
end

% Fill in the short strings of <false> and remove from index lists:
isshortidx = find(falselengths <= 4 * width);
false2delete = zeros(size(isshortidx));
true2delete = zeros(size(isshortidx));
for k = 1:length(isshortidx)
    if firstfalseidx(1) == 1
        y(firstfalseidx(isshortidx(k)) : ...
            firsttrueidx(isshortidx(k)) - 1) = 1;
        false2delete(k) = isshortidx(k);
        true2delete(k) = isshortidx(k);
    else
        if k < length(isshortidx)
            y(firstfalseidx(isshortidx(k)) : ...
                firsttrueidx(isshortidx(k + 1)) - 1) = 1;
            true2delete(k) = isshortidx(k + 1);
        else
            y(firstfalseidx(isshortidx(k)) : end) = 1;
            true2delete(end) = [];
        end
        false2delete(k) = isshortidx(k);
    end
end
firstfalseidx(false2delete) = [];
firsttrueidx(true2delete) = [];

% Now smooth the transitions:
hw = hanning(2*width+1)';
for k = 1:length(firstfalseidx)
    yidx = firstfalseidx(k) : min(length(y), firstfalseidx(k) + width - 1);
    hwidx = width + (2 : length(yidx) + 1);
    y(yidx) = hw(hwidx);
    yidx = max(1, firsttrueidx(k) - width) : firsttrueidx(k) - 1;
    hwidx = width - length(yidx) + 1 : width;
    y(yidx) = hw(hwidx);
end


