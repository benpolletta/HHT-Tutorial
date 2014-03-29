%extracts the indices of extrema
% Writtezn by Gabriel Rilling
function [indmin, indmax, indzer] = memd_extr(x,t)

if(nargin==1)
  t=1:length(x);
end

m = length(x);

% Finding zeros.

if nargout > 2
  x1=x(1:m-1);
  x2=x(2:m);
  indzer = find(x1.*x2<0); % Zeros are places where x changes from positive to negative.

  if any(x == 0)
    iz = find( x==0 );  % Zeros are places where x equals zero.
    indz = [];
    if any(diff(iz)==1) % Unless x equals zero for more than one time point.
      zer = x == 0;
      dz = diff([0 zer 0]); % Zero-one vector; equals one where x equals zero.
      debz = find(dz == 1); % Indices where x becomes zero after being nonzero.
      finz = find(dz == -1)-1; % Indices where x is zero before becoming nonzero.
      indz = round((debz+finz)/2); % Middle index of a stretch of zeros.
    else
      indz = iz;
    end
    indzer = sort([indzer indz]); % Concatenates places where x equals zero and places where x changes sign.
  end
end

% Now maxima & minima...

d = diff(x); % Looking at "derivative".

n = length(d);
d1 = d(1:n-1);
d2 = d(2:n);
indmin = find(d1.*d2<0 & d1<0)+1; % Places where derivative changes sign from negative to positive are minima.
indmax = find(d1.*d2<0 & d1>0)+1; % Places where derivative changes sign from positive to negative are maxima.


% when two or more successive points have the same value we consider only one extremum in the middle of the constant area
% (only works if the signal is uniformly sampled)

if any(d==0)

  imax = [];
  imin = [];

  bad = (d==0);
  dd = diff([0 bad 0]);
  debs = find(dd == 1);
  fins = find(dd == -1);
  if debs(1) == 1
    if length(debs) > 1
      debs = debs(2:end);
      fins = fins(2:end);
    else
      debs = [];
      fins = [];
    end
  end
  if length(debs) > 0
    if fins(end) == m
      if length(debs) > 1
        debs = debs(1:(end-1));
        fins = fins(1:(end-1));

      else
        debs = [];
        fins = [];
      end
    end
  end
  lc = length(debs);
  if lc > 0
    for k = 1:lc
      if d(debs(k)-1) > 0
        if d(fins(k)) < 0
          imax = [imax round((fins(k)+debs(k))/2)];
        end
      else
        if d(fins(k)) > 0
          imin = [imin round((fins(k)+debs(k))/2)];
        end
      end
    end
  end

  if length(imax) > 0
    indmax = sort([indmax imax]);
  end

  if length(imin) > 0
    indmin = sort([indmin imin]);
  end

end
end
