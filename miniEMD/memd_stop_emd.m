% tests if there are enough (3) extrema to continue the decomposition
% Written by Gabriel Rilling
function stop = memd_stop_emd(r)
[indmin,indmax] = memd_extr(r);
ner = length(indmin) + length(indmax);
stop = ner < 3;
end
