function [r, ack] = memd_doNothingPre(r, aux)
ack = fliplr(aux);
disp(aux);
fprintf('That did nothing!\n');

