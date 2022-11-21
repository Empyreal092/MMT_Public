% Action and energy diagnostics
Conserve_Calc
Nsample = [Nsample Ntotal];
Hsample = [Hsample Htotal]; H1sample = [H1sample H1]; H2sample = [H2sample H2];
Psample = [Psample Ptotal];
% Slope calculations and save
LSFit
SlopeSampleI = [SlopeSampleI gI]; % keep track of the slopes
SlopeSampleDs = [SlopeSampleDs gDs]; % keep track of the slopes
SlopeSampleDl = [SlopeSampleDl gDl]; % keep track of the slopes