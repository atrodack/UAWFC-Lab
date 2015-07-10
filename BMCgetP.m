function [ Ppos, masked_dOTF ] = BMCgetP( dOTF, DM, lambda, onAct_locations, ALGO)
%[ Ppos, masked_dOTF ] = BMCgetP( dOTF, DM, onAct_locations, ALGO)
%   Detailed explanation goes here


if nargin < 5
    ALGO = 'gold';
end

k = (2*pi) / lambda;

% Build a grid Object and store dOTF in it
G = AOGrid(1);
G.spacing(DM.spacing);
G.grid(dOTF);

% Use Grid functionality to mask the dOTF
[X, Y] = G.COORDS;
GT = G.copy;
GT.grid(~(Y - tand(150)*X > -0.0004)); 
G * GT;


masked_dOTF = G.grid;
dOTF = -1i * conj(dOTF);

phase = angle(dOTF);
unwrapped_phase = uwrap(phase,ALGO);
OPL = unwrapped_phase / k;
OPL = OPL .* GT.grid;




%     for n = 1:length(DM.OnActs)
%         piston_area = OPL .* ActMasks{n};
%         Ppos2(n,1) = mean(mean(abs(piston_area)>0));
%     end

for n = 1:length(DM.OnActs)
    Ppos(n,1) = OPL(onAct_locations{n}(2),onAct_locations{n}(1));
end

SELECT = Ppos(Ppos~=0);
MEAN = mean(SELECT);
for n = 1:length(Ppos)
    if Ppos(n) ~= 0
        Ppos(n) = Ppos(n) - MEAN;
    end
end



end

