function [ Ppos, masked_dOTF ] = BMCgetP( dOTF, DM, lambda, onAct_locations, ALGO)
%[ Ppos, masked_dOTF ] = BMCgetP( dOTF, DM, onAct_locations, ALGO)
%   Detailed explanation goes here
load dOTF_act_698_mask.mat;

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

phase = angle(dOTF);
unwrapped_phase = uwrap(phase,ALGO);
OPL = unwrapped_phase / k;
OPL = OPL .* dOTF_act_698_mask .* GT.grid;


for n = 1:length(DM.OnActs)
    x = onAct_locations{n}(1);
    y = onAct_locations{n}(2);
    piston_area = OPL(y-7:y+7,x-7:x+7);
    Ppos(n,1) = mean(mean(piston_area(abs(piston_area)>0)));
    if isnan(Ppos(n,1))
        Ppos(n,1) = 0;
    end
end


% for n = 1:length(DM.OnActs)
%     Ppos(n,1) = OPL(onAct_locations{n}(2),onAct_locations{n}(1));
% end

SELECT = Ppos(Ppos~=0);
MEAN = mean(SELECT);
for n = 1:length(Ppos)
    if Ppos(n) ~= 0
        Ppos(n) = Ppos(n) - MEAN;
    end
end



end

