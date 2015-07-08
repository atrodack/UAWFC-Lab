function [ PTTpos ] = floatpiston_IrisAO( PTTpos )
%[ PTTpos ] = floatpiston_IrisAO( PTTpos )
%   centers the pistons around 0 to minimize the amount of stroke used

pistonlist = PTTpos(:,1);

maxpiston = max(pistonlist(abs(pistonlist)>0));
minpiston = min(pistonlist(abs(pistonlist)>0));
piston_float = (abs(maxpiston)+abs(minpiston))/2;
if abs(maxpiston) >= abs(minpiston)
    pistonlist(abs(pistonlist)>0) = pistonlist(abs(pistonlist)>0) - piston_float;
else
    pistonlist(abs(pistonlist)>0) = pistonlist(abs(pistonlist)>0) + piston_float;
end


PTTpos(:,1) = pistonlist;


end

