function [DM] = dmsetcol(DM, colnum, VAL)


DM.flatten;

vec = linspace(1,32,32)';
column = ((colnum-1)*32) + vec;

actuatorlist = DM.actuators(:,3);
actuatorlist(column) = VAL;

DM.setActs(actuatorlist);
DM.render;
end
