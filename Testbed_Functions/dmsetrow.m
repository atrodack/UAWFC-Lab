function [ DM ] = dmsetrow( DM, rownum, VAL )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


DM.flatten;

vec = linspace(1,1024-31,32)';
row = ((rownum-1)) + vec;

actuatorlist = DM.actuators(:,3);
actuatorlist(row) = VAL;

DM.setActs(actuatorlist);
DM.render;
end


