function [ pixel_act_map, Areal_Averaging_radius ] = computeBMCactpixelmap( DM, A, Field, pokeact )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

F = Field.copy;
F.name = 'Actuator Pixel Map Field';

testacts1 = [1,32,993];
% testacts1 = [33,63,962]
testacts2 = randi(1024,250,1);


DM.flatten;
Ppos_flat = DM.actuators(:,3);
Ppos_flat(testacts1) = 1e-5;

dOTF = BMCcomputedOTF( DM, pokeact, Ppos_flat, false, Field, A);
figure;
plotComplex(dOTF,5);
axis xy;
sqar;



fprintf('Pick the Bottom Left Point\n\n');
btleft = pickPoint;
hold on
plot(btleft(2),btleft(1),'r*');
hold off


fprintf('Pick the Top Left Point\n\n');
topleft = pickPoint;
hold on
plot(topleft(2),topleft(1),'r*');
hold off


fprintf('Pick the Bottom Right Point\n\n');
btright = pickPoint;
hold on
plot(btright(2),btright(1),'r*');
hold off

x = linspace(btleft(2),btright(2),32);
y = linspace(btleft(1),topleft(1),32);
[X,Y] = meshgrid(x,y);



DM.flatten;
Ppos_flat = DM.actuators(:,3);
Ppos_flat(testacts2) = 1e-6;


dOTF = BMCcomputedOTF( DM, 698, Ppos_flat, false, Field, A);

clf;
plotComplex(dOTF,5);
axis xy;
sqar;


hold on
plot(ceil(X),ceil(Y),'r*')
hold off








end

