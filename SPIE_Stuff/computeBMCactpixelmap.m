function [ calibrated_BMC_act_locations ] = computeBMCactpixelmap( DM, A, Field, pokeact )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

F = Field.copy;
F.name = 'Actuator Pixel Map Field';

% testacts1 = [1,32,993];



DM.flatten;
DM = dmsetrow(DM,1,1e-6);
Ppos_flat = DM.actuators(:,3);
DM = dmsetcol(DM,1,1e-6);
Ppos_flat = Ppos_flat + DM.actuators(:,3);
Ppos_flat(Ppos_flat~=0) = 1e-6;

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

DM = dmsetcol(DM,16,1e-6);
Ppos_flat = DM.actuators(:,3);


dOTF = BMCcomputedOTF( DM, 698, Ppos_flat, false, Field, A);

clf;
plotComplex(dOTF,5);
axis xy;
sqar;


% hold on
% plot(ceil(X),ceil(Y),'r*')
% hold off

calibrated_BMC_act_locations = cell(DM.nActs,1);
counter = 1;
for n = 1:32
    for m = 1:32
        hold on
        calibrated_BMC_act_locations{counter} = [X(m,n),Y(m,n)];
        calibrated_BMC_act_locations{counter} = round(calibrated_BMC_act_locations{counter});
        plot(calibrated_BMC_act_locations{n}(1),calibrated_BMC_act_locations{n}(2),'*g');
        drawnow;
        hold off
        counter = counter + 1;
    end
end








end

