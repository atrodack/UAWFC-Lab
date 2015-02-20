clear all
clc
close all

load('dOTF_phase.mat')
load('2015_2_18_14_29_dOTFstructure_data_run_1.mat')


x = 225:233;
y = 292:300;

seg8 = dOTF_phase(x,y);
imagesc(seg8)
p1 = pickPoint;
p2 = pickPoint;
p3 = pickPoint;

p1(1,3) = seg8(p1(2),p1(1));
p2(1,3) = seg8(p2(2),p2(1));
p3(1,3) = seg8(p3(2),p3(1));

[fx,fy] = gradient(dOTF_phase,0.1);
