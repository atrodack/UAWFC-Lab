function [tip,tilt] = calctiptiltdOTF(Phase_Mask,Mask_center_point,shift)

if length(shift) == 2
    shiftx = shift(1);
    shifty = shift(2);
else
    shiftx = shift;
    shifty = shift;
end

center_x = Mask_center_point(2);
center_y = Mask_center_point(1);

p_ul = [center_y-shifty,center_x+shiftx];
p_ur = [center_y+shifty,center_x+shiftx];
p_ll = [center_y-shifty,center_x-shiftx];
p_lr = [center_y+shifty,center_x-shiftx];

z1 = Phase_Mask(p_ul(1),p_ul(2));
z2 = Phase_Mask(p_ur(1),p_ur(2));
z3 = Phase_Mask(p_ll(1),p_ll(2));
z4 = Phase_Mask(p_lr(1),p_lr(2));

p1 = [p_ul(1),p_ul(2),z1];
p2 = [p_ur(1),p_ur(2),z2];
p3 = [p_ll(1),p_ll(2),z3];
p4 = [p_lr(1),p_lr(2),z4];

tilt = (p3(3) - p1(3)) / (p3(2) - p1(2));
tip = (p3(3) - p4(3)) / (p1(1) - p4(1));

tip = tip*1e3;
tilt = tilt*1e3;






end