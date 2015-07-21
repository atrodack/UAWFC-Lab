function DM_Pistons = makeDMfitsfile(DM, rownums, colnums, pistonrow, pistoncol)


if length(pistonrow) == 1
    pistonrow = pistonrow .* ones(length(rownums),1);
end
if length(pistoncol) == 1
    pistoncol = pistoncol .* ones(length(colnums),1);
end

if length(rownums) ~= length(colnums) || length(pistonrow) ~= length(rownums) || length(pistonrow) ~= length(pistoncol)
    error('Vector Length Mismatch');
end

DM.flatten;
Ppos_flat = DM.actuators(:,3);

for n = 1:length(rownums)
    DM = dmsetrow(DM,rownums(n),pistonrow(n));
    Ppos_flat = Ppos_flat + DM.actuators(:,3);
    DM = dmsetcol(DM,colnums(n),pistoncol(n));
    Ppos_flat = Ppos_flat + DM.actuators(:,3);
end

DM.setActs(Ppos_flat);
DM_Pistons = reshape(DM.actuators(:,3),[32,32]);
DM_Pistons = single((DM_Pistons));

tempdir = pwd;
cd /home/alex/Desktop/Testbed_fits_files;
% cd /home/lab/src/scripts
fitswrite(DM_Pistons,'DM_Pistons.fits');
cd(tempdir);


figure(6);
imagesc(DM_Pistons); axis xy; sqar;
axis xy;
sqar;
