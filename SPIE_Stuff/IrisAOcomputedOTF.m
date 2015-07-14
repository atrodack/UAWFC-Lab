function [ dOTF, PSF1, PSF2, OTF1, OTF2 ] = IrisAOcomputedOTF( DM, pokeseg, PTTpos_in, Noise_Parameters,Field, A, Phasescreen1, Phasescreen2 )
%[ dOTF, PSF1, PSF2, OTF1, OTF2 ] = IrisAOcomputedOTF( DM, pokeseg, PTTpos_in, Noise_Parameters,Field, A, Phasescreen1, Phasescreen2 )
%   

if nargin < 7
    Phasescreen1 = 1;
    Phasescreen2 = 1;
elseif nargin < 8
    Phasescreen2 = 1;
end

if ~iscell(Noise_Parameters)
    error('Noise_Parameters must be a cell array');
else
    if isempty(Noise_Parameters{5})
        Noise_Parameters{5} = false; %if flag is unset, don't use noise
    end
end

FOV = Field.FOV;
FoV = Field.FoV;
PLATE_SCALE = Field.PLATE_SCALE;
new_spacing = DM.spacing;

F = Field.copy;
F.name = 'IrisAO Field 1';

F2 = Field.copy;
F2.name = 'IrisAO Field 2';

% display('Initializing DM dOTF Positions');

PTTpos_flat = PTTpos_in;
PTT_flat = mapSegments(PTTpos_flat);
PTTpos_poked = PTTpos_in;
PTTpos_poked(pokeseg,1) = (F.lambda) / 4;
PTT_poked = mapSegments(PTTpos_poked);

% display('Flattening DM');

DM.PTT(PTT_flat);
DM.touch;
DM.render;

F.planewave * Phasescreen1 * Phasescreen2 * A * DM;
PSF1 = F.mkPSF(FOV,PLATE_SCALE);

% fprintf('Poking Segment %d \n',pokeseg);

DM.PTT(PTT_poked);
DM.touch;
DM.render;

% display('Making PSF 2');

F2.planewave * Phasescreen1 * Phasescreen2 * A * DM;
PSF2 = F2.mkPSF(FOV,PLATE_SCALE);

if Noise_Parameters{5} == true
    [ PSF1, PSF2 ] = average_noisy_images( PSF1, PSF2, F.grid, Noise_Parameters );
end
    

F.touch; F2.touch;
F.grid(PSF1); F2.grid(PSF2);

% display('Making OTFs');

OTF1 = F.mkOTF2(FoV,new_spacing(1));
OTF2 = F2.mkOTF2(FoV,new_spacing(1));

F.touch; F2.touch;

% display('Computing dOTF');

dOTF = OTF1 - OTF2;
end

