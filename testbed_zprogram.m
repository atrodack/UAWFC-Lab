function RECON = testbed_zprogram(D,Nmax,DM,dOTF,lambda,RECON)
%RECON = testbed_zprogram(D,Nmax,DM,dOTF,lambda,RECON)
%
% Author:
% ATRodack
%
% Date:
%  201554
%
% Function to attempt programming a Reconstructor matrix for dOTF
%
% Inputs:
% D -- Input Diameter to scale Zernikes appropriately
% nmax -- maximum Zernike n order to go to
% lambda -- wavelength of reconstruction
%
% Output:
% RECON -- constructed reconstructor matrix

if(nargin < 6)
    RECON = struct;
elseif(nargin < 5 || nargin > 6)
    fprintf('\nPrinting Help Message:\nRECON = testbed_zprogram(D,Nmax,DM,dOTF,lambda,RECON)\n');
    error('Inputs are incorrect'); 
elseif(nargin == 6)
    isa(RECON,'struct');
end


RECON.TrainingMethod = 'zernike';
RECON.OWD = Nmax;
RECON.DM = DM;
RECON.WFS = dOTF;
RECON.D = D;

NZmodes = (Nmax+1)*(Nmax+2)/2;
RECON.ACTS = zeros(RECON.DM.nActs,NZmodes);




end
    