function [OTF,TT] = alignOTF(OTF,UV,SELECT)

% [OTF,TT] = alignOTF(OTF,UV,SELECT)
%
% JLCodona 20100629.  At SPIE 2010 San Diego!

% OTF = fftshift2d(fft2(circshift(CUBE(:,:,n),1-PEAKS(n,:))));
% OTF = OTF/OTF(CENTER(1),CENTER(2)); % Normalize!
PHASE = angle(OTF(SELECT));
X = UV(SELECT,1);
Y = UV(SELECT,2);
Xo = UV(:,1);
Yo = UV(:,2);

XCOEFS = polyfit(X,PHASE,1);
PHASEx1 = polyval(XCOEFS,X);

YCOEFS = polyfit(Y,PHASE-PHASEx1,1);
PHASEy1 = polyval(YCOEFS,Y);

XCOEFS = polyfit(X,PHASE-PHASEy1,1);
PHASEx = polyval(XCOEFS,X);

YCOEFS = polyfit(Y,PHASE-PHASEx,1);
% PHASEy = polyval(YCOEFS,Y);

TT = [XCOEFS(1),YCOEFS(1)];

dPHASE = polyval(XCOEFS,Xo)+polyval(YCOEFS,Yo);
OTF(:) = OTF(:).*exp(-1i*dPHASE);
