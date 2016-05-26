function [OTF,TT] = vernierAlignOTF(OTF,RADIUS,CEN)

% [OTF,TT] = vernierAlignOTF(OTF,RADIUS,CEN)

SZ = size(OTF);

[X,Y,R] = mkImageCoords(OTF,1,CEN);
UV = [X(:) Y(:)];
MASK = (R(:)<=RADIUS);

[OTF,TT] = alignOTF(OTF,UV,MASK);

%CUBE(:,:,nn) = real(fftshift(ifft2(circshift(OTF,1-CENTER))));

