Zdm = 2e4;
Zconj = Zprop/2;
Zcc = Zconj;

FOV = 2;
dFOV = FOV/100;

P = [144 127];

NPertAct = 1;
N = 12;

TT = 2e-6*[1 0];
dPISTON = F.lambda/8;

% mkIrisAODM;
colormap(gray);

ZCC=0:Zdm/8:2*Zdm;
ZP = (0:1/16:10)*F.lambda

for n=1:length(ZCC)
% for n=1:length(ZP)
    
    Zcc = ZCC(n);
%     dPISTON = ZP(n);
    
    DM.trueUp.setPiston(N,F.lambda/5).setTipTilt(N,[1 1]*1e-6);
    F.planewave*A;F.propagate(Zdm)*DM;
    F.propagate(Zcc-Zdm);
    OTF0 = computeOTF(F.mkPSF(FOV,dFOV));
    
    DM.setTipTilt(NPertAct,TT);
    DM.setPiston(NPertAct,dPISTON);
    
    subplot(2,2,1);
    [x,y] = F.coords;
    imagesc(x,y,F.abs.^2);sqar;
    %title(sprintf('Intensity@CP: Zcc/Zdm=%.2g',Zcc/Zdm));
    title(sprintf('Intensity@CP: piston=%g \pi Zcc/Zdm=%.2g',...
        dPISTON*F.k,Zcc/Zdm));
    axis xy
    axis off
    sqar;
    
    F.planewave*A;F.propagate(Zdm)*DM;
    F.propagate(Zcc-Zdm);
    %F.planewave*A;F.propagate(Zdm)*DM;
    %F.propagate(Zcc-Zdm);
    OTF = computeOTF(F.mkPSF(FOV,dFOV));
    %axis off
    
    subplot(2,2,2);
    F.show;
    title('Pert Field in CP');
    axis xy
    axis off
    
    dOTF = OTF - OTF0;
    
    subplot(2,2,3);
    imagesc(abs(dOTF));sqar;
    colorbar;
    title('abs(dOTF)');
    axis xy
    
    subplot(2,2,4);
    PHASE = uwrap(angle(dOTF));
    %PHASE = PHASE - PHASE(101,101);
    PHASE = PHASE - PHASE(P(1),P(2));
    %imagesc(PHASE);sqar;
    imagesc(PHASE,[-1 1]*2);sqar;
    colorbar;
    title('angle(dOTF)');
    axis xy
    
    drawnow;
    
    %saveJPEG(sprintf('PFRAME_%02d.jpg',n),120);
    saveJPEG(sprintf('FRAME_%02d.jpg',n),120);
    
end


