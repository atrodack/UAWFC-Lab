function [ F ] = gothroughcoronagraph( F,FPMASK,LYOTSTOP )
%[ F ] = gothroughcoronagraph( F,FPMASK,LYOTSTOP )

if ~isa(F,'AOField')
    error('F must be an AOField Object');
end

fprintf('Going through the Coronagraph\n');
F.grid(F.fft/F.nx); % Go to the focal plane.
F*FPMASK; % Pass through the focal plane mask.
F.grid(F.fft/F.nx); % Go to the Lyot pupil plane.
F*LYOTSTOP; % Pass through the Lyot Stop.

end

