function Deconv = Wiener(Field,Blur,noise,gamma)
%Input a field with a noise estimation (variance, data or some other information for the noise).
%Options are to input the field and blur in the Fourier domain to
%immediately start with multiplications (default), or to convert from the spatial
%domain to apply filtering techniques (uncomment the other FIELD and BLUR
%options). If you want to emulate the "Hotdog" filter, just set noise = 1,
%and use gamma as in Hotdog.
%% Fields
% FIELD = fftshift(fft2(fftshift(Field)));
FIELD = Field;

% BLUR = fftshift(fft2(fftshift(Blur)));
BLUR = Blur;

%% Estimate Noise
% Pxx = abs(FIELD).^2; %Power spectrum (including noise effects) of the signal, assuming FIELD input in the Fourier domain
% SNR = (Pxx)./noise;
SNR = 1/noise;

%% Wiener filtering
dWiener = conj(BLUR)./(abs(BLUR).^2 + gamma*(1/SNR));

%Go back to the spatial domain
DECONV = FIELD.*dWiener;
Deconv = fftshift(ifft2(ifftshift(DECONV)));

end