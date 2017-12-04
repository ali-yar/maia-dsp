% Note: This script was written in Octave.
% Some small differences may happen if running in Matlab

clear; clc;

% Exercice 1 - 2D frequency plan

%1.1
N = 64; % img size
T = 1;
Ts = 1/N; Fs = N; df = Fs/N; % sampling
Im (N/8:N/4, (N/4)+1:N/2) = 1;
Im (1:N/4, (N/2)+1:N) = Im;
Im ((N/4)+1:N/2, :) = Im;
Im ((N/2)+1:3*N/4, :) = Im (1:N/4 , :);
Im ((3*N/4)+1:N, :) = Im (1:N/4 , :);

figure(1);
subplot(321);
imagesc(Im);
colormap('gray');
title('Image');

%1.2
dft = fft2(Im)/(N.^2); % fft normalized by img size
dft = fftshift(dft);
mag = abs(dft);
arg = angle(dft) * 180/pi; % phase in degrees

subplot(323);
imagesc(real(dft));
title('2D DFT - Real part');

subplot(324);
imagesc(imag(dft));
title('2D DFT - Imaginary part');

subplot(325);
imagesc(mag);
title('2D DFT - Amplitude');

subplot(326);
imagesc(arg);
title('2D DFT - Phase');

%%

clear; clc;

% Exercice 2 - Reconstruction with amplitude and/or phase

% define common size for both images
height = 256;
width = 256;

%2.1 load image 1
figure(2);
Im1 = imread('cameraman.tif'); % if using matlab
%Im1 = imread('./cameraman.png'); % if using local image
Im1 = imresize(Im1, [height, width]);
subplot(221);
colormap('gray');
imagesc(Im1);
title('Image 1');

%2.2 find dft magnitude and phase of image 1
dft1 = fft2(Im1);
shifted_dft1 = fftshift(dft1);
mag1 = abs(shifted_dft1);
arg1 = angle(shifted_dft1);
log1 = log(mag1+1); % log transformation to display low amplitudes 

subplot(222);
imagesc(mag1);
colormap gray;
title('2D DFT - Amplitude');

subplot(223);
imagesc(arg1);
title('2D DFT - Phase');

subplot(224);
imagesc(log1);
title('Log Transform');

% 2.3
figure(3);
reconstructed_mag1 = mag1 .* exp(i*0);
reconstructed_mag1 = ifft2(reconstructed_mag1);
subplot(121);
colormap('gray');
imagesc(uint8(abs(reconstructed_mag1)));
title('Image 1 reconstructed from mag ');

reconstructed_arg1 = 1 .* exp(i*arg1);
reconstructed_arg1 = ifft2(reconstructed_arg1);
subplot(122);
colormap('gray');
imagesc(uint8(abs(reconstructed_arg1)));
title('Image 1 reconstructed from Phase');

% 2.4 load image 2
figure(4);
Im2 = imread('moon.tif'); % if using matlab
%Im2 = imread('./moon.jpg'); % if using local image
Im2 = imresize(Im2, [height, width]);
subplot(221);
colormap('gray');
imagesc(Im2);
title('Image 2');

% 2.5 find dft magnitude and phase of image 2
dft2 = fft2(Im2);
shifted_dft2 = fftshift(dft2);
mag2 = abs(shifted_dft2);
arg2 = angle(shifted_dft2);
log2 = log(mag2+1);

subplot(222);
imagesc(mag2);
colormap gray;
title('2D DFT - Amplitude');

subplot(223);
imagesc(arg2);
title('2D DFT - Phase');

subplot(224);
imagesc(log2);
title('Log Transform');

% 2.6 Reconstruct image 1 and image 2
figure(5);
reconstructed1 = mag1 .* exp(i*arg1);
reconstructed1 = ifft2(reconstructed1);
subplot(121);
colormap('gray');
imagesc(uint8(abs(reconstructed1)));
title('Image 1 reconstructed');

reconstructed2 = mag2 .* exp(i*arg2);
reconstructed2 = ifft2(reconstructed2);
subplot(122);
colormap('gray');
imagesc(uint8(abs(reconstructed2)));
title('Image 2 reconstructed');

% 2.7 Reconstruct from mixed magnitudes and phases
figure(6);
reconstructed3 = mag2 .* exp(i*arg1);
reconstructed3 = ifft2(reconstructed3);
subplot(121);
colormap('gray');
imagesc(uint8(abs(reconstructed3)));
title('Image 2 Mag & Image 1 Phase');

reconstructed4 = mag1 .* exp(i*arg2);
reconstructed4 = ifft2(reconstructed4);
subplot(122);
colormap('gray');
imagesc(uint8(abs(reconstructed4)));
title('Image 1 Mag & Image 2 Phase');