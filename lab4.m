% Note: This script was written in Octave.
% There might be small differences when running it in Matlab

% Exercise 1 - Noise

clear; clc;

% If using Matlab, move this function to a separate file
function [] = func (s, N, bins, fs, title1, title2, title3, title4, fig)
  figure(fig);
  
  subplot(2,2,1);
  plot(s);
  title(title1);
  
  subplot(2,2,2);
  hist(s);
  title(title2);
  
  [h, xh] = hist(s,bins); % 50 bins histogram
  h = h/N; % normalized histogram => pdf (sum(h)=1)
  subplot(2,2,3);
  plot(xh,h,'-*');
  title(title3);
  
  dft = fft(s); % discrete fourier transform
  dft = fftshift(dft); % shift dft to center
  step = fs/N;
  n = -fs/2 : step : fs/2 - step;
  mag = abs(dft);
  subplot(2,2,4);
  plot(n, mag);
  title(title4);
end


% 1.1, 1.2, 1.3
fig = 1;
N = 10^5;
fs = 10000; % sampling frequency = 10 kHz
bins = 50;
title1 = sprintf('s[n] (Gaussian noise) N=%d', N);
title2 = 'Histogram';
title3 = 'PDF (hsitogram)';
title4 = 'DFT - Magnitude';
s = randn(1,N); % gaussian noise
s = s - mean(s);
func (s, N, bins, fs, title1, title2, title3, title4, fig);
% The length of the DFT is 10000 frequencies

% 1.4 change of spectrum
fig = 2;
sb = s(1:2:end);
N = length(sb);
title1 = sprintf('Sub-sampled signal (factor of 2) N=%d', N);
title2 = 'Histogram';
title3 = 'PDF (hsitogram)';
title4 = 'DFT - Magnitude';
func (sb, N, bins, fs, title1, title2, title3, title4, fig);
% The pdf is similar
% The frequency distribution is different

% 1.5
fig = 3;
sc = sin(s);
N = length(sc);
title1 = sprintf('sin(s[n]) N=%d', N);
title2 = 'Histogram';
title3 = 'PDF (hsitogram)';
title4 = 'DFT - Magnitude';
func (sc, N, bins, fs, title1, title2, title3, title4, fig);
% The pdf is different
% The frequency distribution is same


% 1.6
fig = 4;
m = [1 1;2 2] / 2;
sd = conv2(s,m,'same');
K = max(sd) * max(s);
sd = sd / K; % normalized signal
N = length(sd);
title1 = sprintf('s[d] = [1 1;2 2]/2 * s[n] N=%d', N);
title2 = 'Histogram';
title3 = 'PDF - filtering 0.5*[1 1;2 2]';
title4 = 'DFT - Magnitude';
func (sd, N, bins, fs, title1, title2, title3, title4, fig);
% The pdf is same
% The frequency distribution is different

% Exercise 2 - DFT

clear; clc;

% 2.1
figure(5);
f = 5;
fs = 50;
t = 0:1/fs:1-1/fs;
xn = sin(2*pi*f*t);
N = length(xn);
freq = (-N/2 : N/2 -1) * fs/N;
xf = fftshift(fft(xn));
subplot(221);
plot(t,xn);
title('Sine Signal');
xlabel('Times (sec)');
ylabel('Amplitude');
subplot(222);
plot(freq,abs(xf));
title('Magnitude');
xlabel('Frequency');
ylabel('|X(f)|');
subplot(223);
plot(t,real(xf));
title('Real');
xlabel('Frequency');
ylabel('Re(X(f))');
subplot(224);
plot(t,imag(xf));
title('Imaginary');
xlabel('Frequency');
ylabel('Im(X(f))');

% 2.2
figure(6);
f = 5;
fs = 50;
t = 0:1/fs:1-1/fs;
xn = cos(2*pi*f*t);
N = length(xn);
freq = (-N/2 : N/2 -1) * fs/N;
xf = fftshift(fft(xn));
subplot(221);
plot(t,xn);
title('Cosine Signal');
xlabel('Times (sec)');
ylabel('Amplitude');
subplot(222);
plot(freq,abs(xf));
title('Magnitude');
xlabel('Frequency');
ylabel('|X(f)|');
subplot(223);
plot(t,real(xf));
title('Real');
xlabel('Frequency');
ylabel('Re(X(f))');
subplot(224);
plot(t,imag(xf));
title('Imaginary');
xlabel('Frequency');
ylabel('Im(X(f))');

% The magnitude of the DFT of sine and cosine signals are same.
% The real parts are different.
% The imaginary parts have similar shape but different amplitudes.

% 2.3
figure(7);
f = 5;
fs = 50;
t = 0:1/fs:1-1/fs;
xn = square(2*pi*f*t);
N = length(xn);
freq = (-N/2 : N/2 -1) * fs/N;
xf = fftshift(fft(xn));
subplot(221);
plot(t,xn);
title('Signal');
xlabel('Times (sec)');
ylabel('Amplitude');
subplot(222);
plot(freq,abs(xf));
title('Magnitude');
xlabel('Frequency');
ylabel('|X(f)|');
subplot(223);
plot(t,real(xf));
title('Real');
xlabel('Frequency');
ylabel('Re(X(f))');
subplot(224);
plot(t,imag(xf));
title('Imaginary');
xlabel('Frequency');
ylabel('Im(X(f))');