% Note: Written using Octave

% Exercise 1 - Noise

clear; clc;

function [] = func (s, N, bins, fs, title1, title2, title3, title4, fig)
  figure(fig);
  
  subplot(2,2,1);
  plot(s);
  title(title1);
  
  subplot(2,2,2);
  hist(s);
  title(title2);
  
  [h, xh] = hist(s,bins,1); % 50 bins normalized histogram => pdf (sum(h)=1)
  subplot(2,2,3);
  plot(xh,h,'-*');
  title(title3);
  
  dft = fft(s); % discrete fourier transform
  dft = fftshift(dft); % shift dft to center
  step = fs/N;
  n = -fs/2 : step : fs/2 - step;
  mag = abs(dft);
  subplot(2,2,4);
  stem(n, mag);
  title(title4);
endfunction


% 1.1, 1.2, 1.3
fig = 1;
N = 10**5;
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
% The length of the signal is half of the original
% The histogram has half the counts but the pdf curve is same
% The length of the DFT is the same (10000 frequencies)

% 1.5
fig = 3;
sc = sin(s);
N = length(sc);
title1 = sprintf('sin(s[n]) N=%d', N);
title2 = 'Histogram';
title3 = 'PDF (hsitogram)';
title4 = 'DFT - Magnitude';
func (sc, N, bins, fs, title1, title2, title3, title4, fig);
% Frequency distribution is unchanged
% the spatial distribution pdf is modified


% 1.6
m = [1 1;2 2] / 2;
sd = conv2(s,m,'same');
Nd = length(sd);
K = max(sd) * max(s);
sd = sd / K;
[h3, xh3] = hist(sd,50);
h3 = h3/length(sd);
figure();
plot(xh3,h3,'*-');
title('pdf - filtering 0.5*[1 1;2 2]');

dft = fft(sd);
f = fftshift(dft);
df = fs/Nd;
figure();
plot(-fs/2:df:fs/2-df,abs(f));
% Frequency distribution is modified (filtering)
% the spatial distribution pdf is unchanged

clear; clc; close all;
% Exercise 2 - DFT
% 2.1

figure();
f = 5;
fs = 50;
t = 0:1/fs:1-1/fs;
xn = sin(2*pi*f*t);
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

figure();
% 2.2
f = 5;
fs = 50;
t = 0:1/fs:1-1/fs;
xn = cos(2*pi*f*t);
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


figure();
% 2.3
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


