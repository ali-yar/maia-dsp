% Exercise 1 - 1D COnvolution
% 1.1

figure(1);

x = [1 -2 3 -4 3 2 1];
h = [3 2 1 -2 1 0 -4 0 3];
y = conv(x,h);
N = length(x)+ length(h) - 1;
n = 0:N-1;
subplot(2,1,1);
stem(n, y);
title('Convolution');

x2 = [x zeros(1,length(y)-length(x))];
y1 = filter(h,1,x2);
subplot(2,1,2);
stem(n, y1);
title('Filter');

% The convolution and filter produces the same output signal

%1.2
clear; clc;

figure(2);

x = randi([-5 5],1,10);
h = randi([-5 5],1,15);
y = conv(x,h);
lx = length(x);
lh = length(h);
ly = length(y);

x2 = [x zeros(1,ly-lx)];
y1 = filter(h,1,x2);

n = 0: ly-1;
subplot(3,2,1);
stem((0:lx-1),x);
title('Signal x[n]');

subplot(3,2,2);
stem((0:lh-1),h);
title('Impulse h[n]');

subplot(3,2,3:4);
stem(n,y);
title('Convolution');

subplot(3,2,5:6);
stem(n,y1);
title('Filter');

% Exercice 2 - 2D Convolution

% 2.1
clear; clc;

figure(3);

subplot(1,2,1);
I = imread('cameraman.tif');
imshow(I);
title('Original');

K = [1 4 6 4 1; 4 16 24 16 4; 6 24 36 24 6; 4 16 24 16 4; 1 4 6 4 1] / 256.;
subplot(1,2,2);
I = im2double(I);
I2 = conv2(I, K);
imshow(I2);
title('Smoothing with Gaussian kernel');

% 2.2

figure(4);

S = [1 0 -1; 2 0 -2; 1 0 -1];
subplot(1,2,1);
imshow(conv2(I,S));
title('Sobel Filter (Horizontal)');

subplot(1,2,2);
imshow(conv2(I,S'));
title('Sobel Filter (Vertical)');

figure(5);
imshow(conv2(I,S+S'));