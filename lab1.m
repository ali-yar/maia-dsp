% Exercice 1 - Deterministic signals
clear; clc;
figure(1);
n = -5:1:10;

%1.1.a
subplot(3,1,1)
y = zeros(16,1);
y(6) = 1;
stem(n,y); % Dirac function plotted
xlim([-6 11]);
ylim([0 1.1]);
title('Dirac function');
xlabel('Time (n)');
ylabel('Amplitude');

%1.1.b
subplot(3,1,2);
y = [zeros(5,1);ones(11,1)];
stem(n,y); % Unit step function plotted
xlim([-6 11]);
ylim([0 1.1]);
title('Unit step function');
xlabel('Time (n)');
ylabel('Amplitude');

%1.1.c
subplot(3,1,3);
y = [zeros(5,1);ones(11,1)];
y = n .* y';
stem(n,y); % Ramp function plotted
xlim([-6 11]);
ylim([0 10.1]);
title('Ramp function');
xlabel('Time (n)');
ylabel('Amplitude');

% 1.2
clear; clc;
n = 0:1:100;

figure(2);
% 1.2.a
a = -1/12;
w = j * (pi/6);
y = exp((a + w) * n);

subplot(2,1,1);
stem(n,real(y),'o'); % real part of complex exponential plotted
xlim([-1 101]);
ylim([-1.1 1.1]);
grid on;
title('Real part');
xlabel('Time (n)');
ylabel('Amplitude');

subplot(2,1,2);
stem(n,imag(y),'o'); % imaginary part of complex exponential plotted
xlim([-1 101]);
ylim([-1.1 1.1]);
grid on;
title('Imaginary part');
xlabel('Time (n)');
ylabel('Amplitude');

figure(3);
% 1.2.b
a = -0.4;
w = j * (pi/5);
y = 2.5 * exp((a + w) * n);

subplot(2,1,1);
stem(n,real(y),'o'); % real part of complex exponential plotted
xlim([-1 101]);
ylim([-0.6 2.8]);
grid on;
title('Real part');
xlabel('Time (n)');
ylabel('Amplitude');

subplot(2,1,2);
stem(n,imag(y),'o'); % imaginary part of complex exponential plotted
xlim([-1 101]);
ylim([-0.2 1.2]);
grid on;
title('Imaginary part');
xlabel('Time (n)');
ylabel('Amplitude');

% 1.3
clear; clc;
figure(4);
n = 0:39;
sf = 100; % sampling frequency
ts = 1 / sf; % sampling interval
f = 10; % frequency
y = sin(2*pi*f*n*ts);
stem(n, y); % signle period plotted
grid on;
title('Single period signal f=10 Hertz Ts=1/100');
xlabel('Time (n)');
ylabel('Amplitude');


% Exercice 2 - Random signals

% 2.1
clear; clc;
figure(5);
n = 0:999;
N = length(n);
y = -1 + (1 - (-1)) * rand(1,N); % generate N random numbers from -1 to 1
stem(n,y); % uniformly distributed sequence plotted
title('Uniform random sequence');
xlabel('Time (n)');
ylabel('Amplitude');

y_mean = mean(y);
y_variance = sum((y - y_mean).^2) / N;
% mean by Matlab =  -0.010299, theoritical mean = -0.010299057
% variance by Matlab = 0.343731, theoritical variance = 0.344074664
% The difference between Matlab and the theoritical calculation is due to
% matlab using floating-point arithmetic which may occasionaly involve
% round-Off errors.

% 2.2
clear; clc;
figure(6);
n = 0:199;
N = length(n);
y_mean = 0;
y_variance = 5;
y_std = sqrt(y_variance);
y = normrnd(y_mean,y_std,1,N); % random samples from the normal distribution
stem(n,y); % Gaussian random signal plotted
title('Gaussian random signal');
xlabel('Time (n)');
ylabel('Amplitude');
%hist(y,-10:.02:10) ;