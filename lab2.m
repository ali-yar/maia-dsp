% Reminder 1
clear; clc;
figure(1);

subplot(2,1,1);
t = -20:0.01:20;
f = 1; % frequency
x = sin(2*pi*f*t);
plot(t, x);
title('1Hz Sin function');
xlabel('t');
ylabel('Amplitude');

subplot(2,1,2);
n = -20:20;
fs = 20; % sampling frequency
f = 1; % frequency
x = sin(2*pi*(f/fs)*n);
stem(n, x);
grid on;
title('1Hz Sin function sampled at 20Hz');
xlabel('n');
ylabel('Amplitude');


% Exercise 1 - Causality
clear; clc;
figure(2);

subplot(3,1,1);
s = 0.01;
t = -5:s:10;
x = zeros(size(t));
x(t>=4) = 1; % x(k) = H(k-4)
plot(t, x);
title('x(k) = H(k-4)');
xlim([-5 10]);
ylim([0 1.1]);
xlabel('t');
ylabel('Amplitude');

subplot(3,1,2);
x2 = zeros(size(t));
x2(t>=3) = 1; % x(k+1) = H(k-3)
y = (x + x2) / 2.;
plot(t, y);
title('y(k) = (x(k) + x(k+1)) / 2 (non-causal)');
xlim([-5 10]);
ylim([0 1.1]);
xlabel('t');
ylabel('Amplitude');
% It can be seen that there are output values at k<4, however the input values
% are 0 at k<4. So, the system is non-causal.

subplot(3,1,3);
x3 = zeros(size(t));
x3(t>=5) = 1; % x(k-1) = H(k-5)
y = (x + x3) / 2.;
plot(t, y);
title('y(k) = (x(k) + x(k-1)) / 2 (causal)');
xlim([-5 10]);
ylim([0 1.1]);
xlabel('t');
ylabel('Amplitude');
% By replacing x(k+1) with x(k-n) such that n>=0, we obtain a causal system,
% because now, the output value at time k only depends on input values up to time k.  

% Exercise 2 - Stability
clear; clc;

% Defining the accumulator system
% Input:  - x: list
%         - N: (optional) length of ouptut. Default is length of x
% Output: - y: list
function [y] = prim (x, N)
  if ~exist('x')
		error('The input sequence is missing.')
	end
  if ~exist('N')
    N = length(x);
  end
  y(1) = 0 + x(1); 
  for k = 2:N 
    y(k) = y(k-1) + x(k);
  end
endfunction

% 2.1
figure(3);
n = 0:10;
N = length(n);
x = zeros(1,N);
x(5:end) = 1;
stem(n, prim(x));
title('Accumulator system with step signal H(k-4) (unstable)');
xlabel('n');
ylabel('Amplitude');
% The output of the primitive operator is not bounded.
% Therefore, the primitive operator is not stable.

% 2.2
figure(4);
n = -5:5;
N = length(n);
x = zeros(1,N);
x(6) = 1;
stem(n, prim(x));
title('Impulse response of Accumulator system');
xlabel('n');
ylabel('Amplitude');
% The impulse response of the primitive operator is the unit step signal

% 2.3
figure(5);
n = -5:5;
N = length(n);
x = zeros(1,N);
x(6) = 1;
y = [];
y(1) = 0 + x(1); 
for k = 2:N 
  y(k) = x(k) + 2 * y(k-1);
end
stem(n, y);
title('Impulse response of y[n] = x[n] + 2*y[n-1] (unstable)');
xlabel('n');
ylabel('Amplitude');
% The system y[n] = x[n] + 2*y[n-1] is not stable because the output is unbounded
% for a least one bounded input such as the dirac signal.

% 2.4
figure(6);
n = -5:5;
N = length(n);
x = zeros(1,N);
x(6) = 1;
y = [];
y(1) = 0 + x(1); 
for k = 2:N 
  y(k) = x(k) + ( y(k-1) / 3 );
end
stem(n, y);
title('Impulse response of y[n] = x[n] + y[n-1] / 3 (stable)');
xlabel('n');
ylabel('Amplitude');
% The system y[n] = x[n] + y[n-1] / 3 is stable because its impulse response
% is absolutely summable

% Exercise 3 - Invariance and linearity
clear; clc;
figure(7);

% 3.1
xa = [0 0 0 0 1 2 3 4 5 0 0 0 0 0 0 0 0 0 0];
xb = [0 0 0 0 0 0 0 0 0 4 3 2 1 0 0 0 0 0 0];

ya = 3*[0 xa(1:end-1)] - 2*xa + [xa(2:end) 0];
yb = 3*[0 xb(1:end-1)] - 2*xb + [xb(2:end) 0];

x = xa + xb;
y = 3*[0 x(1:end-1)] - 2*x + [x(2:end) 0];

% 3.2
% If y is equal to (ya + yb) then the system is linear
n = 0:length(y)-1;
stem(n,y,'color','r');
hold on;
stem(n+0.05,ya + yb,'color','b'); % shifting the signal a little bit to see difference
title('y = 3*x[k-1] - 2*x[k] + x[k+1] (linear)');
xlabel('n');
ylabel('Amplitude');
% the plot of y is on top of the plot of (ya + yb), so the system is linear

