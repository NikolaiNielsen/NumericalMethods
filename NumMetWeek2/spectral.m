clear all; close all; clc;

L = 2*pi;
N = 128;
x = linspace(-L/2,L/2,N);

f = cos(2*x);

fk = fft(f);
k = kvalues(N); % MATLAB magic....

% Let's calculate some derivatives
% Calculating the derivatives in k-space
dk = 1i*k.*fk;
dk2 = -k.^2.*fk;

% converting to real space
df = real(ifft(dk));
df2 = real(ifft(dk2));

% stuff for plots
deriv = df+df2;
anal = -2*sin(2*x)-4*cos(2*x);
res = deriv-anal;


% plots
figure
hold on
% plot(x,f)
plot(x,deriv)
plot(x,anal)
hold off
legend('numerical derivative','analytical derivative')

figure
plot(x,res)