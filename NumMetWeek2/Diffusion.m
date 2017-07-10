clear all; close all; clc

% some constants
L = 2*pi;
N = 128;
D = 0.1;

x = linspace(-L/2,L/2);
h = abs(x(1)-x(2));

tend = 10;
dt = h^2/(10*2*D);
nt = ceil(tend/dt);
t = zeros(1,nt);


C = zeros(N,nt); % Rows are spatial, coloumns are temporal

if mod(N,2) == 0
    C(N/2:N/2+1,1) = 1/(2*h);
else
    C(ceil(N/2),1) = 1/h;
end