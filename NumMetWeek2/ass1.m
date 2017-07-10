clear all; close all; clc

% some constants
L = 2*pi;
N = 128;
D = 0.1;

x = linspace(-L/2,L/2,N);
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


% Creating the differentiation matrix
% First we create the ones above the diagonal
diag1 = diag(ones(N-1,1),1);
% Then the main diagonal (-2)
diag2 = diag(-2*ones(N,1));
% And the ones below the diagonal
diag3 = diag(ones(N-1,1),-1);
T = D*dt/(h^2)*(diag1+diag2+diag3);

% C(i,t+dt) = C(i,t)+T*C(:,t)
figure
plot(x,C(:,1))
ax = gca;
for i = 2:nt
    t(i) = t(i-1)+dt;
    C(:,i) = C(:,i-1)+T*C(:,i-1);
    plot(x,C(:,i))
    ax.YLim = [0 1/(2*h)];
    drawnow
end

surf(x,t,C')