clear all; close all; clc

% some constants
L = 2*pi;
N = 128;
D = 0.1;


tend = 1;



% Creating the differentiation matrix
% First we create the ones above the diagonal
diag1 = diag(ones(N-1,1),1);
% Then the main diagonal (-2)
diag2 = diag(-2*ones(N,1));
% And the ones below the diagonal
diag3 = diag(ones(N-1,1),-1);


Ns = 10:110;
Ns = Ns.^2;
dts = linspace(dt0,dt0/50,100);
ts = zeros(size(dts));
errdts = zeros(length(x),length(dts));
for j = 1:length(dts)
    
    x = linspace(-L/2,L/2,N);
    h = abs(x(1)-x(2));

    dt0 = h^2/(2*D);
    
    dt = dts(j);
    nt = ceil(tend/dt);
    t = zeros(1,nt);


    C = zeros(N,nt); % Rows are spatial, coloumns are temporal

    if mod(N,2) == 0
        C(N/2:N/2+1,1) = 1/(2*h);
    else
        C(ceil(N/2),1) = 1/h;
    end

    T = D*dt/(h^2)*(diag1+diag2+diag3);
    

    % C(i,t+dt) = C(i,t)+T*C(:,t)
    AnalDelta = @(x,t) 1/(sqrt(2*D*t*2*pi))*exp(-(x.^2)/(4*D*t));
    CAnal = C;
    res = C-CAnal;


    for i = 2:nt
        t(i) = t(i-1)+dt;
        C(:,i) = C(:,i-1)+T*C(:,i-1);
        CAnal(:,i) = AnalDelta(x,t(i));
        res(:,i) = C(:,i)-CAnal(:,i);
    end
    ts(j) = t(end);
    errdts(:,j) = res(:,end);
end

figure
ax = gca;
surf(x,dts,abs(errdts'))
ax.ZScale = 'log';
ax.XLim = [min(x),max(x)];
ax.YLim = [min(dts),max(dts)];