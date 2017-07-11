clear all; close all; clc

% New differential equaiton:
% dudt = Sdelta(x-x0) + D d^2udx^2-u/tau;

% some constants
L = 2*pi;
N = 129;
D = 0.1;
tend = 10;
dtMult = 2;
S = 1;
tau = 1;


% Creating the differentiation matrix
% First we create the ones above the diagonal
diag1 = diag(ones(N-1,1),1);
% Then the main diagonal (-2)
diag2 = diag(ones(N,1));
% And the ones below the diagonal
diag3 = diag(ones(N-1,1),-1);



x = linspace(-L/2,L/2,N);
h = abs(x(1)-x(2));

dt = h^2/(dtMult*2*D);

nt = ceil(tend/dt);
t = zeros(1,nt);


C = zeros(N,nt); % Rows are spatial, coloumns are temporal
res = C;

% if mod(N,2) == 0
% 	C(N/2:N/2+1,1) = 1/(2*h);
% else
	C(ceil(N/2),1) = 1/h;
% end

% Does all the differentiation. (1-dt/tau)*u+(d^2/dt^2)u.
T = (1-dt/tau-2*D*dt/h^2)*diag2+D*dt/(h^2)*(diag1+diag3);

figure
plot(x,res(:,1))
drawnow

CSteady = 1/h*exp(-abs(x)'/(sqrt(tau*D)));
res(:,1) = C(:,1)-CSteady;

for i = 2:nt
	
	t(i) = t(i-1)+dt;
	C(:,i) =T*C(:,i-1);
	% Let's add the constant source
	if mod(N,2) == 0
		C(N/2:N/2+1,i) = 1/h;
	else
		C(ceil(N/2),i) = 1/h;
	end
	res(:,i) = C(:,i)-CSteady;
	plot(x,res(:,i))
	drawnow
end

figure
hold on
plot(x,CSteady)
plot(x,C(:,end))
legend('Steady-state','Numerical')

 
% figure
% surf(x,t,C');
% xlabel('x')
% ylabel('t')
% shading interp