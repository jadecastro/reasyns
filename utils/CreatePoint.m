
clear all; clc; 
%close all

global m J Jw b bu r rw mu muk mu2 muk2 g 
global xd yd thd ud wd T T1 T2 taun Bsa Pia Ka1 Ka2 cycle
global thd_tau tau_st tau_end
global q locationhit tau
global x_look e

% Constants

m = 0.5;
b = 0.0005;    % wheel damping (rotational units)
bu = 1;     % damping from ground
r = 0.2;
J = 1/2*m*r^2;
Jw = 0.0005;
rw = 0.01;  
mu = 0.9;  % friction coeff
muk = 0.8*mu;  % kinetic friction coeff
g = 9.81;

% mu duty cycle

mu2 = 0.02;
muk2 = 0.8*mu2;
tau_st = -0.25;
tau_end = tau_st - 0.05;

% Reference trajectory - trace out a circle

T = 0.01;
Tfin = 5;
t = 0:T:Tfin;

% Controller
x_look = 1;
e = 0.2;

udlen = length(ud);


% Initial conditions
u0 = 0;
v0 = 0;
w0 = 0;
x0 = 0;
y0 = 3;
th0 = 1*pi/4;

wl0 = 1/rw*u0 - r/2/rw*w0;
wr0 = 1/rw*u0 + r/2/rw*w0;


X0 = [wl0; wr0; u0; v0; w0; x0; y0; th0];


%% Kinematics model with transverse control

[t,Xk] = ode45(@ CreateKinematicsPointCtrl, t, X0(6:8));

X = Xk(:,1); Y = Xk(:,2); Th = Xk(:,3);

% determine the control inputs indirectly from the simulation
for i = 1:length(Xk)-1
    Xdot(i) = 1/(t(i+1)-t(i))*(X(i+1) - X(i));
    Ydot(i) = 1/(t(i+1)-t(i))*(Y(i+1) - Y(i));
    Thdot(i) = 1/(t(i+1)-t(i))*(Th(i+1) - Th(i));
end
uV = sqrt(Xdot.^2 + Ydot.^2);
uW = Thdot;

%% Kinematics model with transverse control - Polynomial approximations

[tpoly,Xkpoly] = ode45(@ CreateKinematicsPointCtrlPoly, t, X0(6:8));

Xpoly = Xkpoly(:,1); Ypoly = Xkpoly(:,2); Thpoly = Xkpoly(:,3);

%% Dynamic Model with transverse control - No Slip

% cycle = 0;
% thd_tau = 0;
% 
% [t,Xk] = ode45(@ CreateDynamicsCircleTverseCtrl, [0:0.01:Tfin], X0([3,5:8]));
% 
% X = Xk(:,3); Y = Xk(:,4); Th = Xk(:,5);


%% Hybrid Dynamic Model with transverse control - With Slip

% taun0 = atan2(y0,x0)/2/pi;
% q0 = 1;  % initially no slip
% q = q0;
% 
% cycle = 0;
% thd_tau = 0;
% locationhit = 1;
% 
% odeset('Mass',[eye(8) zeros(8,2);zeros(2,10)]);  % kludge to get q out of ODE function
% odeset('RelTol',1e-10);
% % odeset('AbsTol',[1e-12*ones(8,1); 1e0]);
% [t,Xk] = ode23s(@ CreateHybridDynamicsCircleTverseCtrl, [0:T/10:Tfin], [X0; q0; taun0]);
% % [t,Xk] = ode23s(@ CreateHybridDynamicsCircleTverseCtrl, [0:T/10:Tfin], [X0; q0]);
% 
% X = Xk(:,6); Y = Xk(:,7); Th = Xk(:,8); Qq = Xk(:,9);% Taun = Xk(:,10);
% td = t;
% 
% save CreateSim


%% Plotting

figure(1), clf

Xlimit = max(X(end)+0.5, 2*Y(1)+0.5);
Ylimit = max(X(end)/2+0.5, Y(1)+0.5);

Xvals = [-0.5 Xlimit];

% plot(Xvals,[0 0],'k--','LineWidth',2);

%axis([-0.5 Xlimit -Ylimit Ylimit])
axis equal
hold on

Xinc = downsample(X,ceil(length(X)/10));
Yinc = downsample(Y,ceil(length(X)/10));
Thinc = downsample(Th,ceil(length(X)/10));
ang = linspace(-pi,pi,100);
P = r*eye(2);
for i = 1:length(Xinc)
    circ = P*[cos(ang); sin(ang)] + [Xinc(i); Yinc(i)]*ones(1,length(ang));
    fill(circ(1,:),circ(2,:),[0.8 0.8 0.2])
    H = line([Xinc(i) Xinc(i) + r*cos(Thinc(i))], [Yinc(i);Yinc(i) + r*sin(Thinc(i))]);
    set(H,'LineWidth',3,'Color',[0.8 0 0])
end

plot(X,Y,'LineWidth',2);
plot(Xpoly,Ypoly,'c','LineWidth',2);

figure(2)
subplot(211)
plot(t,X)
subplot(212)
plot(t,Y)

