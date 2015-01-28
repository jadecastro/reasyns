function [xdot] = uniModel(t,x)
% 
% x = [theta_dot, theta, psi_dot, psi, alpha_dot, alpha]'
% u = [tau];

global psi_dot k1 k2 k3 k4

theta = x(1);
theta_dot = x(2);
phi = x(3);
phi_dot = x(4);
alpha = x(5);
alpha_dot = x(6);


%tau = 0;
tau = k1*theta + k2*alpha + k3*theta_dot + k4*alpha_dot;

K = 307/8 + 3/2*alpha^2 - 303/8*cos(2*theta) + 15*alpha*sin(2*theta) + 3/2*alpha^2*cos(2*theta);

theta_ddot = 1/1.5*sin(theta) + 3/8/1.5*phi_dot^2*sin(2*theta) - 3/4/1.5*phi_dot*psi_dot*cos(theta) + ...
    7.5/1.5*phi_dot^2*alpha*cos(2*theta) - 3/1.5*alpha*cos(theta) - 3/2/1.5*phi_dot^2*alpha^2*sin(2*theta) + ...
    3/1.5*phi_dot*psi_dot*alpha*sin(theta) - 15/2/1.5*phi_dot^2*alpha - 5/1.5*tau;

phi_ddot = 63/4/K*psi_dot*phi_dot*cos(theta) + 3/K*psi_dot*alpha_dot*cos(theta) - 303/4/K*phi_dot*theta_dot*sin(2*theta) - ...
    15/K*phi_dot*alpha_dot*sin(2*theta) - 3/K*phi_dot*alpha_dot*alpha + 3/K*phi_dot*theta_dot*alpha^2*sin(2*theta) - ...
    3/K*psi_dot*theta_dot*alpha*sin(theta) - 30*phi_dot*psi_dot*alpha*cos(2*theta) - 3*phi_dot*alpha_dot*alpha*cos(2*theta);

alpha_ddot = -13/3*sin(theta) + 3.75/3*phi_dot^2*sin(2*theta) + 4.5/3*phi_dot*psi_dot*cos(theta) - ...
    73.5/3*phi_dot^2*alpha*cos(2*theta) + 75.75/3*phi_dot^2*alpha + 30/3*alpha*cos(theta) + ...
    15/3*phi_dot^2*alpha^2*sin(2*theta) - 30/3*phi_dot*psi_dot*alpha*sin(theta) + 51/3*tau;

% Energy shaping
g = 1;
md = 1;
mc = 3;
R = 1;
l = 4;
Ja = md*R^2/2;
Jt = md*R^2/4;

k = 10000000;
lambda = -0.00000001;

Meq = Ja + md*R^2 + k*sin(theta);

%v_des = 0;
v_c = 51/10*(1/Meq*(0.5*k*cos(theta)*theta_dot^2 - g*(md*R+mc*(R+l))*sin(theta) + mc*g*alpha*cos(theta)) - ...
    0.1830*sin(theta) + 0.4951*phi_dot^2*sin(2*theta) - 0.2059*phi_dot*psi_dot*cos(theta) + ...
    0.1961*phi_dot^2*alpha*cos(2*theta) - 0.0392*alpha*cos(theta) - 0.0196*phi_dot^2*alpha^2*sin(2*theta) + ...
    0.0392*phi_dot*psi_dot*alpha*sin(theta) - 0.049*phi_dot^2*alpha);

v_d = -lambda*Meq/(Ja + md*R^2)*(alpha_dot + k*theta_dot*sin(theta));
v = v_c + v_d;

theta_ddot = -10/51*v - 0.1830*sin(theta) + 0.4951*phi_dot^2*sin(2*theta) - 0.2059*phi_dot*psi_dot*cos(theta) + ...
    0.1961*phi_dot^2*alpha*cos(2*theta) - 0.0392*alpha*cos(theta) - 0.0196*phi_dot^2*alpha^2*sin(2*theta) + ...
    0.0392*phi_dot*psi_dot*alpha*sin(theta) - 0.049*phi_dot^2*alpha;

alpha_ddot = v;

% linear model
% theta_ddot = 1/1.5*theta - 5/1.5*tau;
% phi_ddot = 0;
% alpha_ddot = -13/3*theta + 51/3*tau;

xdot = [theta_dot; theta_ddot; phi_dot; phi_ddot; alpha_dot; alpha_ddot];

