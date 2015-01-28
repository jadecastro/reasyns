function [xdot] = uniModelLin(t,x)
% 

global psi_dot k1 k2 k3 k4

theta = x(1);
theta_dot = x(2);
alpha = x(3);
alpha_dot = x(4);

%tau = 0;
tau = k1*theta + k2*alpha + k3*theta_dot + k4*alpha_dot;

% linear model
theta_ddot = 1/1.5*theta - 5/1.5*tau;
alpha_ddot = -13/3*theta + 51/3*tau;

% theta_dot = theta_dot;
% theta_ddot = alpha;
% alpha_dot = alpha_dot;
% alpha_ddot = -0.08*theta - 0.06*theta_dot - 0.61*alpha - 0.2*alpha_dot;

xdot = [theta_dot; theta_ddot; alpha_dot; alpha_ddot];

