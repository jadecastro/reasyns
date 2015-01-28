function [Xdot] = CreateKinematicsNLWayptCtrl1(t,X)

global e gotopt xyPath closeEnough

x = X(1); y = X(2); th = X(3);

% xyPath
% gotopt

[U,gotopt] = CreateKinematicsNLWayptCtrl1_ctrl(X,gotopt,e,xyPath,closeEnough);

uV = U(1); uW = U(2);

% Governing equations
xdot = uV*cos(th);
ydot = uV*sin(th);
thdot = uW;

Xdot = [xdot; ydot; thdot];
% if gotopt < length(xyPath(:,1))+1
%     Xdot = [xdot; ydot; thdot];
% else
% %     xyPath
% %     dist2NextPt
%     Xdot = [0;0;0];
% end
