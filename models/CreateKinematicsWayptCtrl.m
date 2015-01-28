function [Xdot] = CreateKinematicsWayptCtrl(t,X)

global e gotopt xyPath closeEnough

x = X(1); y = X(2); th = X(3);

% xyPath
% gotopt
if gotopt <= size(xyPath,1)
    desXR = xyPath(gotopt,1);  desYR = xyPath(gotopt,2);
else
    desXR = xyPath(end,1);  desYR = xyPath(end,2);
end
    
dist2NextPt = sqrt((x - desXR)^2 + (y - desYR)^2);

% If robot is within acceptance radius, index to next waypoint
if abs(dist2NextPt) < abs(closeEnough)
    gotopt = gotopt + 1;
%     disp(['gotopt = ',num2str(gotopt)])
end

Vx = desXR - x;
Vy = desYR - y;
[uV, uW] = feedbackLin(Vx, Vy, th, e);

% Governing equations
xdot = uV*cos(th);
ydot = uV*sin(th);
thdot = uW;

Xdot = [xdot; ydot; thdot];
