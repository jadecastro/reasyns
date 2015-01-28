function [U,gotopt] = CarKinematicsNLWayptCtrl1_ctrl(X,gotopt,l,e,xyPath,closeEnough)

x = X(1); y = X(2); th = X(3); phi = X(4);

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

Vx = desXR - x;% + 1e-9;
Vy = desYR - y;% + 1e-9;
[uV, uW] = feedbackLinCar(Vx, Vy, th, phi, l, e);

% Limit on velocity
% uV = 2; %sign(uV)*max(abs(uV),0.01);  % Limits velocity on (1,0)
%uW = max(min(uW,3),-3);  % Limits angular velocity on (-1,1)

U = [uV; uW];
