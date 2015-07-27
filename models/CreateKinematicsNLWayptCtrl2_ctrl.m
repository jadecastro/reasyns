function [U,gotopt] = CreateKinematicsNLWayptCtrl2_ctrl(X,params,xyPath,gotopt)

e = params.e;
closeEnough = params.closeEnough;

x = X(1); y = X(2); th = X(3);

if gotopt <= size(xyPath,1)
    desXR = xyPath(gotopt,1);  desYR = xyPath(gotopt,2);
else
    desXR = xyPath(end,1);  desYR = xyPath(end,2);
end
    
euclidDist2NextPt = sqrt((x - desXR)^2 + (y - desYR)^2);

% If robot is within acceptance radius, index to next waypoint
if abs(euclidDist2NextPt) < abs(closeEnough)
    gotopt = gotopt + 1;
%     disp(['gotopt = ',num2str(gotopt)])
end

Vx = desXR - x;
Vy = desYR - y;
[uV, uW] = feedbackLin(Vx, Vy, th, e);

% Limit on velocity
uV = 0.05;%max(uV,0);  % Limits velocity on (1,0)
uW = max(min(uW,0.25),-0.25);  % Limits angular velocity on (-1,1)

U = [uV; uW];
