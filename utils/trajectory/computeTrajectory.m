function [u0,x0] = computeTrajectory(sys,X0,xyPath,options,varargin)

%TODO: most of this lives in two places-move this all to a common class method to call.

n = sys.sysparams.n;

type = 'output';
odeOptions = odeset('Events',@eventsReachedWaypoint);

if ~isempty(varargin)
    if ~isempty(varargin{1})
        type = varargin{1};
    end
    
    if length(varargin) > 1
        if varargin{2} == 'RRT'
            odeOptions = odeset('Events',@eventsPathLength);
            global xNormCum stepSizeGlob yLast %TODO: how to remove this dependency???

            xNormCum = 0;
            stepSizeGlob = varargin{3};
            yLast = sys.state2SEconfig([],X0,[]);
        end
    end
end

trimTraj = true;
forcedEndIndx = Inf;

Tfin = 200;
t = 0:options.TstepTraj:Tfin;

% Model and Controller 

% TODO: move all of this to a 'simulate' method in systemdynamics
global gotopt xyPathGlob sysGlob %TODO: remove this global!!
gotopt = 1;         % Initialize waypoint index
xyPathGlob = xyPath;
sysGlob = sys;

figure(500)
plot(xyPath(end,1),xyPath(end,2),'o')
% [t,Xk] = ode45(@(tt,X) sys.dynamics(tt,X,sys.sysparams,xyPath,gotopt), t, X0, odeOptions);
[t,Xk] = ode45(@(tt,X) sys.dynamicsWaypointSteering(tt,X,xyPath,gotopt), t, X0, odeOptions);

for indx = 1:n
    if sys.sysparams.isCyclic(indx)
        % Wrap theta from -pi to pi
        Xk(:,indx) = mod(Xk(:,indx)+pi,2*pi)-pi;
    end
end

% X = Xk(:,1); Y = Xk(:,2); 
if strmatch(type,'state');
    outputk = Xk;
elseif strmatch(type,'output')
    outputk = sys.state2SEconfig([],Xk(:,1:2),[]);
end

% Get the command input profile
gotopt = 1;         % Initialize waypoint index
Uk(1:size(Xk,1),1:sys.sysparams.m) = 0;
for k = 1:1:size(Xk,1)
    [U,gotopt] = sys.steerToXYWaypoints(Xk(k,:),xyPath,gotopt);
    Uk(k,:) = U;
end

% [~,idxNonRegStates] = getRegNonRegStates(sys);
% nonRegStatek = Xk(:,idxNonRegStates);

% Remove data beyond when the final waypoint has been achieved
% TODO: Instead of acceptance radius, make this dependent on the prescribed sets in the algorithm!
% dist2LastPt = sqrt((X - xyPath(end,1)).^2 + (Y - xyPath(end,2)).^2);

dist2LastPt = sqrt(sum((outputk - repmat(xyPath(end,:),size(outputk,1),1)).^2,2));

%indx = find(abs(dist2LastPt) < abs(distAccept),1,'first');
dist2LastPt(size(Xk,1)+1) = inf;
indx = size(Xk,1);
if trimTraj
    for k = 1:size(Xk,1)
        if (abs(dist2LastPt(k+1)) > abs(dist2LastPt(k)) && dist2LastPt(k) < sys.sysparams.closeEnough) || k == forcedEndIndx
            indx = k;
            break
        end
    end
else
    indx = size(Xk,1);
    for k = 1:size(Xk,1)
        if abs(dist2LastPt(k)) > 1e6
            indx = k-1;
            break
        end
    end
end


if ~isempty(indx)
    t(indx+1:end) = [];
    Xk(indx+1:end,:) = [];
    Uk(indx+1:end,:) = [];
end

u0 = Traject(t',Uk');
x0 = Traject(t',Xk');


end

function [value,isterminal,direction] = eventsReachedWaypoint(t,x)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.

global xyPathGlob sysGlob

y = sysGlob.state2SEconfig([],x,[]);

desXR = xyPathGlob(end,1);  desYR = xyPathGlob(end,2);

euclidDist2NextPt = sqrt((y(1) - desXR)^2 + (y(2) - desYR)^2);

% If robot is within acceptance radius, index to next waypoint
value = double(abs(euclidDist2NextPt) - sysGlob.sysparams.closeEnough);

isterminal = 1; % stop the integration
direction = 0; % negative direction

end

function [value,isterminal,direction] = eventsPathLength(t,x)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.

global xNormCum sysGlob stepSizeGlob yLast

y = sysGlob.state2SEconfig([],x,[]);

xNormCum = norm(y(1:2) - yLast(1:2)) + xNormCum;

% If robot is within acceptance radius, index to next waypoint
value = xNormCum - stepSizeGlob;

isterminal = 1; % stop the integration
direction = 1; % positive direction

yLast = y;

end