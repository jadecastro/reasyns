function [u0,x0] = computeTrajectory(sys,X0,xyPath,varargin)

n = sys.params.n;

if length(sys.H) > 2, error('Workspaces of dimension greater than two not yet supported. Sorry.'); end
if ~isempty(varargin)
    type = varargin{1};
else
    type = 'output';
end

trimTraj = true;
forcedEndIndx = Inf;

TstepTraj = 0.02;
Tfin = 100;
t = 0:TstepTraj:Tfin;

% Model and Controller 

% TODO: move all of this to a 'simulate' method in systemdynamics
global gotopt  %TODO: how to remove this dependency???
gotopt = 1;         % Initialize waypoint index

[t,Xk] = ode45(@(tt,X) sys.mdlFun(tt,X,sys.params,xyPath,gotopt), t, X0);

for indx = 1:n
    if sys.params.isCyclic(indx)
        % Wrap theta from -pi to pi
        Xk(:,indx) = mod(Xk(:,indx)+pi,2*pi)-pi;
    end
end

% X = Xk(:,1); Y = Xk(:,2); 
if strmatch(type,'state');
    outputk = Xk;
elseif strmatch(type,'output')
    outputk = Xk(:,1:length(sys.H))*sys.H;
end
nonRegStatek = Xk(:,length(sys.H)+1:end);

% Remove data beyond when the final waypoint has been achieved
% TODO: Instead of acceptance radius, make this dependent on the prescribed sets in the algorithm!
% dist2LastPt = sqrt((X - xyPath(end,1)).^2 + (Y - xyPath(end,2)).^2);
dist2LastPt = sqrt(sum((outputk - repmat(xyPath(end,:),size(outputk,1),1)).^2,2));

%indx = find(abs(dist2LastPt) < abs(distAccept),1,'first');
dist2LastPt(size(Xk,1)+1) = inf;
indx = size(Xk,1);
if trimTraj
    for k = 1:size(Xk,1)
        if (abs(dist2LastPt(k+1)) > abs(dist2LastPt(k)) && dist2LastPt(k) < sys.params.closeEnough) || k == forcedEndIndx
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

gotopt = 1;         % Initialize waypoint index
for k = 1:1:size(Xk,1)
    [U,gotopt] = sys.ctrlFun(Xk(k,:),sys.params,xyPath,gotopt);
    Uk(k,:) = U;
end

if ~isempty(indx)
    t(indx+1:end) = [];
    Xk(indx+1:end,:) = [];
    Uk(indx+1:end,:) = [];
end

u0 = Traject(t',Uk');
x0 = Traject(t',Xk');


end