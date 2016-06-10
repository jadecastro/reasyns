function x0 = computeOpenLoopTrajectory(sys,X0,U0,stepSize,options)

clearvars -global U0
global U0

n = sys.sysparams.n;

Tfin = 100;
t = 0:options.TstepTraj:Tfin;

% Model and Controller

% TODO: move all of this to a 'simulate' method in systemdynamics
global xNormCum sysGlob stepSizeGlob yLastGlob

xNormCum = 0;
sysGlob = sys;
stepSizeGlob = stepSize;
yLastGlob = sys.state2SEconfig([],X0,[]);

odeOptions = odeset('Events',@events);

[t,Xk] = ode45(@(tt,X) sys.dynamics(tt,X,U0), t, X0, odeOptions);
% Uk = repmat(U0,length(t),1);

for indx = 1:n
    if sys.sysparams.isCyclic(indx)
        % Wrap any cyclic states from -pi to pi
        Xk(:,indx) = mod(Xk(:,indx)+pi,2*pi)-pi;
    end
end

% X = Xk(:,1); Y = Xk(:,2);

%u0 = Traject(t',Uk');
x0 = Traject(t',Xk');

end

function [value,isterminal,direction] = events(t,x)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.

global xNormCum sysGlob stepSizeGlob yLastGlob

y = sysGlob.state2SEconfig([],x,[]);

xNormCum = norm(y(1:2) - yLastGlob(1:2)) + xNormCum;

% If the trajectory is sufficiently long, index to the next waypoint
value = xNormCum - stepSizeGlob;

isterminal = 1; % stop the integration
direction = 1; % positive direction

yLastGlob = y;

end
