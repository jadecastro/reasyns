function [ac,redindx] = computeInvariant(u0,x0,sys,Qf);

x = msspoly('x',3);

[xRed,t] = downsampleUniformly(x0);
[uRed,t] = downsampleUniformly(u0);

xtraj = PPTrajectory(foh(tRed',xRed')); % should we be using xRed, uRed here?
utraj = PPTrajectory(foh(tRed',uRed(:,2)'));

% Declare Dubins car model
p = DubinsPlant();

% Set input limits
p = setInputLimits(p,-Inf,Inf);

% Do tvlqr
Q = eye(3);
R = 1;
% Qf = 40*Q;  % in drake, rho is apparently 1 at end of trajectory -- todo: can adjust this to get desired containment

utraj = setOutputFrame(utraj,p.getInputFrame);
xtraj = setOutputFrame(xtraj,p.getStateFrame);
[c,V] = tvlqr(p,xtraj,utraj,Q,R,Qf);
poly = taylorApprox(feedback(p,c),xtraj,[],3);
ts = xtraj.getBreaks();

% Options for funnel computation
options = struct();
options.rho0_tau = 2; % Determine initial guess for rho
options.max_iterations = 5; % Maximum number of iterations to run for
options.stability = false;

% Do funnel computation
Vtraj = sampledFiniteTimeVerification(poly,xtraj.getBreaks(),Qf,V,options);
disp('done');

% Convert V back to state frame
Vxframe = Vtraj.inFrame(p.getStateFrame());

for i = 1:length(tRed)
    S0 = ppval(V.S.pp,tRed(i));
    S1 = ppval(Vxframe.S.pp,tRed(i));
    tmp = S0./S1;
    rhok(i,:) = tmp(1,1);
end

% upsample the trajectories and populate the data structs (todo: update
% with new reasyns classes)
rhopp = interp1(tRed,rhok,'linear','pp');
rhoup = ppval(rhopp,tk);
xpp = spline(tk,xk');
xup = ppval(xpp,tk);
upp = spline(tk,uk');
uup = ppval(upp,tk);
for i = 1:length(tk)
    Kup(:,:,i) = [0 0 0; ppval(c.D.pp,tk(i))]; 
    Pup(:,:,i) = ppval(Vxframe.S.pp,tk(i));
    PupInvTmp = inv(Pup(:,:,i))/rhoup(i);
    PupInvTmp = (PupInvTmp'+PupInvTmp)/2;  % to ensure it is symmetric
    Pup(:,:,i) = inv(PupInvTmp);
end
funnel.t = tk;
funnel.x = xup';
funnel.u = uup';
funnel.P = Pup;
funnel.K = Kup;
funnel.rho = rhoup;
for i = 1:length(tk)
    Vquad(1,i) = (x - xup(:,i))'*Pup(:,:,i)/rhoup(i)*(x - xup(:,i));
end
funnel.V = Vquad;



ac = resample(ac); % resample back to the original resolution


% plot the projection
figure(20)
hold on
options.inclusion = 'projection';
options.plotdims = [1 2];
plotFunnel(Vxframe,options);
fnplt(xtraj,[1 2]);
axis equal


% Tests to make sure simulated trajectories stay inside computed funnel
doTest = 0;
if doTest
    
    % Create closed loop system with optimized controller
    sysCl = feedback(p,c);
    V0 = Vtraj.getPoly(0); % Get inlet of funnel
    xinit = getLevelSet(decomp(V0),V0,[0;0;0]);
    
    Vsall = [];
    figure(30)
    hold on
    grid on
    for j = 1:5
        Vs = [];
        x0 = 0.95*xinit(:,j) + xtraj.eval(0); % Simulate from 0.95*boundary of funnel
        xsim = sysCl.simulate([0 ts(end)],x0);
        for k = 1:length(ts)
            Vs = [Vs, Vxframe.eval(ts(k),xsim.eval(ts(k)))];
        end
        Vsall = [Vsall, Vs];
        
        plot(ts,Vs)
        plot(ts,ones(1,length(ts)),'ro')
        drawnow;
    end
    
    if ~(all(Vsall < 1))
        success = false;
        disp('A Trajectory left the funnel!!')
    else
        success = true;
        disp('All simulated trajectories stayed inside funnel.')
    end
    
end

