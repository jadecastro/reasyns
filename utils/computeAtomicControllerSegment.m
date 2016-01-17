function [ac,c] = computeAtomicControllerSegment(u0,x0,sys,sampSkip,rho_if,isMaximization)

Q = sys.params.ctrloptions.Q;
R = sys.params.ctrloptions.R;
Qf = sys.params.ctrloptions.Qf;

x = msspoly('x',sys.params.n);

[~,tk] = double(x0);
if nargin > 4
    [xRed,~] = downsampleUniformly(x0,sampSkip);
    [uRed,tRed] = downsampleUniformly(u0,sampSkip);
else
    [xRed,~] = double(x0);
    [uRed,tRed] = double(u0);    
end

xtraj = PPTrajectory(foh(tRed,xRed)); % should we be using xRed, uRed here?
if strfind(func2str(sys.polyMdlFun),'CreateKinematics')
    utraj = PPTrajectory(foh(tRed,uRed(2,:)));
elseif strfind(func2str(sys.polyMdlFun),'Holonomic')
    utraj = PPTrajectory(foh(tRed,uRed));
end

% Declare plant model
p = sys.drakeplant;

% Set input limits
p = setInputLimits(p,-Inf,Inf);

% Do tvlqr
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
options.max = isMaximization;

% Do funnel computation

Vtraj = sampledFiniteTimeVerification(poly,xtraj.getBreaks(),Qf,rho_if,V,options);
% Vtraj = sampledFiniteTimeVerification(poly,xtraj.getBreaks(),Qf,V,options,xMssExt);
disp('done');

% Convert V back to state frame
Vxframe = Vtraj.inFrame(p.getStateFrame());

for i = 1:length(tRed)
    S0 = ppval(V.S.pp,tRed(i));
    S1 = ppval(Vxframe.S.pp,tRed(i));
    tmp = S0./S1;
    rhok(i,:) = tmp(1,1);
end

% upsample the trajectories and populate the data structs 
xup = double(x0);
uup = double(u0);
rhopp = interp1(tRed,rhok,'linear','pp');
rho0 = Traject(rhopp);
rhoup = double(rho0,tk);
%TODO: fix this hack to get Drake data
if verLessThan('matlab','7.15')
    zeromatrix = repmat(zeros(length(utraj),sys.params.n),[1,1,length(tk)]);
else
    zeromatrix = repmat(zeros(length(utraj),sys.params.n),1,1,length(tk));
end
tmp1 = Traject(tk,zeromatrix);
tmp2 = Traject(c.D.pp);
K0 = [tmp1; tmp2];
P0 = Traject(Vxframe.S.pp);
for i = 1:length(tk)
    Kup(:,:,i) = [zeros(1,sys.params.n); ppval(c.D.pp,tk(i))];  % TODO: Need to generalize!! (specific only to the case where v & w are inputs, with v zero'd out)
    Pup(:,:,i) = ppval(Vxframe.S.pp,tk(i));
    PupInvTmp = inv(Pup(:,:,i))/rhoup(i);
    PupInvTmp = (PupInvTmp'+PupInvTmp)/2;  % to ensure it is symmetric
    Pup(:,:,i) = inv(PupInvTmp);
end

for i = 1:length(tk)
    Vquad(1,i) = (x - xup(:,i))'*Pup(:,:,i)/rhoup(i)*(x - xup(:,i));
end
% V0 = Traject(tk,Vquad);

ac = QuadraticAC(x0,u0,K0,P0,rho0,Vquad,sys);

% ac = resample(ac); % resample back to the original resolution


% plot the projection
% figure(5)
% hold on
% options.inclusion = 'projection';
% options.plotdims = [1 2];
% plotFunnel(Vxframe,options);
% fnplt(xtraj,[1 2]);
% axis equal

figure(5)
hold on
plot(ac,sys,5)


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

