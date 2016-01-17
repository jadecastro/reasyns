
function [ac, B] = computeConformingFunnel(u01, x01, sys, regMode, ellToCompose, options)
%
% Main function for computing funnels using the control barrier functions approach.
%
% Elements adapted from Russ Tedrake's Underactuated Robotics Courseware
% https://courses.edx.org/courses/MITx/6.832x

% Initialize SOS program
%close all
clear prog
prog = spotsosprog;

% x = msspoly('x',3);

% load tRed, xRed, and uRed
%load('/home/jon/Dropbox/Research/LowLevelControllerSynthesis/reasyns/scripts/barrierTestData.mat')

Q = sys.params.ctrloptions.Q;
R = sys.params.ctrloptions.R;
Qf = sys.params.ctrloptions.Qf;

[ellq0,ellQ0] = double(ellToCompose);

sampSkip = options.sampSkipFun;
[xRed,~] = downsampleUniformly(x01,10);
[uRed,tRed] = downsampleUniformly(u01,10);

xtraj = PPTrajectory(foh(tRed,xRed)); % should we be using xred, ured here?
utraj = PPTrajectory(foh(tRed,uRed(2,:)));

%[x00Red,~] = downsampleUniformly(x00,sampSkip);
%[u00Red,t00Red] = downsampleUniformly(u00,sampSkip);
[x00,t00] = double(x00);
u00 = double(u00);

x00traj = PPTrajectory(foh(t00,x00)); % should we be using xred, ured here?
u00traj = PPTrajectory(foh(t00,u00(2,:)));

% Declare Dubins car model
if strfind(func2str(sys.polyMdlFun),'CreateKinematics')
    p = DubinsPlant();
    hybridProblem = true;
elseif strfind(func2str(sys.polyMdlFun),'Holonomic')
    p = HolonomicPlant();
    hybridProblem = false;
end

% Set input limits
p = setInputLimits(p,-Inf,Inf);

% Do tvlqr
utraj = setOutputFrame(utraj,p.getInputFrame);
xtraj = setOutputFrame(xtraj,p.getStateFrame);
u00traj = setOutputFrame(u00traj,p.getInputFrame);
x00traj = setOutputFrame(x00traj,p.getStateFrame);
[c,V] = tvlqr(p,x00traj,u00traj,Q,R,Qf);
% poly = taylorApprox(feedback(p,c),xtraj,[],3);
%poly = taylorApprox(feedback(p,c),0,Point(p.getStateFrame(),[finalState';1]),[],6);  % take the Taylor series expansion about the final goal point
poly = taylorApprox(feedback(p,c),x00traj,[],3);
ts00 = x00traj.getBreaks();

ts = xtraj.getBreaks();

num_xc = poly.getNumContStates();
if (isa(V,'Trajectory'))
    V1 = V.inFrame(poly.getStateFrame);
    Q = eye(num_xc);
    V0 = tvlyap(poly,V1,Q,Q);
else
    V0 = V;
end

N = length(ts);
%N = 5;

poly = poly.inStateFrame(V0.getFrame); % convert system to Lyapunov function coordinates

d = 6;


% TODO: the following is specific only to the problem at hand... we will need to generalize this
g_X0 = 1 - (x - ellq0)'*inv(ellQ0)*(x - ellq0);  % intial set
%g_X0 = g_X0.inStateFrame(msubs(V0,x,ts(1)).getFrame);

g_Xdomain = 2^2 - (x(1) + ellq0(1))^2 - (x(2) + ellq0(2))^2; % Conservative domain for finding an initial B

prog = prog.withIndeterminate(x);

%% Find a feasible initial B within some conservative domain
% evaluate dynamics and Vtraj at every ts once (for efficiency/clarity)
for ii = 1:N
    % Initialize barrier function as a free d-th degree polynomial

    [prog,Bnew] = prog.newFreePoly(monomials(x,0:d));
    B{ii} = Bnew;
    
    f{ii} = poly.getPolyDynamics(ts(ii));
    if (poly.getNumInputs>0)   % zero all inputs
        f{ii} = subs(f{ii},poly.getInputFrame.getPoly,zeros(poly.getNumInputs,1));
    end
    g{ii} = [0;0;1];
    
    
    % Generate the safe consraint set as an intersection of the complement
    % of any unsafe sets, including the largest ellipse at this point in
    % the funnel
    [polyH,polyK] = double(regMode.p);
    
    idx00 = find(abs(ts(ii)-ts00) < 1e-5, 1);
    tmpArray = ac.ellipsoid;
    [ellq,ellQ] = double(tmpArray(idx00));
    g_Xu1 = [];
    %g_Xu1 = 1 - (x - ellq)'*inv(ellQ)*(x - ellq);
    if ii == N
        g_Xu1 = 1 - (x - ellq)'*inv(ellQ)*(x - ellq);
    end
    for iplane = 1:length(polyK)%find(isectArray(idx00,:) == 0)
        g_Xu1 = [g_Xu1; polyH(iplane,:)*x(1:2) - polyK(iplane)]; % unsafe set
    end
    %g_Xu1 = g_Xu1.inStateFrame(V0(ii).getFrame);
    
    
    % Compute time derivatives and Lie derivatives of B
    if ii > 1
        dBdt = (B{ii} - B{ii-1})/(ts(ii) - ts(ii-1));  % take an approximate derivative - valid for small delta t's
    else
        dBdt = 0;
    end
    dBdxtimesf{ii} = diff(B{ii},x)*f{ii} + dBdt;
    dBdxtimesg{ii} = diff(B{ii},x)*g{ii};
    
    % SOS constraints
    % C1: B(x) < 0 inside initial condition set
    if ii == 1
        [prog,L1] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
        prog = prog.withSOS(-B{ii} - 0.01*L1*g_X0 - 1e-5); % The 1e-4 is there to make sure the inequality is strict (but you don't have to worry about this for the other constraints)
    else
        tmpArray = ac.ellipsoid;
        [ellq,ellQ] = double(tmpArray(idx00));
        g_X01 = 0.01 - (x - ellq)'*inv(ellQ)*(x - ellq);
        [prog,L1] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
        prog = prog.withSOS(-B{ii} - 1*L1*g_X01 - 1e-5); % The 1e-4 is there to make sure the inequality is strict (but you don't have to worry about this for the other constraints)        
    end
    
    % C2: B(x) >= 0 inside unsafe set
    if ii == N
        weightMultiplier = 0.01;
    else
        weightMultiplier = 0.01;
    end
    [prog,L21] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
    prog = prog.withSOS(B{ii} - weightMultiplier*L21*g_Xu1(1) );
    for iu = 2:length(g_Xu1)
        [prog,L21] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
        prog = prog.withSOS(B{ii} - 0.1*L21*g_Xu1(iu) );
    end
    
    % C3: Bdot <= 0 everywhere
    %     prog = prog.withSOS(-Bdot);
    [prog,Ld0] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
    prog = prog.withSOS(-dBdxtimesf{ii} - 0.01*Ld0*g_Xdomain);
    %prog = prog.withSOS(-Bdot);
    
end

% Solve SOS program (for B with L fixed)
SOSoptions = spot_sdp_default_options();
SOSoptions.verbose = 1;
sol = prog.minimize(0,options.solver,SOSoptions);

% Check if SOS program ran correctly
if ~sol.isPrimalFeasible
    error('The SOS problem is not feasible');
end

% One more check for SeDuMi
if strcmp(options.solver_name,'sedumi')
    if sol.info.solverInfo.feasratio < 0
        error('The SOS problem is not feasible');
    end
end

for ii = 1:N
    B{ii} = sol.eval(B{ii});
    sol.eval(B{ii})
%    prog(B{ii})
end


%%  
if true
    for iter = 1:1
        % TODO: check for convergence
        
        %% Fix B and solve for L
        % evaluate dynamics and Vtraj at every ts once (for efficiency/clarity)
        N = length(ts);
        
        for ii = 1:N
            
            clear prog
            prog = spotsosprog;
            
            prog = prog.withIndeterminate(x);            
            
            % Compute time derivatives and Lie derivatives of B
            if ii > 1
                dBdt = (B{ii} - B{ii-1})/(ts(ii) - ts(ii-1));  % take an approximate derivative - valid for small delta t's
            else
                dBdt = 0;
            end
            dBdxtimesf{ii} = diff(B{ii},x)*f{ii} + dBdt;
            dBdxtimesg{ii} = diff(B{ii},x)*g{ii};
            
%             [prog,gamma] = prog.newFree(1);
            gamma = 0;
            
            f{ii} = poly.getPolyDynamics(ts(ii));
            if (poly.getNumInputs>0)   % zero all inputs
                f{ii} = subs(f{ii},poly.getInputFrame.getPoly,zeros(poly.getNumInputs,1));
            end
            g{ii} = [0;0;1];
            
            % Generate the safe consraint set as an intersection of the complement
            % of any unsafe sets, including the largest ellipse at this point in
            % the funnel
%             idx00 = find(ts(ii) == ts00, 1);
            idx00 = find(abs(ts(ii) - ts00) == min(abs(ts(ii) - ts00)), 1);
            tmpArray = ac.ellipsoid;
            [ellq,ellQ] = double(tmpArray(idx00));
            g_Xu1 = [];
            %g_Xu1 = 1 - (x - ellq)'*inv(ellQ)*(x - ellq);
            if ii == N
                g_Xu1 = 1 - (x - ellq)'*inv(ellQ)*(x - ellq);
            end
            for iplane = 1:length(polyK)%iplane = find(isectArray(idx00,:) == 0)
                g_Xu1 = [g_Xu1; polyH(iplane,:)*x(1:2) - polyK(iplane)]; % unsafe set
            end
            %Vt = msubs(V0,x,ts(ii));
            %g_Xu1 = g_Xu1.inStateFrame(Vt.getFrame);
            
            % SOS constraints
            % C1: B(x) < 0 inside initial condition set
            if ii == 1
                [prog,L1] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
                prog = prog.withSOS(-B{ii} - 0.1*L1*g_X0 - 1e-5); % The 1e-4 is there to make sure the inequality is strict (but you don't have to worry about this for the other constraints)
            else
                tmpArray = ac.ellipsoid;
                [ellq,ellQ] = double(tmpArray(idx00));
                g_X01 = 0.01 - (x - ellq)'*inv(ellQ)*(x - ellq);
                [prog,L1] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
                prog = prog.withSOS(-B{ii} - 1*L1*g_X01 - 1e-5); % The 1e-4 is there to make sure the inequality is strict (but you don't have to worry about this for the other constraints)
            end
            
            % C2: B(x) >= 0 inside unsafe set
            if ii == N
                weightMultiplier = 1;
            else
                weightMultiplier = 0.01;
            end
            [prog,L21] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
            prog = prog.withSOS(B{ii} - weightMultiplier*L21*g_Xu1(1) );
            for iu = 2:length(g_Xu1)
                [prog,L21] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
                prog = prog.withSOS(B{ii} - 0.1*L21*g_Xu1(iu) );
            end
            
            % C3: Bdot <= 0 when dB/dx*g = 0
            [prog,Ld0] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
            [prog,Ld1] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
            [prog,Ld2] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
            prog = prog.withSOS(gamma - dBdxtimesf{ii} - 0.01*Ld1*(dBdxtimesg{ii} + 1e-4) - 0.01*Ld0*g_Xdomain);
            prog = prog.withSOS(gamma - dBdxtimesf{ii} + 0.01*Ld2*(dBdxtimesg{ii} + 1e-4) - 0.01*Ld0*g_Xdomain);
            
%             prog = prog.withSOS(-dBdxtimesf{ii} - 0.1*Ld0*g_Xdomain);
            
            % Solve SOS program (for L with B fixed)
            SOSoptions = spot_sdp_default_options();
            SOSoptions.verbose = 1;
%             sol = prog.minimize(gamma,options.solver,SOSoptions);
            sol = prog.minimize(0,options.solver,SOSoptions);
            
            % Check if SOS program ran correctly
            if ~sol.isPrimalFeasible
                error('The SOS problem is not feasible');
            end
            
            % One more check for SeDuMi
            if strcmp(options.solver_name,'sedumi')
                if sol.info.solverInfo.feasratio < 0
                    error('The SOS problem is not feasible');
                end
            end
            
%             slack{ii} = double(sol.eval(gamma));
            LD1{ii} = sol.eval(Ld1);
            LD2{ii} = sol.eval(Ld2);
            
        end
        
        
        %% Fix L and solve for B
        clear prog
        prog = spotsosprog;
        
        prog = prog.withIndeterminate(x);
        
        % evaluate dynamics and Vtraj at every ts once (for efficiency/clarity)
        for ii = 1:N
            % Initialize barrier function as a free d-th degree polynomial
            [prog,Bnew] = prog.newFreePoly(monomials(x,0:d));
            B{ii} = Bnew;
            
            f{ii} = poly.getPolyDynamics(ts(ii));
            if (poly.getNumInputs>0)   % zero all inputs
                f{ii} = subs(f{ii},poly.getInputFrame.getPoly,zeros(poly.getNumInputs,1));
            end
            g{ii} = [0;0;1];
            
            % Generate the safe consraint set as an intersection of the complement
            % of any unsafe sets, including the largest ellipse at this point in
            % the funnel
%             idx00 = find(ts(ii) == ts00, 1);
            idx00 = find(abs(ts(ii) - ts00) == min(abs(ts(ii) - ts00)), 1);
            tmpArray = ac.ellipsoid;
            [ellq,ellQ] = double(tmpArray(idx00));
            g_Xu1 = [];
            %g_Xu1 = 1 - (x - ellq)'*inv(ellQ)*(x - ellq);
            if ii == N
                g_Xu1 = 1 - (x - ellq)'*inv(ellQ)*(x - ellq);
            end
            for iplane = 1:length(polyK)%iplane = find(isectArray(idx00,:) == 0)
                g_Xu1 = [g_Xu1; polyH(iplane,:)*x(1:2) - polyK(iplane)]; % unsafe set
            end
            %g_Xu1 = g_Xu1.inStateFrame(V0(ii).getFrame);
            
            
            % Compute time derivatives and Lie derivatives of B
            if ii > 1
                dBdt = (B{ii} - B{ii-1})/(ts(ii) - ts(ii-1));  % take an approximate derivative - valid for small delta t's
            else
                dBdt = 0;
            end
            dBdxtimesf{ii} = diff(B{ii},x)*f{ii} + dBdt;
            dBdxtimesg{ii} = diff(B{ii},x)*g{ii};
            
            % SOS constraints
            % C1: B(x) < 0 inside initial condition set
            if ii == 1
                [prog,L1] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
                prog = prog.withSOS(-B{ii} - 0.1*L1*g_X0 - 1e-5); % The 1e-4 is there to make sure the inequality is strict (but you don't have to worry about this for the other constraints)
            else
                tmpArray = ac.ellipsoid;
                [ellq,ellQ] = double(tmpArray(idx00));
                g_X01 = 0.01 - (x - ellq)'*inv(ellQ)*(x - ellq);
                [prog,L1] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
                prog = prog.withSOS(-B{ii} - 1*L1*g_X01 - 1e-5); % The 1e-4 is there to make sure the inequality is strict (but you don't have to worry about this for the other constraints)
            end
           
            % C2: B(x) >= 0 inside unsafe set
            if ii == N
                weightMultiplier = 1;
            else
                weightMultiplier = 0.01;
            end
            [prog,L21] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
            prog = prog.withSOS(B{ii} - weightMultiplier*L21*g_Xu1(1) );
            for iu = 2:length(g_Xu1)
                [prog,L21] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
                prog = prog.withSOS(B{ii} - 0.01*L21*g_Xu1(iu) );
            end
            
            % C3: Bdot <= 0 when dB/dx*g = 0
            [prog,Ld0] = prog.newSOSPoly(monomials(x,0:d-2)); % Quadratic SOS multiplier polynomial
            prog = prog.withSOS(-dBdxtimesf{ii} - 0.01*LD1{ii}*(dBdxtimesg{ii} + 1e-4) - 0.01*Ld0*g_Xdomain);
            prog = prog.withSOS(-dBdxtimesf{ii} + 0.01*LD2{ii}*(dBdxtimesg{ii} + 1e-4) - 0.01*Ld0*g_Xdomain);
            
        end
        
        % Solve SOS program (for B with L fixed)
        SOSoptions = spot_sdp_default_options();
        SOSoptions.verbose = 1;
        sol = prog.minimize(0,options.solver,SOSoptions);
        
        % Check if SOS program ran correctly
        if ~sol.isPrimalFeasible
            error('The SOS problem is not feasible');
        end
        
        % One more check for SeDuMi
        if strcmp(options.solver_name,'sedumi')
            if sol.info.solverInfo.feasratio < 0
                error('The SOS problem is not feasible');
            end
        end
        
        for ii = 1:N
            B{ii} = sol.eval(B{ii});
            sol.eval(B{ii})
            %    prog(B{ii})
        end
        
    end
end
        
%%
% Print out B after zeroing out very small coefficient
xprint = msspoly('x',3); % print in terms of x1 and x2
disp(' ');
disp('Barrier function:');

for ii = 1:N
    color = colorArray(floor(ii*size(colorArray,1)/N),:);
    B_sol{ii} = clean(subs(sol.eval(B{ii}),x,xprint),1e-5)
    
    B{ii} = sol.eval(B{ii});
    sol.eval(B{ii})
%     prog(B{ii})

    % check for sanity- at equilibrium point.  TODO: check more points in
    % the IC set too
    %             idx00 = find(ts(ii) == ts00, 1);
    idx00 = find(abs(ts(ii) - ts00) == min(abs(ts(ii) - ts00)), 1);
    [ellq,ellQ] = double(tmpArray(idx00));
    subs(B{ii},x,ellq)

end


