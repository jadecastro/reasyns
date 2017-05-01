function [ac_result, err] = computeFunnel(sysArray, x0, u0, reg, aut, acTrans, iMode, iModeSuccessor, options)
% Construct funnels for a transition of the state machine.
%

global ME

debug = false;
doPlot = false;

transFlag = true;

maxFunnelsTrans = length(x0);

IsectInAll{iMode} = 0;%zeros(1,size(qCover{iMode},1));

% Avoid regions for starting region
% regSafeS = getReg(reg,regBnd,aut,iMode);

ac_result = [];
err = [];
acLast = [];
xTest = [];

for itrans = 1:length(aut.trans)
%     if iMode == aut.trans{itrans}(1) && iModeSuccessor == aut.trans{itrans}(2), break; end
    if iMode == aut.trans{itrans}(1) && ~isempty(intersect(iModeSuccessor, aut.trans{itrans}(2))), break; end
end

% Initially set to fail - if the trajectory checks out then this will be
% set to false.
funFail = true;

% Get the system model and region pair for this transition
sys = sysArray(aut.f{itrans});
regTrans.init = reg(aut.label{vertcat(aut.state{:}) == repmat(iMode,length(aut.state), 1)});
regTrans.goal = reg(aut.label{vertcat(aut.state{:}) == repmat(iModeSuccessor, length(aut.state), 1)});
regTrans.union = union(regTrans.init, regTrans.goal);

for j = 1:length(acTrans)
    if isempty(acTrans{j}), continue, end
    if ~any(vertcat(aut.label{acTrans{j}.post}) == repmat(aut.label{iMode}, length(acTrans{j}.post), 1))
        continue
    end
    acLast = [acLast; acTrans{j}];
end

%     initState
initOutput = sys.state2SEconfig([], x0.double(0), []);

% Validity-check the supplied trajectory with respect to the regions and supplied parameters.
disp('Validity-checking the trajectory....')
try
    figure(3), hold on, plot(x0, [], 3)
    
    % Check the trajectory for violation of the region invariants.
    disp('Checking the provided trajectory for consistency...')
    ballRandomSamples = downsampleUniformly(x0, options.sampSkipColl)';
    noInvarViolation = isinside([regTrans.init, regTrans.goal], sys, ballRandomSamples);
    ttmp = x0.getTimeVec;
    
    % Check the final point of the trajectory.
    noFinalPointViolation = isinside(regTrans.goal, sys, double(x0,ttmp(end))');
    
    % We skip the check for composition with funnels, for now.
    noFinalEllipsoidFunnelViolation = true;

    % Check the final ellipsoid to determine if it fits within the goal region.
    ballTest = ellipsoid(double(x0, ttmp(end)), options.rhof^2*inv(sys.sysparams.Qf));
    
    ballTestProj = reg.projection(sys, ballTest);
    if ~isempty(ballTestProj)
        figure(3), plot(ballTestProj, 'g', 3)
    end
    
    noFinalEllipsoidRegionViolation = regTrans.goal.regionContainsEllipsoid(sys, ballTest);
    
    if noInvarViolation && noFinalPointViolation && noFinalEllipsoidFunnelViolation && ...
            noFinalEllipsoidRegionViolation && length(x0) > 1,
        funFail = false;
    end
    if ~noInvarViolation,                disp('computeFunnel: Trajectory not contained in the region.'); end
    if ~noFinalPointViolation,           disp('computeFunnel: Final point not contained within the goal region.'); end
    if ~noFinalEllipsoidFunnelViolation, disp('computeFunnel: Final ellipse not contained within the successor funnel.'); end
    if ~noFinalEllipsoidRegionViolation, disp('computeFunnel: Final ellipse not contained within the goal region.'); end
catch ME
    %  rethrow(ME)
    disp('computeFunnel: Trajectory validity check failed for the following reason:');
    disp(ME.message);
    disp('    stack trace:');
    for i = 1:length(ME.stack)
        disp(['line: ', num2str(ME.stack(i).line), '  function: ', ME.stack(i).name]);
    end
end

% Compute the funnel.
if ~funFail
    disp('Computing funnel....')
    
    try
        if ~isempty(acLast)
            last_ellipsoid = acLast.ellipsoid;
            ellToCompose = last_ellipsoid(end);
        else
            ellToCompose = [];
        end
        rhof = options.rhof;   %final rho.
        options.isMaximization = true;
        
        [ac, ~] = computeAtomicController(u0, x0, sys, regTrans, ellToCompose, options, transFlag, rhof);
        ac = ac.setTransition(iMode,iModeSuccessor);
        
        funFail = false;
        
        % TODO: check if inside the polytopes        
    catch ME
        %rethrow(ME)
        disp('computeFunnel: Funnel computation failed for the following reason:')
        disp(ME.message);
        disp('    stack trace:');
        for i = 1:length(ME.stack)
            disp(['line: ', num2str(ME.stack(i).line), '  function: ', ME.stack(i).name]);
        end
        for itrans = 1:length(aut.trans)
            if aut.trans{itrans}(2)==iMode
                err = itrans;
            end
        end
    end
    
    % verify that the funnel we obtained meets the containment criteria
    % TODO: create the isinside method
    if exist('ac')
        if debug
            funFail = false;
        elseif ~isinside(ac, [regTrans.init; regTrans.goal], sys, options.sampSkipValid)
            disp('computeFunnel: Computed funnel does not lie in the union of the designated regions.');
            %keyboard
            funFail = true;
        end
    else
        funFail = true;
    end
end

if funFail
    error('computeFunnel: Error: Funnel computation failed.')
    % Optional: plot the failed result
    % plot(ac.x0, 'k', 5)
else
    ac_result = [ac_result; ac];
    
    if doPlot
        plot(ac,sys,5)
    end
    % plot(ac.x0,'k',3)
    % disp(['Iteration #',num2str(i),' / ',num2str(maxFunnelsTrans(iMode))])
    
end
