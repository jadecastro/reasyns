function [ac_result, err] = computeFunnel(sysArray, x0, u0, reg, regDefl, regBnd, aut, iMode, iModeSuccessor, options)
% Construct funnels for a transition of the state machine.
%

global ME

debug = false;
doPlot = false;

transFlag = true;

maxFunnelsTrans = options.maxFunnelsTrans;
maxFunTrials = options.maxFunTrials;
maxTrajLength = options.maxTrajLength;
maxTrials2 = options.maxTrials2;

IsectInAll{iMode} = 0;%zeros(1,size(qCover{iMode},1));

% Avoid regions for starting region
% regSafeS = getReg(reg,regBnd,aut,iMode);

ac_result = [];
err = [];
acLast = [];
xTest = [];

for itrans = 1:length(aut.trans)
%     if iMode == aut.trans{itrans}(1) && iModeSuccessor == aut.trans{itrans}(2), break; end
    if iMode == aut.trans{itrans}(1) && ~isempty(intersect(iModeSuccessor,aut.trans{itrans}(2))), break; end
end

funFail = true;

for funindx = 1:maxFunTrials
    funindx
    
    i = 1;
    %         while true  % will only increment i if atomic controller generation has succeeded.

    % Get the system model and region pair for this transition
    aut.f{itrans}
    sys = sysArray(aut.f{itrans});
    regTrans.init = reg(aut.label{vertcat(aut.state{:}) == repmat(iMode,length(aut.state),1)});
    regTrans.goal = reg(aut.label{vertcat(aut.state{:}) == repmat(iModeSuccessor,length(aut.state),1)});
    regTrans.union = union(regTrans.init, regTrans.goal);
    
    for j = 1:length(acTrans)
        if ~isempty(acTrans{j})
            if any(vertcat(aut.label{acTrans{j}.post}) == repmat(aut.label{iMode},length(acTrans{j}.post),1))
                acLast = [acLast; acTrans{j}];
            end
        end
    end
    
%     initState
    initOutput = sys.state2SEconfig([],initState,[]);
    initOutput = initOutput(1:2)
    
    trial2 = 1;
    for trial2 = 1:maxTrials2
        
        disp('Computing final point....')
        try
            
            figure(500), hold on, plot(x0,[],500)
            
            % Check the trajectory for invariant violation, final point
            % violation of a successor funnel (if only one), and
            % violation of the successor region.
            disp('Checking the provided trajectory for consistency...')
            noInvarViolation = isinside([regTrans.init, regTrans.goal],sys,downsampleUniformly(x0,options.sampSkipColl)');
            ttmp = x0.getTimeVec;
            
            noFinalPointViolation = isinside(regTrans.goal,sys,double(x0,ttmp(end))');
            
            noFinalEllipsoidFunnelViolation = true;
            ballTest = ellipsoid(double(x0,ttmp(end)),options.rhof^2*inv(sys.sysparams.Qf));
            %                 if ~isempty(acNext) % any transition funnels have been already computed for any of the successors and only one outgoing transition from the successor
            %                     noFinalEllipsoidFunnelViolation = acNext.funnelContainsEllipsoid(sys,ballTest,100);
            %                 end
            
            ballTestProj = reg.projection(sys,ballTest);
            if ~isempty(ballTestProj)
                figure(500), plot(ballTestProj,'g',5)
            end
            
            noFinalEllipsoidRegionViolation = regTrans.goal.regionContainsEllipsoid(sys,ballTest);
            
            if noInvarViolation && noFinalPointViolation && noFinalEllipsoidFunnelViolation && noFinalEllipsoidRegionViolation ...
                    && length(x0) > 1 && length(x0) < maxTrajLength,
                break,
            end
            disp('Trajectory incompatible with constraints; recomputing...')
            if ~noInvarViolation,                disp('... trajectory not contained in the region'); end
            if ~noFinalPointViolation,           disp('... final point not contained within the goal region'); end
            if ~noFinalEllipsoidFunnelViolation, disp('... final ellipse not contained within the successor funnel'); end
            if ~noFinalEllipsoidRegionViolation, disp('... final ellipse not contained within the goal region'); end
        catch ME
            %  rethrow(ME)
            disp('something went wrong with the trajectory computation... recomputing')
        end
    end
    if trial2 == maxTrials2
        funFail = true;
        % error('Cannot find a feasible trajectory for this mode. Consider increasing the tree depth.')
    else
        funFail = false;
    end
    
    if ~funFail
        disp('Computing funnel....')
        
        try
            if ~isempty(acLast)
                tmp = acLast.ellipsoid;
                ellToCompose = tmp(end);
            else
                ellToCompose = [];
            end
            rhof = options.rhof;   %final rho.
            options.isMaximization = true;
            
            [ac, c] = computeAtomicController(u0,x0,sys,regTrans,ellToCompose,options,transFlag,rhof);
            ac = ac.setTransition(iMode,iModeSuccessor);
            
            funFail = false;
            
            % TODO: check if inside the polytopes
            
            plot(ac.x0,'k',5)
            
        catch ME
            %rethrow(ME)
            disp(ME.message)
            disp('something went wrong with the funnel computation...  kicking back out to the main script.')
            for itrans = 1:length(aut.trans)
                if aut.trans{itrans}(2)==iMode
                    err = itrans;
                end
            end
        end
        
        % verify that the funnel we obtained meets the containment criteria
        % TODO: create the isinside method
        if ~isinside(ac,[regTrans.init; regTrans.goal],sys,options.sampSkipValid)
            disp('... computed funnel does not lie in the union of the designated regions'); 
            %keyboard
            if ~debug
                funFail = true;
            end
        end
    end
    
    if funFail     % Need to repeat this iteration
        disp('funnel computation failed.')
        % Optional: plot the failed result
        %             plot(ac.x0,'k',5)
    else
        
        ac_result = [ac_result; ac];
        
        if doPlot
            plot(ac,sys,5)
        end
        %plot(ac.x0,'k',3)
        
%         disp(['Iteration #',num2str(i),' / ',num2str(maxFunnelsTrans(iMode))])

    end
        
    if ~isempty(err)
        break
    end

    if ~funFail || ~isempty(err)
        break
    end
end

if funFail && funindx == maxFunTrials
    error('Cannot find feasible funnels for this mode.')
end


