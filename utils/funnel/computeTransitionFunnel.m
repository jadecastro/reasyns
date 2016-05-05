function [ac_trans, errTrans] = computeTransitionFunnel(sysArray,reg,regDefl,regBnd,aut,acTrans,iModeToPatch,iModeSuccessor,options)
%
% Reach operation -- construct transition funnels
%

global ME

debug = true;

maxFunnelsTrans = options.maxFunnelsTrans;
maxFunTrials = options.maxFunTrials;
maxTrajLength = options.maxTrajLength;
maxTrials2 = options.maxTrials2;

IsectInAll{iModeToPatch} = 0;%zeros(1,size(qCover{iModeToPatch},1));

% Avoid regions for starting region
% regSafeS = getReg(reg,regBnd,aut,iModeToPatch);

ac_trans = [];
errTrans = [];
acLast = [];
xTest = [];

for itrans = 1:length(aut.trans)
    if iModeToPatch == aut.trans{itrans}(1) && iModeSuccessor == aut.trans{itrans}(2), break; end
end

% ==========================
% Compute transition funnels

funFail = true;

for funindx = 1:maxFunTrials
    funindx
    
    i = 1;
    %         while true  % will only increment i if atomic controller generation has succeeded.

    % Get the system model and region pair for this transition
    aut.f{itrans}
    sys = sysArray(aut.f{itrans});
    regTrans.init = reg(aut.label{vertcat(aut.state{:}) == iModeToPatch});
    regTrans.goal = reg(aut.label{vertcat(aut.state{:}) == iModeSuccessor});
    regTrans.union = union(regTrans.init, regTrans.goal);
    
    for j = 1:length(acTrans)
        if ~isempty(acTrans{j})
            if any(acTrans{j}.post == iModeToPatch)
                acLast = [acLast; acTrans{j}];
            end
        end
    end
    
    if length(acLast) == 1  % if precisely one incoming funnel (constructed so far), select the initial point as the final point in that funnel.  TODO: don't join together funnels if another incoming funnel is constructed later on
        ttmp = getTimeVec(acLast.x0);
        qCenter = double(acLast.x0,ttmp(end));
        Qsav = sys.sysparams.Qrand;
        sys.sysparams.Qrand = 1e-4*eye(sys.sysparams.n);
        for trial2 = 1:maxTrials2
            try
                initState = getCenterRand(sys,regDefl(aut.label{vertcat(aut.state{:}) == iModeToPatch}),[],qCenter);
                break
            catch ME
                disp(ME.message)
                disp('something went wrong with the initial state generation... recomputing')
            end
        end
        if trial2 == maxTrials2
            funFail = true;
        end
        sys.sysparams.Qrand = Qsav;
    else
        for trial2 = 1:maxTrials2
            try
                initState = getCenterRand(sys,regDefl(aut.label{vertcat(aut.state{:}) == iModeToPatch}),[]);
                break
            catch ME
                disp(ME.message)
                disp('something went wrong with the initial state generation... recomputing')
            end
        end
    end
    initState
    initOutput = sys.state2SEconfig([],initState,[]);
    initOutput = initOutput(1:2);
    
    regSafeSG = reg.getRegTrans(regBnd,aut,itrans);
    %         regSafeG = getReg(reg,regBnd,aut,iModeSuccessor);
    %                 ellTransS = ellipsoid(ac_trans(itrans,1));
    
    trial2 = 1;
    for trial2 = 1:maxTrials2
        
        disp('Computing final point....')
        try
            finalState = getCenterRand(sys,regDefl(aut.label{vertcat(aut.state{:}) == iModeSuccessor}),[]); %,vReg{aut.label{iModeToPatch}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
            goalOutput = sys.state2SEconfig([],finalState,[]);
            goalOutput = goalOutput(1:2);
            
            finalState
            path = [initState; finalState];
            type = 'state';
            disp('Computing nominal trajectory....')
            
            % TODO: make it an option to fall back on RRT if the
            % optional method 'dynamicsWaypointSteering' is not
            % specified!
            [u0,x0] = computeTrajectory(sys,initState,path,options,type);
            
            figure(500), hold on, plot(x0,[],500)
            
            % Check the trajectory for invariant violation, final point
            % violation of a successor funnel (if only one), and
            % violation of the successor region.
            disp('Possible trajectory found. Checking for consistency...')
            noInvarViolation = isinside([regTrans.init, regTrans.goal],sys,downsampleUniformly(x0,options.sampSkipColl)');
            ttmp = x0.getTimeVec;
            
            noFinalPointViolation = isinside(regTrans.goal,sys,double(x0,ttmp(end))');
            
            noFinalEllipsoidFunnelViolation = true;
            ballTest = ellipsoid(double(x0,ttmp(end)),inv(sys.sysparams.Qf));
            %                 if ~isempty(acNext) % any transition funnels have been already computed for any of the successors and only one outgoing transition from the successor
            %                     noFinalEllipsoidFunnelViolation = acNext.funnelContainsEllipsoid(sys,ballTest,100);
            %                 end
            
            [H,K] = double(regTrans.goal.p);
            hpp = hyperplane(H',K');
            
            [~,~,H] = sys.getRegNonRegStates([],double(x0,ttmp(end)),[]);
            ballTestProj = projection(ballTest,H(1:2,:)');
            noFinalEllipsoidRegionViolation = ~any(intersect(ballTestProj,hpp,'u'));
            
            if noInvarViolation && noFinalPointViolation && noFinalEllipsoidFunnelViolation && noFinalEllipsoidRegionViolation ...
                    && length(x0) > 1 && length(x0) < maxTrajLength,
                break,
            end
            disp('Trajectory incompatible with constraints; recomputing...')
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
            
            [ac, c] = computeAtomicController(u0,x0,sys,regTrans,ellToCompose,options,rhof);
            ac = ac.setTransition(iModeToPatch,iModeSuccessor);

            ac_trans = [ac_trans; ac];
            
            funFail = false;
            
            % TODO: check if inside the polytopes
            
            plot(ac.x0,'k',5)
            
        catch ME
            rethrow(ME)
            disp(ME.message)
            disp('something went wrong with the funnel computation...  kicking back out to the main script.')
            errTrans = trans(:,2)==iModeToPatch;
        end
        
        % verify that the funnel we obtained meets the containment criteria
        % TODO: create the isinside method
        if ~debug
            if ~isinside(ac,regSafeSG,options.sampSkipValid) || max(rho_d) < 1e-6
                funFail = true;
            end
        end
    end
    
    if funFail     % Need to repeat this iteration
        disp('funnel computation failed.')
        % Optional: plot the failed result
        %             plot(ac.x0,'k',5)
    else
        plot(ac.x0,'k',3)
        
        disp(['Iteration #',num2str(i),' / ',num2str(maxFunnelsTrans(iModeToPatch))])

    end
        
    if ~isempty(errTrans)
        break
    end

    if ~funFail || ~isempty(errTrans)
        break
    end
end

if funFail && funindx == maxFunTrials
    error('Cannot find feasible funnels for this mode.')
end


