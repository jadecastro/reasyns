function [ac_inward, errTrans] = computeInwardJoiningFunnel(sysArray,reg,regDefl,regBnd,aut,acTrans,indexState,options)
% Construct inward-facing funnels joining together all transition funnel
% combinations that enter and exit a particular state.
%

global ME

debug = true;

maxFunnelsTrans = options.maxFunnelsTrans;
maxFunTrials = options.maxFunTrials;
maxTrajLength = options.maxTrajLength;
maxTrials2 = options.maxTrials2;
trans = vertcat(aut.trans{:});

IsectInAll{indexState} = 0;%zeros(1,size(qCover{indexState},1));

% Avoid regions for starting region
% regSafeS = getReg(reg,regBnd,aut,indexState);

ac_inward = [];
errTrans = [];

regTrans.init = reg(aut.label{vertcat(aut.state{:}) == indexState});
regTrans.goal = regTrans.init;
regTrans.union = union(regTrans.init, regTrans.goal);

% ==========================
% Compute inward funnels

funFail = true;

for funindx = 1:maxFunTrials
    
    i = 1;
    %         while true  % will only increment i if atomic controller generation has succeeded.
    
    %lastTrans = trans(:,2)==indexState;
    %acLast = [acTrans{lastTrans}];
    
    indexTransPostVect = [];
    candidateIndexTransPreVect = [];
    for indexTrans = 1:length(trans)
        if ~isempty(acTrans{indexTrans})            
            if aut.label{acTrans{indexTrans}.pre} == aut.label{indexState}
                indexTransPostVect = [indexTransPostVect; indexTrans];
            end
        end
    end
    for indexTrans = 1:length(trans)
        if ~isempty(acTrans{indexTrans})            
            for postState = acTrans{indexTrans}.post
                if aut.label{postState} == aut.label{indexState}
                    candidateIndexTransPreVect = [candidateIndexTransPreVect; indexTrans];
                end
            end
        end
    end
    
    % does there exist a path to any of the candidateIndexTransPostVect?
    for j = 1:length(indexTransPostVect)
        indexTransPreVect{j} = [];
        tmpPre = trans(indexTransPostVect(j),1);
        for k = 1:10
            for indexTrans = 1:length(trans)
                if ismember(trans(indexTrans,2),tmpPre) && (aut.label{trans(indexTrans,1)} == aut.label{indexState})  % a region-preserving transition has been found
                    tmpPre = [tmpPre; trans(indexTrans,1)];
                    [foundTransition, i] = intersect(trans(candidateIndexTransPreVect,2),tmpPre);
                    if ~isempty(foundTransition)
                        indexTransPreVect{j} = [indexTransPreVect{j}; candidateIndexTransPreVect(i)];
                    end
                end
            end
        end
        indexTransPreVect{j} = unique(indexTransPreVect{j});
    end
    
    if isempty(indexTransPreVect), funFail = false; return; end

    acNext = [acTrans{indexTransPostVect}];
    
    figure(3)
    clf
    hold on, axis equal
    plot(reg(aut.label{vertcat(aut.state{:}) == indexState}),'r')
    plot(acNext,sys,3)
    
    for jpost = 1:length(indexTransPostVect)
        indexTransPost = indexTransPostVect(jpost);
        
        if isempty(indexTransPreVect{jpost}), funFail = false; return; end
        
        for jpre = 1:length(indexTransPreVect{jpost})
            indexTransPre = indexTransPreVect{jpost}(jpre);
            
            acLast = acTrans{indexTransPre};
            
            sys = sysArray(aut.f{vertcat(aut.state{:}) == aut.trans{indexTransPost}(1)}); % Assign the dynamics according to the transition to be taken.
            
            regSafeSG = getRegTrans(reg,regBnd,aut,indexTransPre);
            %         regSafeG = getReg(reg,regBnd,aut,iModeSuccessor);
            
            %                 ellTransS = ellipsoid(ac_trans(indexTransPre,1));
            
            ttmp = getTimeVec(acTrans{indexTransPre}.x0);
            
            qCenter = double(acTrans{indexTransPre}.x0,ttmp(end));
            Qsav = sys.sysparams.Qrand;
            sys.sysparams.Qrand = 1e-4*eye(length(sys.sysparams.Qrand));
            for trial2 = 1:maxTrials2
                try
                    initState = getCenterRand(sys,regTrans.init,[],qCenter);
                    break
                catch ME
                    disp(ME.message)
                    disp('something went wrong with the initial state generation... recomputing')
                end
            end
            if trial2 == maxTrials2
                funFail = true;
            end
            options.Qrand = Qsav;
            
            try initState; catch, keyboard; end
            initOutput = sys.state2SEconfig([],initState,[]);
            initOutput = initOutput(1:2);
            
            trial2 = 1;
            for trial2 = 1:maxTrials2
                
                disp('Computing final point....')
                try
                    qCenter = double(acNext(jpost).x0,0);  % in the vicinity of the first point of the funnel
                    %                 finalState = getCenterRand_new(sys,regDefl(aut.label{indexState}),acNext,options,qCenter) %,vReg{aut.label{indexState}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                    finalState = double(acNext(jpost).x0,0)'; %,vReg{aut.label{indexState}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                    goalOutput = sys.state2SEconfig([],finalState,[]);
                    goalOutput = goalOutput(1:2);
                    
                    stepSize = options.TstepRRT;
                    [path] = buildReachabilityRRT(regBnd.v,{regTrans.init.v},{regTrans.init.v},[],[],[],[],initState,finalState,stepSize,sys,regTrans.init,acNext(jpost),options);
                                    %                 [path] = buildCtrlSpaceRRT(regBnd.v,{reg(aut.label{indexState}).v},{reg(aut.label{indexState}).v},acTrans{indexTransPre}.ellipsoid,[],[],[],initState,finalState,stepSize,sys,regBnd,acNext,options);
                    disp('Computing nominal trajectory....')
                    
                    x0 = Traject(path.t',path.x');
                    u0 = Traject(path.t',path.u');
                    
                    % Check the trajectory for invariant violation, final point
                    % violation of a successor funnel (if only one), and
                    % violation of the successor region.
                    disp('Possible trajectory found. Checking for consistency...')
                    noInvarViolation = isinside([regTrans.init, regTrans.goal],sys,downsampleUniformly(x0,options.sampSkipColl)');
                    ttmp = x0.getTimeVec;
                    
                    noFinalPointViolation = isinside(regTrans.goal,sys,double(x0,ttmp(end))');
                    
                    noFinalEllipsoidFunnelViolation = true;
                    ballTest = ellipsoid(double(x0,ttmp(end)),options.rhof^2*inv(sys.sysparams.Qf));
                    if ~debug
                        if ~isempty(acNext(jpost)) % any transition funnels have been already computed for any of the successors and only one outgoing transition from the successor
                            noFinalEllipsoidFunnelViolation = acNext(jpost).funnelContainsEllipsoid(sys,ballTest,100);
                        end
                    end
                    
                    [H,K] = double(regTrans.goal.p);
                    hpp = hyperplane(H',K');
                    
                    [~,~,H] = sys.getRegNonRegStates([],double(x0,ttmp(end)),[]);
                    ballTestProj = projection(ballTest,H(1:2,:)');
                    
                    figure(500), plot(ballTestProj,'g',5)
                    
                    noFinalEllipsoidRegionViolation = ~any(intersect(ballTestProj,hpp,'u'));
                    
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
                    % rethrow(ME)
                    %                     keyboard
                    disp(ME.message)
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
                    rhof = options.rhof;   %final rho.  TODO: handle the more general case and get it from containment
                    options.isMaximization = true;
                    
                    [ac, c] = computeAtomicController(u0,x0,sys,regTrans,ellToCompose,options,rhof);
                    
                    funFail = false;
                    
                    % TODO: check if inside the polytopes
                    
                    plot(ac.x0,'k',5)
                    
                catch ME
                    %                 rethrow(ME)
                    disp(ME.message)
                    disp('something went wrong with the funnel computation...  kicking back out to the main script.')
                    errTrans = trans(:,2)==indexState;
                    break
                end
                
                % verify that the funnel we obtained meets the containment criteria
                % TODO: create the isinside method
                if ~isinside(ac,regTrans.init,sys,options.sampSkipValid)
                    disp('... computed funnel does not lie in the designated region');
                    if ~debug
                        funFail = true;
                    end
                end
            end
            
            if funFail     % Need to repeat this iteration
                disp('funnel computation failed.')
                % Optional: plot the failed result
                %             plot(ac.x0,'k',5)
                break
            end

            ac_inward = [ac_inward; ac];
            
            plot(ac,sys,5)
            %plot(ac.x0,'k',3)
            
%             disp(['Iteration #',num2str(i),' / ',num2str(maxFunnelsTrans(indexState))])
        end
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


