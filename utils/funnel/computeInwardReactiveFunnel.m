function [ac_inward, bc_inward, errTrans, existingReg, newRegArray, reg] = computeInwardReactiveFunnel(sysArray,reg,regDefl,regBnd,aut,acTrans,acIn,indexState,options,fileName)
% Construct inward-facing reactive funnels for all outgoing transitions
% from a given state.
%

global ME

debug = false;

maxFunnelsTrans = options.maxFunnelsTrans;
maxFunTrials = options.maxFunTrials;
maxTrajLength = options.maxTrajLength;
maxTrials1 = options.maxTrials1;
maxTrials2 = options.maxTrials2;
trans = vertcat(aut.trans{:});

IsectInAll{indexState} = 0;%zeros(1,size(qCover{indexState},1));

% Avoid regions for starting region
% regSafeS = getReg(reg,regBnd,aut,indexState);

ac_inward = [];
bc_inward = [];
errTrans = [];
newRegArray = [];
ac = [];
bc = [];
acNew = [];
bcNew = [];
xTest = [];

% ==========================
% Compute inward funnels

containsLevelSet = false;

acLast = [];
acNextAll = [];
indexTransPostVect = [];
candidateIndexTransPreVect = [];
for indexTrans = 1:length(trans)
    if ~isempty(acTrans{indexTrans})
        if aut.label{acTrans{indexTrans}.pre} == aut.label{indexState}
            acNextAll = [acNextAll; acTrans{indexTrans}];
            indexTransPostVect = [indexTransPostVect; indexTrans];
        end
    end
end
for indexTrans = 1:length(trans)
    if ~isempty(acTrans{indexTrans})
        for postState = acTrans{indexTrans}.post
            if aut.label{postState} == aut.label{indexState}
                acLast = [acLast; acTrans{indexTrans}];
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

count = 0;

% Loop over all successor modes
for indexTrans = indexTransPostVect'
    count = count+1;
    
    sys = sysArray(aut.f{indexTrans});  % For now, delay changing the dynamics until the transition reach tube is reached.
    
    acNext = [acTrans{indexTrans}];
    regMode.init = reg(aut.label{vertcat(aut.state{:}) == indexState});
    regGoal = reg(aut.label{vertcat(aut.state{:}) == acNext.post});
    
    indexTransPre = indexTransPreVect{count}(1);  % choose the first one if more than one option to choose from
    acLast = acTrans{indexTransPre};
    
    % ==========================
    % Verify invariance of the transition funnels up to the edge of the region
      
    % 1. Verify the region boundary itself
    % TODO: extract all the indices for which there is an overlap with the region boundary; compute a hull containing the intersection with the boundary
    existingReg = [];
    newReg = [];
    
    funStepSize = 130;
    ttmp = acNext.x0.getTimeVec();
    idxToChkReactivity{indexTrans} = 1;
    for j = length(ttmp):-funStepSize:20
        if isinside(reg(aut.label{vertcat(aut.state{:}) == indexState}),sys,double(acNext.x0,ttmp(j))')
            idxToChkReactivity{indexTrans} = j:-funStepSize:20;  % you can adjust the step amount smaller to improve the "tightness" of the resulting new region to the true unverified set.
            break
        end
    end
    
    indexToCompose = idxToChkReactivity{indexTrans}(1);
    
    % DEBUGGING 
    %     if indexTrans == 4
    %         indexToCompose = 20;  % first- good
    %     else
    %         indexToCompose = 20;  % second- good
    %     end
    %indexToCompose = 700;  % first- bad
    %indexToCompose = 400;  % second- bad
    %     if indexTrans == 4
    %         indexToCompose = 650;  % first- good
    %     else
    %         indexToCompose = 80;  % second- good
    %     end
    %indexToCompose = 750;  % first (lab)- more friendly
    %indexToCompose = 340;  % second (lab)- more friendly
    
    figure(3)
    clf
    hold on, axis equal
    plot(reg(aut.label{vertcat(aut.state{:}) == indexState}),'r')
    plot(acNext,sys,3)
%     plot(acLast,sys,3)
    
    % ==========================
    % Now, search for funnels that verify reachability from the initial set which are invariant to the region 
    % We do this until containment is reached, failing upon reaching maxTrials1 iterations.
    trial1 = 0;
    funFail = true;
    while funFail && (trial1 <= maxTrials1)
        trial1 = trial1 + 1;
        
        disp(['  Attempting to compose outgoing transition funnel Index ',num2str(indexToCompose)])
        if indexToCompose <= 0
            funFail = true;
            disp('Cannot generate reactive funnels because I have exhausted all the level sets in the transition funnel without success.')
            break
        end
        ttmp = acNext.x0.getTimeVec();
        qCenter = double(acNext.x0,ttmp(indexToCompose));
        initState = [];
        for trial2 = 1:maxTrials2
            try
                if trial1 == 1  % use a greedy approach: choose a random point "close" to the level set trajectory
                    Qsav = sys.sysparams.Qrand;
                    sys.sysparams.Qrand = 1e-3*eye(sys.sysparams.n);
                    
                    initState = getCenterRand(sys,reg(aut.label{vertcat(aut.state{:}) == indexState}),[],qCenter);
                    
                    sys.sysparams.Qrand = Qsav;
                else % relax the restriction on the random point computation
                    initState = getCenterRand(sys,reg(aut.label{vertcat(aut.state{:}) == indexState}),[],qCenter);
                end
                break
            catch ME
                rethrow(ME)
%                 disp(ME.message)
%                 disp('something went wrong with the initial state generation... recomputing')
            end
        end
        if isempty(initState)
            funFail = true;
        end
        
%         initState
        %initOutput = sys.state2SEconfig([],initState,[]);
        %initOutput = initOutput(1:2);
        
        funFail = false;
        for trial2 = 1:maxTrials2
            
            disp('Computing final point....')
            try
                % Choose the appropriate region for acceptance of the goal states
                if isempty(acIn)  
                    if isempty(acLast)
                        acNextFirst = acNextAll(1);
                        ttmpFirst = acNextFirst.x0.getTimeVec();
                        qCenter = double(acNextFirst.x0,ttmpFirst(1));  % in the vicinity of the first point of the funnel
                        finalState = getCenterRand(sys,regDefl(aut.label{vertcat(aut.state{:}) == indexState}),acNextAll,qCenter); %,vReg{aut.label{indexState}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                        acAcceptCriterion = acNextAll;
                    else
                        ttmpLast = acLast.x0.getTimeVec();
                        qCenter = double(acLast.x0,ttmpLast(end));  % in the vicinity of the last point of the funnel
                        finalState = getCenterRand(sys,regDefl(aut.label{vertcat(aut.state{:}) == indexState}),acLast,qCenter); %,vReg{aut.label{indexState}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                        acAcceptCriterion = acLast;
                    end
                else  % make use of the provided inward funnels as our goal funnel
                    acInward = acIn{indexState}(1);  % get the first one in the set.
                    ttmpLast = acInward.x0.getTimeVec();
                    qCenter = double(acInward.x0,ttmpLast(1));  % in the vicinity of the first point of the funnel
                    finalState = getCenterRand(sys,regDefl(aut.label{vertcat(aut.state{:}) == indexState}),acInward,qCenter); %,vReg{aut.label{indexState}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                    acAcceptCriterion = acInward;
                end
%                 finalState
                goalOutput = sys.state2SEconfig([],finalState,[]);
                goalOutput = goalOutput(1:2);

                disp('Computing nominal trajectory....')
                [path] = buildReachabilityRRT(regBnd.v,{regMode.init.v},{regMode.init.v},[],[],[],[],initState,finalState,options.TstepRRT,sys,regMode.init,acAcceptCriterion,options);
                
                x0 = Traject(path.t',path.x');
                u0 = Traject(path.t',path.u');
                
                % Check the trajectory for invariant violation, and
                % violation of the current funnel with the
                % acAcceptCriterion funnel.
                disp('Possible trajectory found. Checking for consistency...')
                noInvarViolation = isinside(regMode.init,sys,downsampleUniformly(x0,options.sampSkipColl)');
                ttmp = x0.getTimeVec;
                
                noFinalEllipsoidFunnelViolation = true;
                if ~debug
                    ballTest = ellipsoid(double(x0,ttmp(end)),options.rhof^2*inv(sys.sysparams.Qf));
                    if ~isempty(acAcceptCriterion) % any transition funnels have been already computed for any of the successors and only one outgoing transition from the successor
                        noFinalEllipsoidFunnelViolation = acAcceptCriterion.funnelContainsEllipsoid(sys,ballTest,100);
                    end
                end
                
                if noInvarViolation && noFinalEllipsoidFunnelViolation && length(x0) > 1 && length(x0) < maxTrajLength,
                    funFail = false;
                    break,
                end
                disp('Trajectory incompatible with constraints; recomputing...')
                if ~noInvarViolation,                       disp('... trajectory not contained in the region'); end
                if ~noFinalEllipsoidFunnelViolation,        disp('... final ellipse not contained within the successor funnel'); end
            catch ME
                %  rethrow(ME)
                disp(ME)
                disp('something went wrong with the trajectory computation... recomputing')
                funFail = true;
            end
        end
        if funFail && (trial2 == maxTrials2)
            indexToCompose = indexToCompose - funStepSize;
            % error('Cannot find a feasible trajectory for this mode. Consider increasing the tree depth.')
        else
            funFail = false;
        end
        
        if ~funFail
            disp('Computing funnel....')
            % save 'barrier_test'
            
            try
                
                % [ac,existingRegNew,newRegNew] = computePolytopeAtomicControllerDubins(u0,x0,sys,acNext,reg(aut.label{indexState}),options.ctrloptions_trans,options.sampSkipFun,xMssExt);
                x00 = x0;  u00 = u0;
                
                tmp = acNext.ellipsoid;
                ellToCompose = tmp(indexToCompose);
                
                rhof = options.rhof;   %final rho.  TODO: handle the more general case and get it from containment
                options.isMaximization = true;
                [acNew, cNew, isectIdx] = computeAtomicController(u00,x00,sys,regMode,ellToCompose,options,rhof);
                acNew.setTransition(indexState);
                
                plot(acNew.x0,'k',5)
                rhoi = double(acNew.rho,0);
                
                % verify that the funnel we obtained meets the containment criteria
                % TODO: move this into computeAtomicController
                if ~isinside(acNew,regMode.init,sys,options.sampSkipValid)
                    disp('... computed funnel does not lie in the designated region');
                    if ~debug
                        funFail = true;
                    end
                else
                    funFail = false;
                end
                
                if ~isempty(isectIdx)
                    
                    if options.flagBuildConformingFunnel
                        % NB: the following assumes only one contiguous interval where the funnel left the region.
                        % split the funnel into three parts:
                        %   - one from the start to the time of the first intersection
                        %   - another during the interval of intersection (this is created using CBFs)
                        %   - another following the intersection
                        t = acNew.x0.getTimeVec();
                        acPre = [];
                        
                        % TODO: put all this into quadraticAC
                        t1 = t(max(isectIdx)+1:end);
                        x1 = Traject(t1,double(acNew.x0,t1));
                        u1 = Traject(t1,double(acNew.u0,t1));
                        K1 = Traject(t1,double(acNew.K,t1));
                        P1 = Traject(t1,double(acNew.P,t1));
                        rho1 = Traject(t1,double(acNew.rho,t1));
                        Vquad = acNew.V(max(isectIdx)+1:end);
                        acPost = QuadraticAC(x1,u1,K1,P1,rho1,Vquad,sys,acNew.pre);
                        
                        if false %min(isectIdx) > 1
                            t0 = t(1:min(isectIdx)-1);
                            x0 = Traject(t0,double(ac.x0,t0));
                            u0 = Traject(t0,double(ac.u0,t0));
                            
                            % Compute a new atomic controller that minimizes the funnel (in contrast to the original method that maximizes it).
                            % The minimization is used here to find a suitable initial condition for the CBF.
                            options.isMaximization = false;
                            [acPre] = computeAtomicController(u0,x0,sys,regMode,ellToCompose,options,rhoi);
                            acPre.setTransition(indexState);
                            
                        end
                        
                        % Now, construct the barriers
                        t01 = t(min(isectIdx):max(isectIdx));
                        x01 = Traject(t01,double(acNew.x0,t01));
                        u01 = Traject(t01,double(acNew.u0,t01));
                        
                        % define the initial set as the end of the prefix funnel
                        if ~isempty(acPre)
                            tmp = acPre.ellipsoid;
                            ellToCompose = tmp(end);
                        else
                            tmp = acNext.ellipsoid;
                            ellToCompose = tmp(indexToCompose);
                        end
                        
                        [bcNew] = computeConformingFunnel(u00,x00,u01,x01,sys,acNew,regMode.init,ellToCompose,options);
                        bcNew.setTransition(indexState);
                        
                        % Plot stuff
%                         figure(90)
%                         axis equal
%                         hold on
%                         plot(reg(1),'r')
%                         plot(reg(2),'g')
                        
                        bcNew.plot(acNew.x0,5)
                        
                        % pause
                        % plot(acNext,sys,90,[],[0,0,1])
                        % plot(acPost,sys,90,[],[0,1,0])
                        plot(acNext.x0,'k',5)
                        plot(acNew.x0,'k',5)
                        
                        funFail = false;
                        
                    else
                        funfail = true;
                        disp('computation of inward funnels failed.')
                    end
                end
                
                % TODO: check if inside the polytopes
                
            catch ME
                %                     rethrow(ME)
                disp(ME.message)
                disp('something went wrong with the funnel computation... attempting to find a more conservative initial set, and restarting reachability analysis.')
%                 keyboard
                indexToCompose = indexToCompose - funStepSize;
                errTrans = trans(:,2)==indexState;
            end
            
        end
        
        if funFail && (trial1 == maxTrials1)
            indexToCompose = indexToCompose - funStepSize;
            % error('Cannot find a feasible funnel for this mode. Consider increasing the tree depth.')
        end
        
        if ~funFail && isempty(errTrans)
            
            if ~isempty(xTest), 
                [containsLevelSet, xTest] = acNew.funnelContainsEllipsoid(sys,xTest);
            else
                [containsLevelSet, xTest] = acNew.funnelContainsEllipsoid(sys,ellToCompose);
            end
            
            % If the computation hasn't failed up to this point, keep the results
            ac = [ac; acNew];
            bc = [bc; bcNew];
            
            plot(acNew,sys,5)
            
            if ~containsLevelSet
                funFail = true;
            end
        end
    end % while
    
    if isempty(errTrans) && ~funFail  % if a funnel was generated that passes the test, 
        
        % Create a new region based on the last unverified index of the next funnel.
        idxLast = indexToCompose + funStepSize;
        
    else
        disp('No funnel generated for this mode. Updating the input files accordingly.')
        
        % Set the index to the beginning of the next funnel and don't store an atomic controller.
        idxLast = 1;

    end
    
    % Create the new region
    newRegVert = [];
    ttmp = acNext.x0.getTimeVec();
    ellAcNext = projection(acNext,sys);
    for j = length(ttmp):-funStepSize:idxLast
        newRegVert = [buildNewRegion(ellAcNext(j), true); newRegVert];
    end
    [newRegConvHullIdx] = convhull(newRegVert(:,1),newRegVert(:,2));
    newReg = Region([regMode.init.name,'_',regGoal.name],newRegVert(newRegConvHullIdx,:));
    newReg = intersect(newReg,reg(aut.label{vertcat(aut.state{:}) == indexState}));
    
    % Subtract the underapproximated reactive funnel
    %     newRegVert = [];
    %     ellAcReact = projection(ac,sys);
    %     ttmp = ac.x0.getTimeVec();
    %     for j = length(ttmp):-10:1
    %         newRegVert = [buildNewRegion(ellAcReact(j), false); newRegVert];
    %     end
    newRegVert = buildNewRegion(ellAcNext(idxLast), true);
    
    [newRegConvHullIdx] = convhull(newRegVert(:,1),newRegVert(:,2));
    subReg = Region([regMode.init.name,'_',regGoal.name],newRegVert(newRegConvHullIdx,:));
    subReg = intersect(subReg,reg(aut.label{vertcat(aut.state{:}) == indexState}));
    
    newReg = regiondiff(newReg.p,subReg.p);
    
    %reg = [reg; newReg];
    
    % plot it!
    plot(newReg,'m')
    
    %newRegVert = extreme(newReg);
    %newRegArray = [newRegArray; newReg];
    
    % update the region file
    %addNewRegionToFile(newRegVert, aut, reg, trans, indexTrans, fileName)    
    
    % update the abstraction via the specification
    %modifySpecForNewRegion(aut, trans, indexTrans, fileName)
    
    % add to the set of inward funnels
    ac_inward{count} = ac;
    bc_inward{count} = bc;
    
end

