function [ac_inward, bc_inward, errTrans, existingReg, newReg] = reactiveJoinOp_new(sysArray,reg,regDefl,regBnd,aut,acTrans,acIn,iModeToPatch,calibMatrix,options)
%
% Reach operation -- construct transition funnels
%

global ME

maxFunnelsTrans = options.maxFunnelsTrans;
maxFunTrials = options.maxFunTrials;
maxTrajLength = options.maxTrajLength;
maxTrials2 = 5; % options.maxTrials2;
trans = vertcat(aut.trans{:});
options

IsectInAll{iModeToPatch} = 0;%zeros(1,size(qCover{iModeToPatch},1));

% Avoid regions for starting region
% regSafeS = getReg(reg,regBnd,aut,iModeToPatch);

ac_inward = [];
bc_inward = [];
errTrans = [];

% ==========================
% Compute inward funnels

funFail = true;

isectReact = false;

lastTrans = trans(:,2)==iModeToPatch;
acLast = [acTrans{lastTrans}];
acLast = acLast(1); % NB: for simplicity, only support one incoming transition..

for itrans = 5%find(trans(:,1)==iModeToPatch)'
    
    sys = sysArray(aut.f{itrans});  % For now, delay changing the dynamics until the transition reach tube is reached.
    
    acNext = [acTrans{itrans}];
    
    % ==========================
    % Verify invariance of the transition funnels up to the edge of the region
    x = msspoly('x',4); % instantiate the symbolic variables we will be using in a bit
    
    % 1. Verify the region boundary itself
    % TODO: extract all the indices for which there is an overlap with the region boundary; compute a hull containing the intersection with the boundary
    existingReg = [];
    newReg = [];
    
    funStepSize = 130;
    ttmp = acNext.x0.getTimeVec();
    idxToChkReactivity{itrans} = 1;
    for j = length(ttmp):-funStepSize:20
        if isinside(reg(aut.q{trans(itrans,1)}),sys,double(acNext.x0,ttmp(j))')
            idxToChkReactivity{itrans} = j:-funStepSize:20;  % you can adjust the step amount smaller to improve the "tightness" of the resulting new region to the true unverified set.
            break
        end
    end
    
    isRegionVerified = false;
    
    % 2. If the region itself is not verified as an invariant, then find one that is
    
    if ~isRegionVerified
        
        % pick any valid final point
        %             qCenter = double(acLast.x0,ttmp(end));  % in the vicinity of the first point of the funnel
        qCenter = double(acNext.x0,0);  % in the vicinity of the first point of the funnel
        finalState = getCenterRand_new(sys,reg(aut.q{iModeToPatch}),acLast,options,qCenter); %,vReg{aut.q{iModeToPatch}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
        goalOutput = finalState(1:length(sys.H))*sys.H
        g = goalOutput;
        global g   % TODO: a smarter way of doing this
        
        % plot things
        figure(90), hold on, axis equal
        plot(reg(aut.q{iModeToPatch}),'r')
        plot(reg(1),'b')
        plot(reg(2),'m')
        plot(reg(3),'y')
        plot(reg(4),'g')
        plot(acNext,sys,90)
        
        figure(3), hold on, axis equal
        plot(reg(aut.q{iModeToPatch}),'r')
        plot(reg(1),'b')
        plot(reg(2),'m')
        plot(reg(3),'y')
        plot(reg(4),'g')
        plot(acLast,sys,3)
        
        % get the first ellipse just inside of the boundary
        % TODO: using ellipses since we don't have the barrier yet, and cannot do an intersection with it.. perhaps there is an underapproximation?
        
        
        % continue iterating backwards until we have found a barrier
        for idx = idxToChkReactivity{itrans}(1:end)
            barrierSuccess = true;
            %                 idx=idx-500
            idx
            tmpArray = acNext.ellipsoid;
            [ellq,ellQ] = double(tmpArray(idx));
            ellq = [ellq; 1];
            ellQ = blkdiag(ellQ, 0.99);
            % TODO: the following is specific only to the problem at hand... we will need to generalize this
            g_X0 = 1 - (x - ellq)'*inv(ellQ)*(x - ellq);  % intial set
            [polyH,polyK] = double(reg(aut.q{iModeToPatch}).p);
            for iplane = 1:length(polyK)
                g_Xu1 = polyH(iplane,:)*x(1:2) - polyK(iplane); % unsafe set
                %                     g_Xu1 = polyH*x(1:2) - polyK; % unsafe set
                for maxAngularVelocity = -0.15:0.3:0.15
                    try
                        % barrierDubinsBoxPushingExamp
                        barrierDubinsBoxPushingExampSingleMode
                        soln = [0 0 0];
                        for isol = 1:length(Bf_sol), soln(isol) = subs(Bf_sol{isol},x,ellq), end
                    catch
                        disp('invariance check failed.')
                        soln = -Inf;
                    end
                    if any(soln > 0)  % TODO: make this check more complete
                        disp('Either the barrier is not nonpositive inside the initial set or the optimizer failed at finding a solution.')
                        barrierSuccess = false;
                    else
                        barrierSuccess = true;
                        break
                    end
                end
                if ~barrierSuccess, break; end
                b{iplane} = B_sol;
            end
            if barrierSuccess == true
                disp('success!')
                idx
                break
            end
        end
    end
    if barrierSuccess == false
        warning('Did not find any certificates for invariance to the region!! It is likely that the entire transition funnel will inevitably cause entry to the successor region under these dynamics.')
    end
    
    % Visualize the region based on invariance only (for debugging purposes)
    idxLast = idx + funStepSize;
    newRegVert = [];
    ellAcNext = projection(acNext,sys);
    for j = length(ttmp):-funStepSize:idxLast
        ell = ellAcNext(j);
        isOverApprox = true;
        buildNewRegion
        newRegVert = [newVertT; newRegVert];
    end
    [newRegConvHullIdx] = convhull(newRegVert(:,1),newRegVert(:,2));
    newReg = region(newRegVert(newRegConvHullIdx,:));
    newReg = intersect(newReg,reg(aut.q{iModeToPatch}));
    plot(newReg,'g')
    
    %%%%%% For testing only!
    
    %     % update the region file
    %     addNewRegionToFile
    %
    %     % update the abstraction via the specification
    %     modifySpecForNewRegion
    
    %%%%%%
    
    
    % ==========================
    % Now, search for funnels that verify reachability from the initial set verified to be invariant
    while true
        idx
        if idx <= 0
            funFail = true;
            disp('Cannot verify reachability from the transition funnel! No reactive funnel was created.')
            break
        end
        qCenter = double(acNext.x0,ttmp(idx));
        Qsav = options.Qrand;
        options.Qrand = 1e-3;
        for trial2 = 1:maxTrials2
            try
                initState = getCenterRand_new(sys,reg(aut.q{iModeToPatch}),[],options,qCenter);
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
        initState
        initOutput = initState(1:length(sys.H))*sys.H;
        
        trial2 = 1;
        funFail = true;
        for trial2 = 1:maxTrials2
            
            disp('Computing final point....')
            try
                if isempty(acIn)
                    ttmpLast = acLast.x0.getTimeVec();
                    qCenter = double(acLast.x0,ttmpLast(end));  % in the vicinity of the last point of the funnel
                    finalState = getCenterRand_new(sys,regDefl(aut.q{iModeToPatch}),acLast,options,qCenter); %,vReg{aut.q{iModeToPatch}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                    acAcceptCriterion = acLast;
                else  % use the provided inward funnels
                    acInward = acIn{iModeToPatch}(1);  % just get the first one in the set.
                    ttmpLast = acInward.x0.getTimeVec();
                    qCenter = double(acInward.x0,ttmpLast(0));  % in the vicinity of the first point of the funnel
                    finalState = getCenterRand_new(sys,regDefl(aut.q{iModeToPatch}),acInward,options,qCenter); %,vReg{aut.q{iModeToPatch}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                    acAcceptCriterion = acInward;
                end
                finalState
                goalOutput = finalState(1:length(sys.H))*sys.H;
                
                % TODO: make this work for both types of dynamics!
                %                 if strfind(func2str(sys.polyMdlFun),'CreateKinematics')
                %                     path = [initOutput; goalOutput];
                %                     type = 'output';
                %                 elseif strfind(func2str(sys.polyMdlFun),'Holonomic')
                %                     path = [initState; finalState];
                %                     type = 'state';
                %                 end

                stepSize = 0.15;
                %  [path] = buildCtrlSpaceRRT(regBnd.v,{reg(aut.q{trans(itrans,1)}).v},{reg(aut.q{trans(itrans,1)}).v},acTrans{itrans}.ellipsoid,[],[],[],initState,finalState,stepSize,sys,reg(aut.q{iModeToPatch}),acLast,options);
                %  [path] = buildCtrlSpaceRRT(regBnd.v,{reg(aut.q{trans(itrans,1)}).v},{reg(aut.q{trans(itrans,1)}).v},acTrans{itrans}.ellipsoid,[],[],[],initState,finalState,stepSize,sys,[reg(3), reg(4), reg(5)],acLast,options);
                [path] = buildReachabilityRRT(regBnd.v,{reg(aut.q{iModeToPatch}).v},{reg(aut.q{iModeToPatch}).v},[],[],[],[],initState,finalState,stepSize,sys,reg(aut.q{iModeToPatch}),acAcceptCriterion,options);
                disp('Computing nominal trajectory....')
                
                x0 = traject(path.t',path.x');
                u0 = traject(path.t',path.u');
                isect = isinside([reg(aut.q{trans(itrans,1)}),reg(aut.q{trans(itrans,2)})],sys,downsampleUniformly(x0,options.sampSkipColl)');
                if isect && length(x0) > 1 && length(x0) < maxTrajLength,
                    funFail = false;
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
            idx = idx - funStepSize;
            % error('Cannot find a feasible trajectory for this mode. Consider increasing the tree depth.')
        else
            funFail = false;
        end
        
        if ~funFail
            disp('Computing funnel....')
            
            try
                
                % [ac,existingRegNew,newRegNew] = computePolytopeAtomicControllerDubins(u0,x0,sys,acNext,reg(aut.q{iModeToPatch}),options.ctrloptions_trans,options.sampSkipFun,xMssExt);
                [ac] = computeAtomicControllerSegmentDubins(u0,x0,sys,options.ctrloptions_trans,options.sampSkipFun,[]);
                funFail = false;
                
                % TODO: check if inside the polytopes
                
                plot(ac.x0,'k',5)
                
            catch ME
                %                     rethrow(ME)
                disp(ME.message)
                disp('something went wrong with the funnel computation... attempting to find a more conservative initial set, and restarting reachability analysis.')
                keyboard
                idx = idx - funStepSize;
                errTrans = trans(:,2)==iModeToPatch;
            end
            
            % verify that the funnel we obtained meets the containment criteria
            % TODO: create the isinside method
            %                     if ~isinside(ac,regSafeSG,options.sampSkipValid) || max(rho_d) < 1e-6
            %                         funFail = true;
            %                     end
        end
        
        if ~funFail && isempty(errTrans)
            
            plot(ac.x0,'k',3)
            
            %                 % check for coverage
            %                 ellArray = ac.ellipsoid;
            %                 ellBndReact{1} = ellArray;  % for all reactive funnels found so far for this transition
            %                 tmp = inv(reactEllBndInv21{1,1}.P);
            %                 Q = (tmp+tmp')/2;
            %                 %                             [xbar,Q] = double(E);  % ellipse #idx in this transition funnel
            %                 qTest = ellipsoidrand(reactEllBndInv21{1,1}.x,Q,100);
            %                 clear badIndx isectReact
            %                 [isectReact,badIndx] = checkIntersection5(regBnd{1}.v,vReg,ellBndReact,[],qTest,sys);
            %                 isectReact = isinside(ac,regSafeSG,options.sampSkipValid);
            %
            
            % for now, just create a few reach sets and then break the loop
            % if funindx >= 2
            %     isectReact = true;
            % elseif funindx > 1  % build upon the existing regions
            %     existingReg = existingReg; % intersection(existingReg, existingRegNew);
            %     newReg = union(newReg, newRegNew);
            % else  % funindx = 1; initalize the (likely partially-verified) regions
            %     existingReg = existingRegNew;
            %     newReg = newRegNew;
            % end
            isectReact = true;
            
            if isectReact  % Success!
                funFail = false;
                %                 disp(['Iteration #',num2str(i),' / ',num2str(maxFunnelsTrans(iModeToPatch))])
                break
            end
        end
    end
    
    if isempty(errTrans) && ~funFail  % if a funnel was generated that passes the test, 
        
        % Create a new region based on the last unverified index of the next funnel.
        idxLast = idx + funStepSize;
        
        % add to the set of inward funnels
        ac_inward = [ac_inward; ac];

        
    else
        disp('No funnel generated for this mode. Updating the input files accordingly.')
        
        % Set the index to the beginning of the next funnel and don't store an atomic controller.
        idxLast = 1;

    end
    
    % Create the new region
    newRegVert = [];
    ellAcNext = projection(acNext,sys);
    for j = length(ttmp):-funStepSize:idxLast
        ell = ellAcNext(j);
        isOverApprox = true;
        buildNewRegion
        newRegVert = [newVertT; newRegVert];
    end
    [newRegConvHullIdx] = convhull(newRegVert(:,1),newRegVert(:,2));
    newReg = region(newRegVert(newRegConvHullIdx,:));
    newReg = intersect(newReg,reg(aut.q{iModeToPatch}));
    
    % Subtract the underapproximated reactive funnel
    newRegVert = [];
    ellAcReact = projection(ac,sys);
    for j = length(ttmp):-funStepSize:idxLast
        ell = ellAcNext(j);
        isOverApprox = false;
        buildNewRegion
        newRegVert = [newVertT; newRegVert];
    end
    [newRegConvHullIdx] = convhull(newRegVert(:,1),newRegVert(:,2));
    subReg = region(newRegVert(newRegConvHullIdx,:));
    subReg = intersect(subReg,reg(aut.q{iModeToPatch}));
    
    newReg = regiondiff(newReg.p,subReg.p);
    
    % plot it!
    plot(newReg,'m')
    
    % update the region file
    addNewRegionToFile
    
    % update the abstraction via the specification
    modifySpecForNewRegion
    
    % save the barriers
    bc_inward = [bc_inward; b];
    
end

