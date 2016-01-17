function [ac_inward, bc_inward, errTrans, existingReg, newRegArray, reg] = computeInwardFunnel(fileName,sysArray,reg,regDefl,regBnd,aut,acTrans,acIn,iModeToPatch,options)
%
% Compute inward funnels 
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
newRegArray = [];
bc = [];

% ==========================
% Compute inward funnels

funFail = true;

isectReact = false;

lastTrans = trans(:,2)==iModeToPatch;
acLast = [acTrans{lastTrans}];
if ~isempty(acLast)
    acLast = acLast(1); % NB: for simplicity, only support one incoming transition..
end
acNextAll = [];
for itrans = find(trans(:,1)==iModeToPatch)'
    acNextAll = [acNextAll; acTrans{itrans}];
end

for find(trans(:,1)==iModeToPatch)'
    
    sys = sysArray(aut.f{itrans});  % For now, delay changing the dynamics until the transition reach tube is reached.
    
    acNext = [acTrans{itrans}];
    
    % ==========================
    % Verify invariance of the transition funnels up to the edge of the region
    x = msspoly('x',3); % instantiate the symbolic variables we will be using in a bit
    
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
    
    idx = idxToChkReactivity{itrans}(1);
    %idx = 700;  % first- bad
    %idx = 400;  % second- bad
    %idx = 1200;  % first- good
    %idx = 850;  % second- good
    %idx = 750;  % first (lab)- more friendly
    %idx = 340;  % second (lab)- more friendly
    
    % plot things
    figure(90), hold on, axis equal
    %plot(reg)
    plot(acNext,sys,90)
    
    figure(3), hold on, axis equal
    plot(reg(aut.q{iModeToPatch}),'r')
%     plot(acLast,sys,3)
    
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
                    if isempty(acLast)
                        acNextFirst = acNextAll(1);
                        ttmpFirst = acNextFirst.x0.getTimeVec();
                        qCenter = double(acNextFirst.x0,ttmpFirst(1));  % in the vicinity of the last point of the funnel
                        finalState = getCenterRand_new(sys,regDefl(aut.q{iModeToPatch}),acNextAll,options,qCenter); %,vReg{aut.q{iModeToPatch}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                        acAcceptCriterion = acNextAll;
                    else
                        ttmpLast = acLast.x0.getTimeVec();
                        qCenter = double(acLast.x0,ttmpLast(end));  % in the vicinity of the last point of the funnel
                        finalState = getCenterRand_new(sys,regDefl(aut.q{iModeToPatch}),acLast,options,qCenter); %,vReg{aut.q{iModeToPatch}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                        acAcceptCriterion = acLast;
                    end
                else  % use the provided inward funnels
                    acInward = acIn{iModeToPatch}(1);  % just get the first one in the set.
                    ttmpLast = acInward.x0.getTimeVec();
                    qCenter = double(acInward.x0,ttmpLast(1));  % in the vicinity of the first point of the funnel
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
                
                x0 = Traject(path.t',path.x');
                u0 = Traject(path.t',path.u');
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
            save 'barrier_test'
            
            try
                
                % [ac,existingRegNew,newRegNew] = computePolytopeAtomicControllerDubins(u0,x0,sys,acNext,reg(aut.q{iModeToPatch}),options.ctrloptions_trans,options.sampSkipFun,xMssExt);
                x00 = x0;  u00 = u0;
                [ac, c] = computeAtomicController(u00,x00,sys,options);
                plot(ac.x0,'k',5)
                rhoi = double(ac.rho,0);
                [res, isectIdx, isectArray] = isinside(ac,reg(aut.q{trans(itrans,1)}),sys);
                if ~res
                    % NB: the following assumes only one contiguous interval where the funnel left the region. 
                    % split the funnel into three parts: 
                    %   - one from the start to the time of the first intersection
                    %   - another during the interval of intersection (this is created using CBFs)  
                    %   - another following the intersection
                    t = ac.x0.getTimeVec();
                    acPre = [];
                    
                    % TODO: put all this into quadraticAC
                    t1 = t(max(isectIdx)+1:end);
                    x1 = Traject(t1,double(ac.x0,t1));
                    u1 = Traject(t1,double(ac.u0,t1));
                    K1 = Traject(t1,double(ac.K,t1));
                    P1 = Traject(t1,double(ac.P,t1));
                    rho1 = Traject(t1,double(ac.rho,t1));
                    Vquad = ac.V(max(isectIdx)+1:end);
                    acPost = QuadraticAC(x1,u1,K1,P1,rho1,Vquad,sys);
                    
                    if false %min(isectIdx) > 1
                        t0 = t(1:min(isectIdx)-1);
                        x0 = Traject(t0,double(ac.x0,t0));
                        u0 = Traject(t0,double(ac.u0,t0));
                        
                        % Compute a new atomic controller that minimizes the funnel (in contrast to the original method that maximizes it). 
                        % The minimization is used here to find a suitable initial condition for the CBF.
                        [acPre] = computeAtomicController(u0,x0,sys,options,rhoi);

                    end
                    
                    % Now, construct the barriers                    
                    t01 = t(min(isectIdx):max(isectIdx));
                    x01 = Traject(t01,double(ac.x0,t01));
                    u01 = Traject(t01,double(ac.u0,t01));
                                        
                    regMode = reg(aut.q{iModeToPatch});
                    
                    % define the initial set as the end of the prefix funnel
                    if ~isempty(acPre)
                        ellToCompose = acPre.ellipsoid(end);
                    else
                        ellToCompose = acNext.ellipsoid(indexToCompose);
                    end
                    
                    [ac, bc] = computeConformingFunnel(u01,x01,sys,regMode,ellToCompose,options);
                    
                    % Plot stuff
                    figure(90)
                    axis equal
                    hold on
                    plot(reg(1),'r')
                    plot(reg(2),'g')
                    
                    plotBarriers(ac,bc)
                    
                    % pause
                    % plot(acNext,sys,90,[],[0,0,1])
                    % plot(acPost,sys,90,[],[0,1,0])
                    plot(acNext.x0,'k',90)
                    plot(ac.x0,'k',90)

                end
                funFail = false;
                
                % TODO: check if inside the polytopes
                
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
        newRegVert = [buildNewRegion(ellAcNext(j), true); newRegVert];
    end
    [newRegConvHullIdx] = convhull(newRegVert(:,1),newRegVert(:,2));
    newReg = Region(newRegVert(newRegConvHullIdx,:));
    newReg = intersect(newReg,reg(aut.q{iModeToPatch}));
    
    % Subtract the underapproximated reactive funnel
    %     newRegVert = [];
    %     ellAcReact = projection(ac,sys);
    %     ttmp = ac.x0.getTimeVec();
    %     for j = length(ttmp):-10:1
    %         newRegVert = [buildNewRegion(ellAcReact(j), false); newRegVert];
    %     end
    newRegVert = buildNewRegion(ellAcNext(idxLast), true);
    
    [newRegConvHullIdx] = convhull(newRegVert(:,1),newRegVert(:,2));
    subReg = Region(newRegVert(newRegConvHullIdx,:));
    subReg = intersect(subReg,reg(aut.q{iModeToPatch}));
    
    newReg = regiondiff(newReg.p,subReg.p);
    
    %reg = [reg; newReg];
    
    % plot it!
    plot(newReg,'m')
    
    newRegVert = extreme(newReg);
    newRegArray = [newRegArray; newReg];
    
    % update the region file
    addNewRegionToFile(newRegVert, aut, reg, fileName)    
    
    % update the abstraction via the specification
    modifySpecForNewRegion(aut, trans, itrans, fileName)
    
    % save the barriers
    bc_inward = [bc_inward; bc];
    
end

