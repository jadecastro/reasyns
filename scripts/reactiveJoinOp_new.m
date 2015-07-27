function [ac_inward, errTrans, existingReg, newReg] = reactiveJoinOp_new(sys,reg,regDefl,regBnd,aut,acIn,iModeToPatch,options)
%
% Reach operation -- construct transition funnels
%

global ME

maxFunnelsTrans = options.maxFunnelsTrans;
maxFunTrials = options.maxFunTrials;
maxTrajLength = options.maxTrajLength;
maxTrials2 = 5; % options.maxTrials2;
trans = vertcat(aut.trans{:});

IsectInAll{iModeToPatch} = 0;%zeros(1,size(qCover{iModeToPatch},1));

% Avoid regions for starting region
% regSafeS = getReg(reg,regBnd,aut,iModeToPatch);

ac_inward = [];
errTrans = [];

% ==========================
% Compute transition funnels

funFail = true;

isectReact = false;

existingReg = [];
newReg = [];

for funindx = 1:maxFunTrials

    lastTrans = trans(:,2)==iModeToPatch;
    acLast = [acIn{lastTrans}]; 
    acLast = acLast(1); % NB: for simplicity, only support one incoming transition.. 
    
    for itrans = find(trans(:,1)==iModeToPatch)'
        
        acNext = [acIn{itrans}];
        
        % TODO: put this in region class
        [H,K] = double(hull([reg(aut.q{trans(itrans,1)}).p,reg(aut.q{trans(itrans,2)}).p]));
        x = msspoly('x',2);
        mssReg = (H*x(1:2)-K)' + eps*sum(x);
        xMssExt = -mssReg;
        
        regSafeSG = getRegTrans(reg,regBnd,aut,itrans);
%         regSafeG = getReg(reg,regBnd,aut,iModeSuccessor);
        
        %                 ellTransS = ellipsoid(ac_trans(itrans,1));
        
        ttmp = getTimeVec(acNext.x0);
        
        idxToChkReactivity{itrans} = 1;
        for j = length(ttmp):-60:20
            if isinside(reg(aut.q{trans(itrans,1)}),sys,ppval(acNext.x0.pp,ttmp(j))')
                idxToChkReactivity{itrans} = j:-60:20;
                break
            end
        end
        
        for idx = idxToChkReactivity{itrans}(1:end)
            idx
            qCenter = ppval(acNext.x0.pp,ttmp(idx));
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
                    qCenter = ppval(acLast.x0.pp,ttmp(end));  % in the vicinity of the first point of the funnel
                    finalState = getCenterRand_new(sys,regDefl(aut.q{iModeToPatch}),acLast,options,qCenter); %,vReg{aut.q{iModeToPatch}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                    finalState
                    goalOutput = finalState(1:length(sys.H))*sys.H;
                    
                    path = [initOutput; goalOutput];
                    stepSize = 0.1;
%                     [path] = buildCtrlSpaceRRT(regBnd.v,{reg(aut.q{trans(itrans,1)}).v},{reg(aut.q{trans(itrans,1)}).v},acIn{itrans}.ellipsoid,[],[],[],initState,finalState,stepSize,sys,reg(aut.q{iModeToPatch}),acLast,options);
%                     [path] = buildCtrlSpaceRRT(regBnd.v,{reg(aut.q{trans(itrans,1)}).v},{reg(aut.q{trans(itrans,1)}).v},acIn{itrans}.ellipsoid,[],[],[],initState,finalState,stepSize,sys,[reg(3), reg(4), reg(5)],acLast,options);
                    [path] = buildCtrlSpaceRRT(regBnd.v,{reg(aut.q{iModeToPatch}).v},{reg(aut.q{iModeToPatch}).v},acIn{itrans}.ellipsoid,[],[],[],initState,finalState,stepSize,sys,reg(aut.q{iModeToPatch}),acLast,options);
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
                % error('Cannot find a feasible trajectory for this mode. Consider increasing the tree depth.')
            else
                funFail = false;
            end
            
            if ~funFail
                disp('Computing funnel....')
                
                try
                    [ac,existingRegNew,newRegNew] = computePolytopeAtomicControllerDubins(u0,x0,sys,acNext,reg(aut.q{iModeToPatch}),options.ctrloptions_trans,options.sampSkipFun,xMssExt);
                    funFail = false;
                    
                    % TODO: check if inside the polytopes
                    
                    plot(ac.x0,'k',5)
                    
                catch ME
%                     rethrow(ME)
                    disp(ME.message)
                    disp('something went wrong with the funnel computation...  kicking back out to the main script.')
                    keyboard
                    errTrans = trans(:,2)==iModeToPatch;
                    break
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
                if funindx >= 2
                    isectReact = true;
                elseif funindx > 1  % build upon the existing regions 
                    existingReg = intersection(existingReg, existingRegNew);
                    newReg = union(newReg, newRegNew);
                else  % funindx = 1; initalize the (likely partially-verified) regions
                    existingReg = existingRegNew;
                    newReg = newRegNew;
                end
                
                if isectReact
                    disp(['Iteration #',num2str(i),' / ',num2str(maxFunnelsTrans(iModeToPatch))])
                    break
                end
            end
        end
        
    end
    
    if isempty(errTrans) && isectReact
        ac_inward = [ac_inward; ac];
        break
    end
    %             if ~funFail
    %                 if i == maxFunnelsTrans(iModeToPatch)
    %                     break
    %                 end
    %                 i = i+1;
    %             end
    %             disp(['    iteration #: ',num2str(i)]);
    
    %         end
    
    if ~funFail || ~isempty(errTrans)
        break
    end
end

if funindx == maxFunTrials
    error('Cannot find feasible funnels for this mode. Consider increasing the tree depth.')
end


