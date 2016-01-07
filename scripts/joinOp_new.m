function [ac_inward, errTrans] = joinOp_new(sysArray,reg,regDefl,regBnd,aut,acTrans,iModeToJoin,options)
%
% Reach operation -- construct transition funnels
%

global ME

maxFunnelsTrans = options.maxFunnelsTrans;
maxFunTrials = options.maxFunTrials;
maxTrajLength = options.maxTrajLength;
maxTrials2 = options.maxTrials2;
trans = vertcat(aut.trans{:});

IsectInAll{iModeToJoin} = 0;%zeros(1,size(qCover{iModeToJoin},1));

% Avoid regions for starting region
% regSafeS = getReg(reg,regBnd,aut,iModeToJoin);

ac_inward = [];
errTrans = [];

% ==========================
% Compute transition funnels

funFail = true;

for funindx = 1:maxFunTrials
    
    i = 1;
    %         while true  % will only increment i if atomic controller generation has succeeded.
    
    %lastTrans = trans(:,2)==iModeToJoin;
    nextTrans = trans(:,1)==iModeToJoin;  % we are using a deterministic strategy, so the transition out of the state is necessarily a singleton
    %acLast = [acTrans{lastTrans}];
    acNext = [acTrans{nextTrans}];
    
    for itrans = find(trans(:,2)==iModeToJoin)'
        
        sys = sysArray(aut.f{nextTrans}); % Assign the dynamics according to the transition to be taken.
        
        regSafeSG = getRegTrans(reg,regBnd,aut,itrans);
%         regSafeG = getReg(reg,regBnd,aut,iModeSuccessor);
        
        %                 ellTransS = ellipsoid(ac_trans(itrans,1));
        
        ttmp = getTimeVec(acTrans{itrans}.x0);
        
        qCenter = double(acTrans{itrans}.x0,ttmp(end));
        Qsav = options.Qrand;
        options.Qrand = 1e-4;
        for trial2 = 1:maxTrials2
            try
                initState = getCenterRand_new(sys,reg(aut.q{iModeToJoin}),[],options,qCenter);
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
        for trial2 = 1:maxTrials2
            
            disp('Computing final point....')
            try
                qCenter = double(acNext.x0,0);  % in the vicinity of the first point of the funnel
%                 finalState = getCenterRand_new(sys,regDefl(aut.q{iModeToJoin}),acNext,options,qCenter) %,vReg{aut.q{iModeToJoin}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                finalState = double(acNext.x0,0)' %,vReg{aut.q{iModeToJoin}},regAvoidS.vBN,vBnd{1}, [],[],Hout,n,limsNonRegState,'rand',Qrand);
                goalOutput = finalState(1:length(sys.H))*sys.H;
                
                % TODO: make this work for both types of dynamics!
                %                 if strfind(func2str(sys.polyMdlFun),'CreateKinematics')
                %                     path = [initOutput; goalOutput];
                %                     type = 'output';
                %                 elseif strfind(func2str(sys.polyMdlFun),'Holonomic')
                %                     path = [initState; finalState];
                %                     type = 'state';
                %                 end
                
                stepSize = 0.2;
                
                figure(3)
                clf
                plot(acNext,sys,3)
                axis equal
                
                [path] = buildCtrlSpaceRRT(regBnd.v,{reg(aut.q{iModeToJoin}).v},{reg(aut.q{iModeToJoin}).v},[],[],[],[],initState,finalState,stepSize,sys,reg(aut.q{iModeToJoin}),acNext,options);
%                 [path] = buildCtrlSpaceRRT(regBnd.v,{reg(aut.q{iModeToJoin}).v},{reg(aut.q{iModeToJoin}).v},acTrans{itrans}.ellipsoid,[],[],[],initState,finalState,stepSize,sys,regBnd,acNext,options);
                disp('Computing nominal trajectory....')
                
                x0 = traject(path.t',path.x');
                u0 = traject(path.t',path.u');
                isect = isinside([reg(aut.q{trans(itrans,1)}),reg(aut.q{trans(itrans,2)})],sys,downsampleUniformly(x0,options.sampSkipColl)');
                if isect && length(x0) > 1 && length(x0) < maxTrajLength,
                    break,
                end
                disp('Trajectory incompatible with constraints; recomputing...')
            catch ME
                % rethrow(ME)
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
                ac = computeAtomicControllerSegmentDubins(u0,x0,sys,options.ctrloptions_trans,options.sampSkipFun,[]);
                ac_inward = [ac_inward; ac];
                funFail = false;
                
                % TODO: check if inside the polytopes
                
                plot(ac.x0,'k',5)
                
            catch ME
                %                 rethrow(ME)
                disp(ME.message)
                disp('something went wrong with the funnel computation...  kicking back out to the main script.')
                errTrans = trans(:,2)==iModeToJoin;
                break
            end
            
            % verify that the funnel we obtained meets the containment criteria
            % TODO: create the isinside method
            %                     if ~isinside(ac,regSafeSG,options.sampSkipValid) || max(rho_d) < 1e-6
            %                         funFail = true;
            %                     end
        end
        
        if funFail     % Need to repeat this iteration
            disp('funnel computation failed.')
            % Optional: plot the failed result
            %             plot(ac.x0,'k',5)
            break
        end
        
        plot(ac.x0,'k',3)
        %                 plot(ac_red,'r',4)
        
        %                 i1 = i + (k2-1)*maxFunnelsTrans(iModeToJoin);
        %                 ac_trans(i1,itrans) = ac;
        
        disp(['Iteration #',num2str(i),' / ',num2str(maxFunnelsTrans(iModeToJoin))])
    end
    
    if ~isempty(errTrans)
        break
    end
    %             if ~funFail
    %                 if i == maxFunnelsTrans(iModeToJoin)
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


