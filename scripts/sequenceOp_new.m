function [ac_inward] = sequenceOp_new(sysArray,reg,regDefl,regBnd,aut,acTrans,iModeToPatch,options)
%
% Reach operation -- construct transition funnels
%

maxFunnelsInward = options.maxFunnelsInward;
maxFunTrials = options.maxFunTrials;
maxTrajLength = options.maxTrajLength;
maxTrials2 = options.maxTrials2;
trans = vertcat(aut.trans{:});

IsectInAll{iModeToPatch} = 0;%zeros(1,size(qCover{iModeToPatch},1));

ac_inward = [];

for funindx = 1:maxFunTrials
    
    if maxFunnelsInward(iModeToPatch) > 0
        % ==========================
        % Compute inward funnels
        
        funFail = true;
        
        
        i = 1;
        %         while true  % will only increment i if atomic controller generation has succeeded.
        
        %lastTrans = trans(:,2)==iModeToJoin;
        nextTrans = trans(:,1)==iModeToPatch;  % we are using a deterministic strategy, so the transition out of the state is necessarily a singleton
        %acLast = [acTrans{lastTrans}];
        acNext = [acTrans{nextTrans}];
        
        sys = sysArray(aut.f{nextTrans}); % Assign the dynamics according to the transition to be taken.
        
        %regSafeSG = getRegTrans(reg,regBnd,aut,itrans);
        %         regSafeG = getReg(reg,regBnd,aut,iModeSuccessor);
        
        %                 ellTransS = ellipsoid(ac_trans(itrans,1));
        
        Qsav = options.Qrand;
        options.Qrand = 1e-4;
        for trial2 = 1:maxTrials2
            try
                initState = getCenterRand_new(sys,regDefl(aut.q{iModeToPatch}),[],options);
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
                
                path = [initOutput; goalOutput];
                stepSize = 0.05;
                
                figure(3)
                clf
                plot(acNext,sys,3)
                axis equal
                
                [u0,x0] = computeTrajectory(sys,initState,path,type);
                figure(500), hold on, plot(x0,[],500)
                isect = isinside([reg(aut.q{iModeToPatch}),reg(aut.q{iModeToPatch})],sys,downsampleUniformly(x0,options.sampSkipColl)');
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
        
        % Sanity check to determine if the trajectory is a valid one.
        t = x0.getTimeVec;
        if ~isinternal(ac,double(x0,t(end)),'u')
            funFail = true;
        end
        
        if ~funFail
            disp('Computing funnel....')
            
            try
                ac = computeAtomicControllerSegmentDubins(u0,x0,sys,options.sampSkipFun,[]);
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
            
            % Placeholder- check only whether or not the final point is inside the funnel
            if ~isinternal(ac,double(x0,t(end)),'u')
                funFail = true;
            end
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

end

if maxFunnelsInward(iModeToPatch) > 0 && funindx == maxFunTrials
    error('Cannot find feasible funnels for this mode. Consider increasing the tree depth.')
end
