%
% Sequence Operation -- compute inward-facing funnels
%

IsectIn{iModeToPatch} = [];
hindx = 1e8;  % reset Halton seed
hindx1 = 1;  % reset Halton seed

% "Avoid regions" are taken as all but the two in the current transition
[vAvoid,vAvoidB,regAvoid,regSafe] = deal([]);
count = 0;
for j = 1:length(reg)
    if ~(aut.q{iModeToPatch} == j)
        count = count+1;
        vAvoid{count} = reg{j}.v;
        vAvoidB{count} = reg{j}.vB;
        regAvoid{count} = reg{j};
    else
        regSafe{1} = reg{j};
    end
end

% Get transition reach tubes for the current state
ellBndInv11 = [];  ellBnd11 = [];
count = 0;
for indx = 1:length(aut.trans)
    if iModeToPatch == aut.trans{indx}(1),
        count = count+1;
        count1 = 0;
        for ii = 1:size(funnel,1)
            if ~isempty(funnel{ii,indx})
                count1 = count1+1;
                for j = 1:length(funnel{ii,indx}.t)
                    ellBndInv11{count1,count}(j,1).x = funnel{ii,indx}.x(j,:)';
                    ellBndInv11{count1,count}(j,1).P = funnel{ii,indx}.P(:,:,j)/funnel{ii,indx}.rho(j);  % <-- assume rho is constant (here only)
                end
            end
        end
    end
end

for k2 = 1:depthInward
    ellBndInv12 = [];  ellBnd12 = [];
    if k2 == 1  % we're trying to create sequentially-composable reach tubes for the current state;
        %         build funnels for both incoming and outgoing edges for the current region
        count = 0;
        for indx = 1:length(aut.trans)
            if iModeToPatch == aut.trans{indx}(1) || iModeToPatch == aut.trans{indx}(2),
                count = count+1;
                count1 = 0;
                for ii = 1:size(funnel,1)
                    if ~isempty(funnel{ii,indx})
                        count1 = count1+1;
                        for j = 1:length(funnel{ii,indx}.t)
                            ellBndInv12{count1,count}(j,1).x = funnel{ii,indx}.x(j,:)';
                            ellBndInv12{count1,count}(j,1).P = funnel{ii,indx}.P(:,:,j)/funnel{ii,indx}.rho(j);  % <-- assume rho is constant (here only)
                        end
                    end
                end
            end
        end
    else  % build funnels only for the outgoing edges from the current region
        count = 0;
        for indx = 1:length(aut.trans)
            if iModeToPatch == aut.trans{indx}(1),
                count = count+1;
                count1 = 0;
                for ii = 1:size(funnel,1)
                    if ~isempty(funnel{ii,indx})
                        count1 = count1+1;
                        for j = 1:length(funnel{ii,indx}.t)
                            ellBndInv12{count1,count}(j,1).x = funnel{ii,indx}.x(j,:)';
                            ellBndInv12{count1,count}(j,1).P = funnel{ii,indx}.P(:,:,j)/funnel{ii,indx}.rho(j);  % <-- assume rho is constant (here only)
                        end
                    end
                end
            end
        end
    end
    
    clear ellBnd111 ellBndInv111
    if k2 ~= 1  % we're trying to create sequentially-composable reach tubes for the current state
        ellBndInv111 = [];  ellBnd111 = [];
        count1 = 0;
        for ii = 1:size(funnelIn,1)
            if ~isempty(funnelIn{ii,iModeToPatch})
                count1 = count1+1;
                for j = 1:length(funnelIn{ii,iModeToPatch}.t)
                    ellBndInv111{count1,1}(j,1).x = funnelIn{ii,iModeToPatch}.x(j,:)';
                    ellBndInv111{count1,1}(j,1).P = funnelIn{ii,iModeToPatch}.P(:,:,j)/funnelIn{ii,iModeToPatch}.rho(j);  % <-- assume rho is constant (here only)
                end
            end
        end
    end
    
    Isect1{iModeToPatch} = [];
    
    for i = 1:maxFunnelsInward(iModeToPatch)%1:maxFunnelsInward(iModeToPatch)  % Vol(S\cup Ri) < Vol(Ri) and j < Nlimit
        funFail = true;
        
        for funindx = 1:maxFunTrials
            % while funFail
            for trial1 = 1:maxTrials1
                try
                    disp('... Computing final point....')
                    finalPt = [];
                    while isempty(finalPt)
                        if k2 == 1
                            [finalPt,hindx1] = getRandom2(vReg{aut.q{iModeToPatch}},vAvoid,ellBndInv11,[],'discrep','inside',hindx1,Hout,n,limsNonRegState,isCyclic);  % Get a final point inside the existing outgoing funnels
                        else
                            [finalPt,hindx1] = getRandom2(vReg{aut.q{iModeToPatch}},vAvoid,ellBndInv111,[],'discrep','inside',hindx1,Hout,n,limsNonRegState,isCyclic);  % Get final point
                        end
                    end
                    %                     if isempty(finalPt)
                    %                         error('getRandom2: Cannot obtain a point.')
                    %                     end
                    goalOutput = finalPt(1:length(Hout))*Hout;
                    for nonRegTrial = 1:maxNonRegTrials
                        nonRegTrial
                        disp('... Computing initial point....')
                        [initState,hindx] = getRandom2(vReg{aut.q{iModeToPatch}},vAvoid,[],[],'discrep','insidedisjunct',hindx,Hout,n,limsNonRegState,isCyclic);  % Get initial point
                        initState
                        if isempty(initState)
                            disp('getRandom1: Cannot obtain an initial point.')
                            error('getRandom1: Cannot obtain an initial point.')
                        end
                        initOutput = initState(1:length(Hout))*Hout;
                        pathTmp.q = [initOutput; goalOutput];
                        disp('... Simulating trajectory ....')
                        [t,Xk,Uk] = genNominalTrajectory(initState,pathTmp.q,Hout,n,isCyclic,modelType);
                        if k2 == 1
                            [isectNonReg] = checkIntersection4(ellBndInv11,Xk(end,:),Hout,n,isCyclic);
                        else
                            [isectNonReg] = checkIntersection4(ellBndInv111,Xk(end,:),Hout,n,isCyclic);
                        end
                        if ~isectNonReg && length(t) > 1 && length(t) < maxTrajLength, break, end
                    end
                    isect = checkIntersection(vAvoid,vBnd{:},[],Xk,Hout);
                catch error
                    rethrow(error)
                    isect = true;
                end
                if ~isect && ~isectNonReg && length(t) > 1 && length(t) < maxTrajLength,
                    break,
                end
            end
            if trial1 == maxTrials1
                funFail = true;
                %                             error('Cannot find a feasible trajectory for this mode. Consider increasing the tree depth.')
            else
                funFail = false;
            end
            
            tsav = t;
            if ~funFail
                disp('... Computing funnel....')
                chopIfInfeas = true;
                try % to deal with spurious cases where PinvTmp in computeFunnel is not psd.
                    % TODO: separate into two trajectories -- with the split occurring in the "gamma" zone.
                    if k2 == 1
%                         [funnelI,rhoMin,rho_d,info,redindx] = computeFunnel(f0,t,Xk,Uk,S0,Q,R,Hout,n,isCyclic,sampSkipFun,regBnd,reg{aut.q{iModeToPatch}},regAvoid,[],[],[],[],[],ellBndInv11);
                        Qf = 60*eye(n);
                        if chopIfInfeas
                            [funnelI,rhoMin,rho_d,info,redindx] = computeFunnelByDrakeRecursive(t,Uk,Xk,Qf,sampSkipFun);
                        else
                            error('not implemented')
                            %                                         Q = @(t) 1e-4*diag([1 1 1]);
                            %                                         [funnelI,rhoMin,rho_d,info,redindx] = computeFunnel(f0,t,Xk,Uk,S0,Q,R,Hout,n,isCyclic,sampSkipFun,regBnd,reg{aut.q{iModeToReactJoin}},regAvoid1,[],[],[],[],[],ellBndInv11);
                        end
                    else
                        error('not implemented')
                        %                         [funnelI,rhoMin,rho_d,info,redindx] = computeFunnel(f0,t,Xk,Uk,S0,Q,R,Hout,n,isCyclic,sampSkipFun,regBnd,reg{aut.q{iModeToPatch}},regAvoid,[],[],[],[],[],ellBndInv111);
                    end
                    % Check if the funnel is misbehaving
                    if sum(rho_d < 1e9) > 3 && sum(rho_d == 1e9) < 1 && rhoMin ~= 1e9 || max(rho_d) <= 1e-5
                        funFail = false;
                    else
                        funFail = true;
                    end
                catch
                    funFail = true;
                end
                
                if ~funFail
                    isectBnd = chkContainment(funnelI,regAvoid,regSafe,regX{iModeToPatch},Hout,n,sampSkipValid);
                    if ~isectBnd && max(rho_d) > 1e-6
                        funFail = false;
                    else
                        funFail = true;
                    end
                end
            end
            
            if ~funFail,
                break
            end
        end
        
        if funindx == maxFunTrials
            error('Cannot find feasible inward funnels for this mode. Consider increasing the tree depth.')
            %TODO: do another major iteration if k not at maxMajorIter
        end
        
        t = tsav;  % t is being overwritten as an mss symbol (why?!?)
        if ~funFail
            figure(3)
            plot(Xk(:,1),Xk(:,2),'m','LineWidth',2)
            
            figure(4)
            ellArrayCurr = [];
            if debugFlg
                indxSet = 1:1:size(funnelI.x,1);
            else
                clear redindxFun
                for j = 1:length(redindx)
                    redindxFun(j) = find(abs(funnelI.t - t(redindx(j))) == min(abs(funnelI.t - t(redindx(j)))));
                end
                indxSet = redindxFun;
            end
            for j = indxSet
                Psav(:,:,j) = funnelI.P(:,:,j);%Ps{j};
                Ksav(:,:,j) = funnelI.K(:,:,j);%Ks{j};
                tmp = inv(funnelI.P(:,:,j));
                tmp = (tmp+tmp')/2;
                E1 = ellipsoid(funnelI.x(j,:)',tmp*funnelI.rho(j));
                E = projection(E1,[Hout; zeros(n-length(Hout),length(Hout))]);
                ellArrayCurr = [ellArrayCurr; E1];
            end
            clear pAvoid
            %         for mm = 1:length(vAvoid), pAvoid(mm) = polytope(vAvoid{mm}); end
            %         isect = intersect(projection(ellArrayCurr,[Hout; zeros(n-length(Hout),length(Hout))]),pAvoid,'u')
            %         if ~any(isect)
            plot(projection(ellArrayCurr,[Hout; zeros(n-length(Hout),length(Hout))]))
            drawnow
            
            i1 = i + (k2-1)*maxFunnelsInward(iModeToPatch);
            funnelIn{i1,iModeToPatch} = funnelI;
            rhoMinArrayIn{i1,iModeToPatch} = funnelI.rho;
            ellFunnelIn{i1,iModeToPatch} = ellArrayCurr;
            trajFunnelIn{i1,iModeToPatch} = funnelI.x;
            
            % Success for this transition if Init is covered
            disp(['Iteration #',num2str(i),' / ',num2str(maxFunnelsInward(iModeToPatch))])
            IsectIn{iModeToPatch}(i1,:) = isinternal_quick(ellFunnelIn{i1,iModeToPatch},qCover{iModeToPatch}');  % note: any itransReg is fine; all should have the same points
            IsectInAll{iModeToPatch} = any([IsectInAll{iModeToPatch}; IsectIn{iModeToPatch}],1);
            NcoverAct = sum(IsectInAll{iModeToPatch});
            disp([' Percent of region covered: ',num2str(NcoverAct/Ncover),' , ',num2str(NcoverAct),' out of ',num2str(Ncover)]);
            if NcoverAct > coverPct*Ncover
                break
            end
        end
        if funFail
            figure(5)
            plot(Xk(:,1),Xk(:,2),'k','LineWidth',2)
            drawnow
        end
    end
end
