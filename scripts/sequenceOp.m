%
% Sequence Operation -- compute inward-facing funnels
%

IsectIn{imode} = [];
hindx = 1e8;  % reset Halton seed
hindx1 = 1;  % reset Halton seed

% "Avoid regions" are taken as all but the two in the current transition
[vAvoid,vAvoidB,regAvoid,regSafe] = deal([]);
count = 0;
for j = 1:length(reg)
    if ~(aut.q{imode} == j)
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
    if imode == aut.trans{indx}(1),
        count = count+1;
        count1 = 0;
        for ii = 1:size(funnel,1)
            if ~isempty(funnel{ii,indx})
                count1 = count1+1;
                for j = 1:length(funnel{ii,indx}.t)
                    ellBndInv11{count1,count}(j,1).x = funnel{ii,indx}.x(j,:)';
                    ellBndInv11{count1,count}(j,1).P = funnel{ii,indx}.P(:,:,j)/funnel{ii,indx}.rho(1);  % <-- assume rho is constant (here only)
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
            if imode == aut.trans{indx}(1) || imode == aut.trans{indx}(2),
                count = count+1;
                count1 = 0;
                for ii = 1:size(funnel,1)
                    if ~isempty(funnel{ii,indx})
                        count1 = count1+1;
                        for j = 1:length(funnel{ii,indx}.t)
                            ellBndInv12{count1,count}(j,1).x = funnel{ii,indx}.x(j,:)';
                            ellBndInv12{count1,count}(j,1).P = funnel{ii,indx}.P(:,:,j)/funnel{ii,indx}.rho(1);  % <-- assume rho is constant (here only)
                        end
                    end
                end
            end
        end
    else  % build funnels only for the outgoing edges from the current region
        count = 0;
        for indx = 1:length(aut.trans)
            if imode == aut.trans{indx}(1),
                count = count+1;
                count1 = 0;
                for ii = 1:size(funnel,1)
                    if ~isempty(funnel{ii,indx})
                        count1 = count1+1;
                        for j = 1:length(funnel{ii,indx}.t)
                            ellBndInv12{count1,count}(j,1).x = funnel{ii,indx}.x(j,:)';
                            ellBndInv12{count1,count}(j,1).P = funnel{ii,indx}.P(:,:,j)/funnel{ii,indx}.rho(1);  % <-- assume rho is constant (here only)
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
            if ~isempty(funnelIn{ii,imode})
                count1 = count1+1;
                for j = 1:length(funnelIn{ii,imode}.t)
                    ellBndInv111{count1,1}(j,1).x = funnelIn{ii,imode}.x(j,:)';
                    ellBndInv111{count1,1}(j,1).P = funnelIn{ii,imode}.P(:,:,j)/funnelIn{ii,imode}.rho(1);  % <-- assume rho is constant (here only)
                end
            end
        end
    end
    
    Isect1{imode} = [];
    
    for i = 1:maxFunnelsInward(imode)%1:maxFunnelsInward(imode)  % Vol(S\cup Ri) < Vol(Ri) and j < Nlimit
        funFail = true;
        
        for funindx = 1:maxFunTrials
            % while funFail
            for trial1 = 1:maxTrials1
                try
                    disp('... Computing final point....')
                    if k2 == 1
                        [finalPt,hindx1] = getRandom2(vReg{aut.q{imode}},vAvoid,ellBndInv11,[],'discrep','inside',hindx1,Hout,n,limsNonRegState,isCyclic);  % Get a final point inside the existing outgoing funnels
                    else
                        [finalPt,hindx1] = getRandom2(vReg{aut.q{imode}},vAvoid,ellBndInv111,[],'discrep','inside',hindx1,Hout,n,limsNonRegState,isCyclic);  % Get final point
                    end
                    if isempty(finalPt)
                        error('getRandom2: Cannot obtain a point.')
                    end
                    goalOutput = finalPt(1:length(Hout))*Hout;
                    for nonRegTrial = 1:maxNonRegTrials
                        nonRegTrial
                        disp('... Computing initial point....')
                        [initState,hindx] = getRandom2(vReg{aut.q{imode}},vAvoid,[],[],'discrep','insidedisjunct',hindx,Hout,n,limsNonRegState,isCyclic);  % Get initial point
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
                        if ~isectNonReg && length(t) > 1 && length(t) < 300, break, end
                    end
                    isect = checkIntersection(vAvoid,vBnd{:},[],Xk,Hout);
                catch error
                    rethrow(error)
                    isect = true;
                end
                if ~isect && ~isectNonReg && length(t) > 1 && length(t) < 300,
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
                if k2 == 1
                    [funnelI,rhoMin,rho_d,info,redindx] = computeFunnel(f0,t,Xk,Uk,S0,Q,R,Hout,n,isCyclic,sampSkipFun,regBnd,reg{aut.q{imode}},regAvoid,[],[],[],[],[],ellBndInv11);
                else
                    [funnelI,rhoMin,rho_d,info,redindx] = computeFunnel(f0,t,Xk,Uk,S0,Q,R,Hout,n,isCyclic,sampSkipFun,regBnd,reg{aut.q{imode}},regAvoid,[],[],[],[],[],ellBndInv111);
                end
                % Check if the funnel is misbehaving
                if sum(rho_d < 1e9) > 3 && sum(rho_d == 1e9) < 1 && rhoMin ~= 1e9
                    funFail = false;
                else
                    funFail = true;
                end
                
                if ~funFail
                    isectBnd = false;
                    for j = 1:sampSkipValid:size(funnelI.x,1)
                        tmp = inv(funnelI.P(:,:,j));
                        tmp = (tmp+tmp')/2;
                        E1 = ellipsoid(funnelI.x(j,:)',tmp*rhoMin);
                        E = projection(E1,[Hout; zeros(n-length(Hout),length(Hout))]);
                        figure(5)
                        axis([min(regBnd{1}.v(:,1)) max(regBnd{1}.v(:,1)) min(regBnd{1}.v(:,2)) max(regBnd{1}.v(:,2))])
                        plot(E)
                        if any(intersect(E,regX{1}.hExtB))
                            isectBnd = true;
                            break
                        end
                        for ii = 1:length(regAvoid)
                            tmp1(ii) = size(extreme(regAvoid{ii}.p),1) == size(regAvoid{ii}.v,1);
                        end
                        if all(tmp1)  % avoid regions are convex
                            for ii = 1:length(regAvoid)
                                if any(intersect(E,regAvoid{ii}.pB2))
                                    isectBnd = true;
                                    break
                                end
                            end
                        else  % they're not convex; take the safe regions (these should be convex)
                            for ii = 1:length(regSafe)
                                if any(intersect(E,regSafe{ii}.hBN2))
                                    isectBnd = true;
                                    break
                                end
                            end
                        end
                    end
                    if ~isectBnd
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
                E1 = ellipsoid(funnelI.x(j,:)',tmp*rhoMin);
                E = projection(E1,[Hout; zeros(n-length(Hout),length(Hout))]);
                ellArrayCurr = [ellArrayCurr; E1];
            end
            clear pAvoid
            %         for mm = 1:length(vAvoid), pAvoid(mm) = polytope(vAvoid{mm}); end
            %         isect = intersect(projection(ellArrayCurr,[Hout; zeros(n-length(Hout),length(Hout))]),pAvoid,'u')
            %         if ~any(isect)
            plot(projection(ellArrayCurr,[Hout; zeros(n-length(Hout),length(Hout))]))
            drawnow
            
            i1 = i + (k2-1)*maxFunnelsInward(imode);
            funnelIn{i1,imode} = funnelI;
            rhoMinArrayIn{i1,imode} = funnelI.rho;
            ellFunnelIn{i1,imode} = ellArrayCurr;
            trajFunnelIn{i1,imode} = funnelI.x;
            
            % Success for this transition if Init is covered
            disp(['Iteration #',num2str(i),' / ',num2str(maxFunnelsInward(imode))])
            IsectIn{imode}(i1,:) = isinternal_quick(ellFunnelIn{i1,imode},qCover{imode}');  % note: any itransReg is fine; all should have the same points
            IsectInAll{imode} = any([IsectInAll{imode}; IsectIn{imode}],1);
            NcoverAct = sum(IsectInAll{imode});
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
