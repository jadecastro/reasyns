% 
% Join Operation -- join final point of transition funnels to the
% reactively-composable set
% 

joinedModesAll = [];

hindx = 1e8;  % reset Halton seed
hindx1 = 1;  % reset Halton seed

% Avoid regions for starting region
% TODO: helper function
[vAvoid1,vAvoid1B,regAvoid1] = deal([]);
count = 0;
for j = 1:length(reg)
    if ~(aut.q{iModeToJoin} == j)
        count = count+1;
        vAvoid1{count} = reg{j}.v;
        vAvoid1B{count} = reg{j}.vB;
        regAvoid1{count} = reg{j};
    else
        regSafe{1} = reg{j};
    end
end

[Vbnd1,ellBndInv11] = deal([]);
% Get transition reach tubes from initial state
count = 0;
% TODO: helper function
for indx = 1:length(aut.trans)
    if iModeToJoin == aut.trans{indx}(1)  % when isReactive is false, we can ignore both ellCInv11 and ellBndInv11 in funnel construction
        count = count+1;
        count1 = 0;
        for ii = 1:size(funnel,1)
            if ~isempty(funnel{ii,indx})
                count1 = count1+1;
                jj = 0;
                for j = 1:sampSkipColl:length(funnel{ii,indx}.t)
                    jj = jj+1;
                    ellBndInv11{count1,count}(jj,1).x = funnel{ii,indx}.x(j,:)';
                    ellBndInv11{count1,count}(jj,1).P = funnel{ii,indx}.P(:,:,j)/funnel{ii,indx}.rho(j);  % <-- assume rho is constant (here only)
                end
                Vbnd1{count1,count} = ones(size(funnel{ii,indx}.V)) - funnel{ii,indx}.V;
            end
        end
    end
end
% Get inward-directed funnels for the two regions from the
% previous iteration
count1 = 0;
[Vin1,ellCInv1,ellCInv11] = deal([]);
for ii = 1:size(funnelIn,1)
    if ~isempty(funnelIn{ii,iModeToJoin})
        count1 = count1+1;
        jj = 0;
        for j = 1:sampSkipColl:length(funnelIn{ii,iModeToJoin}.t)
            jj = jj+1;
            ellCInv11{count1}(jj,1).x = funnelIn{ii,iModeToJoin}.x(j,:)';
            ellCInv11{count1}(jj,1).P = funnelIn{ii,iModeToJoin}.P(:,:,j)/funnelIn{ii,iModeToJoin}.rho(j);  % <-- assume rho is constant (here only)
        end
        Vin1{count1} = ones(size(funnelIn{ii,iModeToJoin}.V)) - funnelIn{ii,iModeToJoin}.V;
    end
end

for itrans = find(transTmp(:,2)==iModeToJoin)'  % do for all transition funnels entering q(iModeToJoin)
    
    [vAvoid,vAvoid2,vAvoidB,vAvoid2B,regAvoid,regAvoid2] = deal([]);
    % TODO: helper function
    count = 0;
    for j = 1:length(reg)
        if ~(aut.q{aut.trans{itrans}(1)} == j) && ~(aut.q{aut.trans{itrans}(2)} == j)
            count = count+1;
            vAvoid{count} = reg{j}.v;
            vAvoidB{count} = reg{j}.vB;
            regAvoid{count} = reg{j};
        end
    end
    % TODO: helper function
    count = 0;
    for j = 1:length(reg)
        if ~(aut.q{aut.trans{itrans}(2)} == j)
            count = count+1;
            vAvoid2{count} = reg{j}.v;
            vAvoid2B{count} = reg{j}.vB;
            regAvoid2{count} = reg{j};
        end
    end
    
    % Get transition reach tubes from final state
    % TODO: helper function
    [ellBnd1,ellBnd2,Vbnd2,ellC2,Vin2,ellBndInv21,ellCInv21] = deal([]);
    count1 = 0;
    for ii = 1:size(funnel,1)
        if ~isempty(funnel{ii,itrans})
            count1 = count1+1;
            for j = 1:length(funnel{ii,itrans}.t)
                ellBndInv21{count1,1}(j,1).x = funnel{ii,itrans}.x(j,:)';
                ellBndInv21{count1,1}(j,1).P = funnel{ii,itrans}.P(:,:,j)/funnel{ii,itrans}.rho(j);  % <-- assume rho is constant (here only)
            end
            Vbnd2{count1,1} = ones(size(funnel{ii,itrans}.V)) - funnel{ii,itrans}.V;
        end
    end
    
    for i = 1:maxFunnelsTrans(iModeToJoin)
        funnelCond{i} = [];
        disp('... Computing initial point as the outlet of the transition funnel from last iteration....')
        initState = funnel{i,itrans}.x(end,:);  % initial point is the final point of the transition funnel
        initState
        initOutput = initState(1:length(Hout))*Hout;
        
        
        % (not yet used):
        ellBndInv12.x = funnel{i,itrans}.x(end,:)';
        ellBndInv12.P = funnel{i,itrans}.P(:,:,end)/funnel{i,itrans}.rho(end);
        
        funFail = true;
        if ~isempty(ellCInv11)
            chkSuccessful = checkIntersection4(ellBndInv11,initState,Hout,n,isCyclic) || checkIntersection4(ellCInv11,initState,Hout,n,isCyclic);
        else
            chkSuccessful = checkIntersection4(ellBndInv11,initState,Hout,n,isCyclic);
        end
        if chkSuccessful % if an inward funnel happens to succceed, then flag as successful without creating a conduit
            for funindx = 1:maxFunTrials
                funindx
                % while funFail
                count1 = 1;  % TODO: transfer to a different set of transition funnels if paths cannot be found
                for trial3 = 1:maxTrials3
                    try
                        disp('... Computing final point....')
                        finalPt = [];
                        while isempty(finalPt)
                            [finalPt,hindx1] = getRandom2(vReg{aut.q{iModeToJoin}},vAvoid,ellBndInv11,ellCInv11,'discrep','inside',hindx1,Hout,n,limsNonRegState,isCyclic);  % Get a final point inside the existing outgoing funnels
                        end
                        %                             finalPt = ellBndInv11{count1,1}(1,1).x';
                        goalOutput = finalPt(1:length(Hout))*Hout;
                        goalOutput
                        pathTmp.q = [initOutput; goalOutput];
                        disp('... Simulating trajectory ....')
                        %[t,Xk,Uk] = genNominalTrajectory(initState,pathTmp.q,Hout,n,isCyclic,modelType);
                        stepSize = 1.0;
                        [path] = buildCtrlSpaceRRT(vBnd{:},vAvoid1,[],ellBndInv11,[],[],[],initState,finalPt,stepSize,Hout,n,isCyclic,modelType);
                        if ~isempty(path)
                            t = path.t;
                            Xk = path.x;
                            Uk = path.u;
                            if ~isempty(t)
                                isect = checkIntersection(vAvoid,vBnd{:},[],Xk,Hout);
                                [isectFinal] = checkIntersection4(ellBndInv11,Xk(end,:),Hout,n,isCyclic);
                            else
                                isect = true;
                                isectFinal = true;
                            end
                        else
                            isect = true;
                            isectFinal = true;
                        end
                    catch error
                        rethrow(error)
                        isect = true;
                    end
                    if ~isect && ~isectFinal && length(t) < maxTrajLength,
                        break,
                    end
                end
                if trial3 == maxTrials3
                    funFail = true;
                    %                             error('Cannot find a feasible trajectory for this mode. Consider increasing the tree depth.')
                else
                    funFail = false;
                end
                
                if ~funFail
                    disp('... Computing funnel....')
                    
                    chopIfInfeas = true;
                    try % to deal with spurious cases where PinvTmp in computeFunnel is not psd.
                        % TODO: separate into two trajectories -- with the split occurring in the "gamma" zone.
                        Qf = 60*eye(n);
                        if chopIfInfeas
                            [funnelI,rhoMin,rho_d,info,redindx] = computeFunnelByDrakeRecursive(t,Uk,Xk,Qf,sampSkipFun);
                        else
                            error('not implemented')
                            %[funnelI,rhoMin,rho_d,info,redindx] = computeFunnel(f0,t,Xk,Uk,S0,Q,R,Hout,n,isCyclic,sampSkipFun,regBnd,reg{aut.q{iModeToJoin}},regAvoid,[],[],[],[],[],ellBndInv11);
                        end
                        
                        %TODO: check if mouth of join funnel doesn't contain the outlet of the inward funnel
                        
                        % Check if the funnel is misbehaving
                        if sum(rho_d < 1e9) > 1 && sum(rho_d == 1e9) < min(5,length(t)) && rhoMin ~= 1e9 || max(rho_d) <= 1e-6
                            funFail = false;
                        else
                            funFail = true;
                        end
                        
                    catch ERR
%                         rethrow(ERR)
                        funFail = true;
                    end
                    
                    if ~funFail
                        isectBnd = chkContainment(funnelI,regAvoid,regSafe,regX{itrans},Hout,n,sampSkipValid);
                        if ~isectBnd && max(rho_d) > 1e-6  % TODO: 'pass' criteria as properties in a class
%                         if ~isectBnd && all(rho_d > 1e-6)  % TODO: 'pass' criteria as properties in a class
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
                
                funnelJoin{i,itrans} = funnelI;
                rhoMinArrayJoin{i,itrans} = funnelI.rho;
                ellFunnelJoin{i,itrans} = ellArrayCurr;
                trajFunnelJoin{i,itrans} = funnelI.x;
                
                disp(['Iteration #',num2str(i),' / ',num2str(maxFunnelsTrans(iModeToJoin))])
            
                joinedModesAll = [joinedModesAll; true]
            else
                joinedModesAll = [joinedModesAll; false]
                break
            end
        end
    end
    itrans
    joinedModesAll
    if ~all(joinedModesAll)
        break
    end
end

iPreModeToJoin = transTmp(itrans,1);
if all(joinedModesAll)
    joinedModes(iModeToJoin) = true  % succeeded for all funnels for each incoming transition to this mode
else
    joinedModes(iModeToJoin) = false  
end