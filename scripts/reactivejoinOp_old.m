% 
% Reactive Join Operation -- join last point of the "gamma"-zone funnels to
% the reactively-composable set
% 

reactJoinedModesAll = [];

%hindx = 1e8;  % reset Halton seed
%hindx1 = 1;  % reset Halton seed

% Avoid regions for starting region
% TODO: helper function
[vAvoid1,vAvoid1B,regAvoid1,vAvoid1Bcomp] = deal([]);
count = 0;
for j = 1:length(reg)
    if ~(aut.q{iModeToJoin} == j)
        count = count+1;
        vAvoid1{count} = reg{j}.v;
        vAvoid1B{count} = reg{j}.vB;
        regAvoid1{count} = reg{j};
    end
end
count = 0;
for j = 1:length(reg)
    if aut.q{iModeToJoin} == j
        count = count+1;
        vAvoid1Bcomp{count} = reg{j}.vB;
    end
end

% Avoid regions for seeding the initial funnel sample.
% TODO: impose a limit on the size of transition funnels so that at
% least one ellipse fits in between this buffer zone. We'll need to
% check if that ellipse is completely contained by r.c. funnels created
% in this step.
vAvoidReactInit = vAvoid1;
for indx = 1:length(vAvoid1Bcomp)
    vAvoidReactInit{length(vAvoid1)+indx} = vAvoid1Bcomp{indx};
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
                for j = 1:length(funnel{ii,indx}.t)
                    ellBndInv11{count1,count}(j,1).x = funnel{ii,indx}.x(j,:)';
                    ellBndInv11{count1,count}(j,1).P = funnel{ii,indx}.P(:,:,j)/funnel{ii,indx}.rho(j);  % <-- assume rho is constant (here only)
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
        for j = 1:length(funnelIn{ii,iModeToJoin}.t)
            ellCInv11{count1}(j,1).x = funnelIn{ii,iModeToJoin}.x(j,:)';
            ellCInv11{count1}(j,1).P = funnelIn{ii,iModeToJoin}.P(:,:,j)/funnelIn{ii,iModeToJoin}.rho(j);  % <-- assume rho is constant (here only)
        end
        Vin1{count1} = ones(size(funnelIn{ii,iModeToJoin}.V)) - funnelIn{ii,iModeToJoin}.V;
    end
end

for itrans = find(transTmp(:,1)==iModeToJoin)'
    
    [vAvoid,vAvoid2,vAvoidB,vAvoid2B,regAvoid,regAvoid2] = deal([]);
    
    % Get transition reach tubes from final state
    % TODO: helper function
    [ellBnd1,ellBnd2,Vbnd2,ellC2,Vin2,ellBndInv21,ellCInv21,reactEllBndInv21,isectReactZone,isolReactZone] = deal([]);
    count1 = 0;
    for ii = 1:size(funnel,1)
        if ~isempty(funnel{ii,itrans})
            count1 = count1+1;
            for j = 1:length(funnel{ii,itrans}.t)
                ellBndInv21{count1,1}(j,1).x = funnel{ii,itrans}.x(j,:)';
                ellBndInv21{count1,1}(j,1).P = funnel{ii,itrans}.P(:,:,j)/funnel{ii,itrans}.rho(j);  % <-- assume rho is constant (here only)
                isectReactZone{count1,1}(j) = ~checkIntersection(vAvoidReactInit,vBnd{:},[],ellBndInv21{count1,1}(j,1).x',Hout);
            end
            if ~any(isectReactZone{count1,1}), error('No points within the reactvity zone found. Increase either the width of the zone or the sample refinement.'), end
            % go backwards from final point and isolate only the last set of points
            % which are about to exit the region
            isolReactZone{count1,1} = zeros(size(isectReactZone{count1,1}));
            reactEllBndInv21{count1,1} = [];
            for j = length(funnel{ii,itrans}.t):-1:1
                if j < length(funnel{ii,itrans}.t) && isectReactZone{count1,1}(j) < isectReactZone{count1,1}(j+1), break, end
                isolReactZone{count1,1}(j) = isectReactZone{count1,1}(j);
                if isolReactZone{count1,1}(j) ~= 0, reactEllBndInv21{count1,1} = [reactEllBndInv21{count1,1}; ellBndInv21{count1,1}(j,1)]; end
            end
            
            Vbnd2{count1,1} = ones(size(funnel{ii,itrans}.V)) - funnel{ii,itrans}.V;
        end
    end
    
    for i = 1:maxFunnelsReactJoin(iModeToJoin)
        funnelCond{i} = [];
        disp('... Computing initial point inside the transition funnel and within gamma of the boundary ....')
        initState = [];
        while isempty(initState)
            initState = getRandom2(vReg{aut.q{iModeToJoin}},vAvoidReactInit,reactEllBndInv21,[],'discrep','inside',hindx1,Hout,n,limsNonRegState,isCyclic);  % initial point is the final point of the transition funnel
            if initState(1)<1
                initState = [];
            end
        end
        initState
        initOutput = initState(1:length(Hout))*Hout;
        
        % (not yet used):
%         ellBndInv12.x = funnel{i,itrans}.x(end,:)';
%         ellBndInv12.P = funnel{i,itrans}.P(:,:,end)/funnel{i,itrans}.rho(end);
        
        funFail = true;
        if checkIntersection4(ellBndInv11,initState,Hout,n,isCyclic) || checkIntersection4(ellCInv11,initState,Hout,n,isCyclic) % if an inward funnel happens to succceed, then flag as successful without creating a conduit
            for funindx = 1:maxFunTrials
                funindx
                % while funFail
                
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
                        [path] = buildCtrlSpaceRRT(vBnd{:},vAvoid1,[],ellBndInv11,[],ellCInv11,[],initState,finalPt,stepSize,Hout,n,isCyclic,modelType);
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
                    [funnelI,rhoMin,rho_d,info,redindx] = computeFunnel(f0,t,Xk,Uk,S0,Q,R,Hout,n,isCyclic,sampSkipFun,regBnd,reg{aut.q{iModeToJoin}},regAvoid1,[],[],[],[],[],ellBndInv11);
                    
                    %TODO: check if mouth of join funnel doesn't contain the outlet of the inward funnel
                    
                    % Check if the funnel is misbehaving
                    if sum(rho_d < 1e9) > 1 && sum(rho_d == 1e9) < min(5,length(t)) && rhoMin ~= 1e9 || max(rho_d) <= 1e-6
                        funFail = false;
                    else
                        funFail = true;
                    end
                    if ~funFail
                        isectBnd = false;
                        for j = 1:sampSkipValid:size(funnelI.x,1)
                            tmp = inv(funnelI.P(:,:,j));
                            tmp = (tmp+tmp')/2;
                            E1 = ellipsoid(funnelI.x(j,:)',tmp*funnelI.rho(j));
                            E = projection(E1,[Hout; zeros(n-length(Hout),length(Hout))]);
                            if any(intersect(E,regX{1}.hExtB))
                                isectBnd = true;
                                break
                            end
                            for iii = 1:length(regAvoid)
                                tmp1(iii) = size(extreme(regAvoid{iii}.p),1) == size(regAvoid{iii}.v,1);
                            end
                            if all(tmp1)  % avoid regions are convex
                                for iii = 1:length(regAvoid)
                                    if any(intersect(E,regAvoid{iii}.pB2))
                                        isectBnd = true;
                                        break
                                    end
                                end
                            else  % they're not convex; take the safe regions (these should be convex)
                                for iii = 1:length(regSafe)
                                    if any(intersect(E,regSafe{iii}.hBN2))
                                        isectBnd = true;
                                        break
                                    end
                                end
                            end
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
                %         for mm = 1:length(vAvoid1), pAvoid(mm) = polytope(vAvoid1{mm}); end
                %         isect = intersect(projection(ellArrayCurr,[Hout; zeros(n-length(Hout),length(Hout))]),pAvoid,'u')
                %         if ~any(isect)
                plot(projection(ellArrayCurr,[Hout; zeros(n-length(Hout),length(Hout))]))
                drawnow
                
                funnelReactJoin{i,itrans} = funnelI;
                rhoMinArrayReactJoin{i,itrans} = funnelI.rho;
                ellFunnelReactJoin{i,itrans} = ellArrayCurr;
                trajFunnelReactJoin{i,itrans} = funnelI.x;
                
                disp(['Iteration #',num2str(i),' / ',num2str(maxFunnelsReactJoin(iModeToJoin))])
                
                reactJoinedModesAll = [reactJoinedModesAll; true]
            else % we've exceeded the number of trial iterations; failure.
                reactJoinedModesAll = [reactJoinedModesAll; false]
                break
            end
        end
    end
end

if all(reactJoinedModesAll)
    reactJoinedModes(iModeToJoin) = true  % succeeded for all funnels for each incoming transition to this mode
else
    reactJoinedModes(iModeToJoin) = false
end
