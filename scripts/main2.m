
clear ellFunnel ellInit ellFinal trans

format compact
% clear all
close all
global xyPath x debugFlg

n = 3;  % TODO: conflicts with 'n' for RRT

debugFlg = false;

% Get map and transitions
% fourRegionMap
% threeRegionMap
pursuerEvaderMap

Nterm = 10;  % number of consecutive failures before coverage terminates
coverPct = 0.8;

Ntrans = length(aut.trans);
Nmodes = length(aut.q);

% RRT parameters
stepSize = 1;
n = 100;
radius = 0; %0.2;

maxTrials1 = 100; % set to high value
maxTrials2 = 10;

% for state = 1:Ntrans
%     randPt = getRandom(vReg{state},'discrep');
%     regPoint{state} = randPt;
% end
% regPoint{state+1} = regPoint{1};  %close the loop

%%  Step A
clear t Xk
ellFunnel = [];
ellFunnelC = [];
for i = 1:Ntrans, 
    ellFinal{i} = [];
end

maxFunnels(1:Ntrans) = 40;
% maxFunnels([2 3 5]) = 40;
% maxFunnels([2 3 5 6 8]) = 40;

for k2 = 1:10
    if k2 == 2, break, end
    for itrans = 1:Ntrans
        hindx = 1;  % reset Halton seed
        %     ellInit{itrans} = [];
        %     ellFinal{itrans} = [];
        
        % "Avoid regions" are taken as all but the two in the current transition
        % TODO: make v's obsolete
        [vAvoid,vAvoid1,vAvoid2,vAvoidB,vAvoid1B,vAvoid2B,regAvoid,regAvoid1,regAvoid2] = deal([]);
        
        count = 0;
        for j = 1:length(reg)
            if ~(aut.q{aut.trans{itrans}(1)} == j) && ~(aut.q{aut.trans{itrans}(2)} == j)
                count = count+1;
                vAvoid{count} = reg{j}.v;
                vAvoidB{count} = reg{j}.vB;
                regAvoid{count} = reg{j};
            end
        end
        count = 0;
        for j = 1:length(reg)
            if ~(aut.q{aut.trans{itrans}(1)} == j)
                count = count+1;
                vAvoid1{count} = reg{j}.v;
                vAvoid1B{count} = reg{j}.vB;
                regAvoid1{count} = reg{j};
            end
        end
        count = 0;
        for j = 1:length(reg)
            if ~(aut.q{aut.trans{itrans}(2)} == j)
                count = count+1;
                vAvoid2{count} = reg{j}.v;
                vAvoid2B{count} = reg{j}.vB;
                regAvoid2{count} = reg{j};
            end
        end
        
        if k2 ~= 1  % we're trying to create sequentially-composable reach tubes for this transition
            ellBnd111 = [];
            count1 = 0;
            for ii = 1:size(funnel,1)
                if ~isempty(funnel{ii,itrans,1})
                    count1 = count1+1;
                    for j = 1:length(funnel{ii,itrans,1}.t)
                        tmp = inv(funnel{ii,itrans,1}.P(:,:,j));
                        tmp = (tmp+tmp')/2;
                        ellBnd111{count1,1}(j,1) = ellipsoid(funnel{ii,itrans,1}.x(j,:)',tmp*funnel{ii,itrans,1}.rho(j))';
                    end
                end
            end
        end
        
        Ncover = 10000;
        qReg = aut.q{aut.trans{itrans}(1)};
        qCover{itrans} = getCoverPts(vReg,{qReg},1,Ncover);
        Isect{itrans,1} = [];
        
        for i = 1:maxFunnels(itrans)  % will also break if region is considered covered as per a coverage metric
            % TODO: break if fixpoint is reached
            funFail = true;
            while funFail
                if k2 == 1
                    for trial1 = 1:maxTrials1
                        [randPt,hindx] = getRandom2(vReg{aut.q{aut.trans{itrans}(1)}},vAvoid1,[],[],'discrep',[],hindx);  % Get initial point
                        initXYTh = randPt;
                        initXY = randPt(1:2);
                        finalPt = getCenter(vReg{aut.q{aut.trans{itrans}(2)}},vAvoid2,[]);  % Get final point
                        finalPt
                        goalXY = finalPt(1:2);
                        try
                            disp('Building RRT....')
                            [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                            [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                            isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                        catch error
%                             rethrow(error)
                            isect = true;
                        end
                        if ~isect && length(t) > 1 && length(t) < 300, break, end
                        for trial2 = 1:maxTrials2
                            disp('Trajectory incompatible with constraints; recomputing...')
                            
                            [finalPt,hindx] = getRandom2(vReg{aut.q{aut.trans{itrans}(2)}},vAvoid2,[],[],'discrep',[],hindx);  % Give up on trying to have a nice compact set; randomly select the final point
                            finalPt
                            goalXY = finalPt(1:2);
                            try
                                disp('Building RRT....')
                                [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                                [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                                isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                            catch error
%                                 rethrow(error)
                                isect = true;
                            end
                            if ~isect && length(t) > 1 && length(t) < 300, break, end
                        end
                        if ~isect && length(t) > 1 && length(t) < 300, break, end
                    end
                else
                    try
                        for trial1 = 1:maxTrials1
                            [finalPt,hindx] = getRandom2(vReg{aut.q{aut.trans{itrans}(1)}},vAvoid1,ellBnd111,[],'discrep','inside',hindx);  % Get final point
                            goalXY = finalPt(1:2);
                            % TODO: ensure each final funnel ellipse is contained inside the transition reach tube
                            for thetaTrial = 1:100
                                thetaTrial
                                [randPt,hindx] = getRandom2(vReg{aut.q{aut.trans{itrans}(1)}},vAvoid1,[],[],'discrep',[],hindx);  % Get initial point
                                randPt
                                if isempty(randPt)
                                    disp('getRandom1: Cannot obtain an initial point.')
                                    error('getRandom1: Cannot obtain an initial point.')
                                end
                                initXYTh = randPt;
                                initXY = randPt(1:2);
                                [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                                if ~isempty(pathTmp.q)
                                    [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                                    [isectTh] = checkIntersection4(ellBnd111,Xk(end,:));
                                    if ~isectTh && length(t) > 1 && length(t) < 300, break, end
                                end
                            end
                            isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                            if ~isect && length(t) > 1 && length(t) < 300, break, end
                        end
                    catch error
                        rethrow(error)
                        isect = true;
                    end
                end
                
                if trial1 == maxTrials1 && trial2 == maxTrials2
                    error('Cannot construct reach tube for the current transition.')
                end
                
                [funnelI,rhoMin,rho_d] = computeFunnel(t,Xk,Uk,regBnd,regX{itrans},regAvoid);
                
                % Check if the funnel is misbehaving
                isectBnd = false;
                for j = 1:10:size(funnelI.x,1)
                    tmp = inv(funnelI.P(:,:,j));
                    tmp = (tmp+tmp')/2;
                    E1 = ellipsoid(funnelI.x(j,:)',tmp*rhoMin);
                    E = projection(E1,[1 0 ;0 1 ;0 0]);
                    if any(intersect(E,regX{itrans}.hExtB))
                        isectBnd = true;
                        break
                    end
                    if ~isempty(regX{itrans}.pIntB)
                        if any(intersect(E,regX{itrans}.pIntB2))
                            j
                            isectBnd = true;
                            break
                        end
                    end
                end
                if sum(rho_d < 1e9) > 3 && sum(rho_d == 1e9) < 1 && ~isectBnd
                    funFail = false;
                end
            end
            
            figure(3)
            plot(Xk(:,1),Xk(:,2),'m','LineWidth',2)
            
            figure(4)
            ellArrayCurr = [];
            if debugFlg
                indxSet = 1:1:size(funnelI.x,1);
            else
                indxSet = [1:20:size(funnelI.x,1) size(funnelI.x,1)];
            end
            for j = indxSet
                Psav(:,:,j) = funnelI.P(:,:,j);%Ps{j};
                Ksav(:,:,j) = funnelI.K(:,:,j);%Ks{j};
                tmp = inv(funnelI.P(:,:,j));
                tmp = (tmp+tmp')/2;
                E1 = ellipsoid(funnelI.x(j,:)',tmp*rhoMin);
                E = projection(E1,[1 0 ;0 1 ;0 0]);
                ellArrayCurr = [ellArrayCurr; E1];
            end
            clear pAvoid
            plot(projection(ellArrayCurr,[1 0;0 1; 0 0]))
            
            i1 = i + (k2-1)*maxFunnels(itrans);
            funnel{i1,itrans,1} = funnelI;
            % TODO: the rest shoudl be made obsolete by the funnel struct
            rhoMinArray{i1,itrans,1} = funnelI.rho;
            ellFunnel{i1,itrans,1} = ellArrayCurr;
            trajFunnel{i1,itrans,1} = funnelI.x;
            
            % Success for this transition if Init is covered
            Isect{itrans,1}(i1,:) = isinternal(ellFunnel{i1,itrans,1},qCover{itrans}');
            NcoverAct = sum(any(Isect{itrans,1},1));
            
            disp(['Iteration #',num2str(i)])
            disp([' Percent of region covered: ',num2str(NcoverAct/Ncover),' , ',num2str(NcoverAct),' out of ',num2str(Ncover)]);

%             if NcoverAct > coverPct*Ncover
%                 break
%             end
        end
    end
end

%%

save pursuitEvasion10

%%  Step B
% Now, try to compute some inward-directed funnels
maxFunnels(1:Nmodes) = 80;
maxFunnels([2 4 6]) = 140;

for k2 = 1:10
    if k2 == 2, break, end
    for imode = 1:Nmodes
        hindx = 1e8;  % reset Halton seed
        hindx1 = 1;  % reset Halton seed
        %     ellInit{itrans} = [];
        %     ellFinal{itrans} = [];
        
        % "Avoid regions" are taken as all but the two in the current
        % transition
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
        
        % Get outgoing reach tubes from initial state
        ellBnd1 = [];  ellBnd11 = [];
        count = 0;
        for indx = 1:length(aut.trans)
            if imode == aut.trans{indx}(1),
                count = count+1;
                count1 = 0;
                for ii = 1:size(funnel,1)
                    if ~isempty(funnel{ii,indx,1})
                        count1 = count1+1;
                        for j = 1:length(funnel{ii,indx,1}.t)
                            tmp = inv(funnel{ii,indx,1}.P(:,:,j));
                            tmp = (tmp+tmp')/2;
                            ellBnd11{count1,count}(j,1) = ellipsoid(funnel{ii,indx,1}.x(j,:)',tmp*funnel{ii,indx,1}.rho(j))';
                        end
                        %                     ellBnd1{count1,count} = ellFunnel{ii,indx,1};
                    end
                end
            end
        end
        
        if k2 ~= 1  % we're trying to create sequentially-composable reach tubes for this mode
            ellBnd111 = [];
            count1 = 0;
            for ii = 1:size(funnelIn,1)
                if ~isempty(funnelIn{ii,imode,1})
                    count1 = count1+1;
                    for j = 1:length(funnelIn{ii,imode,1}.t)
                        tmp = inv(funnelIn{ii,imode,1}.P(:,:,j));
                        tmp = (tmp+tmp')/2;
                        ellBnd111{count1,1}(j,1) = ellipsoid(funnelIn{ii,imode,1}.x(j,:)',tmp*funnelIn{ii,imode,1}.rho(j))';
                    end
                end
            end
        end
        
        Isect1{imode} = [];
        
        for i = 1:maxFunnels(imode)  % TODO: replace with while Vol(S\cup Ri) < Vol(Ri) and j < Nlimit
            funFail = true;
            while funFail
                for trial1 = 1:maxTrials1
                    try
                        if k2 == 1
                            [finalPt,hindx1] = getRandom2(vReg{aut.q{imode}},vAvoid,ellBnd11,[],'discrep','inside',hindx1);  % Get a final point inside the existing outgoing funnels
                        else
                            [finalPt,hindx1] = getRandom2(vReg{aut.q{imode}},vAvoid,ellBnd111,[],'discrep','inside',hindx1);  % Get final point
                        end
                        if isempty(finalPt)
                            error('getRandom1: Cannot obtain a point.')
                        end
                        goalXY = finalPt(1:2);
                        for thetaTrial = 1:100
                            thetaTrial
                            [randPt,hindx] = getRandom2(vReg{aut.q{imode}},vAvoid,[],[],'discrep','outside',hindx);  % Get initial point
                            randPt
                            if isempty(randPt)
                                disp('getRandom1: Cannot obtain an initial point.')
                                error('getRandom1: Cannot obtain an initial point.')
                            end
                            initXYTh = randPt;
                            initXY = randPt(1:2);
                            [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                            if ~isempty(pathTmp.q)
                                [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                                if k2 == 1
                                    [isectTh] = checkIntersection4(ellBnd11,Xk(end,:));
                                else
                                    [isectTh] = checkIntersection4(ellBnd111,Xk(end,:));
                                end
                                if ~isectTh && length(t) > 1 && length(t) < 300, break, end
                            end
                        end
                        isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                    catch error
                        isect = true;
                                            rethrow(error)
                    end
                    if ~isect && length(t) > 1 && length(t) < 300, break, end
                end
                
                [funnelI,rhoMin,rho_d] = computeFunnel(t,Xk,Uk,regBnd,reg{aut.q{imode}},regAvoid);
                
                % Check if the funnel is misbehaving
                isectBnd = false;
                for j = 1:10:size(funnelI.x,1)
                    tmp = inv(funnelI.P(:,:,j));
                    tmp = (tmp+tmp')/2;
                    E1 = ellipsoid(funnelI.x(j,:)',tmp*rhoMin);
                    E = projection(E1,[1 0 ;0 1 ;0 0]);
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
                if sum(rho_d < 1e9) > 3 && sum(rho_d == 1e9) < 1 && ~isectBnd
                    funFail = false;
                end
            end
            
            figure(3)
            plot(Xk(:,1),Xk(:,2),'m','LineWidth',2)
            
            figure(4)
            ellArrayCurr = [];
            if debugFlg
                indxSet = 1:1:size(funnelI.x,1);
            else
                indxSet = [1:20:size(funnelI.x,1) size(funnelI.x,1)];
            end
            for j = indxSet
                Psav(:,:,j) = funnelI.P(:,:,j);%Ps{j};
                Ksav(:,:,j) = funnelI.K(:,:,j);%Ks{j};
                tmp = inv(funnelI.P(:,:,j));
                tmp = (tmp+tmp')/2;
                E1 = ellipsoid(funnelI.x(j,:)',tmp*rhoMin);
                E = projection(E1,[1 0 ;0 1 ;0 0]);
                ellArrayCurr = [ellArrayCurr; E1];
            end
            clear pAvoid
            %         for mm = 1:length(vAvoid), pAvoid(mm) = polytope(vAvoid{mm}); end
            %         isect = intersect(projection(ellArrayCurr,[1 0;0 1; 0 0]),pAvoid,'u')
            %         if ~any(isect)
            plot(projection(ellArrayCurr,[1 0;0 1; 0 0]))
            %             figure(7)
            %             plot(ellArrayCurr)
            
            i1 = i + (k2-1)*maxFunnels(imode);
            funnelIn{i1,imode,1} = funnelI;
            rhoMinArrayC{i1,imode,1} = funnelI.rho;
            ellFunnelC{i1,imode,1} = ellArrayCurr;
            trajFunnelC{i1,imode,1} = funnelI.x;
            
            % Success for this transition if Init is covered
            try
                icover = [];
                isect = [];
                for itrans = 1:length(aut.trans)
                    if aut.trans{itrans}(1) == imode
                        itransMode = itrans;
                        %                 for mm = 1:size(ellFunnel,1)
                        %                     isect(mm,:) = isinternal(ellFunnel{mm,itrans,1},qCover{itransReg}');
                        %                 end
                        icover = [icover; any(Isect{itrans,1},1)];
                    end
                end
                Isect1{imode}(i1,:) = isinternal(ellFunnelC{i1,imode,1},qCover{itransMode}');  % note: any itransReg is fine; all should have the same points
                
                insideIntersect = all(icover,1);
                insideCentralRT = any(Isect1{imode},1);
                Ncover1 = sum(~insideIntersect);  % Number of points not in the intersection region
                icover1 = insideCentralRT.*(~insideIntersect(1:length(Isect1{imode})));  % Number of points covered by the inward-directed RT
                %         icover1 = [all(icover,1); any(isect,1)];
                NcoverAct = sum(icover1);
            catch
                
            end
            
            disp(['Iteration #',num2str(i)])
            disp([' Percent of region covered: ',num2str(NcoverAct/Ncover1),' , ',num2str(NcoverAct),' out of ',num2str(Ncover1)]);
%             if NcoverAct > coverPct*Ncover1
%                 break
%             end
        end
    end
end

%%
save pursuitEvasion20

%%  Step C

k = 2;

maxFunnels(1:Ntrans) = 20;
maxFunnels([2 3 5]) = 40;
% maxFunnels([2 3 5 6 8]) = 40;

for k2 = 1:10
    if k2 == 2, break, end
    for itrans = 1:Ntrans
        hindx = 1;  % reset Halton seed
        hindx1 = 1;  % reset Halton seed
        
        %     ellInit{itrans} = [];
        %     ellFinal{itrans} = [];
        % "Avoid regions" are taken as all but the two in the current
        % transition + the ellipsoids
        [vAvoid,vAvoid1,vAvoid2,vAvoidB,vAvoid1B,vAvoid2B,regAvoid,regAvoid1,regAvoid2] = deal([]);
        count = 0;
        for j = 1:length(reg)
            if ~(aut.q{aut.trans{itrans}(1)} == j) && ~(aut.q{aut.trans{itrans}(2)} == j)
                count = count+1;
                vAvoid{count} = reg{j}.v;
                vAvoidB{count} = reg{j}.vB;
                regAvoid{count} = reg{j};
            end
        end
        count = 0;
        for j = 1:length(reg)
            if ~(aut.q{aut.trans{itrans}(1)} == j)
                count = count+1;
                vAvoid1{count} = reg{j}.v;
                vAvoid1B{count} = reg{j}.vB;
                regAvoid1{count} = reg{j};
            end
        end
        count = 0;
        for j = 1:length(reg)
            if ~(aut.q{aut.trans{itrans}(2)} == j)
                count = count+1;
                vAvoid2{count} = reg{j}.v;
                vAvoid2B{count} = reg{j}.vB;
                regAvoid2{count} = reg{j};
            end
        end
        
        % NB: reach tube constructs are m x n cell arrays of funnels (arrays of p
        % ellipsoids), where:
        % m: # of funnels in the transition
        % n: # of outgoing transitions for the current region
        % p: # of ellipsoids making up the funnel
        
        % Get outgoing reach tubes from initial state
        count = 0;
        [ellBnd1,ellBnd2,Vbnd1,Vbnd2,ellC1,ellC2,Vin1,Vin2,ellBnd11,ellBnd21,ellC11,ellC21] = deal([]);
        for indx = 1:length(aut.trans)
            if aut.trans{itrans}(1) == aut.trans{indx}(1)
                count = count+1;
                count1 = 0;
                for ii = 1:size(funnel,1)
                    if ~isempty(funnel{ii,indx,k-1})
                        count1 = count1+1;
                        for j = 1:length(funnel{ii,indx,k-1}.t)
                            tmp = inv(funnel{ii,indx,k-1}.P(:,:,j));
                            tmp = (tmp+tmp')/2;
                            ellBnd11{count1,count}(j,1) = ellipsoid(funnel{ii,indx,k-1}.x(j,:)',tmp*funnel{ii,indx,k-1}.rho(j))';
                        end
                        %                     ellBnd1{count1,count} = ellFunnel{ii,indx,k-1};
                        Vbnd1{count1,count} = ones(size(funnel{ii,indx,k-1}.V)) - funnel{ii,indx,k-1}.V;
                    end
                end
            end
        end
        % Get outgoing reach tubes from final state
        count = 0;
        for indx = 1:length(aut.trans)
            if aut.trans{itrans}(2) == aut.trans{indx}(1)
                count = count+1;
                count1 = 0;
                for ii = 1:size(funnel,1)
                    if ~isempty(funnel{ii,indx,k-1})
                        count1 = count1+1;
                        for j = 1:length(funnel{ii,indx,k-1}.t)
                            tmp = inv(funnel{ii,indx,k-1}.P(:,:,j));
                            tmp = (tmp+tmp')/2;
                            ellBnd21{count1,count}(j,1) = ellipsoid(funnel{ii,indx,k-1}.x(j,:)',tmp*funnel{ii,indx,k-1}.rho(j));
                        end
                        %                     ellBnd2{count1,count} = ellFunnel{ii,indx,k-1};
                        Vbnd2{count1,count} = ones(size(funnel{ii,indx,k-1}.V)) - funnel{ii,indx,k-1}.V;
                    end
                end
            end
        end
        % Get inward-directed funnels for the two regions from the
        % previous iteration
        count1 = 0;
        count2 = 0;
        for ii = 1:size(funnelIn,1)
            if ~isempty(funnelIn{ii,aut.trans{itrans}(1),k-1})
                count1 = count1+1;
                for j = 1:length(funnelIn{ii,aut.trans{itrans}(1),k-1}.t)
                    tmp = inv(funnelIn{ii,aut.trans{itrans}(1),k-1}.P(:,:,j));
                    tmp = (tmp+tmp')/2;
                    ellC11{count1}(j,1) = ellipsoid(funnelIn{ii,aut.trans{itrans}(1),k-1}.x(j,:)',tmp*funnelIn{ii,aut.trans{itrans}(1),k-1}.rho(j));
                end
                ellC1{count1} = ellFunnelC{ii,aut.trans{itrans}(1),k-1};
                Vin1{count1} = ones(size(funnelIn{ii,aut.trans{itrans}(1),k-1}.V)) - funnelIn{ii,aut.trans{itrans}(1),k-1}.V;
            end
            if ~isempty(funnelIn{ii,aut.trans{itrans}(2),k-1})
                count2 = count2+1;
                for j = 1:length(funnelIn{ii,aut.trans{itrans}(2),k-1}.t)
                    tmp = inv(funnelIn{ii,aut.trans{itrans}(2),k-1}.P(:,:,j));
                    tmp = (tmp+tmp')/2;
                    ellC21{count2}(j,1) = ellipsoid(funnelIn{ii,aut.trans{itrans}(2),k-1}.x(j,:)',tmp*funnelIn{ii,aut.trans{itrans}(2),k-1}.rho(j));
                end
                ellC2{count2} = ellFunnelC{ii,aut.trans{itrans}(2),k-1};
                Vin2{count2} = ones(size(funnelIn{ii,aut.trans{itrans}(2),k-1}.V)) - funnelIn{ii,aut.trans{itrans}(2),k-1}.V;
            end
        end
        
        if k2 ~= 1  % we're trying to create sequentially-composable reach tubes for this transition
            ellBnd111 = [];
            count1 = 0;
            for ii = 1:size(funnel,1)
                if ~isempty(funnel{ii,itrans,k})
                    count1 = count1+1;
                    for j = 1:length(funnel{ii,itrans,k}.t)
                        tmp = inv(funnel{ii,itrans,k}.P(:,:,j));
                        tmp = (tmp+tmp')/2;
                        ellBnd111{count1,1}(j,1) = ellipsoid(funnel{ii,itrans,k}.x(j,:)',tmp*funnel{ii,itrans,k}.rho(j))';
                    end
                end
            end
        end
        
        Isect{itrans,k} = [];
        
        % clear entries in ellFunnel for current transition
        if k2 == 1
            for mm = 1:size(funnel,1)
                ellFunnel{mm,itrans,k} = [];
                funnel{mm,itrans,k} = [];
            end
        end
        
        goodIndx = [];
        
        for i = 1:maxFunnels(itrans)  % will also break if region is considered covered as per a coverage metric
            % ellFunnel{:,itrans} = [];  % Clear the reach tube for the current transition
            funFail = true;
            while funFail
                if k2 == 1
                    for trial1 = 1:maxTrials1
                        disp('Computing initial point....')
                        try
                            [randPt,hindx] = getRandom2(vReg{aut.q{aut.trans{itrans}(1)}},vAvoid1,ellBnd11,[],'discrep','inside',hindx);  % Get initial point
                            if isempty(randPt)
                                error('getRandom1: Cannot obtain an initial point.')
                            end
                            initXYTh = randPt;
                            initXY = randPt(1:2);
                            disp('Computing final point....')
                            [finalPt,hindx] = getRandom2(vReg{aut.q{aut.trans{itrans}(2)}},vAvoid2,ellBnd21,ellC21,'discrep','inside',hindx);  % Give up on trying to have a nice compact set; randomly select the final point
                            if isempty(finalPt)
                                error('getRandom1: Cannot obtain a point.')
                            end
                            goalXY = finalPt(1:2);
                            disp('Building RRT....')
                            [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                            [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                            xRed = downsample(Xk,5);
                            disp('Simulating trajectory....')
                            [isect,goodIndx] = checkIntersection3(vBnd{:},vAvoid1,vAvoid2,ellBnd11,ellBnd21,ellC11,ellC21,xRed);
                            if ~isempty(goodIndx) && checkIntersection(vAvoid1B,vBnd{:},[],xRed(goodIndx,:))
                                isect = false;
                                goodIndx1 = find(repmat(xRed(goodIndx,:),size(Xk,1),1)==Xk,1,'first');
                                t(goodIndx1+1:end) = [];
                                Xk(goodIndx1+1:end,:) = [];
                                Uk(goodIndx1+1:end,:) = [];
                            end
                        catch error
%                             rethrow(error)
                            isect = true;
                        end
                        if ~isect && length(t) > 1 && length(t) < 300, break, end
                        for trial2 = 1:maxTrials2
                            disp('Trajectory incompatible with constraints; recomputing...')
                            try
                                disp('... Computing final point....')
                                [finalPt,hindx] = getRandom2(vReg{aut.q{aut.trans{itrans}(2)}},vAvoid2,ellBnd21,ellC11,'discrep','inside',hindx);  
                                if isempty(finalPt)
                                    error('getRandom1: Cannot obtain a point.')
                                end
                                goalXY = finalPt(1:2);
                                disp('Building RRT....')
                                [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                                [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                                xRed = downsample(Xk,5);
                                disp('Simulating trajectory....')
                                [isect,goodIndx] = checkIntersection3(vBnd{:},vAvoid1,vAvoid2,ellBnd11,ellBnd21,ellC11,ellC21,xRed);
                                if ~isempty(goodIndx) && checkIntersection(vAvoid1B,vBnd{:},[],xRed(goodIndx,:))
                                    isect = false;
                                    goodIndx1 = find(repmat(xRed(goodIndx,:),size(Xk,1),1)==Xk,1,'first');
                                    t(goodIndx1+1:end) = [];
                                    Xk(goodIndx1+1:end,:) = [];
                                    Uk(goodIndx1+1:end,:) = [];
                                end
                            catch error
%                                 rethrow(error)
                                isect = true;
                            end
                            if ~isect && length(t) > 1 && length(t) < 300, break, end
                        end
                        if ~isect && length(t) > 1 && length(t) < 300, break, end
                    end
                else  % create sequentially-composable reach tubes
                    try
                        for trial1 = 1:maxTrials1
                            [finalPt,hindx] = getRandom2(vReg{aut.q{aut.trans{itrans}(1)}},vAvoid1,ellBnd111,[],'discrep','inside',hindx);  % Get final point
                            goalXY = finalPt(1:2);
                            % TODO: ensure each final funnel ellipse is contained inside the transition reach tube
                            for thetaTrial = 1:100
                                thetaTrial
                                [randPt,hindx] = getRandom2(vReg{aut.q{aut.trans{itrans}(1)}},vAvoid1,ellBnd11,ellC11,'discrep','inside',hindx);  % Get initial point
                                randPt
                                if isempty(randPt)
                                    disp('getRandom1: Cannot obtain an initial point.')
                                    error('getRandom1: Cannot obtain an initial point.')
                                end
                                initXYTh = randPt;
                                initXY = randPt(1:2);
                                [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                                if ~isempty(pathTmp.q)
                                    [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                                    xRed = downsample(Xk,5);
                                    [isectTh] = checkIntersection4(ellBnd111,Xk(end,:));
                                    if ~isectTh && length(t) > 1 && length(t) < 300, break, end
                                end
                            end
                            [isect,goodIndx] = checkIntersection3(vBnd{:},vAvoid1,vAvoid2,ellBnd11,ellBnd21,ellC11,ellC21,xRed);
                            if ~isect && length(t) > 1 && length(t) < 300, break, end
                        end
                    catch error
                        rethrow(error)
                        isect = true;
                    end
                end
                
                if trial1 == maxTrials1 && trial2 == maxTrials2
                    error('Cannot construct reach tube for the current transition.')
                end
                
                disp('Computing funnel....')
                [funnelI,rhoMin,rho_d] = computeFunnel(t,Xk,Uk,regBnd,regX{itrans},regAvoid1,regAvoid2,Vbnd1,Vbnd2,Vin1,Vin2,ellBnd11,ellBnd21,ellC11,ellC21);
                
                % Check if the funnel is misbehaving
                isectBnd = false;
                for j = 1:10:size(funnelI.x,1)
                    tmp = inv(funnelI.P(:,:,j));
                    tmp = (tmp+tmp')/2;
                    E1 = ellipsoid(funnelI.x(j,:)',tmp*rhoMin);
                    E = projection(E1,[1 0 ;0 1 ;0 0]);
                    if any(intersect(E,regX{itrans}.hExtB))
                        isectBnd = true;
                        break
                    end
                    if ~isempty(regX{itrans}.pIntB)
                        if any(intersect(E,regX{itrans}.pIntB2))
                            isectBnd = true;
                            break
                        end
                    end
                end
                if sum(rho_d < 1e9) > 3 && sum(rho_d == 1e9) < 5 && min(rho_d(rho_d < 1e9)) < 1e3 && ~isectBnd
                    funFail = false;
                end
            end
            
            figure(3)
            plot(Xk(:,1),Xk(:,2),'m','LineWidth',2)
            
            figure(4)
            ellArrayCurr = [];
            if debugFlg
                indxSet = 1:1:size(funnelI.x,1);
            else
                indxSet = [1:20:size(funnelI.x,1) size(funnelI.x,1)];
            end
            for j = indxSet
                Psav(:,:,j) = funnelI.P(:,:,j);%Ps{j};
                Ksav(:,:,j) = funnelI.K(:,:,j);%Ks{j};
                tmp = inv(funnelI.P(:,:,j));
                tmp = (tmp+tmp')/2;
                E1 = ellipsoid(funnelI.x(j,:)',tmp*rhoMin);
                E = projection(E1,[1 0 ;0 1 ;0 0]);
                ellArrayCurr = [ellArrayCurr; E1];
            end
            
            plot(projection(ellArrayCurr,[1 0;0 1; 0 0]))
            %         figure(7)
            %         plot(ellArrayCurr)
            
            i1 = i + (k2-1)*maxFunnels(itrans);
            funnel{i1,itrans,k} = funnelI;
            rhoMinArray{i1,itrans,k} = funnelI.rho;
            ellFunnel{i1,itrans,k} = ellArrayCurr;
            trajFunnel{i1,itrans,k} = funnelI.x;
            
            try
                % Success for this transition if Init is covered
                icover = [];
                %             isect = [];
                for itrans1 = 1:length(aut.trans)
                    if aut.trans{itrans1}(1) == aut.trans{itrans}(1)
                        itransMode = itrans1;
                        icover = [icover; any(Isect{itrans,k-1},1)];
                    end
                end
                Isect{itrans,k}(i1,:) = isinternal(ellFunnel{i1,itrans,k},qCover{itransMode}');  % covereage of current reach tube
                
                insideIntersect = all(icover,1);
                insideCurrRT = any(Isect{itrans,k},1);
                Ncover1 = sum(insideIntersect);  % Number of points within the intersection region
                icover1 = insideCurrRT.*(insideIntersect(1:length(Isect{itransMode})));  % Number of points covered by the inward-directed RT
                NcoverAct = sum(icover1);
            catch
                
            end
            
            disp(['Iteration #',num2str(i)])
            disp([' Percent of region covered: ',num2str(NcoverAct/Ncover1),' , ',num2str(NcoverAct),' out of ',num2str(Ncover1)]);
%             if NcoverAct > coverPct*Ncover1
%                 break
%             end
        end
        save threeRegions30_2
        % save pursuitEvasion30
    end
end

%%
save pursuitEvasion30

%%  Step D
% Now, try to compute some inward-directed funnels
maxFunnels(1:Nmodes) = 80;
% maxFunnels([2 4]) = 100;

for k2 = 1:10
    if k2 == 2, break, end
    for imode = 1:Nmodes
        hindx = 1e8;  % reset Halton seed
        hindx1 = 1;  % reset Halton seed
        %     ellInit{itrans} = [];
        %     ellFinal{itrans} = [];
        
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
        
        % Get outgoing reach tubes from initial state
        ellBnd1 = [];  ellBnd11 = [];
        count = 0;
        for indx = 1:length(aut.trans)
            if imode == aut.trans{indx}(1),
                count = count+1;
                count1 = 0;
                for ii = 1:size(funnel,1)
                    if ~isempty(funnel{ii,indx,k})
                        count1 = count1+1;
                        for j = 1:length(funnel{ii,indx,k}.t)
                            tmp = inv(funnel{ii,indx,k}.P(:,:,j));
                            tmp = (tmp+tmp')/2;
                            ellBnd11{count1,count}(j,1) = ellipsoid(funnel{ii,indx,k}.x(j,:)',tmp*funnel{ii,indx,k}.rho(j))';
                        end
                        ellBnd1{count1,count} = ellFunnel{ii,indx,k};
                    end
                end
            end
        end
        
        clear ellBnd111
        if k2 ~= 1  % we're trying to create sequentially-composable reach tubes for this mode
            ellBnd111 = [];
            count1 = 0;
            for ii = 1:size(funnelIn,1)
                if ~isempty(funnelIn{ii,imode,k})
                    count1 = count1+1;
                    for j = 1:length(funnelIn{ii,imode,k}.t)
                        tmp = inv(funnelIn{ii,imode,k}.P(:,:,j));
                        tmp = (tmp+tmp')/2;
                        ellBnd111{count1,1}(j,1) = ellipsoid(funnelIn{ii,imode,k}.x(j,:)',tmp*funnelIn{ii,imode,k}.rho(j))';
                    end
                end
            end
        end
        
        Isect1{imode} = [];
        
        for i = 1:maxFunnels(imode)  % TODO: replace with while Vol(S\cup Ri) < Vol(Ri) and j < Nlimit
            funFail = true;
            while funFail
                for trial1 = 1:maxTrials1
                    try
                        if k2 == 1
                            [finalPt,hindx1] = getRandom2(vReg{aut.q{imode}},vAvoid,ellBnd11,[],'discrep','inside',hindx1);  % Get a final point inside the existing outgoing funnels
                        else
                            [finalPt,hindx1] = getRandom2(vReg{aut.q{imode}},vAvoid,ellBnd111,[],'discrep','inside',hindx1);  % Get final point
                            
                        end
                        if isempty(finalPt)
                            error('getRandom1: Cannot obtain a point.')
                        end
                        goalXY = finalPt(1:2);
                        for thetaTrial = 1:100
                            thetaTrial
                            [randPt,hindx] = getRandom2(vReg{aut.q{imode}},vAvoid,[],[],'discrep','outside',hindx);  % Get initial point
                            randPt
                            if isempty(randPt)
                                disp('getRandom1: Cannot obtain an initial point.')
                                error('getRandom1: Cannot obtain an initial point.')
                            end
                            initXYTh = randPt;
                            initXY = randPt(1:2);
                            [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                            if ~isempty(pathTmp.q)
                                [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                                if k2 == 1
                                    [isectTh] = checkIntersection4(ellBnd11,Xk(end,:));
                                else
                                    [isectTh] = checkIntersection4(ellBnd111,Xk(end,:));
                                end
                                if ~isectTh && length(t) > 1 && length(t) < 300, break, end
                            end
                        end
                        isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                    catch error
                        %                     rethrow(error)
                        isect = true;
                    end
                    if ~isect && ~isectTh && length(t) > 1 && length(t) < 300, break, end
                end
                
                [funnelI,rhoMin,rho_d,info] = computeFunnel(t,Xk,Uk,regBnd,reg{aut.q{imode}},regAvoid);
                
                % Check if the funnel is misbehaving
                isectBnd = false;
                for j = 1:10:size(funnelI.x,1)
                    tmp = inv(funnelI.P(:,:,j));
                    tmp = (tmp+tmp')/2;
                    E1 = ellipsoid(funnelI.x(j,:)',tmp*rhoMin);
                    E = projection(E1,[1 0 ;0 1 ;0 0]);
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
                if sum(rho_d < 1e9) > 3 && sum(rho_d == 1e9) < 5 && ~isectBnd
                    funFail = false;
                end
            end
            
            figure(3)
            plot(Xk(:,1),Xk(:,2),'m','LineWidth',2)
            
            figure(4)
            ellArrayCurr = [];
            if debugFlg
                indxSet = 1:1:size(funnelI.x,1);
            else
                indxSet = [1:20:size(funnelI.x,1) size(funnelI.x,1)];
            end
            for j = indxSet
                Psav(:,:,j) = funnelI.P(:,:,j);%Ps{j};
                Ksav(:,:,j) = funnelI.K(:,:,j);%Ks{j};
                tmp = inv(funnelI.P(:,:,j));
                tmp = (tmp+tmp')/2;
                E1 = ellipsoid(funnelI.x(j,:)',tmp*rhoMin);
                E = projection(E1,[1 0 ;0 1 ;0 0]);
                ellArrayCurr = [ellArrayCurr; E1];
            end
            clear pAvoid
            %         for mm = 1:length(vAvoid), pAvoid(mm) = polytope(vAvoid{mm}); end
            %         isect = intersect(projection(ellArrayCurr,[1 0;0 1; 0 0]),pAvoid,'u')
            %         if ~any(isect)
            plot(projection(ellArrayCurr,[1 0;0 1; 0 0]))
            %             figure(7)
            %             plot(ellArrayCurr)
            
            i1 = i + (k2-1)*maxFunnels(imode);
            funnelIn{i1,imode,k} = funnelI;
            rhoMinArrayC{i1,imode,k} = funnelI.rho;
            ellFunnelC{i1,imode,k} = ellArrayCurr;
            trajFunnelC{i1,imode,k} = funnelI.x;
            
            % Success for this transition if Init is covered
            icover = [];
            %             isect = [];
            for itrans = 1:length(aut.trans)
                if aut.trans{itrans}(1) == aut.q{imode}
                    itransMode = itrans;
                    %                 for mm = 1:size(ellFunnel,1)
                    %                     isect(mm,:) = isinternal(ellFunnel{mm,itrans,1},qCover{itransMode}');
                    %                 end
                    icover = [icover; any(Isect{itrans,k},1)];
                end
            end
            Isect1{imode}(i1,:) = isinternal(ellFunnelC{i1,imode,k},qCover{itransMode}');  % note: any itransReg is fine; all should have the same points
            
            insideIntersect = all(icover,1);
            insideCentralRT = any(Isect1{imode},1);
            Ncover1 = sum(~insideIntersect);  % Number of points not in the intersection region
            icover1 = insideCentralRT(1:length(insideIntersect)).*(~insideIntersect);  % Number of points covered by the inward-directed RT
            %         icover1 = [all(icover,1); any(isect,1)];
            NcoverAct = sum(icover1);
            
            disp(['Iteration #',num2str(i)])
            disp([' Percent of region covered: ',num2str(NcoverAct/Ncover1),' , ',num2str(NcoverAct),' out of ',num2str(Ncover1)]);
            if NcoverAct > coverPct*Ncover1
                break
            end
        end
    end
end

%%
save pursuitEvasion40

