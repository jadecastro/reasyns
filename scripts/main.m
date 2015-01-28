
clear ellFunnel ellInit ellFinal trans

% clear all
close all
global xyPath x debugFlg

n = 3;  % TODO: conflicts with 'n' for RRT

debugFlg = false;

% Get map and transitions
% fourRegionMap
threeRegionMap

Nterm = 10;  % number of consecutive failures before coverage terminates
coverPct = 0.8;

Ntrans = length(trans);

figure(7)
clf
hold on

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

count = 0;
for itrans = 1:Ntrans
    hindx = 1;  % reset Halton seed
%     ellInit{itrans} = [];
%     ellFinal{itrans} = [];
    
    % "Avoid regions" are taken as all but the two in the current transition
    % TODO: make v's obsolete
    [vAvoid,vAvoid1,vAvoid2,vAvoidB,vAvoid1B,vAvoid2B,regAvoid,regAvoid1,regAvoid2] = deal([]);

    count = 0;
    for j = 1:length(reg)
        if ~(trans{itrans}(1) == j) && ~(trans{itrans}(2) == j)
            count = count+1;
            vAvoid{count} = reg{j}.v;
            vAvoidB{count} = reg{j}.vB;
            regAvoid{count} = reg{j};
        end
    end
    count = 0;
    for j = 1:length(reg)
        if ~(trans{itrans}(1) == j)
            count = count+1;
            vAvoid1{count} = reg{j}.v;
            vAvoid1B{count} = reg{j}.vB;
            regAvoid1{count} = reg{j};
        end
    end
    count = 0;
    for j = 1:length(reg)
        if ~(trans{itrans}(2) == j)
            count = count+1;
            vAvoid2{count} = reg{j}.v;
            vAvoid2B{count} = reg{j}.vB;
            regAvoid2{count} = reg{j};
        end
    end
    vAvoid
    
    Ncover = 10000;
    qCover{itrans} = getCoverPts(vReg,trans,itrans,Ncover);
    Isect{itrans,1} = [];
    
    if itrans == 2 || itrans == 4 || itrans == 6
        maxFunnels = 80;
    else
        maxFunnels = 20;
    end
    for i = 1:maxFunnels  % will also break if region is considered covered as per a coverage metric
        % TODO: break if fixpoint is reached
        funFail = true;
        while funFail
            for trial1 = 1:maxTrials1
                [randPt,hindx] = getRandom1(vReg{trans{itrans}(1)},vAvoid1,[],'discrep',[],hindx);  % Get initial point
                initXYTh = randPt;
                initXY = randPt(1:2);
                
                finalPt = getCenter(vReg{trans{itrans}(2)},vAvoid2,[]);  % Get final point
                finalPt
                goalXY = finalPt(1:2);
                try
                    disp('Building RRT....')
                    [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                    [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
%                     figure(3)
%                     plot(Xk(:,1),Xk(:,2),'b','LineWidth',2)
                    isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                catch error
%                     rethrow(error)
                    isect = true;
                end
                if ~isect && length(t) > 1 && length(t) < 30, break, end
                for trial2 = 1:maxTrials2
                    disp('Trajectory incompatible with constraints; recomputing...')
                    
                    [finalPt,hindx] = getRandom1(vReg{trans{itrans}(2)},vAvoid2,[],'discrep',[],hindx);  % Give up on trying to have a nice compact set; randomly select the final point
                    finalPt
                    goalXY = finalPt(1:2);
                    try
                        disp('Building RRT....')
                        [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                        [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
%                         figure(3)
%                         plot(Xk(:,1),Xk(:,2),'b','LineWidth',2)
                        isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                    catch error
%                         rethrow(error)
                        isect = true;
                    end
                    if ~isect && length(t) > 1 && length(t) < 30, break, end
                end
                if ~isect && length(t) > 1 && length(t) < 30, break, end
            end
            
            if trial1 == maxTrials1 && trial2 == maxTrials2
                error('Cannot construct reach tube for the current transition.')
            end
            
            [funnelI,rhoMin,rho_d] = computeFunnel(t,Xk,Uk,regBnd,regX{itrans},regAvoid);
            if sum(rho_d < 1e9) > 3 && sum(rho_d == 1e9) < 5
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
        
        funnel{i,itrans,1} = funnelI;
        % TODO: the rest shoudl be made obsolete by the funnel struct
        rhoMinArray{i,itrans,1} = funnelI.rho;
        ellFunnel{i,itrans,1} = ellArrayCurr;
        trajFunnel{i,itrans,1} = funnelI.x;
        
        %NB: in actuality, want the intersection of the funnels and the
        %initial set
        %             ellInit{itrans} = [ellInit{itrans}; ellFunnel{i,itrans,1}(1)];
        %             ellFinal{itrans} = [ellFinal{itrans}; ellFunnel{i,itrans,1}(end)];
        %         else
        %             i = i-1;
        %         end
        %         if i > 1
        %             if isinsideunion(ellArrayCurr(1),ellInit{itrans})  % if funnel is completely inside the reach tube, it will not be added
        %                 ellFunnel{i,itrans} = [];
        %                 i = i-1;
        %                 count = count+1;
        %                 if count > Nterm-1
        %                     break
        %                 end
        %             end
        %         end

        % Success for this transition if Init is covered
        Isect{itrans,1}(i,:) = isinternal(ellFunnel{i,itrans,1},qCover{itrans}');
        NcoverAct = sum(any(Isect{itrans,1},1));
        
        disp(['Iteration #',num2str(i)])
        disp([' Percent of region covered: ',num2str(NcoverAct/Ncover),' , ',num2str(NcoverAct),' out of ',num2str(Ncover)]);
        
%         if i > 1 && isinsideunion(ellArrayCurr(1),ellInit{itrans})  % if funnel is completely inside the reach tube, it will not be added
%                 ellFunnel{i,itrans} = [];
%                 i = i-1;
%         end
        
        if NcoverAct > coverPct*Ncover
            break
        end

    end
end

%%
save threeRegions10

%% Remove trivial funnels
% funnel0 = funnel;
% ellFunnel0 = ellFunnel;
% trajFunnel0 = trajFunnel;
% clear elim funnel ellFunnel trajFunnel
% for j = 1:size(funnel0,2)
%     for i = 1:size(funnel0,1)
%         if ~isempty(funnel0{i,j,1})
% %             elim(i) = max(maxeig(ellFunnel0{i,j,1})) < 0.05;
%             elim(i,j) = max(funnel0{i,j,1}.rho) < 1;
%         end
%     end
% end
% % TODO: make the loop below obsolete
% for j = 1:size(funnel0,2)    
%     count = 0;
%     for i = 1:size(funnel0,1)
%         if ~elim(i,j),
%             count = count+1;
%             funnel{count,j,1} = funnel0{i,j,1};
%             ellFunnel{count,j,1} = ellFunnel0{i,j,1};
%             trajFunnel{count,j,1} = trajFunnel0{i,j,1};
%         end
%     end
% end
% 
% %%
% figure(4)
% options.fill = 1;
% options.color = [1 0 0];
% for j = 1:size(funnel,2)
%     for i = 1:size(funnel,1)
%         if ~isempty(ellFunnel{i,j,1}), plot(projection(ellFunnel{i,j,1},[1 0;0 1;0 0]),options); end
%     end
% end
% for j = 1:size(funnel,2)
%     for i = 1:size(funnel,1)
%         if ~isempty(ellFunnel{i,j,1}), plot(trajFunnel{i,j,1}(:,1),trajFunnel{i,j,1}(:,2),'k','LineWidth',2); end
%     end
% end

%%  Step B
% Now, try to compute some inward-directed funnels
count = 0;
for iregion = sort(unique(cell2mat(trans)))
    hindx = 1e8;  % reset Halton seed
    hindx1 = 1;  % reset Halton seed
%     ellInit{itrans} = [];
%     ellFinal{itrans} = [];
    
    % "Avoid regions" are taken as all but the two in the current
    % transition
    [vAvoid,vAvoidB,regAvoid] = deal([]);
    count = 0;
    for j = 1:length(reg)
        if ~(iregion == j)
            count = count+1;
            vAvoid{count} = reg{j}.v;
            vAvoidB{count} = reg{j}.vB;
            regAvoid{count} = reg{j};
        end
    end
    
    % Get outgoing reach tubes from initial state
    ellBnd1 = [];  ellBnd11 = [];
    count = 0;
    for indx = 1:length(trans)
        if iregion == trans{indx}(1),
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
    
    Isect1{iregion} = [];
    
    if iregion == 2
        maxFunnels = 120;
    else
        maxFunnels = 40;
    end
    for i = 1:maxFunnels  % TODO: replace with while Vol(S\cup Ri) < Vol(Ri) and j < Nlimit
        funFail = true;
        while funFail
            for trial1 = 1:maxTrials1
                try
                    [finalPt,hindx1] = getRandom1(vReg{iregion},vAvoid,ellBnd11,'discrep','inside',hindx1);  % Get a final point inside the existing outgoing funnels
                    if isempty(finalPt)
                        error('getRandom1: Cannot obtain a point.')
                    end
                    goalXY = finalPt(1:2);
                    for thetaTrial = 1:100
                        thetaTrial
                        [randPt,hindx] = getRandom1(vReg{iregion},vAvoid,[],'discrep','outside',hindx);  % Get initial point
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
                            [isectTh] = checkIntersection4(ellBnd11,Xk(end,:));
                            if ~isectTh && length(t) > 1 && length(t) < 30, break, end
                        end
                    end
                    isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                catch error
                    isect = true;
%                     rethrow(error)
                end
                if ~isect && length(t) > 1 && length(t) < 30, break, end
            end
            
            [funnelI,rhoMin,rho_d] = computeFunnel(t,Xk,Uk,regBnd,reg{iregion},regAvoid);
            if sum(rho_d < 1e9) > 3 && sum(rho_d == 1e9) < 5
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
            funnelIn{i,iregion,1} = funnelI;
            rhoMinArrayC{i,iregion,1} = funnelI.rho;
            ellFunnelC{i,iregion,1} = ellArrayCurr;
            trajFunnelC{i,iregion,1} = funnelI.x;
            
            %NB: in actuality, want the intersection of the funnels and the
            %initial set
%             ellInit{itrans} = [ellInit{itrans}; ellFunnel{i,itrans}(1)];
%             ellFinal{itrans} = [ellFinal{itrans}; ellFunnel{i,itrans}(end)];
%         else
%             i = i-1;
%         end
        % Success for this transition if Init is covered (TODO: coverage metric)
%         if isinsideunion(ellInit{itrans},vReg{trans{itrans}(1)})
%             break
%         end
%         if i > 1
%             if isinsideunion(ellArrayCurr(1),ellInit{itrans})  % if funnel is completely inside the reach tube, it will not be added
%                 ellFunnel{i,itrans} = [];
%                 i = i-1;
%                 count = count+1;
%                 if count > Nterm-1
%                     break
%                 end
%             end
%         end

        % Success for this transition if Init is covered
        icover = [];
        isect = [];
        for itrans = 1:length(trans)
            if trans{itrans}(1) == iregion
                itransReg = trans{itrans}(1);
%                 for mm = 1:size(ellFunnel,1)
%                     isect(mm,:) = isinternal(ellFunnel{mm,itrans,1},qCover{itransReg}');
%                 end
                icover = [icover; any(Isect{itrans,1},1)];
            end
        end
        Isect1{iregion}(i,:) = isinternal(ellFunnelC{i,iregion,1},qCover{itransReg}');  % note: any itransReg is fine; all should have the same points

        insideIntersect = all(icover,1);
        insideCentralRT = any(Isect1{iregion},1);
        Ncover1 = sum(~insideIntersect);  % Number of points not in the intersection region
        icover1 = insideCentralRT.*(~insideIntersect(1:length(Isect1{iregion})));  % Number of points covered by the inward-directed RT
%         icover1 = [all(icover,1); any(isect,1)];
        NcoverAct = sum(icover1);
        
        disp(['Iteration #',num2str(i)])
        disp([' Percent of region covered: ',num2str(NcoverAct/Ncover1),' , ',num2str(NcoverAct),' out of ',num2str(Ncover1)]);
        if NcoverAct > coverPct*Ncover1
            break
        end
    end
end

%%
save threeRegions20

%%
% funnelIn0 = funnelIn;
% ellFunnelC0 = ellFunnelC;
% trajFunnelC0 = trajFunnelC;
% clear elim funnelIn ellFunnelC trajFunnelC
% for j = 1:size(funnelIn0,2)
%     for i = 1:size(funnelIn0,1)
%         if ~isempty(funnelIn0{i,j,1})
% %             elim(i,j) = max(maxeig(ellFunnelC0{i,j,1})) < 0.05;
%             elim(i,j) = max(funnelIn0{i,j,1}.rho) < 1;
%         end
%     end
% end
% % elim([10 16 18 21 25 34 40],2) = 1;
% 
% for j = 1:size(funnelIn0,2)    
%     count = 0;
%     for i = 1:size(funnelIn0,1)
%         if ~elim(i,j),
%             count = count+1;
%             funnelIn{count,j,1} = funnelIn0{i,j,1};
%             ellFunnelC{count,j,1} = ellFunnelC0{i,j,1};
%             trajFunnelC{count,j,1} = trajFunnelC0{i,j,1};
%         end
%     end
% end
% 
% %%
% figure(4)
% options.fill = 1;
% options.color = [1 0 0];
% for j = 1:size(funnelIn,2)
%     for i = 1:size(funnelIn,1)
%         if ~isempty(ellFunnelC{i,j,1}), plot(projection(ellFunnelC{i,j,1},[1 0;0 1;0 0]),options); end
%     end
% end
% for j = 1:size(funnelIn,2)
%     for i = 1:size(funnelIn,1)
%         if ~isempty(ellFunnelC{i,j,1}), plot(trajFunnelC{i,j,1}(:,1),trajFunnelC{i,j,1}(:,2),'k','LineWidth',2); end
%     end
% end
% 
% % plot initial funnels
% figure(7)
% hold on
% for j = 1:size(funnel,2)
%     for i = 1:size(funnel,1)
%         if ~isempty(ellFunnel{i,j,1}), plot(ellFunnel{i,j}(1)); end
%     end
% end
% 
% for j = 1:size(funnelIn,2)
%     for i = 1:size(funnelIn,1)
%         if ~isempty(ellFunnelC{i,j,1}), plot(ellFunnelC{i,j}(1),'g'); end
%     end
% end
% figure(8)
% hold on
% for j = 1:size(funnel,1)
%     for i = 1:size(funnel,2)
%         if ~isempty(ellFunnel{i,j,1}), plot(ellFunnel{i,j}(end)); end
%     end
% end

%%  Step C
% load testcase2

% reduce the number of ellipses if necessary
% ellFunnelNew = [];
% if length(ellFunnel{1,1,1}) > 50
%     disp('Downsampling the funnels...')
%     for i = 1:size(ellFunnel,1)
%         for j = 1:size(ellFunnel,2)
%             ellFunnelNew{i,j,1} = ellFunnel{i,j,1}(1:10:length(ellFunnel{i,j,1}));
%         end
%     end
%     ellFunnel = ellFunnelNew;
% end

k = 2;
% for k = 2:2  % TODO: coverage
%     hindx = 1;  % reset Halton seed

% Initially set all reach tubes equal to their values from the previous iteration
% for mm = 1:size(funnel,1)
%     for nn = 1:size(funnel,2)
%         ellFunnel{mm,nn,k} = ellFunnel{mm,nn,k-1};
%     end
% end

for itrans = 4%1:Ntrans
    hindx = 1;  % reset Halton seed
    hindx1 = 1;  % reset Halton seed
    
    %     ellInit{itrans} = [];
    %     ellFinal{itrans} = [];
    % "Avoid regions" are taken as all but the two in the current
    % transition + the ellipsoids
    [vAvoid,vAvoid1,vAvoid2,vAvoidB,vAvoid1B,vAvoid2B,regAvoid,regAvoid1,regAvoid2] = deal([]);
    count = 0;
    for j = 1:length(reg)
        if ~(trans{itrans}(1) == j) && ~(trans{itrans}(2) == j)
            count = count+1;
            vAvoid{count} = reg{j}.v;
            vAvoidB{count} = reg{j}.vB;
            regAvoid{count} = reg{j};
        end
    end
    count = 0;
    for j = 1:length(reg)
        if ~(trans{itrans}(1) == j)
            count = count+1;
            vAvoid1{count} = reg{j}.v;
            vAvoid1B{count} = reg{j}.vB;
            regAvoid1{count} = reg{j};
        end
    end
    count = 0;
    for j = 1:length(reg)
        if ~(trans{itrans}(2) == j)
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
    for indx = 1:length(trans)
        if trans{itrans}(1) == trans{indx}(1),
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
    for indx = 1:length(trans)
        if trans{itrans}(2) == trans{indx}(1),
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
        if ~isempty(funnelIn{ii,trans{itrans}(1),k-1})
            count1 = count1+1;
            for j = 1:length(funnelIn{ii,trans{itrans}(1),k-1}.t)
                tmp = inv(funnelIn{ii,trans{itrans}(1),k-1}.P(:,:,j));
                tmp = (tmp+tmp')/2;
                ellC11{count1}(j,1) = ellipsoid(funnelIn{ii,trans{itrans}(1),k-1}.x(j,:)',tmp*funnelIn{ii,trans{itrans}(1),k-1}.rho(j));
            end
            ellC1{count1} = ellFunnelC{ii,trans{itrans}(1),k-1};
            Vin1{count1} = ones(size(funnelIn{ii,trans{itrans}(1),k-1}.V)) - funnelIn{ii,trans{itrans}(1),k-1}.V;
        end
        if ~isempty(funnelIn{ii,trans{itrans}(2),k-1})
            count2 = count2+1;
            for j = 1:length(funnelIn{ii,trans{itrans}(2),k-1}.t)
                tmp = inv(funnelIn{ii,trans{itrans}(2),k-1}.P(:,:,j));
                tmp = (tmp+tmp')/2;
                ellC21{count2}(j,1) = ellipsoid(funnelIn{ii,trans{itrans}(2),k-1}.x(j,:)',tmp*funnelIn{ii,trans{itrans}(2),k-1}.rho(j));
            end
            ellC2{count2} = ellFunnelC{ii,trans{itrans}(2),k-1};
            Vin2{count2} = ones(size(funnelIn{ii,trans{itrans}(2),k-1}.V)) - funnelIn{ii,trans{itrans}(2),k-1}.V;
        end
    end
    Isect{itrans,k} = [];
    
    % clear entries in ellFunnel for current transition
    for mm = 1:size(funnel,1)
        ellFunnel{mm,itrans,k} = [];
        funnel{mm,itrans,k} = [];
    end
    
    goodIndx = [];
    
    if itrans == 2 || itrans == 4 || itrans == 6
        maxFunnels = 40;
    else
        maxFunnels = 20;
    end
    for i = 1:maxFunnels  % will also break if region is considered covered as per a coverage metric
        % ellFunnel{:,itrans} = [];  % Clear the reach tube for the current transition
        funFail = true;
        while funFail
            for trial1 = 1:maxTrials1
                disp('Computing initial point....')
                try
                    [randPt,hindx] = getRandom1(vReg{trans{itrans}(1)},vAvoid1,ellBnd11,'discrep','inside',hindx);  % Get initial point
                    if isempty(randPt)
                        error('getRandom1: Cannot obtain an initial point.')
                    end
                    initXYTh = randPt;
                    initXY = randPt(1:2);
                    disp('Computing final point....')
                    %             finalPt = getCenter(vReg{trans{itrans}(2)},vAvoid2,ellBnd2)  % Get final point
                    [finalPt,hindx] = getRandom1(vReg{trans{itrans}(2)},vAvoid2,ellBnd21,'discrep','inside',hindx);  % Give up on trying to have a nice compact set; randomly select the final point
                    if isempty(finalPt)
                        error('getRandom1: Cannot obtain a point.')
                    end
                    goalXY = finalPt(1:2);
                    disp('Building RRT....')
                    [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                    [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                    disp('Simulating trajectory....')
                    [isect,goodIndx] = checkIntersection3(vBnd{:},vAvoid1,vAvoid2,ellBnd11,ellBnd21,ellC11,ellC21,Xk);
                    if ~isempty(goodIndx)
                        isect = false;
                        t(goodIndx+1:end) = [];
                        Xk(goodIndx+1:end,:) = [];
                        Uk(goodIndx+1:end,:) = [];
                    end
                catch error
%                     rethrow(error)
                    isect = true;
                end
                if ~isect && length(t) > 1 && length(t) < 30, break, end
                for trial2 = 1:maxTrials2
                    disp('Trajectory incompatible with constraints; recomputing...')
                    
                    try
                        disp('... Computing final point....')
                        [finalPt,hindx] = getRandom1(vReg{trans{itrans}(2)},vAvoid2,ellBnd21,'discrep','inside',hindx);  % Give up on trying to have a nice compact set; randomly select the final point
                        if isempty(finalPt)
                            error('getRandom1: Cannot obtain a point.')
                        end
                        goalXY = finalPt(1:2);
                        disp('Building RRT....')
                        [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
                        [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                        disp('Simulating trajectory....')
                        [isect,goodIndx] = checkIntersection3(vBnd{:},vAvoid1,vAvoid2,ellBnd11,ellBnd21,ellC11,ellC21,Xk);
                        if ~isempty(goodIndx)
                            isect = false;
                            t(goodIndx+1:end) = [];
                            Xk(goodIndx+1:end,:) = [];
                            Uk(goodIndx+1:end,:) = [];
                        end
                    catch error
%                         rethrow(error)
                        isect = true;
                    end
                    if ~isect && length(t) > 1 && length(t) < 30, break, end
                end
                if ~isect && length(t) > 1 && length(t) < 30, break, end
            end
            
            if trial1 == maxTrials1 && trial2 == maxTrials2
                error('Cannot construct reach tube for the current transition.')
            end
            
            disp('Computing funnel....')
            [funnelI,rhoMin,rho_d] = computeFunnel(t,Xk,Uk,regBnd,regX{itrans},regAvoid1,regAvoid2,Vbnd1,Vbnd2,Vin1,Vin2,ellBnd11,ellBnd21,ellC11,ellC21);
            if sum(rho_d < 1e9) > 3 && sum(rho_d == 1e9) < 5 && min(rho_d(rho_d < 1e9)) < 1e2
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
        
        funnel{i,itrans,k} = funnelI;
        rhoMinArray{i,itrans,k} = funnelI.rho;
        ellFunnel{i,itrans,k} = ellArrayCurr;
        trajFunnel{i,itrans,k} = funnelI.x;
        
        %NB: in actuality, want the intersection of the funnels and the
        %initial set
        %         ellInit{itrans} = [ellInit{itrans}; ellFunnel{i,itrans}(1)];
        %         ellFinal{itrans} = [ellFinal{itrans}; ellFunnel{i,itrans}(end)];
        
        % Success for this transition if Init is covered (TODO: coverage metric)
        %         if isinsideunion(ellInit{itrans},vReg{trans{itrans}(1)})
        %             break
        %         end
        %         if i > 1
        %             if isinsideunion(ellArrayCurr(1),ellInit{itrans})  % if funnel is completely inside the reach tube, it will not be added
        %                 ellFunnel{i,itrans} = [];
        %                 i = i-1;
        %                 count = count+1;
        %                 if count > Nterm-1
        %                     break
        %                 end
        %             end
        %         end
        
        % Success for this transition if Init is covered
        icover = [];
        %             isect = [];
        for itrans1 = 1:length(trans)
            if trans{itrans1}(1) == iregion
                itransReg = trans{itrans1}(1);
                %                     for mm = 1:size(ellFunnel,1)
                %                         isect(mm,:) = isinternal(ellFunnel{mm,itrans1,k-1},qCover{itransReg}');  % take reach tubes from PREVIOUS iteration
                %                     end
                icover = [icover; any(Isect{itrans,k-1},1)];
            end
        end
        Isect{itrans,k}(i,:) = isinternal(ellFunnel{i,itrans,k},qCover{itransReg}');  % covereage of current reach tube
        
        insideIntersect = all(icover,1);
        insideCurrRT = any(Isect{itrans,k},1);
        Ncover1 = sum(insideIntersect);  % Number of points within the intersection region
        icover1 = insideCurrRT.*(insideIntersect(1:length(Isect{itransReg})));  % Number of points covered by the inward-directed RT
        NcoverAct = sum(icover1);
        
        disp(['Iteration #',num2str(i)])
        disp([' Percent of region covered: ',num2str(NcoverAct/Ncover1),' , ',num2str(NcoverAct),' out of ',num2str(Ncover1)]);
        if NcoverAct > coverPct*Ncover1
            break
        end
    end
end

%%
save threeRegions30

% %%  Remove trivial funnels
% funnel0 = funnel;
% ellFunnel0 = ellFunnel;
% trajFunnel0 = trajFunnel;
% clear ellFunnel trajFunnel
% for j = 1:size(funnel0,2)
%     for i = 1:size(funnel0,k)
%         if ~isempty(funnel0{i,j,k})
%             %             elim(i) = max(maxeig(ellFunnel0{i,j,1})) < 0.05;
%             elim(i,j) = max(funnel{i,j,k}.rho) < 1;
%         end
%     end
% end
% for j = 1:size(funnel0,2)
%     count = 0;
%     for i = 1:size(funnel0,1)
%         for k1 = 1:k-1
%             ellFunnelC{i,j,k1} = ellFunnelC0{i,j,k1};
%             trajFunnelC{i,j,k1} = trajFunnelC0{i,j,k1};
%         end
%         if ~elim(i,j),
%             count = count+1;
%             ellFunnel{count,j,k} = ellFunnel0{i,j,k};
%             trajFunnel{count,j,k} = trajFunnel0{i,j,k};
%         end
%     end
% end
% 
% %%
% figure(4)
% options.fill = 1;
% options.color = [1 0 0];
% for j = 1:size(funnel,2)
%     for i = 1:size(funnel,1)
%         if ~isempty(ellFunnel{i,j,k}), plot(projection(ellFunnel{i,j,k},[1 0;0 1;0 0]),options); end
%     end
% end
% for j = 1:size(funnel,2)
%     for i = 1:size(funnel,1)
%         if ~isempty(ellFunnel{i,j,k}), plot(trajFunnel{i,j,k}(:,1),trajFunnel{i,j,k}(:,2),'k','LineWidth',2); end
%     end
% end


%%  Step D
% Now, try to compute some inward-directed funnels
count = 0;
for iregion = sort(unique(cell2mat(trans)))
    hindx = 1e8;  % reset Halton seed
    hindx1 = 1;  % reset Halton seed
    %     ellInit{itrans} = [];
    %     ellFinal{itrans} = [];
    
    % "Avoid regions" are taken as all but the two in the current transition
    [vAvoid,vAvoidB,regAvoid] = deal([]);
    count = 0;
    for j = 1:length(reg)
        if ~(iregion == j)
            count = count+1;
            vAvoid{count} = reg{j}.v;
            vAvoidB{count} = reg{j}.vB;
            regAvoid{count} = reg{j};
        end
    end
    
    % Get outgoing reach tubes from initial state
    ellBnd1 = [];  ellBnd11 = [];
    count = 0;
    for indx = 1:length(trans)
        if iregion == trans{indx}(1),
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
    
    Isect1{iregion} = [];
    
    if iregion == 2
        maxFunnels = 80;
    else
        maxFunnels = 40;
    end
    for i = 1:maxFunnels  % TODO: replace with while Vol(S\cup Ri) < Vol(Ri) and j < Nlimit
        funFail = true;
        while funFail
            for trial1 = 1:maxTrials1
                try
                    [finalPt,hindx1] = getRandom1(vReg{iregion},vAvoid,ellBnd11,'discrep','inside',hindx1);  % Get a final point inside the existing outgoing funnels
                    if isempty(finalPt)
                        error('getRandom1: Cannot obtain a point.')
                    end
                    goalXY = finalPt(1:2);
                    for thetaTrial = 1:100
                        thetaTrial
                        [randPt,hindx] = getRandom1(vReg{iregion},vAvoid,[],'discrep','outside',hindx);  % Get initial point
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
                            [isectTh] = checkIntersection4(ellBnd11,Xk(end,:));
                            if ~isectTh && length(t) > 1 && length(t) < 30, break, end
                        end
                    end
                    isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                catch error
%                     rethrow(error)
                    isect = true;
                end
                if ~isect && ~isectTh && length(t) > 1 && length(t) < 30, break, end
            end
            
            [funnelI,rhoMin,rho_d,info] = computeFunnel(t,Xk,Uk,regBnd,reg{iregion},regAvoid);
            if sum(rho_d < 1e9) > 3 && sum(rho_d == 1e9) < 5
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
        
        funnelIn{i,iregion,k} = funnelI;
        rhoMinArrayC{i,iregion,k} = funnelI.rho;
        ellFunnelC{i,iregion,k} = ellArrayCurr;
        trajFunnelC{i,iregion,k} = funnelI.x;
        
        %NB: in actuality, want the intersection of the funnels and the
        %initial set
        %             ellInit{itrans} = [ellInit{itrans}; ellFunnel{i,itrans}(1)];
        %             ellFinal{itrans} = [ellFinal{itrans}; ellFunnel{i,itrans}(end)];
        %         else
        %             i = i-1;
        %         end
        % Success for this transition if Init is covered (TODO: coverage metric)
        %         if isinsideunion(ellInit{itrans},vReg{trans{itrans}(1)})
        %             break
        %         end
        %         if i > 1
        %             if isinsideunion(ellArrayCurr(1),ellInit{itrans})  % if funnel is completely inside the reach tube, it will not be added
        %                 ellFunnel{i,itrans} = [];
        %                 i = i-1;
        %                 count = count+1;
        %                 if count > Nterm-1
        %                     break
        %                 end
        %             end
        %         end
        
        % Success for this transition if Init is covered
        icover = [];
        %             isect = [];
        for itrans = 1:length(trans)
            if trans{itrans}(1) == iregion
                itransReg = trans{itrans}(1);
                %                 for mm = 1:size(ellFunnel,1)
                %                     isect(mm,:) = isinternal(ellFunnel{mm,itrans,1},qCover{itransReg}');
                %                 end
                icover = [icover; any(Isect{itrans,k},1)];
            end
        end
        Isect1{iregion}(i,:) = isinternal(ellFunnelC{i,iregion,k},qCover{itransReg}');  % note: any itransReg is fine; all should have the same points
        
        insideIntersect = all(icover,1);
        insideCentralRT = any(Isect1{iregion},1);
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
% end

%%
save threeRegions40

%%
% figure(10)
% clf
% hold on
% axis equal
% plot(pReg{2},'g')
% plot(pReg{1},'c')
% plot(pReg{3},'b')
% 
% funnelIn0 = funnelIn;
% ellFunnelC0 = ellFunnelC;
% trajFunnelC0 = trajFunnelC;
% clear ellFunnelC trajFunnelC
% for j = 1:size(funnelIn0,2)
%     for i = 1:size(funnelIn0,1)
%         if ~isempty(funnelIn0{i,j,k})
% %             elim(i,j) = max(maxeig(ellFunnelC0{i,j,1})) < 0.05;
%             elim(i,j) = max(funnelIn{i,j,k}.rho) < 1;
%         end
%     end
% end
% for j = 1:size(funnelIn0,2)    
%     count = 0;
%     for i = 1:size(funnelIn0,1)
%         for k1 = 1:k-1
%             ellFunnelC{i,j,k1} = ellFunnelC0{i,j,k1};
%             trajFunnelC{i,j,k1} = trajFunnelC0{i,j,k1};
%         end
%         if ~elim(i,j),
%             count = count+1;
%             ellFunnelC{count,j,k} = ellFunnelC0{i,j,k};
%             trajFunnelC{count,j,k} = trajFunnelC0{i,j,k};
%         end
%     end
% end
% options.fill = 1;
% options.color = [1 0 0];
% for j = 1:size(funnelIn,2)
%     for i = 1:size(funnelIn,1)
%         if ~isempty(ellFunnelC{i,j,k}), plot(projection(ellFunnelC{i,j,k},[1 0;0 1;0 0]),options); end
%     end
% end
% for j = 1:size(funnelIn,2)
%     for i = 1:size(funnelIn,1)
%         if ~isempty(ellFunnelC{i,j,k}), plot(trajFunnelC{i,j,k}(:,1),trajFunnelC{i,j,k}(:,2),'k','LineWidth',2); end
%     end
% end
% 
