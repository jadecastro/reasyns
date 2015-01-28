% Now, try to compute some inward-directed funnels
count = 0;
for iregion = sort(unique(cell2mat(trans)))
    hindx = 1;  % reset Halton seed
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
    ellBnd1 = [];
    count = 0;
    for indx = 1:length(trans)
        if iregion == trans{indx}(1),
            count = count+1;
            for ii = 1:size(funnel,1)
                ellBnd1{ii,count} = ellFunnel{ii,indx,1};
            end
        end
    end
    
    Isect1{iregion} = [];
    
    for i = 1:20  % TODO: replace with while Vol(S\cup Ri) < Vol(Ri) and j < Nlimit
        
        for trial1 = 1:maxTrials1
            [randPt,hindx] = getRandom1(vReg{iregion},vAvoid,ellBnd1,'discrep','outside',hindx);  % Get initial point
            %         [randPt,hindx] = getRandom1(vReg{iregion},vAvoid,[],'discrep',[],hindx);  % Get initial point
            hindx
            initXYTh = randPt;
            initXY = randPt(1:2);
            
            [finalPt,hindx1] = getRandom1(vReg{iregion},vAvoid,ellBnd1,'discrep','inside',hindx1);  % Get a final point inside the existing outgoing funnels
            goalXY = finalPt(1:2);
            
            [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);  %TODO: Fix check on bounds
            
            [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);  % TODO: Need to make sure
            isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
            if ~isect && length(t) > 1, break, end
            for trial2 = 1:maxTrials2
                disp('Trajectory incompatible with constraints; recomputing...')
                
                [finalPt,hindx1] = getRandom1(vReg{iregion},vAvoid,ellBnd1,'discrep','inside',hindx1);  % Give up on trying to have a nice compact set; randomly select the final point
                goalXY = finalPt(1:2);
                
                [pathTmp,rrtIndx] = buildRRT(vAvoid,vBnd{:},initXY,goalXY,stepSize,n,radius);
                
                [t,Xk,Uk] = genNominalTrajectory(initXYTh,pathTmp.q);
                isect = checkIntersection(vAvoid,vBnd{:},[],Xk);
                if ~isect && length(t) > 1, break, end
            end
        end
        
        figure(3)
        plot(Xk(:,1),Xk(:,2),'k','LineWidth',2)

        [funnelI,rhoMin,rho_d] = computeFunnel(t,Xk,Uk,regBnd,reg{iregion},regAvoid);
        
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
        icover1 = insideCentralRT.*(~insideIntersect);  % Number of points covered by the inward-directed RT
%         icover1 = [all(icover,1); any(isect,1)];
        NcoverAct = sum(icover1);
        
        disp(['Iteration #',num2str(i)])
        disp([' Percent of region covered: ',num2str(NcoverAct/Ncover1),' , ',num2str(NcoverAct),' out of ',num2str(Ncover1)]);
        if NcoverAct > coverPct*Ncover1
            break
        end
    end
end
