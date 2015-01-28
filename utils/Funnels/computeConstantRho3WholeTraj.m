function [rhoMin,rho_d,Info] = computeConstantRho3WholeTraj(f0,ts,Ac,K,xMu,uk,regBnd,X,Eprior,P1,reg1,Xbnd1,Xin1,P2,reg2,Xbnd2,Xin2,ellBndInv1,ellBndInv2,ellInInv1,ellInInv2,H,n,isCyclic,succRegChk,bloatFlag)
% This is a routine to get max rho assuming it is constant throughout the
% trajectory.  

if sum(isCyclic) > 1, error('number of cyclic dimensions cannot exceed 1'), end
if length(isCyclic) ~= n, error('number of entries in isCylic must be equal to n!'), end
if length(ts) ~= size(xMu,1), error('length of t must match the number the rows in x.'), end
if length(ts) ~= size(uk,1), error('length of t must match the number the rows in u.'), end

global t x debugFlg bloatFlag

t = msspoly('t');

rho_d = [];
rho_runningMin = inf;
sclr = 1;%2000;
N1 = 100;
scalar = 1000;

trialMat = repmat([0 2*pi -2*pi],n,1).*repmat(isCyclic,1,3);

for i = 1:size(xMu,1)
    E = shape(Eprior(i),1/sclr);
%     xbar = EpriorInv(i).x;
%     Qinv = EpriorInv(i).P*sclr;
    [xbar,Q] = double(E);
    V(i) = (x - xbar)'*inv(Q)*(x - xbar);
    Ki = K{i};
    Aci = Ac{i};
    ubar = uk(i,:)';
    
    xdot(:,i) = f0(t,x,ubar - Ki*(x - xbar)) - f0(t,xbar,ubar);
    % xdot = Aci*(x - xbar);  % <--- for testing
    % xdot = f0(t,x,[2;ubar(2)] - Ki*(x - xbar));  % attempt at treating limits
    
    % [L1,Lui,Lue] = computeL(V(i),xdot,regBnd,reg1);
%     [rho_sos,Info{i},posDefMultiplier] = computeRhoFeas(V(i),xdot);
end

tau_sos = ts(2:end) - ts(1:end-1);
% [rho_sos,Info,posDefMultiplier] = computeRhoFeasBilinearAlternationWhole(V,xdot,tau_sos,xMu);
[rho_sos,Info,posDefMultiplier] = computeRhoBilinearAlternationWhole(V,xdot,tau_sos,X);
rho_ratio = ones(size(rho_sos));

% Compute rho based on distance to polytopes
for i = 1:size(xMu,1)
    i;
    E = shape(Eprior(i),1/sclr);
%     xbar = EpriorInv(i).x;
%     Qinv = EpriorInv(i).P*sclr;
    [xbar,Q] = double(E);
    V(i) = (x - xbar)'*inv(Q)*(x - xbar);
    Ki = K{i};
    Aci = Ac{i};
    ubar = uk(i,:)';
    
    xdot = f0(t,x,ubar - Ki*(x - xbar)) - f0(t,xbar,ubar);
    % xdot = Aci*(x - xbar);  % <--- for testing
    % xdot = f0(t,x,[2;ubar(2)] - Ki*(x - xbar));  % attempt at treating limits
    
    if isempty(Xbnd1) && isempty(reg2) && isempty(Xbnd2)  % % we have to construct an inward-facing funnel   
        
        % [L1,Lui,Lue] = computeL(V(i),xdot,regBnd,reg1);
%         [rho_sos,Info{i}] = computeRho(V(i),xdot,X);
%         Info{i};
        
        if ~isempty(ellBndInv1) && i == size(xMu,1)  % be sure to check if ellBnd1 fully contains the funnel terminus
            
            reg = reg1;
            %             Xbnd = Xbnd1;
            %             Xin = Xin1;
            ellBndInv = ellBndInv1;
            ellInInv = ellInInv1;
            
            % sanity check to see if final point is indeed inside the required set
            for nn = 1:size(ellBndInv,2)  % For all outgoing transitions
                isLastPointInside(nn) = false;
                for m = 1:size(ellBndInv,1) % for all transition funnels
                    if ~isempty(ellBndInv{m,nn})
                        for j = 1:length(ellBndInv{m,nn})
                            if any(isCyclic) && any(isinternal_quickInv(ellBndInv{m,nn}(j),repmat(xbar,1,3)+trialMat,'u'))
                                isLastPointInside(nn) = true;
                                break
                            elseif ~any(isCyclic) && isinternal_quickInv(ellBndInv{m,nn}(j),xbar,'u');
                                isLastPointInside(nn) = true;
                                break
                            end
                        end
                    end
                    if isLastPointInside(nn), break; end
                end
            end
            disp(['   last point inside?  ',num2str(all(isLastPointInside))])
            if ~all(isLastPointInside)
                rho = 1e-9;
                rho_ratio(i) = rho/rho_sos(i);
%                 rho_d = [rho_d; rho];
%                 rho_runningMin = min(rho_d);
                break
            end
            
            for j = 1:length(reg)
                if bloatFlag
                    vReg{j} = reg{j}.vB;
                else
                    vReg{j} = reg{j}.v;
                end
            end
                
            for nn = 1:size(ellBndInv,2)  % For all outgoing transitions
                for m = 1:size(ellBndInv,1) % for all transition funnels
                    %                 d = distance(ellBnd{m,nn},ball);
                    clear d
%                     XbndR{m,nn} = [];
                    ellBndInvR{m,nn} = [];
                    if ~isempty(ellBndInv{m,nn})
                        for j = 1:length(ellBndInv{m,nn})
                            if any(isCyclic)
                                d(j) = any(isinternal_quickInv(ellBndInv{m,nn}(j),repmat(xbar,1,3)+trialMat,'u'));
                            else
                                d(j) = isinternal_quickInv(ellBndInv{m,nn}(j),xbar,'u');
                            end
                        end
                        idxSet = find(d);
                        for p = idxSet
%                             XbndR{m,nn} = [XbndR{{m,nn}; Xbnd{m,nn}(p)];
                            ellBndInvR{m,nn} = [ellBndInvR{m,nn}; ellBndInv{m,nn}(p)];
                        end
                    end
                end
            end
            
            [xbar,Q] = double(E);
%             E = ellipsoid(xbar,Q*rho_runningMin);
            E = ellipsoid(xbar,Q*rho_sos(i));
            [Ei,Ep] = deal(E);
            [xbar,Q] = double(E);
            
            % starting with the trial rho, iteratively check containment
            % while decreasing rho
            for k = 1:500
                [xbar,Q] = double(E);
                qTest = ellipsoidrand(xbar,Q,N1);
                clear badIndx isect
                [isect,badIndx] = checkIntersection5(regBnd{1}.v,vReg,ellBndInvR,[],qTest,H,n,isCyclic);
                if ~isect
                    k
                    break,
                end
                if isect
                    bad1 = badIndx;
                    if ~isempty(bad1)
                        k
                        clear isectBad
                        isectBad = false;
                        isectBad = checkIntersection5(regBnd{1}.v,vReg,ellBndInv,[],qTest(bad1,:),H,n,isCyclic);
                        % [isectBad,~,ellInR,ellBndR] = checkIntersection5(regBnd{1}.v,vReg,ellBnd,[],qTest(bad1(j),:),H,n,isCyclic);
                        if ~isectBad
                            k
                            break,
                        end
                    end
                end
                Ep = E;
                E = shape(E,0.9);  % decrease by 10% and try again..
            end
            if debugFlg
                plot(projection(E,[H; zeros(n-length(H),length(H))]),'g');
                drawnow
                %             keyboard
            end
            tmp = double(E)./double(Ei);
%             rho_ell = tmp(1,1)*rho_runningMin;
            rho_ell = tmp(1,1)*rho_sos(i);
        else
            rho_ell = inf;
        end
        
        rho = min(rho_sos(i),rho_ell);
        rho_ratio(i) = rho/rho_sos(i);
%         rho_d = [rho_d; rho];
%         rho_runningMin = min(rho_d);
        
    elseif ~isempty(Xbnd1) && ~isempty(reg2) && ~isempty(Xbnd2)  % we have to construct a transition funnel
        
        tmpV = regBnd{1}.v;
        for j = 1:length(reg1)
            if bloatFlag
                tmpV1{j} = reg1{j}.vB;
            else
                tmpV1{j} = reg1{j}.v;
            end
        end
        for j = 1:length(reg2)
            if bloatFlag
                tmpV2{j} = reg2{j}.vB;
            else
                tmpV2{j} = reg2{j}.v;
            end
        end
        [isect1] = checkIntersection(tmpV1,tmpV,[],xbar',H);
        [isect2] = checkIntersection(tmpV2,tmpV,[],xbar',H);
        % TODO: fix intersections!
        
        clear reg Xbnd Xin ellBnd ellIn vReg XbndR XinR ellBndR ellInR
        
        if ~isect2  % in R2; we can treat it like an inward funnel
            % [L1,Lui,Lue] = computeL(V(i),xdot,regBnd,reg1);
%             [rho_sos,Info{i}] = computeRhoBloat(V(i),xdot,X);
%             Info{i}
            
            if ~isempty(ellBndInv2) && i == size(xMu,1)  % be sure to check if ellBnd1 fully contains the funnel terminus
                
                reg = reg2;
                ellBndInv = ellBndInv2;
                ellInInv = ellInInv2;
                
                for j = 1:length(reg)
                    if bloatFlag
                        vReg{j} = reg{j}.vB;
                    else
                        vReg{j} = reg{j}.v;
                    end
                end
                
                for nn = 1:size(ellBndInv,2)  % For all outgoing transitions
                    for m = 1:size(ellBndInv,1) % for all transition funnels
                        %                 d = distance(ellBnd{m,nn},ball);
                        clear d
                        %                     XbndR{m,nn} = [];
                        ellBndInvR{m,nn} = [];
                        if ~isempty(ellBndInv{m,nn})
                            for j = 1:length(ellBndInv{m,nn})
                                if any(isCyclic)
                                    d(j) = any(isinternal_quickInv(ellBndInv{m,nn}(j),repmat(xbar,1,3)+trialMat,'u'));
                                else
                                    d(j) = isinternal_quickInv(ellBndInv{m,nn}(j),xbar,'u');
                                end
                            end
                            idxSet = find(d);
                            for p = idxSet
                                %                             XbndR{m,nn} = [XbndR{{m,nn}; Xbnd{m,nn}(p)];
                                ellBndInvR{m,nn} = [ellBndInvR{m,nn}; ellBndInv{m,nn}(p)];
                            end
                        end
                    end
                end
                
                for m = 1:length(ellInInv) % for all transition funnels
                    %             d = distance(ellIn{m},ball);
                    clear d
                    %                 XinR{l}{m} = [];
                    ellInInvR{m} = [];
                    if ~isempty(ellInInv{m})
                        for j = 1:length(ellInInv{m})
                            if any(isCyclic)
                                d(j) = any(isinternal_quickInv(ellInInv{m}(j),repmat(xbar,1,3)+trialMat,'u'));
                            else
                                d(j) = isinternal_quickInv(ellInInv{m}(j),xbar,'u');
                            end
                        end
                        idxSet = find(d);
                        for p = idxSet
                            %                         XinR{l}{m} = [XinR{l}{m}; Xin{l}{m}(p)];
                            ellInInvR{m} = [ellInInvR{m}; ellInInv{m}(p)];
                        end
                    end
                end
                
                [xbar,Q] = double(E);
                E = ellipsoid(xbar,Q*rho_runningMin);
                [Ei,Ep] = deal(E);
                [xbar,Q] = double(E);
                
                % starting with the trial rho, iteratively check containment
                % while decreasing rho
                for k = 1:500
                    [xbar,Q] = double(E);
                    qTest = ellipsoidrand(xbar,Q,N1);
                    clear badIndx isect
                    [isect,badIndx] = checkIntersection5(regBnd{1}.v,vReg,ellBndInvR,ellInInvR,qTest,H,n,isCyclic);
                    if ~isect
                        k
                        break,
                    end
                    if isect
                        bad1 = badIndx;
                        if ~isempty(bad1)
                            k
                            clear isectBad
                            isectBad = false;
                            isectBad = checkIntersection5(regBnd{1}.v,vReg,ellBndInv,ellInInv,qTest(bad1,:),H,n,isCyclic);
                            % [isectBad,~,ellInR,ellBndR] = checkIntersection5(regBnd{1}.v,vReg,ellBnd,[],qTest(bad1(j),:),H,n,isCyclic);
                            if ~isectBad
                                k
                                break,
                            end
                        end
                    end
                    Ep = E;
                    E = shape(E,0.9);  % decrease by 10% and try again..
                end
                if debugFlg
                    plot(projection(E,[H; zeros(n-length(H),length(H))]),'g');
                    drawnow
                    %             keyboard
                end
                tmp = double(E)./double(Ei);
                rho_ell = tmp(1,1)*rho_runningMin;
            else
                rho_ell = inf;
            end
            
            rho = min(rho_sos(i),rho_ell);
            rho_ratio(i) = rho/rho_sos(i);
%             rho_d = [rho_d; rho];
%             rho_runningMin = min(rho_d);
            
        else % still inside R1; enforce reactive composition if given transition funnels; else do sos (ellBndInv1 is empty)
            
            clear reg Xbnd Xin ellBnd ellIn vReg XbndR XinR ellBndR ellInR
            reg{1} = reg1;
%             Xbnd{1} = Xbnd1;
%             Xin{1} = Xin1;
            ellBndInv{1} = ellBndInv1;
            ellInInv{1} = ellInInv1;
            
            if ~isempty(ellBndInv{1})
                for l = 1:length(reg)
                    for j = 1:length(reg{l})
                        if bloatFlag
                            vReg{l}{j} = reg{l}{j}.vB;
                        else
                            vReg{l}{j} = reg{l}{j}.v;
                        end
                    end
                    
                    for nn = 1:size(ellBndInv{l},2)  % For all outgoing transitions
                        for m = 1:size(ellBndInv{l},1) % for all transition funnels
                            %                 d = distance(ellBnd{m,nn},ball);
                            clear d
                            %                     XbndR{l}{m,nn} = [];
                            ellBndInvR{l}{m,nn} = [];
                            if ~isempty(ellBndInv{l}{m,nn})
                                for j = 1:length(ellBndInv{l}{m,nn})
                                    if any(isCyclic)
                                        d(j) = any(isinternal_quickInv(ellBndInv{l}{m,nn}(j),repmat(xbar,1,3)+trialMat,'u'));
                                    else
                                        d(j) = isinternal_quickInv(ellBndInv{l}{m,nn}(j),xbar,'u');
                                    end
                                end
                                idxSet = find(d);
                                for p = idxSet
                                    %                             XbndR{l}{m,nn} = [XbndR{l}{m,nn}; Xbnd{l}{m,nn}(p)];
                                    ellBndInvR{l}{m,nn} = [ellBndInvR{l}{m,nn}; ellBndInv{l}{m,nn}(p)];
                                end
                            end
                        end
                    end
                    
                    for m = 1:length(ellInInv{l}) % for all transition funnels
                        %             d = distance(ellIn{m},ball);
                        clear d
                        %                 XinR{l}{m} = [];
                        ellInInvR{l}{m} = [];
                        if ~isempty(ellInInv{l}{m})
                            for j = 1:length(ellInInv{l}{m})
                                if any(isCyclic)
                                    d(j) = any(isinternal_quickInv(ellInInv{l}{m}(j),repmat(xbar,1,3)+trialMat,'u'));
                                else
                                    d(j) = isinternal_quickInv(ellInInv{l}{m}(j),xbar,'u');
                                end
                            end
                            idxSet = find(d);
                            for p = idxSet
                                %                         XinR{l}{m} = [XinR{l}{m}; Xin{l}{m}(p)];
                                ellInInvR{l}{m} = [ellInInvR{l}{m}; ellInInv{l}{m}(p)];
                            end
                        end
                    end
                    %ellBndR{l}
                    %ellInR{l}
                end            
                % [L1,Lui,Lue] = computeL(V(i),xdot,regBnd,reg1);
                %         [rho,info] = computeRho(V(i),xdot,X,XbndR,XinR);
                %         info
                %         Info{i} = info;
                
                if debugFlg
                    figure(1)
                    clf
                    hold on
                    %             for nn = 1:size(Xbnd,2)  % For all outgoing transitions
                    %                 for m = 1:size(Xbnd,1) % for all transition funnels
                    %                     if ~isempty(ellBnd{m,nn})
                    %                         plot(projection(ellBnd{m,nn},[H; zeros(n-length(H),length(H))]),'c')
                    %                     end
                    %                 end
                    %             end
                    %
                    %             for m = 1:length(Xin) % for all transition funnels
                    %                 if ~isempty(ellIn{m})
                    %                     plot(projection(ellIn{m},[H; zeros(n-length(H),length(H))]),'c')
                    %                 end
                    %             end
                    
%                     for l = 1:length(reg)
%                         for nn = 1:size(ellBndInv{l},2)  % For all outgoing transitions
%                             for m = 1:size(ellBndInv{l},1) % for all transition funnels
%                                 if ~isempty(ellBndInvR{l}{m,nn})
%                                     plot(projection(ellBndInvR{l}{m,nn},[H; zeros(n-length(H),length(H))]),'r')
%                                 end
%                             end
%                         end
%                         
%                         for m = 1:length(ellInInv{l}) % for all transition funnels
%                             if ~isempty(ellInInvR{l}{m})
%                                 plot(projection(ellInInvR{l}{m},[H; zeros(n-length(H),length(H))]),'r')
%                             end
%                         end
%                     end
                end
                
                E = shape(Eprior(i),1/scalar);
                [Ei,Ep] = deal(E);
                [xbar,Q] = double(E);
                
                % starting with a tiny rho, iteratively check containment while
                % increating rho
                for k = 1:10000
                    tmp = double(Ep)./double(Ei);
                    rho_test = tmp(1,1)/scalar^2;
                    if rho_test >= rho_runningMin  % we're wasting our time computing more iterations, so break the loop
                        k
                        rho_test
                        break,
                    end
                    [xbar,Q] = double(E);
                    qTest = ellipsoidrand(xbar,Q,N1);
                    clear badIndx isect
                    for l = 1:length(reg)
                        [isect(l),badIndx{l}] = checkIntersection5(regBnd{1}.v,vReg{l},ellBndInvR{l},ellInInvR{l},qTest,H,n,isCyclic);
                    end
                    if all(isect)
                        bad1 = badIndx{1};
                        for l = 1:length(reg)-1
                            bad1 = intersect(badIndx{l+1},bad1);
                        end
                        if ~isempty(bad1)
                            k
                            clear isectBad
                            isectBad = false;
                            for j = 1:length(bad1)
                                for l = 1:length(reg)
                                    isectBad(l) = checkIntersection5(regBnd{1}.v,vReg{l},ellBndInv{l},ellInInv{l},qTest(bad1(j),:),H,n,isCyclic);
                                    % [isectBad(l),~,ellInR{l},ellBndR{l}] = checkIntersection5(regBnd{1}.v,vReg{l},ellBnd{l},ellIn{l},qTest(bad1(j),:),H,n,isCyclic);
                                end
                                if all(isectBad)
                                    break,
                                end
                            end
                            if all(isectBad)
                                k
                                break,
                            end
                        end
                    end
                    Ep = E;
                    E = shape(E,1.1);  % increase by 10% and try again..
                end
                if debugFlg
                    plot(projection(Ep,[H; zeros(n-length(H),length(H))]),'g');
                    drawnow
                    %             keyboard
                end
                tmp = double(Ep)./double(Ei);
                rho_ell = tmp(1,1)/scalar^2;
                
                % perform feasibility check to determine maximal legitimate rho
%                 [rho_sos,Info{i}] = computeRhoFeas(V(i),xdot);
%                 Info{i}
%                 if rho_sos == 1  % sedumi usually returns a 1 if unbounded
%                     rho_sos = 1e9;
%                 end

            else
                rho_ell = inf;
                % [L1,Lui,Lue] = computeL(V(i),xdot,regBnd,reg1);
%                 [rho_sos,Info{i}] = computeRho(V(i),xdot,X);
%                 Info{i}
            end
            
            rho = min(rho_sos(i),rho_ell);
            rho_ratio(i) = rho/rho_sos(i);
%             rho_d = [rho_d; rho];
%             rho_runningMin = min(rho_d);
        end
    end
end

% if debugFlg
%     plotV
% end

% TODO: verify sos constraint satisfaction for this min(rho_ratio)!

rho_d = rho_sos/sclr;
rho_d = min(rho_ratio)*rho_d;
rho_d
rhoMin = min(rho_d)


