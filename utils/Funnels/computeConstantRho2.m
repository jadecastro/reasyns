function [rhoMin,rho_d,Info] = computeConstantRho2(f0,ts,Ac,K,xMu,uk,regBnd,X,Eprior,P1,reg1,Xbnd1,Xin1,P2,reg2,Xbnd2,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2,H,n,isCyclic,succRegChk)
% This is a routine to get max rho assuming it is constant throughout the
% trajectory.  

if sum(isCyclic) > 1, error('number of cyclic dimensions cannot exceed 1'), end
if length(isCyclic) ~= n, error('number of entries in isCylic must be equal to n!'), end
if length(ts) ~= size(xMu,1), error('length of t must match the number the rows in x.'), end
if length(ts) ~= size(uk,1), error('length of t must match the number the rows in u.'), end

global t x debugFlg

t = msspoly('t');

rho_d = [];
rho_runningMin = inf;
sclr = 1;%2000;
N1 = 100;
scalar = 1000;

trialMat = repmat([0 2*pi -2*pi],n,1).*repmat(isCyclic,1,3);

% Compute rho based on distance to polytopes
for i = 1:size(xMu,1)
    i
    E = shape(Eprior(i),1/sclr);
    [xbar,Q] = double(E);
    V(i) = (x - xbar)'*inv(Q)*(x - xbar);
    Ki = K{i};
    Aci = Ac{i};
    ubar = uk(i,:)';
    
    xdot = f0(t,x,ubar - Ki*(x - xbar));
    % xdot = Aci*(x - xbar);  % <--- for testing
    % xdot = f0(t,x,[2;ubar(2)] - Ki*(x - xbar));  % attempt at treating limits
    
    if isempty(Xbnd1) && isempty(reg2) && isempty(Xbnd2)  % % we have to construct an inward-facing funnel   
        
        % [L1,Lui,Lue] = computeL(V(i),xdot,regBnd,reg1);
        [rho_sos,Info{i}] = computeRho(V(i),xdot,X);
        Info{i}
        
        if ~isempty(ellBnd1) && i == size(xMu,1)  % be sure to check if ellBnd1 fully contains the funnel terminus
            
            reg = reg1;
            %             Xbnd = Xbnd1;
            %             Xin = Xin1;
            ellBnd = ellBnd1;
            ellIn = ellIn1;
            
            for j = 1:length(reg)
                vReg{j} = reg{j}.vB;
            end
                
            for nn = 1:size(ellBnd,2)  % For all outgoing transitions
                for m = 1:size(ellBnd,1) % for all transition funnels
                    %                 d = distance(ellBnd{m,nn},ball);
                    clear d
%                     XbndR{m,nn} = [];
                    ellBndR{m,nn} = [];
                    if ~isempty(ellBnd{m,nn})
                        for j = 1:length(ellBnd{m,nn})
                            if any(isCyclic)
                                d(j) = any(isinternal_quick(ellBnd{m,nn}(j),repmat(xbar,1,3)+trialMat,'u'));
                            else
                                d(j) = isinternal_quick(ellBnd{m,nn}(j),xbar,'u');
                            end
                        end
                        idxSet = find(d);
                        for p = idxSet
%                             XbndR{m,nn} = [XbndR{{m,nn}; Xbnd{m,nn}(p)];
                            ellBndR{m,nn} = [ellBndR{m,nn}; ellBnd{m,nn}(p)];
                        end
                    end
                end
            end
            
            [xbar,Q] = double(E);
            E = ellipsoid(xbar,Q*rho_runningMin);
            [Ei,Ep] = deal(E);
            [xbar,Q] = double(E);
            
            % starting with the trial rho, iteratively check containment
            % while decreasing rho
            for k = 1:10000
                [xbar,Q] = double(E);
                qTest = ellipsoidrand(xbar,Q,N1);
                clear badIndx isect
                [isect,badIndx] = checkIntersection5(regBnd{1}.v,vReg,ellBndR,[],qTest,H,n,isCyclic);
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
                        isectBad = checkIntersection5(regBnd{1}.v,vReg,ellBnd,[],qTest(bad1,:),H,n,isCyclic);
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
        
        rho = min(rho_sos,rho_ell);
        rho_d = [rho_d; rho];
        rho_runningMin = min(rho_d);
        
    elseif ~isempty(Xbnd1) && ~isempty(reg2) && ~isempty(Xbnd2)  % we have to construct a transition funnel
        
        tmpV = regBnd{1}.v;
        for j = 1:length(reg1)
            tmpV1{j} = reg1{j}.vB;
        end
        for j = 1:length(reg2)
            tmpV2{j} = reg2{j}.vB;
        end
        [isect1] = checkIntersection(tmpV1,tmpV,[],xbar',H);
        [isect2] = checkIntersection(tmpV2,tmpV,[],xbar',H);
        % TODO: fix intersections!
        
        clear reg Xbnd Xin ellBnd ellIn vReg XbndR XinR ellBndR ellInR
        if ~isect1 && isect2  % in R1
            reg{1} = reg1;
%             Xbnd{1} = Xbnd1;
%             Xin{1} = Xin1;
            ellBnd{1} = ellBnd1;
            ellIn{1} = ellIn1;
        elseif isect1 && ~isect2  % in R2
            reg{1} = reg2;
%             Xbnd{1} = Xbnd2;
%             Xin{1} = Xin2;
            ellBnd{1} = ellBnd2;
            ellIn{1} = ellIn2;
        else  % point is in the neutral zone
            reg{1} = reg1;  reg{2} = reg2;
%             Xbnd{1} = Xbnd1;  Xbnd{2} = Xbnd2;
%             Xin{1} = Xin1;  Xin{2} = Xin2;
            ellBnd{1} = ellBnd1;  ellBnd{2} = ellBnd2;
            ellIn{1} = ellIn1;  ellIn{2} = ellIn2;
        end
        
        if ~isect2
            [rho_sos,Info{i}] = computeRho(V(i),xdot,X);

        
        for l = 1:length(reg)
            for j = 1:length(reg{l})
                vReg{l}{j} = reg{l}{j}.vB;
            end
            
            for nn = 1:size(ellBnd{l},2)  % For all outgoing transitions
                for m = 1:size(ellBnd{l},1) % for all transition funnels
                    %                 d = distance(ellBnd{m,nn},ball);
                    clear d
%                     XbndR{l}{m,nn} = [];
                    ellBndR{l}{m,nn} = [];
                    if ~isempty(ellBnd{l}{m,nn})
                        for j = 1:length(ellBnd{l}{m,nn})
                            if any(isCyclic)
                                d(j) = any(isinternal_quick(ellBnd{l}{m,nn}(j),repmat(xbar,1,3)+trialMat,'u'));
                            else
                                d(j) = isinternal_quick(ellBnd{l}{m,nn}(j),xbar,'u');
                            end
                        end
                        idxSet = find(d);
                        for p = idxSet
%                             XbndR{l}{m,nn} = [XbndR{l}{m,nn}; Xbnd{l}{m,nn}(p)];
                            ellBndR{l}{m,nn} = [ellBndR{l}{m,nn}; ellBnd{l}{m,nn}(p)];
                        end
                    end
                end
            end
            
            for m = 1:length(ellIn{l}) % for all transition funnels
                %             d = distance(ellIn{m},ball);
                clear d
%                 XinR{l}{m} = [];
                ellInR{l}{m} = [];
                if ~isempty(ellIn{l}{m})
                    for j = 1:length(ellIn{l}{m})
                        if any(isCyclic)
                            d(j) = any(isinternal_quick(ellIn{l}{m}(j),repmat(xbar,1,3)+trialMat,'u'));
                        else
                            d(j) = isinternal_quick(ellIn{l}{m}(j),xbar,'u');
                        end
                    end
                    idxSet = find(d);
                    for p = idxSet
%                         XinR{l}{m} = [XinR{l}{m}; Xin{l}{m}(p)];
                        ellInR{l}{m} = [ellInR{l}{m}; ellIn{l}{m}(p)];
                    end
                end
            end
            ellBndR{l}
            ellInR{l}
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

            for l = 1:length(reg)
                for nn = 1:size(ellBnd{l},2)  % For all outgoing transitions
                    for m = 1:size(ellBnd{l},1) % for all transition funnels
                        if ~isempty(ellBndR{l}{m,nn})
                            plot(projection(ellBndR{l}{m,nn},[H; zeros(n-length(H),length(H))]),'r')
                        end
                    end
                end
                
                for m = 1:length(ellIn{l}) % for all transition funnels
                    if ~isempty(ellInR{l}{m})
                        plot(projection(ellInR{l}{m},[H; zeros(n-length(H),length(H))]),'r')
                    end
                end
            end
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
                [isect(l),badIndx{l}] = checkIntersection5(regBnd{1}.v,vReg{l},ellBndR{l},ellInR{l},qTest,H,n,isCyclic);
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
                            isectBad(l) = checkIntersection5(regBnd{1}.v,vReg{l},ellBnd{l},ellIn{l},qTest(bad1(j),:),H,n,isCyclic);
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
        [rho_sos,Info{i}] = computeRhoFeas(V(i),xdot);
        Info{i}
        if rho_sos == 1  % sedumi usually returns a 1 if unbounded
            rho_sos = 1e9;
        end
        
        rho = min(rho_sos,rho_ell);
        rho_d = [rho_d; rho];
        rho_runningMin = min(rho_d);
    end
    
end

% if debugFlg
%     plotV
% end

rho_d = rho_d/sclr;
rho_d
rhoMin = min(rho_d)

if sclr ~= 1
    indx = find(rho_d == rhoMin, 1 );
    E = shape(Eprior(indx),1/sclr);
    rhoMin = rhoMin/sclr
end

