function [rhoMin,rho_d] = computeConstantRho(t,Ac,xMu,regBnd,Eprior,P1,reg1,ellBound1,ellC1,P2,reg2,ellBound2,ellC2)
% This is a routine to get max rho assuming it is constant throughout the
% trajectory.  Warning: rho must separately be checked for compatibility
% as a funnel (i.e. satisfies Lyapunov inequalities).

rho_d = [];
sclr = 2000;

% Compute rho based on distance to polytopes
for i = 1:length(xMu)
    i
    E = shape(Eprior(i),1/sclr);
    
    if isempty(ellBound1) && isempty(vInv2) && isempty(ellBound2)  % Consider vInv1 as the invariant for BOTH regions
%     % METHOD 1:
%     for j = 1:length(vInv)
%         [dist,status,ellPoint,polyPoint] = distanceEllPoly(E,P{j});  %<-- Fix this... something's not correct
%         if dist == 0, break, end
%         d1 = norm(xMu(i,1:2)' - ellPoint);
%         % d2 = norm(xMu(i,1:2)' - polyPoint);
%         d2 = norm((xMu(i,1:2)' - ellPoint) + dist);
%         rhoTmp = (d2/d1)^2;
%         rho_d = [rho_d; rhoTmp];
%     end
%     if dist == 0, break, end
    
    % METHOD 2:
        rhoTmp = getMinRhoPoly(E,vBound,vInv1,xMu(i,:),10);
        
        rho_d = [rho_d; rhoTmp];
    end
    
    if ~isempty(ellBound1) && ~isempty(vInv2) && ~isempty(ellBound2)  
        % Compute rho based on distance to bounding ellipsoids (NB: method only works for 3D state space)
        if size(xMu,2) > 3, error, end
        
        % get trial ellipsoid
        E1{i} = ellipsoid(xMu(i,:)',PinvTmp);
        
        % Create hyperplanes at each time step and resulting rank-2 ellipsoids
        xDot(i,:) = Ac(t(i))*xMu(i,:)';
        phi = atan2(xDot(i,2),xDot(i,1));
        psi = atan2(xDot(i,3),sqrt(xDot(i,2)^2+xDot(i,2)^2));
        c = xDot(i,:)*xMu(i,:)';
        H{i} = hyperplane(xDot(i,:)',c);
        ellHull{i} = hpintersection(E1{i},H{i});
        T{i} = [-sin(phi) -cos(phi)*sin(psi); cos(phi) -sin(phi)*sin(psi); 0 cos(psi)];
    
        % Project all ellipsoids onto hyperplane
        ellTmp = projection(ellHull{i},T{i});
%         [c,Q] = double(ellTmp);
        ellBndTmp1 = [];
        for nn = 1:size(ellBound1,2)
            count = 0;
            for mm = 1:size(ellBound1,1),
                tmp = [];
                for ii = 1:length(ellBound1{mm,nn})
                    if isinternal_quick(ellBound1{mm,nn}(ii),xMu(i,:)')
%                         i
%                         nn
%                         mm
%                         ii
                        tmp = [tmp; ellBound1{mm,nn}(ii)];
%                         tmp = hpintersection(ellBound1{mm,nn}(ii),H{i});
%                         if ~isempty(tmp), ellBndTmp1 = [ellBndTmp1; tmp]; end
%                         if ~isempty(tmp), ellBndTmp1 = [ellBndTmp1; projection(tmp,[-sin(phi) -cos(phi)*sin(psi); cos(phi) -sin(phi)*sin(psi); 0 cos(psi)])]; end
                    end
                end
                if ~isempty(tmp)
                    count = count+1;
                    ellBndTmp1{count,nn} = tmp;
                end
            end
        end
        ellBndTmp2 = [];
        for nn = 1:size(ellBound2,2)
            count = 0;
            for mm = 1:size(ellBound2,1),
                tmp = [];
                for ii = 1:length(ellBound2{mm,nn})
                    if isinternal_quick(ellBound2{mm,nn}(ii),xMu(i,:)')
%                         i
%                         nn
%                         mm
%                         ii
                        tmp = [tmp; ellBound2{mm,nn}(ii)];
%                         tmp = hpintersection(ellBound2{mm,nn}(ii),H{i});
%                         if ~isempty(tmp), ellBndTmp2 = [ellBndTmp2; tmp]; end
%                         if ~isempty(tmp), ellBndTmp2 = [ellBndTmp2; projection(tmp,[-sin(phi) -cos(phi)*sin(psi); cos(phi) -sin(phi)*sin(psi); 0 cos(psi)])]; end
                    end
                end
                if ~isempty(tmp)
                    count = count+1;
                    ellBndTmp2{count,nn} = tmp;
                end
            end
        end
        count = 0;
        ellCTmp1 = [];
        for mm = 1:size(ellC1,1),
            tmp = [];
            for ii = 1:length(ellC1{mm})
                if isinternal_quick(ellC1{mm}(ii),xMu(i,:)')
%                         i
%                         nn
%                         mm
%                         ii
                    tmp = [tmp; ellC1{mm}(ii)];
                end
            end
            if ~isempty(tmp)
                count = count+1;
                ellCTmp1{count} = tmp;
            end
        end
        count = 0;
        ellCTmp2 = [];
        for mm = 1:size(ellC2,1),
            tmp = [];
            for ii = 1:length(ellC2{mm})
                if isinternal_quick(ellC2{mm}(ii),xMu(i,:)')
%                         i
%                         nn
%                         mm
%                         ii
                    tmp = [tmp; ellC2{mm}(ii)];
                end
            end
            if ~isempty(tmp)
                count = count+1;
                ellCTmp2{count} = tmp;
            end
        end
        
        rhoTmp2 = getMinRhoEll(T{i},ellTmp,E,vBound,vInv1,vInv2,ellBndTmp1,ellBndTmp2,ellCTmp1,ellCTmp2,xMu(i,:),10,1);
%         rhoTmp2 = getMinRhoEll(T{i},ellTmp,E,vBound,vInv1,vInv2,ellBound1,ellBound2,ellC1,ellC2,xMu(i,:),10,1);
%         rhoTmp2 = getMinRhoEll(E,vBound,vInv1,vInv2,ellBndTmp1,ellBndTmp2,c');
        
        % rho_d = [rho_d; min(rhoTmp1,rhoTmp2)];
        rho_d = [rho_d; rhoTmp2];
        
    end
    
end
% if dist == 0,
%     rho_d = NaN;
%     rhoMin = 1;
%     return
% end

rho_d = rho_d/sclr;
rhoMin = min(rho_d)

% Refine estimate of limiting rho
indx = find(rho_d == rhoMin, 1 );
E = shape(Eprior(indx),1/sclr);

if isempty(ellBound1) && isempty(vInv2) && isempty(ellBound2) 
    rhoMin = getMinRhoPoly(E,vBound,vInv1,xMu(indx,:),30);
end
if ~isempty(ellBound1) && ~isempty(vInv2) && ~isempty(ellBound2)    
    ellTmp = projection(ellHull{indx},T{indx});
    rhoMin = getMinRhoEll(T{indx},ellTmp,E,vBound,vInv1,vInv2,ellBndTmp1,ellBndTmp2,ellCTmp1,ellCTmp2,xMu(indx,:),10,10);
end
rhoMin = rhoMin/sclr

end

function rho = getMinRhoPoly(Einit,vBound,vInv,q,N)

    Eprior = Einit;
    Epost = Eprior;
    isect = false;
    for k = 1:1000        
        [xE,yE] = Ellipse_Plot(inv(double(Epost)),q(1:2)',1,N);
        isect = double(any(any([xE' yE'] < repmat(min(vBound),length(xE),1))) || any(any([xE' yE'] > repmat(max(vBound),length(xE),1))));
        for j = 1:length(vInv)
            isect(j+1) = any(inpolygon(xE,yE,vInv{j}(:,1),vInv{j}(:,2)));
        end
        if any(isect), break, end
        Eprior = Epost;
        Epost = shape(Epost,1.1);  % heuristically increase by 10% and try again..
    end
    tmp = double(Eprior)./double(Einit);
    rho = tmp(1,1);
end

function rho = getMinRhoEll(T,Einit,Einit1,vBound,vObs1,vObs2,ellBnd1,ellBnd2,ellCTmp1,ellCTmp2,q,N1,N2)

    Eprior = Einit;
    Epost = Eprior;     % Einit and Epost are 2-D ellipses projected on the hyperplane
    Eprior1 = Einit1; 
    Epost1 = Eprior1;    % Einit1 and Epost1 are 2-D ellipses projected on the x-y axis
    isect = false;
    for k = 1:1000 
        tmp = inv(double(Epost1));
        [xE,yE] = Ellipse_Plot((tmp'+tmp)/2,q(1:2)',1,N1);
%         figure(3), plot(xE,yE,'g')
        isect = double(any(any([xE' yE'] < repmat(min(vBound),length(xE),1))) || any(any([xE' yE'] > repmat(max(vBound),length(xE),1))));
%         isectTmp1 = []; isectTmp2 = [];
%         for mm = 1:length(vObs1),
%             isectTmp1 = [isectTmp1 inpolygon(xE,yE,vObs1{mm}(:,1),vObs1{mm}(:,2))];
%         end
%         for mm = 1:length(vObs2),
%             isectTmp2 = [isectTmp2 inpolygon(xE,yE,vObs2{mm}(:,1),vObs2{mm}(:,2))];
%         end
%         if any([isect ~any(isectTmp1)&&~any(isectTmp2)]), 
%             figure(3), plot(xE,yE,'r')
%             break, 
%         end
        isect = [isect checkIntersection3(vBound,vObs1,vObs2,[],[],[],[],[xE;yE]')];
        if any(isect),
            figure(3), plot(qTest(1,:),qTest(2,:),'r')
            break, 
        end
        % TODO: use isinside before moving to the pointwise criteria?
        [c,Q] = double(Epost);
        [xE,yE] = Ellipse_Plot(inv(Q),[0;0],1,N2);
        qTest = T*[xE;yE] + repmat(q',1,length(xE));
        isect = checkIntersection3(vBound,vObs1,vObs2,ellBnd1,ellBnd2,ellCTmp1,ellCTmp2,qTest');
%         figure(3), plot(qTest(1,:),qTest(2,:))
        if any(isect),
            figure(3), plot(qTest(1,:),qTest(2,:),'r')
            break, 
        end
        Eprior = Epost;
        Eprior1 = Epost1;
        Epost = shape(Epost,1.1);  % heuristically increase by 10% and try again..
        Epost1 = shape(Epost1,1.1);
    end
    k
    tmp = double(Eprior)./double(Einit);
    rho = tmp(1,1);
end
    