function [isect,regID] = checkIntersection3a(regBnd,regObs1,regObs2,Xbnd1,Xbnd2,Xin1,Xin2,q)
% Polytopes are assumed to be 2-D, ellipses are n-D
            
for ii = 1:size(q,1)
%     ii
    for mm = 1:length(vObs1), 
        isectTmp1(mm) = inpolygon(q(ii,1),q(ii,2),vObs1{mm}(:,1),vObs1{mm}(:,2));
    end
    for mm = 1:length(vObs2), 
        isectTmp2(mm) = inpolygon(q(ii,1),q(ii,2),vObs2{mm}(:,1),vObs2{mm}(:,2));
    end
    
    Isect(ii) = double(any(q(ii,1:2) < min(vBound)) || any(q(ii,1:2) > max(vBound)));
    if Isect(ii)
        0
    end
    
    if ~any(isectTmp1) && double(all(q(ii,1:2) > min(vBound)) && all(q(ii,1:2) < max(vBound)))  % point is inside R1.  Now check if inside corresponding reach tube
        isect1 = false;
        isectC = false;
        if ~isempty(ellC1)
            for mm = 1:length(ellC1)
                if isempty(ellC1{mm}), break, end
                isectC(mm) = ~isinternal_quick(ellC1{mm},q(ii,:)','u');
            end
        end
        if ~isempty(ellBound1)
            for nn = 1:size(ellBound1,2)
                for mm = 1:size(ellBound1,1),
                    if isstruct(ellBound1) || iscell(ellBound1)
                        if isempty(ellBound1{mm,nn}), break, end
                        isect1(mm) = ~isinternal_quick(ellBound1{mm,nn},q(ii,:)','u');
                    else % we're dealing with an array
                        isect1(mm) = ~isinternal_quick(ellBound1(mm,nn),q(ii,:)','u');
                    end
                end
                isect2(nn) = all(isect1);  % only if any given point lies outside all of the funnels, then flag it
            end
        else
            isect2 = false;
        end
        if any(isect2) && all(isectC), 
            1
            Isect(ii) = 1;
        end  % if any of the points lie outside the intersection of multiple outgoing reach tubes AND outside all inward funnels, flag it
    elseif ~any(isectTmp2) && double(all(q(ii,1:2) > min(vBound)) && all(q(ii,1:2) < max(vBound)))  % point is inside R2.  Now check if inside corresponding reach tube
        isect1 = false;
        isectC = false;
        if ~isempty(ellC2)
            for mm = 1:length(ellC2)
                if isempty(ellC2{mm}), break, end
                isectC(mm) = ~isinternal_quick(ellC2{mm},q(ii,:)','u');
            end
        end
        if ~isempty(ellBound2)
            for nn = 1:size(ellBound2,2)
                for mm = 1:size(ellBound2,1),
                    if isstruct(ellBound2) || iscell(ellBound2)
                        if isempty(ellBound2{mm,nn}), break, end
                        isect1(mm) = ~isinternal_quick(ellBound2{mm,nn},q(ii,:)','u');
                    else % we're dealing with an array
                        isect1(mm) = ~isinternal_quick(ellBound2(mm,nn),q(ii,:)','u');
                    end
                end
                isect2(nn) = all(isect1);  % only if any given point lies outside all of the funnels, then flag it
            end
        else
            isect2 = false;
        end
        if any(isect2) && all(isectC), 
            2
            Isect(ii) = 1;
        end  % if any of the points lie outside the intersection of multiple outgoing reach tubes AND outside all inward funnels, flag it
    else % point is outside both regions or boundary
        Isect(ii) = 1
    end
end

isect = Isect > 0;
