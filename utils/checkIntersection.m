function [isect] = checkIntersection(vObs,vBound,ellBound,q,H)
% Polytopes are assumed to be 2-D, ellipses are n-D

if length(H) ~= 2, error('Workspaces of dimension other than two not yet supported. Sorry.'); end

isect = 0;
for mm = 1:length(vObs), 
    isect = isect + any(inpolygon(q(:,1),q(:,2),vObs{mm}(:,1),vObs{mm}(:,2))); 
end
for ii = 1:size(q,1)
    isect = isect + double(any(q(ii,1:length(H)) < min(vBound)) || any(q(ii,1:length(H)) > max(vBound)));
    if ~isempty(ellBound)
        for nn = 1:size(ellBound,2)
            isect1 = false;
            for mm = 1:size(ellBound,1),
                if isempty(ellBound{mm,nn}), break, end
                isect1(mm) = any(~isinternal_quick(ellBound{mm,nn},q(ii,:),'u'));
            end
            isect2(nn) = all(isect1);  % if any given point lies outside all of the funnels, then flag it
        end
        if any(isect2), isect = 1; end  % if any of the points lie outside the intersection of multiple reach tubes, flag it
    end
end

isect = isect > 0;
