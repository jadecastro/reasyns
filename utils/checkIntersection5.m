function [isect,qBadIndx,ellInInvR,ellBoundInvR] = checkIntersection5(vBound,vObs1,ellBoundInv1,ellInInv1,q,H,n,isCyclic)
% Polytopes are assumed to be 2-D, ellipses are n-D

global debugFlg

if sum(isCyclic) > 1, error('number of cyclic dimensions cannot exceed 1'), end
if length(isCyclic) ~= n, error('number of entries in isCylic must be equal to n!'), end
if length(H) ~= 2, error('Workspaces of dimension other than two not yet supported. Sorry.'); end

if ~isempty(ellInInv1)
    ellInInv1notEmpty = true;
else
    ellInInv1notEmpty = false;
end

trialMat = repmat([0 2*pi -2*pi],n,1).*repmat(isCyclic,1,3);

% if debugFlg
%     size(q,1)
% end

Isect = false;

qBadIndx = [];
q2Indx = [];
ellInInvR = [];
ellBoundInvR = [];
for ii = 1:size(q,1)
    %     if debugFlg
    %         ii
    %         q(ii,:)
    %     end
    
    if nargout == 4
        Isect(ii) = true;
        isect2(1:size(ellBoundInv1,2)) = true;
        for nn = 1:size(ellBoundInv1,2)  % For all outgoing transitions
            for m = 1:size(ellBoundInv1,1) % for all transition funnels
                clear d
                ellBoundInvR{m,nn} = [];
                if ~isempty(ellBoundInv1{m,nn})
                    for j = 1:length(ellBoundInv1{m,nn})
                        if any(isCyclic)
                            d(j) = any(isinternal_quickInv(ellBoundInv1{m,nn}(j),repmat(q(ii,:)',1,3)+trialMat,'u'));
                        else
                            d(j) = isinternal_quickInv(ellBoundInv1{m,nn}(j),q(ii,:)','u');
                        end
                    end
                    idxSet = find(d);
                    for p = idxSet
                        isect2(nn) = false;
                        ellBoundInvR{m,nn} = [ellBoundInvR{m,nn}; ellBoundInv1{m,nn}(p)];
                    end
                end
            end
        end
        Isect(ii) = any(isect2);
        
        if ellInInv1notEmpty
            for m = 1:length(ellInInv1) % for all transition funnels
                %             d = distance(ellIn{m},ball);
                clear d
                ellInInvR{m} = [];
                if ~isempty(ellInInv1{m})
                    for j = 1:length(ellInInv1{m})
                        if any(isCyclic)
                            d(j) = any(isinternal_quickInv(ellInInv1{m}(j),repmat(q(ii,:)',1,3)+trialMat,'u'));
                        else
                            d(j) = isinternal_quickInv(ellInInv1{m}(j),q(ii,:)','u');
                        end
                    end
                    idxSet = find(d);
                    for p = idxSet
                        Isect(ii) = false;
                        ellInInvR{m} = [ellInInvR{m}; ellInInv1{m}(p)];
                    end
                end
            end
        end
        
    else  % This should be faster ...
        if ~isempty(vObs1)
            for mm = 1:length(vObs1),
                isectTmp1(mm) = inpolygon(q(ii,1),q(ii,2),vObs1{mm}(:,1),vObs1{mm}(:,2));
            end
        else
            isectTmp1 = [];
        end
        
        if ~isempty(vBound)
            Isect(ii) = double(any(q(ii,1:length(H)) < min(vBound)) || any(q(ii,1:length(H)) > max(vBound)));
        else
            Isect(ii) = false;
        end
        
        isectC = false;
        if ellInInv1notEmpty
            for mm = 1:length(ellInInv1)
                %             if isempty(ellInInv1{mm}), break, end
                if ~isempty(ellInInv1{mm})
                    if any(isCyclic)
                        isectC1(mm) = ~any(isinternal_quickInv(ellInInv1{mm},repmat(q(ii,:)',1,3)+trialMat,'u'));  % check +/- 2*pi also
                    else
                        isectC1(mm) = ~any(isinternal_quickInv(ellInInv1{mm},q(ii,:)','u'));
                    end
                else
                    isectC(mm) = true;
                end
            end
        else
            isectC = [];
        end
        
        %         if double(all(q(ii,1:length(H)) > min(vBound)) && all(q(ii,1:length(H)) < max(vBound)))
        if ~any(isectTmp1)  % point is inside R1.  Now check if inside corresponding reach tube
            isect1 = false;
            if ~isempty(ellBoundInv1)
                for nn = 1:size(ellBoundInv1,2)
                    for mm = 1:size(ellBoundInv1,1),
                        if isstruct(ellBoundInv1) || iscell(ellBoundInv1)
                            %                         if isempty(ellBoundInv1{mm,nn}), break, end
                            if ~isempty(ellBoundInv1{mm,nn})
                                if any(isCyclic)
                                    isect1(mm) = ~any(isinternal_quickInv(ellBoundInv1{mm,nn},repmat(q(ii,:)',1,3)+trialMat,'u'));
                                else
                                    isect1(mm) = ~any(isinternal_quickInv(ellBoundInv1{mm,nn},q(ii,:)','u'));
                                end
                            else
                                isect1(mm) = true;
                            end
                        else % we're dealing with an array
                            if ~isempty(ellBoundInv1(mm,nn))
                                if any(isCyclic)
                                    isect1(mm) = ~any(isinternal_quickInv(ellBoundInv1(mm,nn),repmat(q(ii,:)',1,3)+trialMat,'u'));
                                else
                                    isect1(mm) = ~any(isinternal_quickInv(ellBoundInv1(mm,nn),q(ii,:)','u'));
                                end
                            else
                                isect1(mm) = true;
                            end
                        end
                    end
                    isect2(nn) = all(isect1);  % only if any given point lies outside all of the funnels, then flag it
                end
            else
                isect2 = false;
            end
            if any(isect2) && all(isectC),
                Isect(ii) = true;
            end  % if any of the points lie outside the intersection of multiple outgoing reach tubes AND outside all inward funnels, flag it
        else % point is outside the region or boundary
            Isect(ii) = true;
        end
        %         end
    end
    if Isect(ii)
        if nargout == 2 || nargout == 4  % we're expecting to run through the points till completion
            qBadIndx = [qBadIndx; ii];
        else
            break,
        end
    end
end

if exist('Isect')
    isect = any(Isect);
else
    isect = true;
end
