function [isect] = checkIntersection4(ellBoundInv1,q,H,n,isCyclic)
% Polytopes are assumed to be 2-D, ellipses are n-D
            
if length(H) ~= 2, error('Workspaces of dimension other than two not yet supported.'); end
if n == length(H), error('Empty non-region state vector not yet supported.'); end

trialMat = repmat([0 2*pi -2*pi],n,1).*repmat(isCyclic,1,3);

for ii = 1:size(q,1)
    %     ii
    Isect(ii) = false;
    for nn = 1:size(ellBoundInv1,2)
        isect1 = false;
        for mm = 1:size(ellBoundInv1,1),
            if isstruct(ellBoundInv1) || iscell(ellBoundInv1)
                if isempty(ellBoundInv1{mm,nn}), break, end
                if any(isCyclic)
                    isect1(mm) = ~any(isinternal_quickInv(ellBoundInv1{mm,nn},repmat(q(ii,:)',1,3)+trialMat,'u'));
                else
                    isect1(mm) = ~any(isinternal_quickInv(ellBoundInv1{mm,nn},q(ii,:)','u'));
                end
            else % we're dealing with an array
                if any(isCyclic)
                    isect1(mm) = ~any(isinternal_quickInv(ellBoundInv1(mm,nn),repmat(q(ii,:)',1,3)+trialMat,'u'));
                else
                    isect1(mm) = ~any(isinternal_quickInv(ellBoundInv1(mm,nn),q(ii,:)','u'));
                end
            end
        end
        isect2(nn) = all(isect1);  % only if any given point lies outside all of the funnels, then flag it
    end
    if any(isect2),
        Isect(ii) = true;
    end  % if any of the points lie outside the intersection of multiple outgoing reach tubes, flag it
    if Isect(ii)
        break
    end
end

if exist('Isect') 
    isect = any(Isect > 0);
else
    isect = true;
end
