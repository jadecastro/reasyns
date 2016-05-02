function [isect,qGoodIndx] = checkIntersection(vBound,vObs1,vObs2,ellBoundInv1,ellBoundInv2,ellInInv1,ellInInv2,q,succRegChk,ac,reg,sys)
% Polytopes are assumed to be 2-D, ellipses are n-D
            
global debugFlg

n = sys.sysparams.n;
isCyclic = sys.sysparams.isCyclic;

if ~(strcmp(succRegChk,'finalPt') || strcmp(succRegChk,'allPts')), error('succRegChk variable must either be ''finalPt'' or ''allPts'''); end

trialMat = repmat([0 2*pi -2*pi],n,1).*repmat(isCyclic,1,3);

% if debugFlg
%     size(q,1)
% end

regID = [];
qGoodIndx = [];
qGoodIndxVec = [];
q2Indx = [];
for ii = 1:size(q,1)
    %     if debugFlg
    %         ii
    %         q(ii,:)
    %     end
    
    if ~isempty(reg)
        isectTmp1 = ~isinside(reg,sys,q);
    elseif ~isempty(vObs1)
        for mm = 1:length(vObs1),
            isectTmp1(mm) = inpolygon(q(ii,1),q(ii,2),vObs1{mm}(:,1),vObs1{mm}(:,2));
        end
    else
        isectTmp1 = [];
    end
    
    if ~isempty(reg)
        isectTmp2 = ~isinside(reg,sys,q);
    elseif ~isempty(vObs2)
        for mm = 1:length(vObs2),
            isectTmp2(mm) = inpolygon(q(ii,1),q(ii,2),vObs2{mm}(:,1),vObs2{mm}(:,2));
        end
    else
        isectTmp2 = [];
    end
    
    output = sys.state2SEconfig([],q(ii,:),[])';
    Isect(ii) = double(any(output(1:2) < min(vBound)) || any(output(1:2) > max(vBound)));
    Isect1(ii) = Isect(ii);
    
    if Isect(ii)
%         0
    end
    
    isectC1 = false;
    if ~isempty(ellInInv1)
        for mm = 1:length(ellInInv1)
            %             if isempty(ellInInv1{mm}), break, end
            if ~isempty(ellInInv1{mm})
                if any(isCyclic)
                    isectC1(mm) = ~any(isinternal_quickInv(ellInInv1{mm},repmat(q(ii,:)',1,3)+trialMat,'u'));  % check +/- 2*pi also
                else
                    isectC1(mm) = ~any(isinternal_quickInv(ellInInv1{mm},q(ii,:)','u'));
                end
            else
                isectC1(mm) = 1;
            end
        end
    else
        isectC1 = [];
    end
    isectC2 = false;
    if ~isempty(ellInInv2)
        for mm = 1:length(ellInInv2)
            %             if isempty(ellInInv2{mm}), break, end
            if ~isempty(ellInInv2{mm})
                if any(isCyclic)
                    isectC2(mm) = ~any(isinternal_quickInv(ellInInv2{mm},repmat(q(ii,:)',1,3)+trialMat,'u'));
                else
                    isectC2(mm) = ~any(isinternal_quickInv(ellInInv2{mm},q(ii,:)','u'));
                end
            else
                isectC2(mm) = 1;
            end
        end
    else
        isectC2 = [];
    end
    isectC = [isectC1 isectC2];
    
    if ~any(isectTmp1) && double(all(output(1:2) > min(vBound)) && all(output(1:2) < max(vBound)))  % point is inside R1.  Now check if inside corresponding reach tube
        regID(ii) = 0;
        
        isect1 = false;
        if ~isempty(ac)
            for nn = 1:numel(ac)
                isect2(nn) = ~any(isinternal(ac(nn),repmat(q(ii,:)',1,3)+trialMat,'u'));
            end
        elseif ~isempty(ellBoundInv1)
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
                            isect1(mm) = 1;
                        end
                    else % we're dealing with an array
                        if ~isempty(ellBoundInv1(mm,nn))
                            if any(isCyclic)
                                isect1(mm) = ~any(isinternal_quickInv(ellBoundInv1(mm,nn),repmat(q(ii,:)',1,3)+trialMat,'u'));   
                            else
                                isect1(mm) = ~any(isinternal_quickInv(ellBoundInv1(mm,nn),q(ii,:)','u'));
                            end
                        else
                            isect1(mm) = 1;
                        end
                    end
                end
                isect2(nn) = all(isect1);  % only if any given point lies outside all of the funnels, then flag it
            end
        else
            isect2 = false;
        end
        if any(isect2) && all(isectC),
%             1
            Isect(ii) = 1;
            IsectR1(ii) = 1;
            Isect1(ii) = 1;
        end  % if any of the points lie outside the intersection of multiple outgoing reach tubes AND outside all inward funnels, flag it
    elseif ~any(isectTmp2) && double(all(output(1:2) > min(vBound)) && all(output(1:2) < max(vBound)))  % point is inside R2.  Now check if inside corresponding reach tube
        regID(ii) = 1;
        
        q2Indx = [q2Indx; ii];  % Store a running list of indices which are in R2
        
        isect1 = false;
        if ~isempty(ellBoundInv2)
            for nn = 1:size(ellBoundInv2,2)
                for mm = 1:size(ellBoundInv2,1),
                    if isstruct(ellBoundInv2) || iscell(ellBoundInv2)
                        %                         if isempty(ellBoundInv2{mm,nn}), break, end
                        if ~isempty(ellBoundInv2{mm,nn})
                            if any(isCyclic)
                                isect1(mm) = ~any(isinternal_quickInv(ellBoundInv2{mm,nn},repmat(q(ii,:)',1,3)+trialMat,'u'));
                            else
                                isect1(mm) = ~any(isinternal_quickInv(ellBoundInv2{mm,nn},q(ii,:)','u'));
                            end
                        else
                            isect1(mm) = 1;
                        end
                    else % we're dealing with an array
                        if ~isempty(ellBoundInv2(mm,nn))
                            if any(isCyclic)
                                isect1(mm) = ~any(isinternal_quickInv(ellBoundInv2(mm,nn),repmat(q(ii,:)',1,3)+trialMat,'u'));
                            else
                                isect1(mm) = ~any(isinternal_quickInv(ellBoundInv2(mm,nn),q(ii,:)','u'));
                            end
                        else
                            isect1(mm) = 1;
                        end
                    end
                end
                isect2(nn) = all(isect1);  % only if any given point lies outside all of the funnels, then flag it
            end
        else
            isect2 = false;
        end
        Isect1(ii) = 0;
        if q2Indx(end) ~= ii
            Isect1(ii) = 1;
        end
        if any(isect2) && all(isectC2),
%             2
%             q2Indx
            Isect(ii) = 1;
            if length(q2Indx) >= 2   % our criterion for ensuring that the trajectory is sufficiently inside R2
                qGoodIndx = q2Indx(end-1);
                qGoodIndxVec = [qGoodIndxVec qGoodIndx];
            end
        end  % if any of the points lie outside the intersection of multiple outgoing reach tubes AND outside all inward funnels, flag it
    else % point is outside both regions or boundary
        regID(ii) = 0;
        Isect(ii) = 1;
        Isect1(ii) = 1;
    end
    if Isect(ii) && strcmp(succRegChk,'allPts')
%         Isect
        break
    end
    if Isect1(ii) && strcmp(succRegChk,'finalPt')  % when in R2, break only when region is violated
%         Isect
%         Isect1
        break
    end
end

if strcmp(succRegChk,'allPts')
    if exist('Isect')
        isect = any(Isect > 0);
    else
        isect = true;
    end
elseif strcmp(succRegChk,'finalPt')
    if ~any(isnan(regID(1:end-1)))
        indxGoodFinal = ~Isect & regID;
        qGoodIndx = find(indxGoodFinal,1,'last');
    else
        qGoodIndx = [];
    end
    if exist('Isect1')
        isect = any(Isect1 > 0);
    else
        isect = true;
    end
end
