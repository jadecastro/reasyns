function [randPt,hindx] = getRandom2(vReg,vAvoid,ellBoundInv,ellInInv,type,inout,hindx,H,n,limsNonRegState,isCyclic)
% create a point contained inside a polytope with optional check for
% containment inside or outside the union of funnels.

if length(H) ~= 2, error('Workspaces of dimension other than two not yet supported. Sorry.'); end
if n == length(H), error('Empty non-region state vector not yet supported. Sorry.'); end
if size(limsNonRegState,2) ~= n - length(H), error('Dimension mismatch in limsNonRegState.'); end

% if isempty(hindx), hindx = 1; end

trialMat = repmat([0 2*pi -2*pi],n,1).*repmat(isCyclic,1,3);

if isempty(ellBoundInv), ellArray = []; end
for j = 1:size(ellBoundInv,2)
    ellArray{j} = [];
    for i = 1:size(ellBoundInv,1)
        ellArray{j} = [ellArray{j}; ellBoundInv{i,j}];
    end
end

minNonRegStates = limsNonRegState(1,:);
maxNonRegStates = limsNonRegState(2,:);

maxIter = 10000;

switch type
    case 'random'
        for indx = 1:1000
            randRegStates = rand(1,length(H));
            randNonRegStates = rand(1,n-length(H));
            tmpRegStates = (H\(min(vReg) + (max(vReg) - min(vReg)).*randRegStates*H)')';
            tmpNonRegStates = minNonRegStates + (maxNonRegStates - minNonRegStates).*randNonRegStates;
            isect = false;
            for mm = 1:length(vAvoid)
                if inpolygon(tmpRegStates(1),tmpRegStates(2),vAvoid{mm}(:,1),vAvoid{mm}(:,2))
                    isect = true;
                end
            end
            if ~isect
                if inpolygon(tmpRegStates(1),tmpRegStates(2),vReg(:,1),vReg(:,2))
                    if ~isempty(ellArray)
                        for j = 1:length(ellArray)
                            if any(isCyclic)
                                isect(j) = any(isinternal_quickInv(ellArray{j},repmat([tmpRegStates tmpNonRegStates]',1,3) + trialMat));
                            else
                                isect(j) = any(isinternal_quickInv(ellArray{j},[tmpRegStates tmpNonRegStates]'));
                            end
                        end
                        if strcmp(inout,'inside'),
                            if all(isect)
                                break
                            end
                        else  % then we're in collision if outside ellipses
                            if ~all(isect)
                                break
                            end
                        end
                    else
                        break
                    end
                end
            end
        end
        
    case 'discrep'
        haltonTmp = haltonset(n);
        for indx = hindx+1:hindx+maxIter
            randRegStates = haltonTmp(indx,1:length(H));
            randNonRegStates = haltonTmp(indx,length(H)+1:n);
            tmpRegStates = (H\(min(vReg) + (max(vReg) - min(vReg)).*randRegStates*H)')';
            tmpNonRegStates = minNonRegStates + (maxNonRegStates - minNonRegStates).*randNonRegStates;
            isect = false;
            for mm = 1:length(vAvoid)
                if inpolygon(tmpRegStates(1),tmpRegStates(2),vAvoid{mm}(:,1),vAvoid{mm}(:,2))
                    isect = true;
                end
            end
            if ~isect
                if inpolygon(tmpRegStates(1),tmpRegStates(2),vReg(:,1),vReg(:,2))
                    if ~isempty(ellArray)
                        for j = 1:length(ellArray)
%                             if length(ellArray{j}) > 1000
%                                 tmp = downsample(ellArray{j},10);
%                                 isect(j) = any(isinternal_quick(tmp,[tmpRegStates tmpNonRegStates]'));
%                             else
                            if any(isCyclic)
                                isect(j) = any(isinternal_quickInv(ellArray{j},repmat([tmpRegStates tmpNonRegStates]',1,3) + trialMat));
                            else
                                isect(j) = any(isinternal_quickInv(ellArray{j},[tmpRegStates tmpNonRegStates]'));
                            end
%                             end
                        end
                        if strcmp(inout,'inside'),  % accept if funnel conjunction holds
                            if all(isect)
                                break
                            end
                        elseif strcmp(inout,'insidedisjunct'),  % accept if funnel disjunction holds
                            if any(isect)
                                break
                            end
                        end
                        if ~isempty(ellInInv)
                            for k = 1:length(ellInInv)
                                if any(isCyclic)
                                    isect(j) = any(isinternal_quickInv(ellInInv{k},repmat([tmpRegStates tmpNonRegStates]',1,3) + trialMat));
                                else
                                    isect(j) = any(isinternal_quickInv(ellInInv{k},[tmpRegStates tmpNonRegStates]'));
                                end
                            end
                            if any(isect)
                                break
                            end
                        end
                    else
                        break
                    end
                end
            end
        end
end

if indx == hindx+maxIter
    hindx = indx;
    disp('maxIter reached')
    randPt = [];
else
    hindx = indx;
    randPt = [tmpRegStates tmpNonRegStates];
end
