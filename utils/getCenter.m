function [randPt] = getCenter(vReg,vAvoid,ellBound,H,n,limsNonRegState,nonRegType)

if length(H) ~= 2, error('Workspaces of dimension other than two not yet supported. Sorry.'); end
if n == length(H), error('Empty non-region state vector not yet supported. Sorry.'); end
if size(limsNonRegState,2) ~= n - length(H), error('Dimension mismatch in limsNonRegState.'); end

ellArray = [];
for i = 1:size(ellBound,1)
    for j = 1:size(ellBound,2)
        ellArray = [ellArray; ellBound{i,j}];
    end
end

minNonRegStates = limsNonRegState(1,:);
maxNonRegStates = limsNonRegState(2,:);

% Place initial point at center of region
tmpRegStates = (max(vReg) + min(vReg))./2;
if strmatch(nonRegType,'zero')
    tmpNonRegStates = zeros(1,n - length(H));
elseif strmatch(nonRegType,'rand')
    tmpNonRegStates = minNonRegStates + (maxNonRegStates - minNonRegStates).*rand(1,n-length(H));
else
    error('Invalid nonRegType entry.')
end

% If in obstacle or outside invariant, recompute its location in a deterministic fashion
haltonTmp = haltonset(3);
for indx = 2:100
    isect = false;
    for mm = 1:length(vAvoid)
        if inpolygon(tmpRegStates(1),tmpRegStates(2),vAvoid{mm}(:,1),vAvoid{mm}(:,2))
            isect = true;
        end
    end
    if ~isect
        if ~isempty(ellBound)
            if isinternal_quick(ellArray,[tmpRegStates tmpNonRegStates])
                break
            end
        else
            break
        end
    end
    tmpRegStates = (H\(min(vReg) + (max(vReg) - min(vReg)).*haltonTmp(indx,1:length(H))*H)')';
end

if indx == 100
    error('Can''t find a non-colliding point.')
end

if strmatch(nonRegType,'zero')
    randPt = [tmpRegStates];
elseif strmatch(nonRegType,'rand')
    randPt = [tmpRegStates tmpNonRegStates];
else
    error('Invalid nonRegType entry.')
end
