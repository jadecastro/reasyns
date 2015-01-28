function [randPt] = getCenterWithRandOrientation(vReg,vAvoid,ellBound)

ellArray = [];
for i = 1:size(ellBound,1)
    for j = 1:size(ellBound,2)
        ellArray = [ellArray; ellBound{i,j}];
    end
end

% Place initial point at center of region
tmpx = (max(vReg(:,1)) + min(vReg(:,1)))/2;
tmpy = (max(vReg(:,2)) + min(vReg(:,2)))/2;
tmpth = -pi + 2*pi*rand(1);

% If in obstacle or outside invariant, recompute its location in a deterministic fashion
haltonTmp = haltonset(3);
for indx = 2:100
    isect = false;
    for mm = 1:length(vAvoid)
        if inpolygon(tmpx,tmpy,vAvoid{mm}(:,1),vAvoid{mm}(:,2))
            isect = true;
        end
    end
    if ~isect
        if ~isempty(ellBound)
            if isinternal_quick(ellArray,[tmpx;tmpy;tmpth])
                break
            end
        else
            break
        end
    end
    tmpx = min(vReg(:,1)) + (max(vReg(:,1)) - min(vReg(:,1)))*haltonTmp(indx,1);
    %tmpy = min(vReg(:,2)) + (max(vReg(:,2)) - min(vReg(:,2)))*haltonTmp(indx,2);
end

randPt = [tmpx tmpy tmpth];
