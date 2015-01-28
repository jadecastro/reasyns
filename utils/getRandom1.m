function [randPt,hindx] = getRandom1(vReg,vAvoid,ellBound,type,inout,hindx)
% create a point contained either within or outside of an array of
% ellipsoids and also inside a polytope.

% if isempty(hindx), hindx = 1; end

if isempty(ellBound), ellArray = []; end
for j = 1:size(ellBound,2)
    ellArray{j} = [];
    for i = 1:size(ellBound,1)
        ellArray{j} = [ellArray{j}; ellBound{i,j}];
    end
end

maxIter = 1000;

switch type
    case 'random'
        for indx = 1:1000
            tmpx = min(vReg(:,1)) + (max(vReg(:,1)) - min(vReg(:,1)))*rand(1);
            tmpy = min(vReg(:,2)) + (max(vReg(:,2)) - min(vReg(:,2)))*rand(1);
            tmpth = -pi + 2*pi*rand(1);
            isect = false;
            for mm = 1:length(vAvoid)
                if inpolygon(tmpx,tmpy,vAvoid{mm}(:,1),vAvoid{mm}(:,2))
                    isect = true;
                end
            end
            if ~isect
                if inpolygon(tmpx,tmpy,vReg(:,1),vReg(:,2))
                    if ~isempty(ellArray)
                        for j = 1:length(ellArray)
                            isect(j) = any(isinternal_quick(ellArray{j},[tmpx;tmpy;tmpth]));
                        end
                        if strcmp(inout,'inside'),
                            if all(isect)
                                break
                            end
                        else
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
        haltonTmp = haltonset(3);
        for indx = hindx+1:hindx+maxIter
            tmpx = min(vReg(:,1)) + (max(vReg(:,1)) - min(vReg(:,1)))*haltonTmp(indx,1);
            tmpy = min(vReg(:,2)) + (max(vReg(:,2)) - min(vReg(:,2)))*haltonTmp(indx,2);
            tmpth = -pi + 2*pi*haltonTmp(indx,3);
            isect = false;
            for mm = 1:length(vAvoid)
                if inpolygon(tmpx,tmpy,vAvoid{mm}(:,1),vAvoid{mm}(:,2))
                    isect = true;
                end
            end
            if ~isect
                if inpolygon(tmpx,tmpy,vReg(:,1),vReg(:,2))
                    if ~isempty(ellArray)
                        for j = 1:length(ellArray)
%                             if length(ellArray{j}) > 1000
%                                 tmp = downsample(ellArray{j},10);
%                                 isect(j) = any(isinternal_quick(tmp,[tmpx;tmpy;tmpth]));
%                             else
                                isect(j) = any(isinternal_quick(ellArray{j},[tmpx;tmpy;tmpth]));
%                             end
                        end
                        if strcmp(inout,'inside'),
                            if all(isect)
                                break
                            end
                        else
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
end

if indx == hindx+maxIter
    hindx = indx;
    disp('maxIter reached')
    randPt = [];
else
    hindx = indx;
    randPt = [tmpx tmpy tmpth];
end
