function [randPt] = getRandomOutside(vReg,vAvoid,ellAvoid,type)
% create a point contained outside the intersection of an array of
% ellipsoids with a polytope.

global hindx

% if isempty(hindx), hindx = 1; end

for j = 1:size(ellBound,2)
    ellArray{j} = [];
    for i = 1:size(ellBound,1)
        ellArray{j} = [ellArray; ellBound{i,j}];
    end
end

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
                        if all(isect)
                            break
                        end
                    else
                        break
                    end
                end
            end
        end
        
    case 'discrep'
        haltonTmp = haltonset(3);
        for indx = hindx+1:hindx+100
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
                            isect(j) = any(isinternal_quick(ellArray{j},[tmpx;tmpy;tmpth]));
                        end
                        if all(isect)
                            break
                        end
                    else
                        hindx = indx;
                        break
                    end
                end
            end
        end
end

randPt = [tmpx tmpy tmpth]
