function [randPt] = getRandom(ellOrPoly,vAvoid,ellBound,type)
% create a point contained within the intersection of an array of
% ellipsoids with either an array of ellipsoids or a polytope.

global hindx

% if isempty(hindx), hindx = 1; end

for j = 1:size(ellBound,2)
    ellArray{j} = [];
    for i = 1:size(ellBound,1)
        ellArray{j} = [ellArray; ellBound{i,j}];
    end
end
    
if strcmp(class(ellOrPoly),'ellipsoid')
    ellInitArray = ellOrPoly;
    ellEA = ellunion_ea(ellInitArray);
    [qx,px] = double(projection(ellEA,[1;0;0]));
    [qy,py] = double(projection(ellEA,[0;1;0]));
    [qth,pth] = double(projection(ellEA,[0;0;1]));
    qth = mod(qth+pi,2*pi)-pi;
    xMin = qx - sqrt(px);  xMax = qx + sqrt(px);
    yMin = qy - sqrt(py);  yMax = qy + sqrt(py);
    thMin = qth - sqrt(pth);  thMax = qth + sqrt(pth);
    
    switch type
        case 'random'
            for indx = 1:100
                tmpx = xMin + (xMax - xMin)*rand(1);
                tmpy = yMin + (yMax - yMin)*rand(1);
                tmpth = thMin + (thMax - thMin)*rand(1);
                isect = false;
                for mm = 1:length(vAvoid)
                    if inpolygon(tmpx,tmpy,vAvoid{mm}(:,1),vAvoid{mm}(:,2))
                        isect = true;
                    end
                end
                if ~isect
                    if isinternal_quick(ellInitArray,[tmpx;tmpy;tmpth])
                        if ~isempty(ellArray)
                            if isinternal_quick(ellArray,[tmpx;tmpy;tmpth])
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
                tmpx = xMin + (xMax - xMin)*haltonTmp(indx,1);
                tmpy = yMin + (yMax - yMin)*haltonTmp(indx,2);
                tmpth = thMin + (thMax - thMin)*haltonTmp(indx,3);
                isect = false;
                for mm = 1:length(vAvoid)
                    if inpolygon(tmpx,tmpy,vAvoid{mm}(:,1),vAvoid{mm}(:,2))
                        isect = true;
                    end
                end
                if ~isect
                    if isinternal_quick(ellInitArray,[tmpx;tmpy;tmpth])
                        if ~isempty(ellArray)
                            if isinternal_quick(ellArray,[tmpx;tmpy;tmpth])
                                hindx = indx;
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
    
else % a polytope
    vReg = ellOrPoly;
    
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
                            if isinternal_quick(ellArray,[tmpx;tmpy;tmpth])
                                break
                            end
                        else
                            break
                        end
                        %                 inside = 0;
                        %                 for i = 1:size(ellArray,1)
                        %                     [q,P] = double(ellArray{i,1});  % only need to check initial ellipse
                        %                     x = [tmpx; tmpy];
                        %                     if (q-x)'*P(1:2,1:2)*(q-x) <= 1
                        %                         inside = inside+1;
                        %                     end
                        %                 end
                        %                 if ~inside,
                        %                     break
                        %                 end
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
                            if isinternal_quick(ellArray,[tmpx;tmpy;tmpth])
                                hindx = indx;
                                break
                            end
                        else
                            hindx = indx;
                            break
                        end
                        %                 inside = 0;
                        %                 for i = 1:size(ellArray,1)
                        %                     [q,P] = double(ellArray{i,1});  % only need to check initial ellipse
                        %                     x = [tmpx; tmpy];
                        %                     if (q-x)'*P(1:2,1:2)*(q-x) <= 1
                        %                         inside = inside+1;
                        %                     end
                        %                 end
                        %                 if ~inside,
                        %                     hindx = indx;
                        %                     break
                        %                 end
                    end
                end
            end
    end
    
end
    
randPt = [tmpx tmpy tmpth]
