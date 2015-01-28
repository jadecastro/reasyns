function [XbndR,XinR,tmpV,tmpV1,tmpV2] = getCurrentFunnel(X0,regBnd,reg1,reg2,ellBnd1,ellBnd2,Xbnd1,Xbnd2,ellIn1,ellIn2,Xin1,Xin2)

[XbndR,XinR,tmpV,tmpV1,tmpV2] = deal([]);

tmpV = regBnd{1}.v;
for j = 1:length(reg1)
    tmpV1{j} = reg1{j}.v;
end
for j = 1:length(reg2)
    tmpV2{j} = reg2{j}.v;
end
[isect1] = checkIntersection(tmpV1,tmpV,[],X0',eye(2));
[isect2] = checkIntersection(tmpV2,tmpV,[],X0',eye(2));
% TODO: fix intersections!

if ~isect1
    disp('in reg 1')
    reg = reg1;
    Xbnd = Xbnd1;
    Xin = Xin1;
    ellBnd = ellBnd1;
    ellIn = ellIn1;
else
    disp('in reg 2')
    reg = reg2;
    Xbnd = Xbnd2;
    Xin = Xin2;
    ellBnd = ellBnd2;
    ellIn = ellIn2;
end

for n = 1:size(Xbnd,2)  % For all outgoing transitions
    for m = 1:size(Xbnd,1) % for all transition funnels
        %                 d = distance(ellBnd{m,n},ball);
        clear d
        XbndR{m,n} = [];
        if ~isempty(ellBnd{m,n})
            for j = 1:length(ellBnd{m,n})
                d(j) = isinternal_quick(ellBnd{m,n}(j),X0);
            end
            idxSet = find(d);
            for p = idxSet
                XbndR{m,n} = [XbndR{m,n}; Xbnd{m,n}(p)];
            end
        end
    end
end

for m = 1:length(Xin) % for all transition funnels
    %             d = distance(ellIn{m},ball);
    clear d
    XinR{m} = [];
    if ~isempty(ellIn{m})
        for j = 1:length(ellIn{m})
            d(j) = isinternal_quick(ellIn{m}(j),X0);
        end
        idxSet = find(d);
        for p = idxSet
            XinR{m} = [XinR{m}; Xin{m}(p)];
        end
    end
end