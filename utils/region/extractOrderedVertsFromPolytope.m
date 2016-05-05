function verts = extractOrderedVertsFromPolytope(p)

v = [];
P1 = [];
for i = 1:length(p)
    
    vNew = extreme(p(i));
    
    % inflate the regions a very slight amount
    regTmp = Region('tmp',vNew,[],0.01);
    vNew = extreme(regTmp.p);
    
    amean = mean(vNew(:,1));  bmean = mean(vNew(:,2));
    [~, indx] = sort(atan2(vNew(:,1)-amean,vNew(:,2)-bmean));
    vNew = flipud(vNew(indx,:));
    if i > 1
        if true %verLessThan('matlab','7.15')
            P2.x = vNew(:,1)-0.01;  P2.y = vNew(:,2);  P2.hole = 0;
            
            P3 = PolygonClip(P1,P2,3);
            
            P1 = [P1, P3];
        else
            [a,b] = polybool('union',v(:,1),v(:,2),vNew(:,1),vNew(:,2));
            v = [a b];
        end
    else
        if true %verLessThan('matlab','7.15')
            P1.x = vNew(:,1);  P1.y = vNew(:,2);  P1.hole = 0;
        else
            v = vNew;
        end
    end
end
if true %verLessThan('matlab','7.15')
    %                 figure(1000)
    %                 for j = 1:length(P3)
    %                     hold on
    %                     obj = patch(P3(j).x,P3(j).y,1);
    %                 end
    %                 tmp = get(gco,'Vertices');
    %                 v = tmp;
    if i > 1
        v = [P3.x, P3.y];
    else
        if ~isempty(P1)
            v = [P1.x, P1.y];
        else
            v = [];
        end
    end
    
    if ~isempty(v)
        amean = mean(v(:,1));  bmean = mean(v(:,2));
        [~, indx] = sort(atan2(v(:,1)-amean,v(:,2)-bmean));
        a = v(indx,1);  b = v(indx,2);
    else
        a = [];  b = [];
    end
else
    [a,b] = poly2ccw(v(:,1),v(:,2));
end

verts = [a,b];
end