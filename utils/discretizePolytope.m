function v1 = discretizePolytope(v,dX)

if iscell(v)
    vv = v;
else
    vv{1} = v;
end

for i = 1:length(vv)
    V = vv{i};
    indxm = size(V,1);
    V1 = [];
    for j = 1:size(V,1)
        indx = j;
        dist = sqrt((V(indx,1) - V(indxm,1))^2 + (V(indx,2) - V(indxm,2))^2);
        deltax = (V(indx,1) - V(indxm,1));
        deltay = (V(indx,2) - V(indxm,2));
        k = 1;
        V1 = [V1; V(indxm,:)];
        while k*dX < dist
            x = k*dX*deltax/dist + V(indxm,1);
            y = k*dX*deltay/dist + V(indxm,2);
            V1 = [V1; x y];
            k = k+1;
        end
        indxm = indx;
    end
    vv1{i} = V1;
end

v1 = vv1{1};
% if iscell(v)
%     v1 = vv1;
% elseif isnumeric(v)
%     v1 = vv1{1};
% end
