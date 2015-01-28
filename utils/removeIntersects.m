function [Vx,Vy] = removeIntersects(X,Vx,Vy)

global x

% Remove violations on the external facets
if iscell(X)
    ve = X.extVert';
    vi = X.intVert';
    pe = X.ext';
    pi = X.int';
else
    ve{1} = X.extVert';
    vi{1} = X.intVert';
    pe{1} = X.ext';
    pi{1} = X.int';
end

for i = 1:length(ve)
    if ~isempty(ve{i})
        for k = 1:size(Vx,2)
            isect1 = any(double(subs(pe{i},x,[Vx(1,k);Vy(1,k);0;0])) < 0);
            isect2 = any(double(subs(pe{i},x,[Vx(1,k);Vy(1,k);0;0])) < 0);
            if isect1 || isect2,
                Vx(:,k) = [NaN; NaN];
                Vy(:,k) = [NaN; NaN];
            end
        end
        for j = 1:size(ve{i},1)
            if j == 1, jm = size(ve{i},1); end
            for k = 1:size(Vx,2)
                [isect,tmp1,tmp2,tmp3] = intersectPoint(Vx(1,k),Vy(1,k),Vx(2,k),Vy(2,k),ve{i}(j,1),ve{i}(j,2),ve{i}(jm,1),ve{i}(jm,2));
                if isect,
                    Vx(:,k) = [NaN; NaN];
                    Vy(:,k) = [NaN; NaN];
                end
            end
            jm = j;
        end
    end
end

rmindx = isnan(Vx(i,:));
Vx(:,rmindx) = [];
Vy(:,rmindx) = [];

% Now do for internal facets
for i = 1:length(vi)
    if ~isempty(vi{i})
        for k = 1:size(Vx,2)
            isect1 = all(double(subs(pi{i},x,[Vx(1,k);Vy(1,k);0;0])) < 0);
            isect2 = all(double(subs(pi{i},x,[Vx(1,k);Vy(1,k);0;0])) < 0);
            if isect1 || isect2,
                Vx(:,k) = [NaN; NaN];
                Vy(:,k) = [NaN; NaN];
            end
        end
        for j = 1:size(vi{i},1)
            if j == 1, jm = size(vi{i},1); end
            for k = 1:size(Vx,2)
                [isect,tmp1,tmp2,tmp3] = intersectPoint(Vx(1,k),Vy(1,k),Vx(2,k),Vy(2,k),vi{i}(j,1),vi{i}(j,2),vi{i}(jm,1),vi{i}(jm,2));
                if isect,
                    Vx(:,k) = [NaN; NaN];
                    Vy(:,k) = [NaN; NaN];
                end
            end
            jm = j;
        end
    end
end

rmindx = isnan(Vx(i,:));
Vx(:,rmindx) = [];
Vy(:,rmindx) = [];
