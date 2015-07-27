function d = dist2hyperplane(p,H)
% computes the distance from an arbirary point to a hyperplane (specified as a
% hyperplane object)

[h,k] = double(H);

if length(h) ~= length(p)
    error('the point and hyperplane dimensions must agree!')
end

if all(size(p) ~= [length(p) 1])
    p = p';
end

xp = k - h'*p;
d = sqrt(sum((xp.*h).^2));
