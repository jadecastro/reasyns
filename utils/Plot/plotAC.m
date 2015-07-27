function plotAC(ac,options)

t = getTimeVec(ac.x0);
x = ppval(ac.x0.pp,t);
E = ac.ellipsoid;
E1 = [];
for i = 1:20:length(t)
    [q,Q] = double(E(i));
    E(i) = ellipsoid(q+[0,0,100]',Q);
    E1 = [E1; E(i)];
end
plot(E1,options);
% plot3(x(:,1),x(:,2),x(:,3))

