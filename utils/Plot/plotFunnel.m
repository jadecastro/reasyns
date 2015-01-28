function plotFunnel(funnel,options)

for k = 1:length(funnel.t)
    tmp = inv(funnel.P(:,:,k));
    tmp = (tmp+tmp')/2;
    E(k) = ellipsoid(funnel.x(k,:)',tmp*funnel.rho(k));
end
plot(E,options);
plot3(funnel.x(:,1),funnel.x(:,2),funnel.x(:,3))

