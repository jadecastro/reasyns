function plotFunnelProj(funnel,Hout,options)

n = length(funnel.x(1,:));
for k = 1:length(funnel.t)
    tmp = inv(funnel.P(:,:,k));
    tmp = (tmp+tmp')/2;
    E1(k) = ellipsoid(funnel.x(k,:)',tmp*funnel.rho(k));
    E(k) = projection(E1(k),[Hout; zeros(n-length(Hout),length(Hout))]);
    plotEllipse1(E(k),options.color,'2d');
    plotEllipse1(E(k),options.color,'2d_outline');
end
% plot(E,options);

%plot3(funnel.x(:,1),funnel.x(:,2),funnel.x(:,3))

