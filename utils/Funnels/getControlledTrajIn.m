function [t1,Xk1,t,Xk,funIndx] = getControlledTrajIn(X0,iReg,kiter,funnelIn,funIndx,modelType)

global Xk Uk K t Xss Uss Kss

for k = 1:length(funnelIn{funIndx,iReg,kiter}.t)
    tmp = inv(funnelIn{funIndx,iReg,kiter}.P(:,:,k));
    tmp = (tmp+tmp')/2;
    E(k) = ellipsoid(funnelIn{funIndx,iReg,kiter}.x(k,:)',tmp*funnelIn{funIndx,iReg,kiter}.rho(k));
end
ellTmp1 = E;
Xk = funnelIn{funIndx,iReg,kiter}.x;
Uk = funnelIn{funIndx,iReg,kiter}.u;
K = funnelIn{funIndx,iReg,kiter}.K;
t = funnelIn{funIndx,iReg,kiter}.t;

clear nrm
for i = 1:size(Xk,1), nrm(i) = norm(X0' - Xk(i,:)); end
startIdx = min(find(min(nrm) == nrm), length(nrm)-1)
Xk(1:startIdx-1,:) = [];
Uk(1:startIdx-1,:) = [];
K(:,:,1:startIdx-1) = [];
t(1:startIdx-1) = [];
%         isinternal_quick(ellTmp1(1),X0)

xpoly = spline(t,Xk');
upoly = spline(t,Uk');
Kpoly = spline(t,K);
Xss = @(t) ppval(xpoly,t);
Uss = @(t) ppval(upoly,t);
Kss = @(t) ppval(Kpoly,t);

switch modelType
    case 'unicycle'
        [t1,Xk1] = ode45(@ CreateKinematicsNL, [0 t(end)], X0);
        % Wrap theta from -pi to pi
        Xk1(:,3) = mod(Xk1(:,3)+pi,2*pi)-pi;
    case 'car'
        [t1,Xk1] = ode23s(@ CarKinematicsNL, [0 t(end)], X0);
        % Wrap theta from -pi to pi
        Xk1(:,3) = mod(Xk1(:,3)+pi,2*pi)-pi;
end    

