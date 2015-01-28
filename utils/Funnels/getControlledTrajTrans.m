function [t1,Xk1,t,Xk,funIndx] = getControlledTrajTrans(X0,itrans,kiter,funnel,funIndx,modelType)

global Xk Uk K t Xss Uss Kss

for k = 1:length(funnel{funIndx,itrans,kiter}.t)
    tmp = inv(funnel{funIndx,itrans,kiter}.P(:,:,k));
    tmp = (tmp+tmp')/2;
    E(k) = ellipsoid(funnel{funIndx,itrans,kiter}.x(k,:)',tmp*funnel{funIndx,itrans,kiter}.rho(k));
end
ellTmp1 = E;
Xk = funnel{funIndx,itrans,kiter}.x;
Uk = funnel{funIndx,itrans,kiter}.u;
K = funnel{funIndx,itrans,kiter}.K;
t = funnel{funIndx,itrans,kiter}.t;

% clear nrm
% for i = 1:size(Xk,1), nrm(i) = norm(X0' - Xk(i,:)); end
% startIdx = find(min(nrm) == nrm)
% Xk(1:startIdx-1,:) = [];
% Uk(1:startIdx-1,:) = [];
% K(:,:,1:startIdx-1) = [];
% t(1:startIdx-1) = [];
%         isinternal_quick(ellTmp1(1),X0)

xpoly = spline(t,Xk');
upoly = spline(t,Uk');
Kpoly = spline(t,K);
Xss = @(t) ppval(xpoly,t);
Uss = @(t) ppval(upoly,t);
Kss = @(t) ppval(Kpoly,t);

switch modelType
    case 'unicycle'
        [t1,Xk1] = ode45(@ CreateKinematicsNL, t, X0);
        % Wrap theta from -pi to pi
        Xk1(:,3) = mod(Xk1(:,3)+pi,2*pi)-pi;
    case 'car'
        [t1,Xk1] = ode23s(@ CarKinematicsNL, t, X0);
        % Wrap theta from -pi to pi
        Xk1(:,3) = mod(Xk1(:,3)+pi,2*pi)-pi;
end    

