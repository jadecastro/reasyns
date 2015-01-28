function [funnel,rhoMin,rho_d,info] = computeROA(t,x0,u0,regBnd,X,reg1,Xbound,Xin)

global x

Nb = 3;     % Order of the lyapunov function

% warning('off','MATLAB:nearlySingularMatrix');
% rmpath('D:\Cornell\Research\CreatePath\ParametricVerification\ellipsoids');

n = 3;   % Dimension of state vector

f0 = @(t,x,u) CreateKinematicsPoly(t,x,u); % Populate this with your dynamics.

% FOR 3-STATE DIFF. DRIVE ROBOT ONLY!!
T = [cos(-Xk(1,3)) -sin(-Xk(1,3)); sin(-Xk(1,3)) cos(-Xk(1,3))];  % transform from x-y to x'-y' rotated by initial theta
for i = 1:Nsteps;
    x0s_trans(:,i) = [T*Xk(i,1:2)' ; (Xk(i,3)-Xk(1,3))];
end

AB = double(subs(diff(msspoly(f(x,u)),[x;u]),[x;u],[x0;u0]));
A = AB(1:n,1:n);
B = AB(1:n,n+(1:m));

% LQR parameters
Q = 1e0*diag([100 100 0.0001]);
% R = 1e2*eye(2);
R = diag([1e4 1e2]);

[K,S] = lqr(A,B,Q,R);

xdot = f0(t,x,ubar - K*(x - xbar));

% figure(5), clf
% hold on
% for i = 1:50:length(xMu)
%     xMu(i,1:2)
%     Ptmp = inv(Ps{i});
%     Ellipse_plot(inv(Ptmp(1:2,1:2)),xMu(i,1:2));
% end
% figure(6), clf
% hold on
% Ptmp = inv(Ps{1});
% Ellipse_plot(inv(Ptmp(1:2,1:2)),xMu(1,1:2));
% Ptmp = inv(Ps{end});
% Ellipse_plot(inv(Ptmp(1:2,1:2)),xMu(end,1:2));

%% Find certificates respecting stability and automaton invariants

% Set up a grid whose points will be minimized to approximately minimize
% the lyapunov function
[tmp1,tmp2,tmp3] = meshgrid(linspace(min(X.extVert(1,:)),max(X.extVert(1,:)),10),linspace(min(X.extVert(2,:)),max(X.extVert(2,:)),10),-pi:pi/5:pi);
xMinGrid = [tmp1(:) tmp2(:) tmp3(:)];
elim = [];
for i = 1:size(xMinGrid,1)
    if ~isempty(X.ext)
        if double(subs(X.ext,x,xMinGrid(i,1:3)')) > 0
%         if double(subs(X.ext,x,[xyzMinGrid(i,1:2)';sin(xyzMinGrid(i,3));cos(xyzMinGrid(i,3))])) > 0
            elim = [elim; i];
        end
    end
    if ~isempty(X.int)
        if all(double(subs(X.int,x,xMinGrid(i,1:3)')) < 0)
%         if all(double(subs(X.int,x,[xyzMinGrid(i,1:2)';sin(xyzMinGrid(i,3));cos(xyzMinGrid(i,3))])) < 0)
            elim = [elim; i];
        end
    end
    if any(double(subs(X.ext,x,xMinGrid(i,1:3)')) < 0)
%     if any(double(subs(X.ext,x,[xyzMinGrid(i,1:2)';sin(xyzMinGrid(i,3));cos(xyzMinGrid(i,3))])) < 0)
        elim = [elim; i];
    end
end
xMinGrid(unique(elim),:) = [];

[lyapfun,info] = computeVasROA(x,xdot,X,Xbound,Xin,xMinGrid,Nb);
[rhoMin,rho_d,info] = computeConstantRho1(tRed,Ac,K,xRed,Ured,regBnd,X,E1prior,[],reg1,Xbound1,Xin1,[],reg2,Xbound2,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2);
rho_d

rhot = rhoMin*ones(size(taus));
rhopp = interp1(taus,rhot,'linear','pp');

% Upsample the ellipsoids
ts = flipud(ts);
% ts = t;
xup = ppval(xpp,ts);
uup = ppval(upp,ts);
Pup = ppval(Spp,ts);
Kup = ppval(Kpp,ts);
rhoup = ppval(rhopp,ts);

funnel.t = ts;
funnel.x = xup';
funnel.u = uup';
for i = 1:length(ts)
    PupInvTmp = inv(Pup(:,:,i));
    PupInvTmp = (PupInvTmp'+PupInvTmp)/2;  % to ensure it is symmetric
    Pup(:,:,i) = inv(PupInvTmp);
end
funnel.P = Pup;
funnel.K = Kup;
funnel.rho = rhoup;
for i = 1:length(ts)
    V(i) = (x - xup(:,i))'*Pup(:,:,i)/rhoup(i)*(x - xup(:,i));
end
funnel.V = V;


