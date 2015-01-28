function [ts,Ps,rhoMin,rho_d] = getLevelSetTraj(t,Xk,Uk,vObs)

% warning('off','MATLAB:nearlySingularMatrix');
% rmpath('D:\Cornell\Research\CreatePath\ParametricVerification\ellipsoids');

n = 3;   % Dimension of state vector

% Generate a funciton f0 which uses simple arithmetic.
f0 = @(t,x,u) CreateKinematicsPoly(t,x,u); % Populate this with your dynamics.

Nsteps = length(Xk);%2000;
reachDir = 'forward';

tspan = [0 t(Nsteps)]; % duration of input.

x00s = Xk(1:Nsteps,:)';
u00s = Uk(1:Nsteps,:)';

t0s = flipud(t(1:Nsteps));
x0s = flipud(x00s);
u0s = flipud(u00s);

xpp = spline(t0s,x0s);
upp = spline(t0s,u0s);
% xpp = spline(t,x00s);

% Find a quadratic lyapunov function based on application of an LQR control
% law about the trajectory
[A,B] = tv_poly_linearize(f0,@(t) ppval(xpp,t),@(t) ppval(upp,t));

S0 = 10*eye(3);   % Qf in terminal cost

% LQR parameters
Q = @(t) diag([1 1 1]);
R = @(t) 1*eye(2);

[ts,Ss] = tv_lqr_riccati_Abias(tspan,A,B,Q,R,S0);
Spp = spline(ts,Ss);
S = @(t) ppval(Spp,t);
K = @(t) inv(R(t))*B(t)'*S(t);
Ac = @(t) A(t) - B(t)*K(t);
Q0  = @(t) (Q(t) + S(t)*B(t)*inv(R(t))*B(t)'*S(t));
taus = ts;

N = length(taus);
Ppp = interp1(taus,reshape(permute(Ss,[3 1 2]),N,n*n),'linear','pp');

% interpolate to determine Ps
for i = 1:Nsteps;%length(t)
    indx = find(min(abs(t(i) - ts)) == abs(t(i) - ts));
    Ps{i} = Ss(:,:,indx);
    xMu(i,:) = Xk(i,:);
end

figure(5), clf
hold on
for i = 1:50:length(xMu)
%     xMu(i,1:2)
    Ptmp = inv(Ps{i});
    Ellipse_plot(inv(Ptmp(1:2,1:2)),xMu(i,1:2));
end

figure(6), clf
hold on
Ptmp = inv(Ps{1});
Ellipse_plot(inv(Ptmp(1:2,1:2)),xMu(1,1:2));
Ptmp = inv(Ps{end});
Ellipse_plot(inv(Ptmp(1:2,1:2)),xMu(end,1:2));


%% Find certificates respecting stability and automaton invariants

% finding an initial rho to get a feasible lagrange multiplier
c = 40;
% rhot = 50 - 49*exp(c*(taus-max(taus))/(max(taus)-min(taus)));
rhot = exp(c*(taus-max(taus))/(max(taus)-min(taus)));
rhopp = interp1(taus,rhot,'linear','pp');

for indx = 1:length(vObs)
    P{indx} = polytope(vObs{indx});
end

rho_d = [];
% % NB: Need to fix this stuff....
% try
%     % Find most permissive rho which certifies invariance about the trajectory
%     for iter=1:1
%         lagrange_multipliers;
%         all([gammas{:}] < 0)
%         if ~all([gammas{:}] < 0),
%             if iter == 1, error('Initial rho(t) infeasible, increase c.');
%             end
%         end
%         %     plot_vanderpol
%         improve_rho;
%         
%         figure,
%         plot(ts,ppval(rhopp,ts),'-o')
%     end
%     % plot_vanderpol
%     
%     rho_d = flipud(ppval(rhopp,ts));
%     
% catch errmsg
%     disp('no feasible solution found for trajectory certificate!')
%     disp('... finding trajectory barrier based on a constant rho -- no guarantees exist!!')

    for i = 1:length(xMu)
        sclr = 2000;
        PinvTmp = inv(Ps{i})/sclr;
        PinvTmp = (PinvTmp'+PinvTmp)/2;  % to ensure it is symmetric
        E1prior = ellipsoid(xMu(i,:)',PinvTmp);
        Eprior = projection(E1,[1 0 ;0 1 ;0 0]);
        
        for j = 1:length(vObs)
            [dist,status,ellPoint,polyPoint] = distanceEllPoly(Eprior,P{j});
%             dist
%             ellPoint
%             polyPoint
            if dist == 0, break, end
            d1 = norm(xMu(i,1:2)' - ellPoint);
            d2 = norm(xMu(i,1:2)' - polyPoint);
            rhoTmp = (d2/d1)^2;
            rho_d = [rho_d; rhoTmp];
        end
        if dist == 0, break, end
    end
    if dist == 0,
        rho_d = NaN;
        rhoMin = 1;
        return
    end
    
    rho_d = rho_d/sclr;
    rhoMin = min(rho_d);
    
    for i = 1:length(xMu)
        % get final ellipsoid
        E1{i} = ellipsoid(Xk(j,:)',tmp*rhoMin);
        
        % Create hyperplanes at each time step to get the inveriant
        xDot(i,:) = Ac*Xk(i,:);
        c = xDot(i,:)*Xk(i,:)';
        H{i} = hyperplane(xDot(i,:)',c);
        Ehull{i} = hpintersection(E1{i},H{i});
    end
    
% end



