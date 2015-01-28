


% warning('off','MATLAB:nearlySingularMatrix');
% rmpath('D:\Cornell\Research\CreatePath\ParametricVerification\ellipsoids');

global x_look e

n = 3;   % Dimension of state vector

% Generate a funciton f0 which uses simple arithmetic.

f0 = @(t,x,u) CreateKinematicsPoly(t,x,u); % Populate this with your dynamics.


%%

clear Ps xMu

Nsteps = length(Xk);%2000;
reachDir = 'forward';

tspan = [0 t(Nsteps)]; % duration of input.

% [ts,xs] = ode45(@(t,x)f0(t,x,[ud(t);wd(t);xd(t);yd(t);thd(t)]),fliplr(tspan),xT);

x00s = Xk(1:Nsteps,:)';
u00s = [0 uV; 0 uW];

t0s = flipud(t(1:Nsteps));
x0s = flipud(x00s);
u0s = flipud(u00s);

xpp = spline(t0s,x0s);
upp = spline(t0s,u0s);
% xpp = spline(t,x00s);

[A,B] = tv_poly_linearize(f0,@(t) ppval(xpp,t),@(t) ppval(upp,t));
%A = tv_poly_linearize(f0,@(t) ppval(xpp,t));

n = 3;
m = 2;
% xT = Xk(end,:)';
% u0 = zeros(m,1);
% Q = eye(n);
% R = 10*eye(m);
% [K0,S0,rho0] = ti_poly_lqr_roa_Abias(@(x,u) f0(0,x,u),xT,u0,Q,R);
% S0 = 1.01*S0/rho0;

S0 = 0.001*eye(3);   % Qf in terminal cost

% 
% if strcmp(reachDir,'backward')
%     % S0 = eye(3);%P0{1};
%     S0 = lyap(A(t(Nsteps))'+1e1*eps*[1 0 0;0 0 0;0 0 0],Q(0));
%     % S0 = Ss(:,:,end);
% elseif strcmp(reachDir,'forward') 
%     tspan = fliplr(tspan); % duration of input.
%     %t0s = t(1:Nsteps);
%     %x0s = x00s;
%     S0 = Ss(:,:,end);
%     S0 = lyap(A(t(1))'+1e1*eps*[1 0 0;0 0 0;0 0 0],Q(0));
% end

%[ts,Ss] = ti_lyapunov_test([10 0],A(t(Nsteps))'+1e10*eps*[1 0 0;0 0 0;0 0 0],Q,S0);
%[ts,Ss] = tv_lyapunov_bias(tspan,A,Q,S0);

Q = @(t) diag([1 1 1]);
R = @(t) 10000*eye(2);

[ts,Ss] = tv_lqr_riccati_Abias(tspan,...
                              A,B,...
                              Q,R,S0);
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

% figure(7), clf
% hold on
% for i = 1:5:length(ts)
%     Ellipse_plot(Ss,[0 0 0]);
% end

% Spp = spline(ts,Ss);
% S = @(t) ppval(Spp,t);
% K = @(t) inv(R(t))*B(t)'*S(t);
% Ac = @(t) A(t);
% Q0  = @(t) (Q(t) + S(t)*B(t)*inv(R(t))*B(t)'*S(t));
% 
% [taus,Ps] = tv_lyapunov(tspan,@(t) Ac(t),Q0,S0);

% N = length(taus);
% Ppp = interp1(taus,reshape(permute(Ps,[3 1 2]),N,n*n),'linear', ...
%               'pp');
% upp = spline(taus,u0(taus'));

figure(5), clf
hold on
for i = 1:50:length(xMu)
    xMu(i,1:2)
    Ptmp = inv(Ps{i});
    Ellipse_plot(inv(Ptmp(1:2,1:2)),xMu(i,1:2));
end

figure(6), clf
hold on
Ptmp = inv(Ps{1});
Ellipse_plot(inv(Ptmp(1:2,1:2)),xMu(1,1:2));
Ptmp = inv(Ps{end});
Ellipse_plot(inv(Ptmp(1:2,1:2)),xMu(end,1:2));


%% Find level sets

c = 40;
rhot = 50 - 49*exp(c*(taus-max(taus))/(max(taus)-min(taus)));
rhopp = interp1(taus,rhot,'linear','pp');

% Repeat to improve rho
for iter=1:1
    lagrange_multipliers;
    all([gammas{:}] < 0)
%     if ~all([gammas{:}] < 0), 
%         if iter == 1, error('Initial rho(t) infeasible, increase c.'); 
%         end
%     end
%     plot_vanderpol
    improve_rho;
    
    figure,
    plot(ts,ppval(rhopp,ts),'-o')
end
% plot_vanderpol

