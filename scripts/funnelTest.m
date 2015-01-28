
global Xk Uk K t l Xss Uss Kss
global mu

format compact

% Map
threeRegionMap

% Downsampling
sampSkipFun = 100;  % skipped samples in funnel computation

% Model parameters
l = 1.0;

% Model parameters
Hout = eye(2);

modelType = 'unicycle';

switch modelType

    case 'linear'
        f0 = @(t,x,u) LinearModelWithInput(t,x,u);
        fLQR = [];  
        n = 2;
        n1 = n;
        Hout = diag([1 1]);
        
        limsNonRegState = [];
        isCyclic = [0 0]';
        
        % LQR controller parameters
        S0 = NaN*eye(n);%1e4*eye(n);   % Qf in terminal cost; enter NaN for infinite-horizon LQR
        Q = @(t) 1e2*diag([1 1]);
        %R = @(t) 1e2*diag([1 1]);
        R = @(t) 1e2*diag([1 1]);
        
        initState = [0;0];
        
        initState1 = initState;
        % initial valus of test traj
%         X0 = initState1;
        X0 = initState1 + 1*[10 0]';
        
    case 'vdp'
        mu = -2;
        
        f0 = @(t,x,u) VDPModelWithInput(t,x,u);
        fLQR = [];  
        n = 2;
        n1 = n;
        Hout = diag([1 1]);
        
        limsNonRegState = [];
        isCyclic = [0 0]';
        
        % LQR controller parameters
        S0 = NaN*eye(n);%1e4*eye(n);   % Qf in terminal cost; enter NaN for infinite-horizon LQR
        Q = @(t) 1e2*diag([1 1]);
        %R = @(t) 1e2*diag([1 1]);
        R = @(t) 1e2*diag([1 1]);
        
        initState = [0;2];
        
        initState1 = initState;
        % initial valus of test traj
%         X0 = initState1;
        X0 = [0 3.2]';
        
    case 'unicycle'
        f0 = @(t,x,u) CreateKinematicsPoly(t,x,u);
        n = 3;
        n1 = n;

        limsNonRegState = [-pi; pi];
        isCyclic = [0 0 1]';
        
        % LQR controller parameters
        S0 = 1e-4*eye(n);   % Qf in terminal cost
        Q = @(t) 1e-4*diag([1 1 1]);
        % R = @(t) 1e2*eye(2);
        R = @(t) diag([1e-5 1e-6]);% @(t) diag([1e4 1e-6]);
        
        initState = [4;4;0];
%         initState = [0;0;pi/2];
%         initState = [0;0;pi];
        
        initState1 = initState;
        % initial valus of test traj
%         X0 = initState1;
        X0 = initState1 + 1*[1 1 pi/2]';
%         X0 = initState1 + -1*[-1 1 1]';
%         X0 = initState1 + -1*[-1 -1 1]';
                
    case 'unicycle4state'
        f0 = @(t,x,u) CreateKinematicsPolyAlt(t,x,u);
        n = 4;
        n1 = n;

        %limsNonRegState = [-pi; pi];
        isCyclic = [0 0 0 0]';
        
        % LQR controller parameters
        S0 = NaN*eye(n);%1e4*eye(n);   % Qf in terminal cost
        Q = @(t) 1e-4*diag([1 1 1 1]);
        % R = @(t) 1e2*eye(2);
        R = @(t) diag([1e-4 1e-6]);% @(t) diag([1e4 1e-6]);
        
        initState = [0;0;0;0];
%         initState = [0;0;pi/2];
%         initState = [0;0;pi];
        
        initState1 = initState;
        % initial valus of test traj
%         X0 = initState1;
        X0 = initState1 + 1*[1 1 1 1]';
%         X0 = initState1 + -1*[-1 1 1]';
%         X0 = initState1 + -1*[-1 -1 1]';
        
    case 'car'
        f0 = @(t,x,u) CarKinematicsPoly(t,x,u);
        n = 4;
        n1 = n;

        limsNonRegState = [-pi -pi/2; pi pi/2]; % for a car-like robot
        isCyclic = [0 0 1 0]';
        
        % LQR controller parameters
        S0 = 1e4*eye(n);   % Qf in terminal cost
        Q = @(t) 1e0*diag([1000 1000 10 10]);
        % R = @(t) 1e2*eye(2);
        R = @(t) diag([1e3 1e2]);
        
        initState = [0;0;0;0];
        initState1 = initState;
        X0 = initState1 + [0.1 0.1 0.1 0.1]';
end

isCyclic1 = isCyclic;
cycIndx = find(isCyclic);

% if n == 4 && ~any(isCyclic)
%     tmp = initState(3);
%     initState1(3) = sin(tmp);
%     initState1(4) = cos(tmp);
% else
%     initState1 = initState;
% end

initOutput = initState(1:2);

% 2-3 state models
%goalState = [0;0];  
% goalState = initState + [-5;1;3];  % "west"
% goalState = initState + [-1;-3;3];  % "south"
goalState = initState + [10;-1;3];  % "east"

% 4 state models:
%goalState = [-1;-10;0;0];  % "south"
goalState1 = goalState;
% if n == 4 && ~any(isCyclic)
%     tmp = goalState(3);
%     goalState1(3) = sin(tmp);
%     goalState1(4) = cos(tmp);
% else
%     goalState1 = goalState;
% end
goalOutput = goalState(1:2);

[t,Xk,Uk] = genNominalTrajectory(initState,[initOutput';goalOutput'],Hout,n1,isCyclic1,modelType);

% kludge to get Uk to behave for us
% Uk(abs(Uk(:,1))<0.1,1) = 0.1;

if any(isCyclic)
    % Transform from x-y to x'-y' rotated by initial theta
    invTmotion = [cos(-Xk(1,cycIndx)) -sin(-Xk(1,cycIndx)); sin(-Xk(1,cycIndx)) cos(-Xk(1,cycIndx))];
    T = blkdiag(invTmotion,eye(n-length(Hout)));
else
    T = eye(n);
end

Xkur = Xk;
% rotate the whole trajectory
X0 = T*X0 - [0 0 Xk(1,cycIndx)]';
Xk = Xk*T' - repmat([0 0 Xk(1,cycIndx)],size(Xk,1),1);

figure(10)
clf
hold on
axis equal
plot(Xk(:,1),Xk(:,2),'r',Xkur(:,1),Xkur(:,2),'b')

% if n == 4 && ~any(isCyclic)
%     Xk1 = Xk(:,1:2);
%     Xk1(:,3) = sin(Xk(:,3));
%     Xk1(:,4) = cos(Xk(:,3));
% else
    Xk1 = Xk;
% end

[traject] = computeController(f0,t,Xk1,Uk,S0,Q,R,Hout,n,isCyclic,sampSkipFun);

Xk0 = traject.x;
Uk0 = traject.u;
K0 = traject.K;
t0 = traject.t;
Xk = Xk0; Uk = Uk0; K = K0; t = t0;

%%
xpoly = spline(t,Xk');
upoly = spline(t,Uk');
Kpoly = spline(t,K);
Xss = @(t) ppval(xpoly,t);
Uss = @(t) ppval(upoly,t);
Kss = @(t) ppval(Kpoly,t);

%%
%         isinternal(ellTmp1(1),X0)
switch modelType
    case 'vdp'
        [t3,Xk3] = ode45(@ VDPModelLQR, t, X0);    
    case 'linear'
        [t3,Xk3] = ode45(@ LinearModelLQR, t, X0);
    case 'unicycle'
        [t3,Xk3] = ode45(@ CreateKinematicsNL, t, X0);
    case 'unicycle4state'
        [t3,Xk3] = ode45(@ CreateKinematicsNLAlt, t, X0);
    case 'car'
        [t3,Xk3] = ode23s(@ CarKinematicsNL, t, X0);
end  

figure(10)
plot(Xk3(:,1),Xk3(:,2),'g')
drawnow

%%
[traject,rhoMin,rho_d] = computeFunnel(f0,t,Xk1,Uk,S0,Q,R,Hout,n,isCyclic,sampSkipFun);
% [traject,rhoMin,rho_d] = computeFunnel(f0,t,Xk1,Uk,S0,Q,R,Hout,n,isCyclic,sampSkipFun,regBnd,regX{1},[]);

Xk0 = traject.x;
Uk0 = traject.u;
K0 = traject.K;
t0 = traject.t;

% figure(10)
% plot(Xk(:,1),Xk(:,2),'k')

%%

Xk = Xk0; Uk = Uk0; K = K0; t = t0;

% clear nrm K1 K2
% for i = 1:size(Xk,1), nrm(i) = norm(X0' - Xk(i,:)); end
% startIdx = max(find(min(nrm) == nrm), 1);
% Xk(1:startIdx-1,:) = [];
% Uk(1:startIdx-1,:) = [];
% K(:,:,1:startIdx-1) = [];
% t(1:startIdx-1) = [];
% K1 = squeeze(K(1,:,:))';
% K2 = squeeze(K(2,:,:))';
% 
% % reset time
% t = t - t(1);

xpoly = spline(t,Xk');
upoly = spline(t,Uk');
Kpoly = spline(t,K);
Xss = @(t) ppval(xpoly,t);
Uss = @(t) ppval(upoly,t);
Kss = @(t) ppval(Kpoly,t);

%%
%         isinternal(ellTmp1(1),X0)
switch modelType
    case 'vdp'
        [t2,Xk2] = ode45(@ VDPModelLQR, t, X0);    
    case 'linear'
        [t2,Xk2] = ode45(@ LinearModelLQR, t, X0);
    case 'unicycle'
        [t2,Xk2] = ode45(@ CreateKinematicsNL, t, X0);
    case 'car'
        [t2,Xk2] = ode23s(@ CarKinematicsNL, t, X0);
end    
% sim('simCreateKinematicsNL');
% Wrap theta from -pi to pi
Xk1(:,cycIndx) = mod(Xk1(:,cycIndx)+pi,2*pi)-pi;

%%
if n == 4 && any(isCyclic)
    figure(11)
    clf
    hold on
    plot(t,Xk(:,4))
    plot(t,pi/2*ones(1,length(t)),'r:',t,-pi/2*ones(1,length(t)),'r:')
end

figure(12)
clf
hold on 
axis equal
% plot(Xk1(:,1),Xk1(:,2))
plot(Xk2(:,1),Xk2(:,2),'g')
plot(Xk(:,1),Xk(:,2),'r')
drawnow


ellArrayCurr = [];
for j = 1:length(traject.t)
    Psav(:,:,j) = traject.P(:,:,j);%Ps{j};
    Ksav(:,:,j) = traject.K(:,:,j);%Ks{j};
    tmp = inv(traject.P(:,:,j));
    tmp = (tmp+tmp')/2;
    E1 = ellipsoid(traject.x(j,:)',tmp*traject.rho(j));
    E = projection(E1,[Hout; zeros(n-length(Hout),length(Hout))]);
    ellArrayCurr = [ellArrayCurr; E1];
end

plot(projection(ellArrayCurr,[Hout; zeros(n-length(Hout),length(Hout))]))
drawnow