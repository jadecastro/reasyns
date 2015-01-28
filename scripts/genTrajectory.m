function [t, Xk, Uk] = genTrajectory(initPoint,path,sys,options,tin)

if length(H) > 2, error('Workspaces of dimension greater than two not yet supported. Sorry.'); end

e = options.e;
closeEnough = options.closeEnough;

gotopt
H = sys.H;
l 
TstepTraj

trimTraj = true;
forcedEndIndx = Inf;

%clear xyPath
xyPath = path;

T = TstepTraj;
Tfin = 10;
t = 0:T:Tfin;

if nargin > 6
    t = tin;
end

% Model and Controller
l = 1.0;
x_look = 1;
e = 0.4;
closeEnough = 0.5;
distAccept = 0.1;

%% Kinematics model with transverse control

gotopt = 1;         % Initialize waypoint index
switch modelType
    case 'vdp'
        [t,Xk] = ode45(@ VDPModel, t, initPoint);
    case 'linear'
        [t,Xk] = ode45(@ LinearModel, t, initPoint);
    case 'unicycle'
        [t,Xk] = ode45(@ CreateKinematicsNLWayptCtrl1, t, initPoint);
    case 'unicycle4state'
        [t,Xk] = ode45(@ CreateKinematicsNLWayptCtrlAlt, t, initPoint);
    case 'car'
        [t,Xk] = ode45(@ CarKinematicsNLWayptCtrl1, t, initPoint);
end

for indx = 1:n
    if isCyclic(indx)
        % Wrap theta from -pi to pi
        Xk(:,indx) = mod(Xk(:,indx)+pi,2*pi)-pi;
    end
end

% X = Xk(:,1); Y = Xk(:,2); 
outputk = Xk(:,1:length(H))*H;
nonRegStatek = Xk(:,length(H)+1:end);

% Remove data beyond when the final waypoint has been achieved
% TODO: Instead of acceptance radius, make this dependent on the prescribed sets in the algorithm!
% dist2LastPt = sqrt((X - xyPath(end,1)).^2 + (Y - xyPath(end,2)).^2);
if strmatch(modelType,'linear') 
    Hdist = [0 1];
    dist2LastPt = sqrt((sum(repmat(Hdist,size(outputk,1),1).*outputk,2) - repmat(xyPath(end,1),size(outputk,1),1)).^2);
elseif strmatch(modelType,'vdp')
    Hdist = [1 0];
    dist2LastPt = sqrt((sum(repmat(Hdist,size(outputk,1),1).*outputk,2) - repmat(xyPath(end,1),size(outputk,1),1)).^2);    
else 
    dist2LastPt = sqrt(sum((outputk - repmat(xyPath(end,:),size(outputk,1),1)).^2,2));
end
%indx = find(abs(dist2LastPt) < abs(distAccept),1,'first');
dist2LastPt(size(Xk,1)+1) = inf;
if trimTraj
    for k = 1:size(Xk,1)
        if (abs(dist2LastPt(k+1)) > abs(dist2LastPt(k)) && dist2LastPt(k) < closeEnough) || k == forcedEndIndx
            indx = k;
            break
        end
    end
else
    indx = size(Xk,1);
    for k = 1:size(Xk,1)
        if abs(dist2LastPt(k)) > 1e6
            indx = k-1;
            break
        end
    end
end
indx

% Determine the control inputs indirectly from the simulation
% for i = 1:size(Xk,1)-1
%     regStateDot(i,:) = 1/(t(i+1)-t(i))*(regState(i+1,:) - regState(i,:));
%     nonRegStateDot(i,:) = 1/(t(i+1)-t(i))*(nonRegState(i+1,:) - nonRegState(i,:));
% end
% uV = sqrt(Xdot.^2 + Ydot.^2);
% uW = Thdot;
% Uk = [0 0; uV' uW'];

gotopt = 1;         % Initialize waypoint index
for k = 1:1:size(Xk,1)
    switch modelType
        case 'vdp'
            vy = xyPath(2,1) - Xk(k,1);
            U = [0 vy];
        case 'linear'
            vy = xyPath(2,1) - Xk(k,1);
            U = [0 vy];
        case 'unicycle'
            [U,gotopt] = CreateKinematicsNLWayptCtrl1_ctrl(Xk(k,:),gotopt,e,xyPath,closeEnough);
        case 'unicycle4state'
            [U,gotopt] = CreateKinematicsNLWayptCtrl1_ctrl(Xk(k,:),gotopt,e,xyPath,closeEnough);
        case 'car'
            [U,gotopt] = CarKinematicsNLWayptCtrl1_ctrl(Xk(k,:),gotopt,l,e,xyPath,closeEnough);
    end
    Uk(k,:) = U;
end

if ~isempty(indx)
    t(indx+1:end) = [];
%     X(indx+1:end) = [];
%     Y(indx+1:end) = [];
%     Th(indx+1:end) = [];
    Xk(indx+1:end,:) = [];
    Uk(indx+1:end,:) = [];
end

% Uk = [0; uW'];

%% Plotting

% figure(1), clf
% 
% Xlimit = max(X(end)+0.5, 2*Y(1)+0.5);
% Ylimit = max(X(end)/2+0.5, Y(1)+0.5);
% 
% Xvals = [-0.5 Xlimit];
% 
% % plot(Xvals,[0 0],'k--','LineWidth',2);
% 
% %axis([-0.5 Xlimit -Ylimit Ylimit])
% axis equal
% hold on
% 
% Xinc = downsample(X,ceil(length(X)/10));
% Yinc = downsample(Y,ceil(length(X)/10));
% Thinc = downsample(Th,ceil(length(X)/10));
% ang = linspace(-pi,pi,100);
% P = r*eye(2);
% for i = 1:length(Xinc)
%     circ = P*[cos(ang); sin(ang)] + [Xinc(i); Yinc(i)]*ones(1,length(ang));
%     fill(circ(1,:),circ(2,:),[0.8 0.8 0.2])
%     H = line([Xinc(i) Xinc(i) + r*cos(Thinc(i))], [Yinc(i);Yinc(i) + r*sin(Thinc(i))]);
%     set(H,'LineWidth',3,'Color',[0.8 0 0])
% end
% 
% plot(X,Y,'LineWidth',2);
% plot(Xpoly,Ypoly,'c','LineWidth',2);
% 
% figure(2)
% subplot(211)
% plot(t,X)
% subplot(212)
% plot(t,Y)

