

global e gotopt xyPath closeEnough xyzG t x
xyzG = [-3.2923; -5.0189; -0.2948]';%[1;-1;0];
xyzInit = [-4.6161; -4.6169; 0.7830]';%[0;-2;0;1];  % xyz-- slight abuse of naming convention
% xyzInit = [0.5;0.5];

xyz2Init = [xyzInit(1:2)'; sin(xyzInit(3)); cos(xyzInit(3))];

plt3d = false;

addpath('BarrierCertificates')

Nb = 3;     % Order of the barrier function

Nloc = 3;

testSetup


%%
% loc 1: Nominal system
f0 = @(t,x,u) CreateKinematicsNLPointCtrlAlt(t,x);
% f0 = @(t,x,u) testSys(t,x);
xdot{1} = f0(t,x);

% loc 2: System active when w = 1 & v > 0
f0 = @(t,x,u) CreateKinematicsNLPointCtrlAltWplus(t,x);
xdot{2} = f0(t,x);
% loc 3: System active when w = -1 & v > 0
f0 = @(t,x,u) CreateKinematicsNLPointCtrlAltWminus(t,x);
xdot{3} = f0(t,x);
% % loc 4: System active when v = 0
% f0 = @(t,x,u) CreateKinematicsNLPointCtrlPolyVzero(t,x);
% xdot{4} = f0(t,x);

% Domain set
beta = 30;
xbar = [xyzG(1:2)';sin(xyzG(3));cos(xyzG(3))];  % ball centered at xyGoal and theta = 0
[X{1:3}] = deal((x-xbar)'*eye(n)*(x-xbar)/beta - 1);

% Initial set (X0 negative)
beta0 = 0.02;  % Hopefully small enough that non-overlap of the hybrid system modes is not a problem
xbar = [xyzInit(1:2)';sin(xyzInit(3));cos(xyzInit(3))];  % ball centered at xyInit
% determine which mode the initial set is in 
Vx = xyzG(1) - xbar(1);
Vy = xyzG(2) - xbar(2);
[X0{1:3}] = deal([]);
if -(1/e*(-Vx*xbar(3) + Vy*xbar(4)) - wlim) < 0
    X0{2} = (x-xbar)'*eye(n)*(x-xbar)/beta0 - 1;
elseif (1/e*(-Vx*xbar(3) + Vy*xbar(4)) - wlim) < 0
    X0{3} = (x-xbar)'*eye(n)*(x-xbar)/beta0 - 1;
else
    X0{1} = (x-xbar)'*eye(n)*(x-xbar)/beta0 - 1;
end

% Unsafe set (Xu negative)
betau = 0.1;
xbar = [2;0;0;1];

[Xui{1}{1:3}] = deal(reg.int);
[Xue{1:3}] = deal(reg.ext);
% [Xue{1:3}] = deal([]);

% Guard, invariant, reset for w = +/-1, v = 0
xd = xyzG(1);
yd = xyzG(2);
Vx = xd - x(1);
Vy = yd - x(2);
Xg{1}(1) = -(1/e*(-Vx*x(3) + Vy*x(4)) - wlim) + eps*x(1)+eps*x(2);    % guard for w becoming = 1
Xg{1}(2) = (1/e*(-Vx*x(3) + Vy*x(4)) + wlim) + eps*x(1)+eps*x(2);   % guard for w becoming = -1
% Xg{1}(3) = (Vx*cosapprox(x(3),n) + Vy*sinapprox(x(3),n));   % guard for v becoming = 0
Xg{2}(1) = (1/e*(-Vx*x(3) + Vy*x(4)) - wlim) + eps*x(1)+eps*x(2);    % guard for w becoming < 1
Xg{3}(1) = -(1/e*(-Vx*x(3) + Vy*x(4)) + wlim) + eps*x(1)+eps*x(2);   % guard for w becoming > -1
% Xg{4}(1) = -(Vx*cosapprox(x(3),n) + Vy*sinapprox(x(3),n));   % guard for v becoming > 0
Idx{1} = [2 3];   % ordered index set of locations with incoming transition to location 1
Idx{2} = [1];       % ordered index set of locations with incoming transition to location 2
Idx{3} = [1];       % ordered index set of locations with incoming transition to location 3
% Idx{4} = [1];       % ordered index set of locations with incoming transition to location 4

% Set up a grid whose points will be minimized to approximately minimize the barrier
[tmp1,tmp2,tmp3] = meshgrid(linspace(min(reg.extVert(1,:)),max(reg.extVert(1,:)),10),linspace(min(reg.extVert(2,:)),max(reg.extVert(2,:)),10),-pi:pi/5:pi);
xyzMinGrid = [tmp1(:) tmp2(:) tmp3(:)];
elim = [];
for i = 1:size(xyzMinGrid,1)
    if ~isempty(X{1})
        if double(subs(X{1},x,[xyzMinGrid(i,1:2)';sin(xyzMinGrid(i,3));cos(xyzMinGrid(i,3))])) > 0
            elim = [elim; i];
        end
    end
    if ~isempty(reg.int)
        if all(double(subs(reg.int,x,[xyzMinGrid(i,1:2)';sin(xyzMinGrid(i,3));cos(xyzMinGrid(i,3))])) < 0)
            elim = [elim; i];
        end
    end
    if any(double(subs(reg.ext,x,[xyzMinGrid(i,1:2)';sin(xyzMinGrid(i,3));cos(xyzMinGrid(i,3))])) < 0)
        elim = [elim; i];
    end
end
xyzMinGrid(unique(elim),:) = [];

% Compute the barrier function
tic
[L, barrier, prog] = findBarrierHybridIterate(t, x, xdot, X, X0, Xui, Xue, Xg, Idx, Nb, Nloc, xyzMinGrid, xyzInit, xyzG);
% [L{1}, barrier{1}] = findBarrier(t, x, xdot{1}, X{1}, X0{1}, Xu{1}, Nb);

toc

disp('Barrier found. Plotting...')

%% Plotting

triangleRegion
% triangleRegion1

if n == 2
    plot2states
elseif n == 3,
    plot3states
elseif n == 4
    plot4states
end


