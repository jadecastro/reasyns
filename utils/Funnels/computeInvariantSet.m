function [sys,acNom] = computeInvariantSet(sys,acNom,skipFun,indices,map,acBnd,acIn)

global x

if sum(sys.isCyclic) > 1, error('number of cyclic dimensions cannot exceed 1'), end
if length(sys.isCyclic) ~= sys.n, error('number of entries in isCylic must be equal to n!'), end

bloatFlag = false;

x = msspoly('x',sys.n);

nn = length(sys.H)+1; % ??
cycIndx = find(sys.isCyclic);

sosProblemType = 'unbounded';
% warning('off','MATLAB:nearlySingularMatrix');

if nargin <= 3
    [indices,map,acBnd,acIn] = deal([]);
    sosProblemType = 'unbounded';
else
    if length(indices) > 2
        error('Only single states or pairs of states supported.')
    end
    sosProblemType = 'bounded';
end

Nsteps = size(acNom.x,1);%2000;
tspan = [0 acNom.t(Nsteps)]; % duration of input.

t0s = acNom.t(1:Nsteps);
x0s = acNom(1:Nsteps,:)';
xpp = spline(t0s,x0s);
u0s = acNom.u(1:Nsteps,:)';
upp = spline(t0s,u0s);

if any(sys.isCyclic)
    ithSet = 1:3;
    thSgn = [0 -1 1];
else
    ithSet = 1;
    thSgn = 0;
end

%% Find certificates respecting the workspace invariants

[tRed,redindx] = downsampleUniformly(t,skipFun);
xRed = xMu(redindx,:);
Ured = Uk(redindx,:);

E1priorInvRed = []; E1priorRed = []; EpriorRed = [];
count = 0;
for i = redindx
    count = count + 1;
%     E1priorInvRed = [E1priorInvRed E1priorInv(i)];
    E1priorRed = [E1priorRed E1prior(i)];
    EpriorRed = [EpriorRed Eprior(i)];
    AcsRed{count} = Acs{i};
    KsRed{count} = Ks{i};
end

switch sosProblemType
    %     case 'bounded'
    %         [rhoMin,rho_d,info] = computeConstantRho3(f0,tRed,AcsRed,KsRed,xRed,Ured,regBnd,X,E1priorRed,[],reg1,Xbound1,Xin1,[],reg2,Xbound2,Xin2,ellBndInv1,ellBndInv2,ellInInv1,ellInInv2,H,n,isCyclic,[],bloatFlag);
    case 'unbounded'
        %         [rhoMin,rho_d,info] = computeConstantRhoUnbdd(f0,tRed,AcsRed,KsRed,xRed,Ured,E1priorRed,H,n,isCyclic);
        X = [];
end
% [rhoMin,rho_d,info] = computeConstantRho3(f0,tRed,AcsRed,KsRed,xRed,Ured,regBnd,X,E1priorRed,[],reg1,Xbound1,Xin1,[],reg2,Xbound2,Xin2,ellBndInv1,ellBndInv2,ellInInv1,ellInInv2,H,n,isCyclic,[],bloatFlag);
[rhoMin,rho_d,info] = computeConstantRho3WholeTraj(f0,tRed,AcsRed,KsRed,xRed,Ured,regBnd,X,E1priorRed,[],reg1,Xbound1,Xin1,[],reg2,Xbound2,Xin2,ellBndInv1,ellBndInv2,ellInInv1,ellInInv2,H,n,isCyclic,[],bloatFlag);
% [rhoMin,rho_d,info] = computeConstantRhoWholetraj(f0,tRed,AcsRed,KsRed,xRed,Ured,regBnd,X,E1priorRed,H,n,isCyclic);  
rho_d

% rhot = rhoMin*ones(size(taus));
% rhopp = interp1(taus,rhot,'linear','pp');
rhot = rho_d;
rhopp = interp1(tRed,rhot,'linear','pp');

% Upsample the ellipsoids
% taus = t;
uup = ppval(upp,taus);
rhoup = ppval(rhopp,taus);

funnel.t = taus;
funnel.x = xup';
funnel.u = uup';
for i = 1:length(taus)
    PupInvTmp = inv(Pup(:,:,i));
    PupInvTmp = (PupInvTmp'+PupInvTmp)/2;  % to ensure it is symmetric
    Pup(:,:,i) = inv(PupInvTmp);
end
funnel.P = Pup;
funnel.K = Kup;
funnel.rho = rhoup;
for i = 1:length(taus)
    V(i) = (x - xup(:,i))'*Pup(:,:,i)/rhoup(i)*(x - xup(:,i));
end
funnel.V = V;


