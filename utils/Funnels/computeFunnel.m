function [funnel,rhoMin,rho_d,info,redindx] = computeFunnel(f0,t,Xk,Uk,S0,Q,R,H,n,isCyclic,skipFun,regBnd,X,reg1,reg2,Xbound1,Xbound2,Xin1,Xin2,ellBndInv1,ellBndInv2,ellInInv1,ellInInv2)

global x

if sum(isCyclic) > 1, error('number of cyclic dimensions cannot exceed 1'), end
if length(isCyclic) ~= n, error('number of entries in isCylic must be equal to n!'), end

bloatFlag = false;

x = msspoly('x',n);

nn = length(H)+1;
cycIndx = find(isCyclic);

sosProblemType = 'unbounded';
% warning('off','MATLAB:nearlySingularMatrix');
% rmpath('D:\Cornell\Research\CreatePath\ParametricVerification\ellipsoids');
if nargin <= 11
    [regBnd,X,reg1,reg2,Xbound1,Xbound2,Xin1,Xin2,ellBndInv1,ellBndInv2,ellInInv1,ellInInv2] = deal([]);
    sosProblemType = 'unbounded';
elseif nargin <= 14
    [reg2,Xbound1,Xbound2,Xin1,Xin2,ellBndInv1,ellBndInv2,ellInInv1,ellInInv2] = deal([]);
    sosProblemType = 'bounded';
elseif nargin <= 20
    [ellBndInv2,ellInInv1,ellInInv2] = deal([]);
    sosProblemType = 'bounded';
end

Nsteps = size(Xk,1);%2000;
tspan = [0 t(Nsteps)]; % duration of input.

t0s = t(1:Nsteps);
x0s = Xk(1:Nsteps,:)';
xpp = spline(t0s,x0s);
u0s = Uk(1:Nsteps,:)';
upp = spline(t0s,u0s);

if any(isCyclic)
    ithSet = 1:3;
    thSgn = [0 -1 1];
else
    ithSet = 1;
    thSgn = 0;
end

% looping needed for "unwrapping" of theta
for ith = ithSet
    %TODO: generalize the following to arb dimensions
    x0s = Xk(1:Nsteps,:)' + repmat(2*pi*thSgn(ith)*isCyclic',Nsteps,1)';
    xpp1 = spline(t0s,x0s);

    if any(isCyclic)
        % Transform from x-y to x'-y' rotated by initial theta
        invTmotion = [cos(-Xk(1,cycIndx)) -sin(-Xk(1,cycIndx)); sin(-Xk(1,cycIndx)) cos(-Xk(1,cycIndx))];
        T = blkdiag(invTmotion,eye(n-length(H)));
    else
        T = eye(n);
    end
    
    % Find a quadratic lyapunov function based on application of an LQR control
    % law about the trajectory
    [A,B] = tv_poly_linearize_new(f0,@(t) t,@(t) ppval(xpp1,t),@(t) ppval(upp,t));
    
    
    
    tspan = linspace(tspan(1),tspan(2),5000);
    [ts{ith},Ss{ith}] = tv_lqr_riccati_Abias(tspan,A,B,Q,R,S0,T);
    Spp{ith} = spline(ts{ith},Ss{ith});
    S = @(t) ppval(Spp{ith},t);
    K{ith} = @(t) inv(R(t))*B(t)'*S(t);
    Ac{ith} = @(t) A(t) - B(t)*K{ith}(t);
    Q0{ith} = @(t) (Q(t) + S(t)*B(t)*inv(R(t))*B(t)'*S(t));
    clear Ks1
    for ii = 1:length(ts{ith})
        Ks1(:,:,ii) = K{ith}(ts{ith}(ii));
    end
    Kpp{ith} = spline(ts{ith},Ks1);
end

% selecting taus is somewhat arbitrary.. we'll use the fine-grained result produced by tv_lqr_riccati_Abias
%taus = ts{1};
% choose taus to be the time vector from the trajectory generation
taus = flipud(t);

N = length(taus);
    
% interpolate to determine Ps, etc
for i = 1:Nsteps;%length(t)
    if any(isCyclic)
        ith(i) = 1;
        if Xk(i,nn) > pi
            ith(i) = 2;
        elseif Xk(i,nn) < -pi
            ith(i) = 3;
        end
    else
        ith(i) = 1;
    end
        
    indx = find(min(abs(t(i) - ts{ith(i)})) == abs(t(i) - ts{ith(i)}),1,'first');
    Ps{i} = Ss{ith(i)}(:,:,indx);
    Ks{i} = K{ith(i)}(t(i));
    Acs{i} = Ac{ith(i)}(t(i));
    xMu(i,:) = Xk(i,:) + isCyclic'*thSgn(ith(i))*2*pi;
    
    try 
    PinvTmp = inv(Ps{i});
    catch
        keyboard
    end
    PinvTmp = (PinvTmp'+PinvTmp)/2;  % to ensure it is symmetric
%     E1priorInv(i).x = xMu(i,:)';
%     E1priorInv(i).P = Ps{i};
    E1prior(i) = ellipsoid(xMu(i,:)',PinvTmp);
    Eprior(i) = projection(E1prior(i),[H; zeros(n-length(H),length(H))]);
end

% we want to allow theta to wrap so that it stays within limits
% if isCyclic
%     for indx = nn:n
%         % Wrap theta from -pi to pi
%         xup(indx,:) = mod(xup(indx,:)+pi,2*pi)-pi;
%     end
% end

taus = flipud(taus);
xup = ppval(xpp,taus);
for i = 1:N
    if any(isCyclic)
        ith(i) = 1;
        if xup(nn,i) > pi
            ith(i) = 2;
        elseif xup(nn,i) < -pi
            ith(i) = 3;
        end
    else
        ith(i) = 1;
    end
    Pup(:,:,i) = ppval(Spp{ith(i)},taus(i));
    Kup(:,:,i) = ppval(Kpp{ith(i)},taus(i));
end

% taus = taus(1:Nsteps);

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


%% Find certificates respecting system and workspace invariants

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


