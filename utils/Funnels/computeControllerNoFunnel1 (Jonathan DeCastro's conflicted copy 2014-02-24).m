function [funnel] = computeControllerNoFunnel1(t,Xk,Uk,S0,Q,R,H,n,isCyclic)

global x

if sum(isCyclic) > 1, error('number of cyclic dimensions cannot exceed 1'), end
if length(isCyclic) ~= n, error('number of entries in isCylic must be equal to n!'), end

x = msspoly('x',n);

nn = length(H)+1;
cycIndx = find(isCyclic);

% warning('off','MATLAB:nearlySingularMatrix');
% rmpath('D:\Cornell\Research\CreatePath\ParametricVerification\ellipsoids');
if nargin <= 13
    [reg2,Xbound1,Xbound2,Xin1,Xin2,ellBnd1,ellBnd2,ellIn1,ellIn2] = deal([]);
elseif nargin <= 19
    [ellBnd2,ellIn1,ellIn2] = deal([]);
end

f0 = @(t,x,u) CreateKinematicsPoly(t,x,u); % Populate this with your dynamics.

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
    [A,B] = tv_poly_linearize(f0,@(t) ppval(xpp1,t),@(t) ppval(upp,t));
    
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

% selecting taus is somewhat arbitrary.. we'll use the fine-grained result produced by the linearization
taus = ts{1};
N = length(taus);
    
% interpolate to determine Ps, etc
for i = 1:Nsteps;%length(t)
    if Xk(i,nn) > pi && any(isCyclic)
        ith(i) = 2;
    elseif Xk(i,nn) < -pi && any(isCyclic)
        ith(i) = 3;
    else
        ith(i) = 1;
    end
        
    indx = min(abs(t(i) - ts{ith(i)})) == abs(t(i) - ts{ith(i)});
    Ps{i} = Ss{ith(i)}(:,:,indx);
    Ks{i} = K{ith(i)}(t(i));
    Acs{i} = Ac{ith(i)}(t(i));
    xMu(i,:) = Xk(i,:) + [zeros(1,length(H)) thSgn(ith(i))*2*pi];
    
    PinvTmp = inv(Ps{i});
    PinvTmp = (PinvTmp'+PinvTmp)/2;  % to ensure it is symmetric
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
    if xup(nn,i) > pi && any(isCyclic)
        ith(i) = 2;
    elseif xup(nn,i) < -pi && any(isCyclic)
        ith(i) = 3;
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


uup = ppval(upp,taus);

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

