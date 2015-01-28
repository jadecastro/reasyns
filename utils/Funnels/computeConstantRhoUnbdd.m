function [rhoMin,rho_d,Info] = computeConstantRhoUnbdd(f0,ts,Ac,K,xMu,uk,Eprior,H,n,isCyclic,succRegChk)
% This is a routine to get max rho assuming it is constant throughout the
% trajectory.  

if sum(isCyclic) > 1, error('number of cyclic dimensions cannot exceed 1'), end
if length(isCyclic) ~= n, error('number of entries in isCylic must be equal to n!'), end
if length(ts) ~= size(xMu,1), error('length of t must match the number the rows in x.'), end
if length(ts) ~= size(uk,1), error('length of t must match the number the rows in u.'), end

global t x debugFlg

t = msspoly('t');

rho_d = [];
rho_runningMin = inf;
sclr = 1;%2000;
N1 = 100;
scalar = 1000;

trialMat = repmat([0 2*pi -2*pi],n,1).*repmat(isCyclic,1,3);

% Compute rho based on distance to polytopes
for i = 1:size(xMu,1)
    i
    E = shape(Eprior(i),1/sclr);
%     xbar = EpriorInv(i).x;
%     Qinv = EpriorInv(i).P*sclr;
    [xbar,Q] = double(E);
    V(i) = (x - xbar)'*inv(Q)*(x - xbar);
    Ki = K{i};
    Aci = Ac{i};
    ubar = uk(i,:)';
    
    xdot = f0(t,x,ubar - Ki*(x - xbar)) - f0(t,xbar,ubar);
    % xdot = Aci*(x - xbar);  % <--- for testing
    % xdot = f0(t,x,[2;ubar(2)] - Ki*(x - xbar));  % attempt at treating limits
    
    % [L1,Lui,Lue] = computeL(V(i),xdot,regBnd,reg1);
%     [rho_sos,Info{i},posDefMultiplier] = computeRhoFeas(V(i),xdot);
    [rho_sos,Info{i},posDefMultiplier] = computeRhoFeasBilinearAlternation(V(i),xdot);
    Info{i}
    
    rho = rho_sos;
    rho_d = [rho_d; rho];
    rho_runningMin = min(rho_d);
    
end

% if debugFlg
%     plotV
% end

rho_d = rho_d/sclr;
rho_d
rhoMin = min(rho_d)

if sclr ~= 1
    indx = find(rho_d == rhoMin, 1 );
    E = shape(Eprior(indx),1/sclr);
    rhoMin = rhoMin/sclr
end

