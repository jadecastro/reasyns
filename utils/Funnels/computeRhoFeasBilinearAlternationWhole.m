function [rho,info,L1] = computeRhoFeasBilinearAlternationWhole(V,xdot,tau,x_test)

global t x

N = size(xdot,2);
rho_runningMax = 1e-9*ones(N,1);
for k = 1:N
    Vdot(k) = diff(V(k),t) + diff(V(k),x)*xdot(:,k);
end

% Instantiate rho: start small
rho = 1e-5*ones(N,1).*exp(-(0:N-1))';
rhodot = [0; (rho(2:N)-rho(1:N-1))./tau];  % piecewise "derivative"

for i = 1:20
    lastRho = rho;
    
    % Find a feasible multiplier
    for k = 1:N
        prog = mssprog;
        neggamma = msspoly('g',1);
        prog.psd = neggamma;
        
        Lxmonom = monomials(x,0:1*(deg(Vdot(k),x)-deg(V(k),x)));  % order is a design choice
        [prog,l] = new(prog,length(Lxmonom),'free');
        L1t = l'*Lxmonom;
%         prog.sos = L1t;
        prog.sos = -neggamma - (Vdot(k) - rhodot(k) + L1t*(rho(k) - V(k)));
        
        opt = -1;  % +1: minimize the cost function; -1: maximize the cost function; 0: solve a feasibility problem (no cost function)
        
        pars.fid = 0; % 0 to suppress screen output; 1 to show screen output
        [prog,info] = sedumi(prog,1e5*neggamma,0,pars,opt);
        if info.pinf  % if this step fails, then use the previous iteration's rho or the initial rho
            disp('fail.')
            break
            %                 L1(k) = NaN;
            %                 rho(1:N) = 1e9
        else
            disp(num2str(-double(prog([neggamma]))))
            L1(k) = prog([L1t]);
            
%             % sanity check to see if the solver is behaving...
%             PosSemiDef(k) = -(Vdot(k) - rhodot(k) + L1(k)*(rho(k) - V(k)));
%             PSDtest(k) = double(msubs(PosSemiDef(k),x,x_test(k,:)'));
%             Vtest(k) = double(msubs(V(k),x,x_test(k,:)'));
%             L1test(k) = double(msubs(L1(k),x,x_test(k,:)'));
%             disp([' psd test, k = ',num2str(k),', V(x_test) = ',num2str(Vtest(k)),', L1(k) = ',num2str(L1test(k)),':  ',num2str(PSDtest(k))])
        end
    end
    if info.pinf % if this step fails, then use the previous iteration's rho or the initial rho
        if ~exist('L1'), L1 = []; end
        disp('stopping.')
        break
    end
    
    % Next, find a feasible rho
    prog = mssprog;
    [prog,rhot] = new(prog,N,'pos');
    rhodott = [0; rhot(2:N)-rhot(1:N-1)];  % piecewise "derivative"
    for k=1:N-1, rhodott(k) = rhodott(k)/tau(k); end
    prog.sos = -(Vdot' - rhodott + L1'.*(rhot - V'));
    
    %     prog.psd = -rhodott;  % TODO: fix and check out this case
%     prog.psd = 10-rhot(N);  % try bounding final rho
    
%    objFn = sum(rhot) - 10000*log(1:N)*rhot;
    objFn = sum(rhot) - 1e2*exp(1:N)/exp(N)*rhot;
    opt = -1;  % +1: minimize the cost function; -1: maximize the cost function; 0: solve a feasibility problem (no cost function)
    
    pars.fid = 0; % 0 to suppress screen output; 1 to show screen output
    [prog,info] = sedumi(prog,objFn,0,pars,opt);
    info
    if info.pinf % if this step fails, then use the previous iteration's rho or the initial rho
        disp('fail.')
        break
        %             clear rho
        %            rho = 1e9*ones(N,1)
    else
        rho = double(prog([rhot]))
        rhodot = [0; rho(2:N)-rho(1:N-1)];  % piecewise "derivative"
        
%         % sanity check to see if the solver is behaving...
%         PosSemiDef = -(Vdot' - rhodot + L1'.*(rho - V'));
%         for k = 1:N,
%             PSDtest(k) = double(msubs(PosSemiDef(k),x,x_test(k,:)'));
%             disp([' psd test, k = ',num2str(k),':  ',num2str(PSDtest(k))])
%         end
        
        if sum(rho) > sum(rho_runningMax)
            disp('updating rho_runningMax...')
            rho_runningMax = rho;
        end
    end
    
    abs((sum(rho) - sum(lastRho))/sum(rho))
    if abs((sum(rho) - sum(lastRho))/sum(rho)) < 0.01 || max(rho) < 1e-10 || (sum(lastRho)/sum(rho) > 5 && i > 5)
        break; % break if it appears that we are at a fixed point or if we're no longer improving
    end
end

end
