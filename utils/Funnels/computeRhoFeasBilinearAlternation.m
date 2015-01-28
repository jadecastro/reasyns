function [rho,info,L1] = computeRhoFeasBilinearAlternation(V,xdot)

global t x

Vdot = diff(V,t) + diff(V,x)*xdot;

% % Compute initial rho
% prog = mssprog;
% rho = msspoly('rho',1);
% prog.psd = rho;
% 
% Lxmonom = monomials(x,0:deg(Vdot,x)-deg(V,x));  % order is a design choice
% [prog,l] = new(prog,length(Lxmonom),'free');
% L1 = l'*Lxmonom;
% prog.sos = L1;
% prog.sos = -(Vdot + L1*(rho - V));
% 
% opt = 0;  % +1: minimize the cost function; -1: maximize the cost function; 0: solve a feasibility problem (no cost function)
%     
% pars.fid = 0; % 0 to suppress screen output; 1 to show screen output
% [prog,info] = sedumi(prog,rho,0,pars,opt);
% if info.pinf
%     disp('fail.')
%     rho = 1e9;
% else
%     rho = double(prog([rho]))
% end

% Instantiate rho: start small
rho = 0.01;

if true % do only if previous program was successful
    for i = 1:100
        lastRho = rho;
        
        % Find a feasible multiplier
        prog = mssprog;
        gamma = msspoly('g',1);
        prog.free = gamma;
        
        Lxmonom = monomials(x,0:deg(Vdot,x)-deg(V,x));  % order is a design choice
        [prog,l] = new(prog,length(Lxmonom),'free');
        L1 = l'*Lxmonom;
        prog.sos = L1;
        prog.sos = gamma - (Vdot + L1*(rho - V));
        
        opt = 1;  % +1: minimize the cost function; -1: maximize the cost function; 0: solve a feasibility problem (no cost function)
        
        pars.fid = 0; % 0 to suppress screen output; 1 to show screen output
        [prog,info] = sedumi(prog,gamma,0,pars,opt);
        if info.pinf
            disp('fail.')
            L1 = [];
            rho = 1e9
            break
        else
            L1 = prog([L1]);
        end
        
        % Next, find a feasible rho
        prog = mssprog;
%         rho = msspoly('rho',1);
%         prog.psd = rho;
        [prog,rho] = new(prog,1,'pos');
        prog.sos = -(Vdot + L1*(rho - V));
        
        opt = -1;  % +1: minimize the cost function; -1: maximize the cost function; 0: solve a feasibility problem (no cost function)
        
        pars.fid = 0; % 0 to suppress screen output; 1 to show screen output
        [prog,info] = sedumi(prog,rho,0,pars,opt);
        info
        if info.pinf
            disp('fail.')
            rho = 1e9
        else
            rho = double(prog([rho]))
        end
        
        if abs((rho - lastRho)/rho) < 0.01, break; end  % break if it appears that we are at a fixed point
    end
end

end
