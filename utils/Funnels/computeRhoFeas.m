function [rho,info,L1] = computeRhoFeas(V,xdot)

global t x

prog = mssprog;
 
rho = msspoly('rho',1);
prog.psd = rho;
% rho = 0.0001;

Vdot = diff(V,t) + diff(V,x)*xdot;
Lxmonom = monomials(x,0:deg(Vdot,x)-deg(V,x));  % order is a design choice
[prog,l] = new(prog,length(Lxmonom),'free');
L1 = l'*Lxmonom;
prog.sos = L1;
prog.sos = -(Vdot + L1*(rho - V));

opt = 0;  % +1: minimize the cost function; -1: maximize the cost function; 0: solve a feasibility problem (no cost function)
    
% NB: usually can get a reasonable solution by attempting both maximization & minimization 
[prog,info] = sedumi(prog,rho,0,[],opt);
if info.pinf
    disp('fail.')
    rho = 1e9;
else
    rho = double(prog([rho]))
end

end
