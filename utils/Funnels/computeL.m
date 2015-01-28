function [L1,Lui,Lue] = computeL(V,xdot,regBnd,reg)

global t x

prog = mssprog;

gamma = msspoly('g',1);
prog.free = gamma;

% rho = msspoly('rho',1);
% prog.free = rho;
rho = 1e-7;

Vdot = diff(V,t) + diff(V,x)*xdot;
Lxmonom = monomials(x,0:deg(Vdot,x)-deg(V,x));
[prog,l] = new(prog,length(Lxmonom),'free');
L1 = l'*Lxmonom;
% prog.sos = L1;
%     prog.sos = -(Vdot - rhodot + L1*(V - rho) - gamma);
prog.sos = gamma-(Vdot + L1*(rho - V));

Lui = [];
%     for k = 1:length(reg)
%         if ~isempty(reg{k})
%             % Unsafe set constraints for interior sets
%             for i = 1:length(reg{k})
%                 degL = 2;
%                 Lxmonom = monomials(x,0:degL);
%                 [prog,lu] = new(prog,length(Lxmonom),'free');
%                 Lui(i) = lu'*Lxmonom;
%                 prog.sos = Lui(i);
%             end
%             prog.sos = (V - rho) + Lui*reg{k}.mss';
%         end
%     end

% Lue = [];
% Unsafe set constraints for exterior sets
for k = 1:length(regBnd)
    if ~isempty(regBnd{k}.mss)
        for i = 1:length(regBnd{k}.mss)
            degL = 1;
            Lxmonom = monomials(x,0:degL);
            [prog,lu] = new(prog,length(Lxmonom),'free');
            Lue(i) = lu'*Lxmonom;
%             prog.sos = Lue(i);
            prog.sos = (V - rho) + Lue(i)*regBnd{k}.mss(i);
        end
    end
end

prog = sedumi(prog,gamma,0,[],1);  

%     gamma = double(prog(gamma));
L1 = prog([L1]);
Lue = prog([Lue]);

end
