function [rho,info] = computeRhoRelaxedbounds(V,xdot,X,Xbnd,Xin,Xinv)

global t x

prog = mssprog;
 
% gamma = msspoly('g',1);
% prog.free = gamma;

rho = msspoly('rho',1);
prog.psd = rho;
% rho = 0.0001;

Vdot = diff(V,t) + diff(V,x)*xdot;
Lxmonom = monomials(x,0:deg(Vdot,x)-deg(V,x));
[prog,l] = new(prog,length(Lxmonom),'free');
L1 = l'*Lxmonom;
% prog.sos = L1;
%     prog.sos = -(Vdot - rhodot + L1*(V - rho) - gamma);
prog.sos = -(Vdot + L1*(rho - V));

if nargin == 3
    
    opt = 1;
    
    % Invariance for exterior sets
    if ~isempty(X.mssExtB)
        clear Lue
        for i = 1:length(X.mssIntB)
            Lxmonom = monomials(x,0:deg(Vdot,x)-deg(V,x));
            [prog,l] = new(prog,length(Lxmonom),'free');
            L1{i} = l'*Lxmonom;
            
            degL = 1;
            Lxmonom = monomials(x,0:degL);
            [prog,lu] = new(prog,length(Lxmonom),'free');
            Lui(i) = lu'*Lxmonom;
            
            for j = 1:length(X.mssExtB)
                Lxmonom = monomials(x,0:degL);
                [prog,lu] = new(prog,length(Lxmonom),'free');
                Lue(j) = lu'*Lxmonom;
            end
            
            prog.sos = Lue(i);
            prog.sos = -Vdot - L1{i}*(rho - V) + Lue.*X.mssExtB' + Lui;
        end
    end
end

if nargin > 3
    
    opt = 0;
    
    for n = 1:size(Xbnd,2)  % For all outgoing transitions
        Lui = [];
        XbndN = [];
        for m = 1:size(Xbnd,1) % for all transition funnels
            if ~isempty(Xbnd{m,n})
                % Unsafe set constraints for interior sets
                for i = 1:length(Xbnd{m,n})
                    degL = 1;
                    Lxmonom = monomials(x,0:degL);
                    [prog,lu] = new(prog,length(Lxmonom),'free');
                    Lui = [Lui lu'*Lxmonom];
                    prog.sos = Lui(i);
                end
                XbndN = [XbndN Xbnd{m,n}'];
            end
        end
        LuIn = [];
        XinN = [];
        for m = 1:length(Xin) % for all inward funnels
            if ~isempty(Xin{m})
                % Unsafe set constraints for interior sets
                for i = 1:length(Xin{m})
                    degL = 1;
                    Lxmonom = monomials(x,0:degL);
                    [prog,lu] = new(prog,length(Lxmonom),'free');
                    LuIn = [LuIn lu'*Lxmonom];
                    prog.sos = LuIn(i);
                end
                XinN = [XinN Xin{m}'];
            end            
        end
        if isempty(XinN)
            XinN = 0;
            LuIn = 0;
        end
        if isempty(XbndN)
            XbndN = 0;
            Lui = 0;
        end
        prog.sos = (V - rho) + 1*Lui*XbndN' + LuIn*XinN';
    end
end

% TODO: limit final ellipse to within some acceptable radius for seq. comp... or use a ROA?

% NB: usually can get a reasonable solution by attempting both maximization & minimization 
[prog,info] = sedumi(prog,-rho,0,[],opt);
if info.pinf
    disp('fail.')
    rho = 1e9;
else
    rho = double(prog([rho]))
end

end
