function [lyapfun,info] = computeVasROA(x,xdot,X,Xbnd,Xin,xMinGrid,Nb)

prog = mssprog;

Vmonom = monomials(x,0:Nb);
[prog,b] = new(prog,length(Vmonom),'free');
V = b'*Vmonom;

Vdot = diff(V,t) + diff(V,x)*xdot;
Lxmonom = monomials(x,0:deg(Vdot,x)-deg(V,x));
[prog,l] = new(prog,length(Lxmonom),'free');
L1 = l'*Lxmonom;
% prog.sos = L1;
prog.sos = -(Vdot - L1*V);

if nargin == 3
    
    opt = 1;
    
    % Unsafe set constraints for exterior sets
    if ~isempty(X.mssExtB)
        clear Lue
        for i = 1:length(X.mssExtB)
            degL = 1;
            Lxmonom = monomials(x,0:degL);
            [prog,lu] = new(prog,length(Lxmonom),'free');
            Lue(i) = lu'*Lxmonom;
            prog.sos = Lue(i);
            prog.sos = V + Lue(i)*X.mssExtB(i);
        end
    end
    
    % Unsafe set constraints for interior sets
    for j = 1:length(X.mssIntB)
        if ~isempty(X.mssIntB{j})
            clear Lui
            for i = 1:length(X.mssIntB{j})
                degL = 2;
                Lxmonom = monomials(x,0:degL);
                [prog,lu] = new(prog,length(Lxmonom),'free');
                Lui(i) = lu'*Lxmonom;
                prog.sos = Lui(i);
            end
            prog.sos = V + Lui*X.mssIntB{j}';
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
        prog.sos = V + 1*Lui*XbndN' + LuIn*XinN';
    end
end

% NB: usually can get a reasonable solution by attempting both maximization & minimization 
pars.fid = 1;
pars.theta = 0.25;

% Define an objective (optional)
% obj = subs(V,x,zeros(length(x),1));
obj1 = 0;
for i = 1:size(xMinGrid,1)
    obj1 = obj1 + 1*subs(V,x,[xMinGrid(i,1);xMinGrid(i,2);xMinGrid(i,3)]);
%     obj1 = obj1 + 1*subs(V,x,[xMinGrid(i,1);xMinGrid(i,2);sin(xMinGrid(i,3));cos(xMinGrid(i,3))]);
end

% obj = msspoly('g',1);
% prog.free = obj;
prog = sedumi(prog,obj1,0,pars,1);

Lgb = prog([Lgb]);
lyapfun = prog([V]);

