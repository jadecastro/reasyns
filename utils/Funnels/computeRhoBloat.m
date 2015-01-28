function [rho,info] = computeRho(V,xdot,X,Xinv,Xbnd,Xin)

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
    
    % Unsafe set constraints for exterior sets
    if bloatFlag 
        XmssExt = X.mssExtB;
    else
        XmssExt = X.mssExt;
    end
    if ~isempty(XmssExt)
        clear Lue
        for i = 1:length(XmssExt)
            degL = 1;
            Lxmonom = monomials(x,0:degL);
            [prog,lu] = new(prog,length(Lxmonom),'free');
            Lue(i) = lu'*Lxmonom;
            prog.sos = Lue(i);
            prog.sos = (V - rho) + Lue(i)*XmssExt(i);
        end
    end
    
    % Unsafe set constraints for interior sets
    if bloatFlag
        XmssInt = X.mssIntB;
    else
        XmssInt = X.mssInt;
    end
    for j = 1:length(XmssInt)
        if ~isempty(XmssInt{j})
            clear Lui
            for i = 1:length(XmssInt{j})
                degL = 1;
                Lxmonom = monomials(x,0:degL);
                [prog,lu] = new(prog,length(Lxmonom),'free');
                Lui(i) = lu'*Lxmonom;
                prog.sos = Lui(i);
            end
            prog.sos = (V - rho) + Lui*XmssInt{j}';
        end
    end
end

if nargin > 4
    
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
