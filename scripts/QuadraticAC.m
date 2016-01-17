
classdef QuadraticAC < PolynomialAC
    
    properties
        % quadratically-parameterized funnel properties
        P;
    end
    
    properties(SetAccess = private) 
        Einv;
    end
    
    methods
        function obj = QuadraticAC(x0,u0,K,P,rho,V,sys)
            % Constructor
            
            %             if nargin ~= 0
            %                 if length(t) ~= size(x,1), error('length of x must match that of t'), end
            %                 if length(t) ~= size(u,1), error('length of u must match that of t'), end
            %                 if length(t) ~= size(P,3), error('length of P must match that of t'), end
            %                 if length(t) ~= size(K,3), error('length of K must match that of t'), end
            %                 if length(t) ~= size(rho,1), error('length of rho must match that of t'), end
            %                 if length(t) ~= size(V,1), error('length of V must match that of t'), end
            %
            %                 this(length(t)) = atomiccontroller;
            %                 for i = 1:length(t)
            %                     this(i).t = t(i);
            %                     this(i).x = x(i,:);
            %                     this(i).u = u(i,:);
            %                     this(i).P = P(:,:,i);
            %                     this(i).K = K(:,:,i);
            %                     this(i).rho = rho(i);
            %                     this(i).V = V(i);
            %                 end
            %
            %             end
            
            obj = obj@PolynomialAC(x0,u0,K,rho,V,sys);
            
            if ~isa(P,'Traject')
                error('One or more of the arguments is not of type Traject!')
            end
            t = getTimeVec(obj.x0);
            tlen = length(t);
            if tlen ~= length(getTimeVec(P))
                warning('One or more time vectors have dimensions that don''t match')
            end
            
            obj.P = P;
            for i = 1:length(t)
                obj.Einv(i).P = ppval(obj.P.pp,t(i));
                obj.Einv(i).x = ppval(obj.x0.pp,t(i));
            end
        end
        
%         function display(obj)
%             % 
%             
% %             obj.display(obj);
%             
%             fprintf('P:\n'); disp(obj.P);
%             fprintf('rho:\n'); disp(obj.rho);
%             fprintf('V:\n'); disp(obj.V);
%             
%             if isempty(obj)
%                 fprintf('Empty atomic controller.\n');
%                 return;
%             end
%             
%             return;
%         end
        
        function E = ellipsoid(ac)
            % Output an array of ellipsoid objects based on parameters from
            % the atomic controller
            for i = 1:size(ac,1)
                for j = 1:size(ac,2)
                    t = getTimeVec(ac(i,j).x0);
                    x = ppval(ac(i,j).x0.pp,t);
                    P = ppval(ac(i,j).P.pp,t);
                    for k = 1:size(P,3)
                        Pinv = inv(P(:,:,k));
                        Pinv_sym = (Pinv'+Pinv)/2;  % to ensure it is symmetric
                        try 
                            E(k) = ellipsoid(x(:,k),Pinv_sym);
                        catch ME
                            Pinv = inv(double(E(k-1)));
                            Pinv_sym = (Pinv'+Pinv)/2;  % to ensure it is symmetric
                            E(k) = ellipsoid(x(:,k),Pinv_sym);
                            %rethrow(ME)
                        end
                    end
                end
            end
        end
        
        % TODO: generalize to arbitrary polynomials
        function ac = normalizerho(ac)
            % TODO: create a flag to indicate if already normalized
            rho = ppval(ac.rho,getTimeVec(ac));
%             P = ppval(ac.P,getTimeVec(ac));
            if all(rho) == 1.0, return; end
            ac.P = ac.P./ac.rho;
            ac.rho = spline(getTimeVec(ac),ones(length(getTimeVec(ac)),1));
        end
        
        % TODO: generalize to arbitrary polynomials
        function Eproj = projection(ac,sys,idx)
            %
            n = ac.x0.pp.dim;
            if ~isempty(sys.H)
                if nargin > 2
                    % makes sense to assume only one ac in this case
                    ell = ellipsoid(ac);
                    Eproj = projection(ell(idx),[sys.H; zeros(n-length(sys.H),length(sys.H))]);
                else
                    for i = 1:length(ac)
                        try
                            ell = ellipsoid(ac(i));
                        catch
                            keyboard
                        end
                        for j = 1:length(ell)
                            Eproj(j,i) = projection(ell(j),[sys.H; zeros(n-length(sys.H),length(sys.H))]);
                        end
                    end
                end
            else
                error('H not defined. Only linear state-output relationships supported.')
            end
        end
        
%         function newobj = interp(ac,newobj,ts)
%             
%             newobj = ac.interp(ac,newobj,ts);
%             
%             for i = 1:Nsteps;
%                 if any(ac.isCyclic)
%                     ith(i) = 1;
%                     if ac(i).x(nn) > pi
%                         ith(i) = 2;
%                     elseif ac(i).x(nn) < -pi
%                         ith(i) = 3;
%                     end
%                 else
%                     ith(i) = 1;
%                 end
%                 
%                 indx = find(min(abs(ac(i).t - ts{ith(i)})) == abs(ac(i).t - ts{ith(i)}),1,'first');
%                 ac(i).P = ac(i).P{ith(i)};
%             end
%         end
        
        function union
        end
        
        function intersect
        end
        
        function res = isinternal(acobj, X, s, varargin)
            %
            res = isinternal_quickInv(acobj.Einv, X, s);

            if ~isempty(varargin)
                sys = varargin{1};
                %NB: not complete.
                if ~res
                    for i = find(sys.params.isCyclic)'
                        X(i) = X(i) + 2*pi;
                        res = isinternal_quickInv(acobj.Einv, X, s);
                        X(i) = X(i) - 2*pi;
                    end
                end
                if ~res
                    for i = find(sys.params.isCyclic)'
                        X(i) = X(i) - 2*pi;
                        res = isinternal_quickInv(acobj.Einv, X, s);
                        X(i) = X(i) + 2*pi;
                    end
                end
            end
                
        end
        
        function [res, idx, isectArray] = isinside(acobj,regobj,sys)
            % checks for containment by simply checking intersections with
            % any of the boundary hyperplanes.
            [H,K] = double(regobj.p);
            hpp = hyperplane(H',K');
            isectArray = [];
            if nargout > 1
                res = true;
                idx = [];
                for i = 1:length(acobj.ellipsoid)
                    tmp = ~intersect(projection(acobj,sys,i),hpp,'u');
                    isectArray = [isectArray; tmp];
                    if ~all(tmp)
                        res = false;
                        idx = [idx; i];
                    end
                end
            else
                tmp = ~intersect(projection(acobj,sys),hpp,'u');
                isectArray = [isectArray; tmp];
                res = ~all(tmp);
            end
        end
        
%         function sample
%         end
        
        function ac = append(ac)
            error('wip')
        end
        
        function acres = merge(acobj1,acobj2)
            %
            actmp = merge@PolynomialAC(acobj1,acobj2);
            Pn = timecat(acobj1.P,acobj2.P);
            acres = QuadraticAC(actmp.x0,actmp.u0,actmp.K,Pn,actmp.rho,actmp.V,actmp.sys);
        end
        
        function plot(ac,sys,fignum,threedflag,color)
            %
            %TODO: treat both 2d and 3d cases
            
            if ~exist('threedflag')
                threedflag = [];
            end
            
            isColorAsOption = false;
            if nargin < 5
                color = 'b';
            elseif length(color) == 3 && isnumeric(color)
                Options.color = color;
                isColorAsOption = true;
            end
            Options.width = 3;
            
            if isempty(threedflag)
                Eproj = projection(ac,sys);
                figure(fignum)
                if isColorAsOption
                    plot(Eproj,Options)
                else
                    plot(Eproj,color)
                end
            else
                if threedflag
                    tmp = downsample(ac.ellipsoid,1);
                    figure(fignum)
                    hold on
                    plot(tmp,'r')
                end
            end
            drawnow
        end
        
        function disp(ac)
            fprintf('Quadratic Atomic Controller with parameters\n');
        end
        
    end
    
end