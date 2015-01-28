
classdef quadraticAC < polynomialAC
    
    properties
        % quadratically-parameterized funnel properties
        P;
    end
    
    methods
        function obj = quadraticAC(x0,u0,K,P,rho,V,sys)
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
            
            obj = obj@polynomialAC(x0,u0,K,rho,V,sys);
            
            if ~isa(P,'traject')
                error('One or more of the arguments is not of type traject!')
            end
            t = getTimeVec(obj.x0);
            tlen = length(t);
            if tlen ~= length(getTimeVec(P))
                warning('One or more time vectors have dimensions that don''t match')
            end
            
            obj.P = P;
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
                    t = getTimeVec(ac);
                    P = ppval(ac.pp,t);
                    Pinv = inv(P);
                    Pinv_sym = (Pinv'+Pinv)/2;  % to ensure it is symmetric
                    E{i,j} = ellipsoids.ellipsoid(ac.x0',Pinv_sym);
                end
            end
        end
        
        % TODO: generalize to arbitrary polynomials
        function ac = normalizerho(ac)
            %
            rho = ppval(ac.rho,getTimeVec(ac));
%             P = ppval(ac.P,getTimeVec(ac));
            if all(rho) == 1.0, return; end
            ac.P = ac.P./ac.rho;
            ac.rho = spline(getTimeVec(ac),ones(length(getTimeVec(ac)),1));
        end
        
        % TODO: generalize to arbitrary polynomials
        function Eproj = projection(ac,sys)
            %
            if ~isempty(sys.H)
                for i = 1:length(ac)
                    Eproj(i) = projection(ac(i).E,[sys.H; zeros(sys.n-length(sys.H),length(sys.H))]);
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
        
        function res = isinternal(acobj, X, s)
            %
            E = ellipsoid(acobj);
            res = isinternal_quickInv(E, X, s);
        end
        
        function res = isinside(acobj,regobj)
            %
            
        end
        
%         function sample
%         end
        
        function ac = append(ac)
            error('wip')
        end
        
        function acres = merge(acobj1,acobj2)
            %
            actmp = merge@polynomialAC(acobj1,acobj2);
            Pn = timecat(acobj1.P,acobj2.P);
            acres = quadraticAC(actmp.x0,actmp.u0,actmp.K,Pn,actmp.rho,actmp.V,actmp.sys);
        end
        
        function plot(ac,sys,fignum,color)
            %
            %TODO: treat both 2d and 3d cases
            Eproj = projection(ac,sys);
            figure(fignum)
            plot(Eproj)
        end
        
        function disp
        end
        
    end
    
end