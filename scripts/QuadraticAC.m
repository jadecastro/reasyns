
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
        
        function E = ellipsoid(obj)
            % Output an array of ellipsoid objects based on parameters from
            % the atomic controller
            for i = 1:size(obj,1)
                for j = 1:size(obj,2)
                    t = getTimeVec(obj(i,j).x0);
                    x = ppval(obj(i,j).x0.pp,t);
                    P = ppval(obj(i,j).P.pp,t);
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
        function obj = normalizerho(obj)
            % TODO: create a flag to indicate if already normalized
            rho = ppval(obj.rho,getTimeVec(obj));
%             P = ppval(obj.P,getTimeVec(obj));
            if all(rho) == 1.0, return; end
            obj.P = obj.P./obj.rho;
            obj.rho = spline(getTimeVec(obj),ones(length(getTimeVec(obj)),1));
        end
        
        % TODO: generalize to arbitrary polynomials
        function Eproj = projection(obj,sys,idx)
            %
            n = obj.x0.pp.dim;
            if ~isempty(sys.H)
                if nargin > 2
                    % makes sense to assume only one ac in this case
                    ell = ellipsoid(obj);
                    Eproj = projection(ell(idx),[sys.H; zeros(n-length(sys.H),length(sys.H))]);
                else
                    for i = 1:length(obj)
                        try
                            ell = ellipsoid(obj(i));
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
        
%         function newobj = interp(obj,newobj,ts)
%             
%             newobj = obj.interp(obj,newobj,ts);
%             
%             for i = 1:Nsteps;
%                 if any(obj.isCyclic)
%                     ith(i) = 1;
%                     if obj(i).x(nn) > pi
%                         ith(i) = 2;
%                     elseif obj(i).x(nn) < -pi
%                         ith(i) = 3;
%                     end
%                 else
%                     ith(i) = 1;
%                 end
%                 
%                 indx = find(min(abs(obj(i).t - ts{ith(i)})) == abs(obj(i).t - ts{ith(i)}),1,'first');
%                 obj(i).P = obj(i).P{ith(i)};
%             end
%         end
        
        function res = isinternal(obj, X, s, varargin)
            %
            res = isinternal_quickInv(obj.Einv, X, s);

            if ~isempty(varargin)
                sys = varargin{1};
                %NB: not complete.
                if ~res
                    for i = find(sys.params.isCyclic)'
                        X(i) = X(i) + 2*pi;
                        res = isinternal_quickInv(obj.Einv, X, s);
                        X(i) = X(i) - 2*pi;
                    end
                end
                if ~res
                    for i = find(sys.params.isCyclic)'
                        X(i) = X(i) - 2*pi;
                        res = isinternal_quickInv(obj.Einv, X, s);
                        X(i) = X(i) + 2*pi;
                    end
                end
            end
                
        end
        
        function [res, idx, isectArray] = isinside(obj,regobj,sys)
            % checks for containment by simply checking intersections with
            % any of the boundary hyperplanes.
            [H,K] = double(regobj.p);
            hpp = hyperplane(H',K');
            isectArray = [];
            if nargout > 1
                res = true;
                idx = [];
                for i = 1:length(obj.ellipsoid)
                    tmp = ~intersect(projection(obj,sys,i),hpp,'u');
                    isectArray = [isectArray; tmp];
                    if ~all(tmp)
                        res = false;
                        idx = [idx; i];
                    end
                end
            else
                tmp = ~intersect(projection(obj,sys),hpp,'u');
                isectArray = [isectArray; tmp];
                res = ~all(tmp);
            end
        end
        
%         function sample
%         end
        
        function ac = append(ac)
            error('wip')
        end
        
        function obj = merge(obj1,obj2)
            %
            actmp = merge@PolynomialAC(obj1,obj2);
            Pn = timecat(obj1.P,obj2.P);
            obj = QuadraticAC(actmp.x0,actmp.u0,actmp.K,Pn,actmp.rho,actmp.V,actmp.sys);
        end
        
        function plot(obj,sys,fignum,threedflag,color)
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
                Eproj = projection(obj,sys);
                figure(fignum)
                if isColorAsOption
                    plot(Eproj,Options)
                else
                    plot(Eproj,color)
                end
            else
                if threedflag
                    tmp = downsample(obj.ellipsoid,1);
                    figure(fignum)
                    hold on
                    plot(tmp,'r')
                end
            end
            drawnow
        end
        
        function u = execute(obj,x)
            %
            % look up the command input given the state
            
            teval = 0.02;
            
            % Use a weighted Euclidean distance to resolve the current state of the TVLQR controller.
            weights = 1 * (~obj.sys.params.isCyclic) + 0.2 * obj.sys.params.isCyclic;
            
            ell = ellipsoid(obj);
            t_trials = getTimeVec(obj.x0);
            
            minDelta = inf;
            
            for i = 2:1:length(ell)
                xtmp = double(obj.x0, t_trials(i));
                % unwrap the heading cyclic coordinate
                [testMinDist, minIdx] = min([
                    norm(weights.*(xtmp - (x+obj.sys.params.isCyclic*obj.sys.params.limsNonRegState(1)))) 
                    norm(weights.*(xtmp - x)) 
                    norm(weights.*(xtmp - (x+obj.sys.params.isCyclic*obj.sys.params.limsNonRegState(2))))
                    ]);
                if testMinDist < minDelta
                    teval = t_trials(i);
                    minDelta = testMinDist;
                end
            end
            
            % select the best-fit parameters
            K = double(obj.K,teval);
            x0 = double(obj.x0,teval);
            u0 = double(obj.u0,teval);
            
            K = K(end-length(u0)+1:end,:);
            
            u_test = [
                (K*(x+obj.sys.params.isCyclic*obj.sys.params.limsNonRegState(1) - x0))'
                (K*(x - x0))'
                (K*(x+obj.sys.params.isCyclic*obj.sys.params.limsNonRegState(2) - x0))'
                ]';
            % u_idx = abs(u_test(end,:)) == min(abs(u_test(end,:)));  
            u_ctrl = u_test(:,minIdx);
            u = u0 + u_ctrl(:,1);
        end
        
        function disp(obj)
            fprintf('Quadratic Atomic Controller with parameters\n');
        end
        
    end
    
end