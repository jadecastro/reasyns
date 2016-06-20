
classdef QuadraticAC < PolynomialAC
    
    properties
        % quadratically-parameterized funnel properties
        P;
        ell;
    end
    
    properties(SetAccess = private) 
        Einv;
    end
    
    methods
        function obj = QuadraticAC(x0,u0,K,P,rho,V,sys,varargin)
            % Constructor
            
            obj = obj@PolynomialAC(x0,u0,K,rho,V,sys,varargin);
            
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
            obj.ell = obj.ellipsoid();
        end
        
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
        function obj1 = normalizerho(obj)
            % TODO: create a flag to indicate if already normalized
            t = obj.x0.getTimeVec;
            rho = double(obj.rho,t);
%             P = ppval(obj.P,getTimeVec(obj));
            if all(rho == 1.0), obj1 = obj; return; end
            for i = 1:length(t)
                newP(:,:,i) = double(obj.P,t(i))/double(obj.rho,t(i));
            end
            P = Traject(t, newP);
            rho = Traject(t, ones(1,length(t)));
            obj1 = QuadraticAC(obj.x0,obj.u0,obj.K,P,rho,obj.V,obj.sys);
        end
        
        % TODO: generalize to arbitrary polynomials
        function Eproj = projection(obj,sys,idx)
            %
            n = obj.x0.pp.dim;
            
            if nargin > 2
                % makes sense to assume only one ac in this case
                [c,~] = double(obj.ell(idx));
                [~,~,H] = sys.getRegNonRegStates([],c,[]);
                Eproj = projection(obj.ell(idx),H(1:2,:)');
            else
                for i = 1:length(obj)
                    ell = obj(i).ell;
                    for j = 1:length(ell)
                        [c,~] = double(ell(j));
                        [~,~,H] = sys.getRegNonRegStates([],c,[]);
                        Eproj(j,i) = projection(ell(j),H(1:2,:)');
                    end
                end
            end
        end
        
        function res = isinternal(obj, X, s, varargin)
            % Checks if a point X lies within the atomiccontroller
            % funnel.  Returns 'true' if X belongs to the funnel interior
            % and 'false' otherwise.
            
            res = isinternal_quickInv(obj.Einv, X, s);

            if ~isempty(varargin)
                sys = varargin{1};
                %NB: not complete.
                if ~res
                    for i = find(sys.sysparams.isCyclic)'
                        X(i) = X(i) + 2*pi;
                        res = isinternal_quickInv(obj.Einv, X, s);
                        X(i) = X(i) - 2*pi;
                    end
                end
                if ~res
                    for i = find(sys.sysparams.isCyclic)'
                        X(i) = X(i) - 2*pi;
                        res = isinternal_quickInv(obj.Einv, X, s);
                        X(i) = X(i) + 2*pi;
                    end
                end
            end
                
        end
        
        function [res, idx, isectArray] = isinside(obj,regobj,sys,varargin)  % (TODO: 'isinside' is a misnomer- either fix or rename it)
            % Checks for containment of a funnel within a polytope by simply checking intersections with
            % any of the boundary hyperplanes.  Returns 'true' if none of
            % the ellipsoids in the atomiccontroller intersect with the
            % boundary of the region and 'false' otherwise.
            
            if length(regobj) > 2, error('Region arrays greater than two not supported.'); end
            if nargin < 4, sampSkip = 1; else, sampSkip = varargin{1}; end
                        
            isectArray = [];
            for k = 1:length(regobj)
                [H,K] = double(regobj(k).p);
                hpp{k} = hyperplane(H',K');
            end
            
            res = true;
            if nargout > 1 || length(regobj) > 1
                idx = [];
                for i = 1:sampSkip:length(obj.ell)
                    for k = 1:length(regobj)
                        nonIntersectingPlanes{k} = ~intersect(obj.projection(sys,i),hpp{k},'u');
                    end
                    
                    if length(regobj) == 1, nonIntersectingPlanes = nonIntersectingPlanes{1}; end
                    
                    isectArray = [isectArray; nonIntersectingPlanes];
                    
                    if length(regobj) == 1
                        if ~all(nonIntersectingPlanes)
                            res = false;
                            idx = [idx; i];
                        end
                    elseif length(regobj) == 2
                        % Given two regions which are adjacent, and each convex, and an ellipse centered in one of the two regions,
                        % an intersection of either zero or precisely one hyperplane for each region implies that the ellipse is inside the union of the two regions.
%                         ~nonIntersectingPlanes{1}
%                         ~nonIntersectingPlanes{2}
                        res1 = sum(vertcat(~nonIntersectingPlanes{1}));
                        res2 = sum(vertcat(~nonIntersectingPlanes{2}));
                        if ~((res1 == 0 || res2 == 0) || (res1 == 1 && res2 == 1)),
                            res = false;
                            if nargout < 2, break; end
                            idx = [idx; i];
                        end
                    end
                end
            else
                for k = 1:length(regobj)
                    nonIntersectingPlanes{k} = ~intersect(obj.projection(sys),hpp{k},'u');
                end
                
                if length(regobj) == 1, nonIntersectingPlanes = nonIntersectingPlanes{1}; end
                
%                 nonIntersectingPlanes
                isectArray = [isectArray; nonIntersectingPlanes];
                
                res = all(nonIntersectingPlanes);
            end
            
        end
        
        function [res, xFail] = funnelContainsEllipsoid(obj,sys,varargin)
            % Checks whether or not the funnel contains a given ellipsoid.
            % Returns 'true' if contained and 'false' otherwise.
            
            d = 100;
            xFail = [];
            
            if nargin > 4, d = varargin{2}; end
            
            if nargin > 3 && ~isa(varargin{1},'ellipsoid'), 
                x = varargin{1}; 
            else
                ell = varargin{1};
                
                [c,Q] = double(ell);
                x = sampleEllipsoidBoundary(c,Q,d);
            end
            
            res = true;
            
            % Loop over all randomly-selected points
            for i = 1:size(x,1)
                if ~obj.isinternal(x(i,:),'u',sys)
                    res = false;
                    xFail = [xFail; x(i,:)];
                    if nargout < 2
                        return;
                    end
                end
            end
            
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
                    tmp = downsample(obj.ell,1);
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
                    idxEval = minIdx;
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
            u_ctrl = u_test(:,idxEval);
            u = u0 + u_ctrl(:,1);
        end
        
        function disp(obj)
            fprintf('Quadratic Atomic Controller with parameters\n');
        end
        
    end
    
end