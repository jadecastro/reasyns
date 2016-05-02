% BarrierFunctionAC  A class for processing atomic controllers assuming the
% form of barrier functions.
%
% SYNOPSIS
%
%    BarrierFunctionAC(x0,u0,K,B,c,V0,sys)
%        Instatiates an atomic controller, with nominal state trajectory x0
%        (Traject), control input trajectory u0 (Traject), control gain
%        matrix trajectory K (Traject), barrier functions B (cell array),
%        controller c, Lyapunov function V0, and system model sys
%        (SystemDynamics).
%
% METHODS
%
%    result = obj.isinternal(x, s)
%    obj    = obj.merge(obj2)
%    obj.plot(command)
%    u = obj.execute(x, poly)
%    obj.disp()


classdef BarrierFunctionAC < PolynomialAC

    properties
        % Barrier function-based funnel properties
        c;
        V0;
    end
    
    properties(SetAccess = private) 
        Einv;
    end
    
    methods
        function obj = BarrierFunctionAC(x0,u0,K,B,c,V0,sys,varargin)
            % Constructor
            
            t = x0.getTimeVec();
            rho0 = Traject(t,zeros(length(t),1));
            
            obj = obj@PolynomialAC(x0,u0,K,rho0,B,sys,varargin);
            obj.c = c;
            obj.V0 = V0;
        end
        
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
        
        function obj = merge(obj,obj2)
            %
            actmp = merge@PolynomialAC(obj1,obj2);
            error('not implemented')
        end
        
        function plot(obj,x0,fignum)
            %
            %TODO: treat both 2d and 3d cases
            
            figure(fignum)
            N = length(obj.V);
            
            [uRed,tRed] = downsampleUniformly(x0,10);
            colorArray = colormap('cool');
            x = msspoly('x',obj.sys.params.n);
            
            for ii = 1:N
                color = colorArray(floor(ii*size(colorArray,1)/N),:);
                
                %     prog(B{ii})
                
                % check for sanity- check equilibrium point.  TODO: check more points in
                % the IC set too
                
                x00 = x0.double(tRed(ii));
                
                % 0-level set of B
                [X,Y] = meshgrid(linspace(-3,2,100),linspace(-2,2,100));
                Th = x00(obj.sys.params.n)*ones(size(X));
                gs_B = msubs(obj.V{ii},x,[X(:),Y(:),Th(:)]');
                [~,H] = contour(X,Y,reshape(double(gs_B),100,100),[0 0],'LineColor',color,'LineWidth',3);
                %set(H,'LineColor',color)
                hold on
            end
            drawnow
        end
        
        function u = execute(obj,x,poly)
            %
            % look up the command input given the state
            
            gamma = 1;
            
            teval = 0.02;
            xs = msspoly('xs',obj.sys.params.n);
            
            % Use a weighted Euclidean distance to resolve the current state of the TVLQR controller.
            weights = 1 * (~obj.sys.params.isCyclic) + 0.2 * obj.sys.params.isCyclic;
            
            t_trials = getTimeVec(obj.x0);
            
            minDelta = inf;
            
            for i = 1:1:length(t_trials)
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
            
            % select the best-fit parameters and compute the TVLQR feedback controller 
            K = double(obj.K,teval);
            [x0,t0] = double(obj.x0,teval);
            u0 = double(obj.u0,teval);
            
            K = K(end-length(u0)+1:end,:);
            
            u_test = [
                (K*(x+obj.sys.params.isCyclic*obj.sys.params.limsNonRegState(1) - x0))'
                (K*(x - x0))'
                (K*(x+obj.sys.params.isCyclic*obj.sys.params.limsNonRegState(2) - x0))'
                ]';
            % u_idx = abs(u_test(end,:)) == min(abs(u_test(end,:)));  
            u_ctrl = u_test(:,idxEval);
            u_pre = u0 + u_ctrl(:,1);
            
            fs = poly.getPolyDynamics(t_trials(idxEval));
            if (poly.getNumInputs > 0)   % zero all inputs
                fs = subs(fs,poly.getInputFrame.getPoly,zeros(poly.getNumInputs,1));
            end
            gs = [0;0;1];
            
            % Compute time derivatives and Lie derivatives of B
            Bs  = obj.V{idxEval};
            Bsm = obj.V{idxEval-1};
            if idxEval > 1
                dBdt = (Bs - Bsm)/(t_trials(idxEval) - t_trials(idxEval-1));  % take an approximate derivative - valid for small delta t's
            else
                dBdt = 0;
            end
            a = diff(Bs,xs)*fs + dBdt;
            b = diff(Bs,xs)*gs;

            % Compute the TV-CBF control law
            if double(b) == 0
                u_cbf = 0;
            else
                u_cbf = -1/(b'*b)*(a + sqrt(a^2 + gamma^2*b'*b))*b;
            end
            
            u_cbf
            u = u_pre + u_cbf;
            
        end
        
        function disp(obj)
            fprintf('Barrier Function Atomic Controller with parameters\n');
        end
           
%         function sample
%         end
        
    end
    
end