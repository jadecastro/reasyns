
classdef BarrierFunctionAC < PolynomialAC
    
    properties
        % quadratically-parameterized funnel properties
        c;
        V;
    end
    
    properties(SetAccess = private) 
        Einv;
    end
    
    methods
        function obj = BarrierFunctionAC(x0,u0,K,B,c,V,sys)
            % Constructor
            
            obj = obj@PolynomialAC(x0,u0,K,0,B,c,V,sys);
            
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
        
        function plot(obj,x0,sys,fignum)
            %
            %TODO: treat both 2d and 3d cases
            
            figure(fignum)
            N = length(obj.B);
            
            [uRed,tRed] = downsampleUniformly(x0,10);
            colorArray = colormap('cool');
            x = msspoly('x',sys.params.n);
            
            for ii = 1:N
                color = colorArray(floor(ii*size(colorArray,1)/N),:);
                
                %     prog(B{ii})
                
                % check for sanity- check equilibrium point.  TODO: check more points in
                % the IC set too
                
                x0 = double(x0,tRed(ii));
                
                % 0-level set of B
                [X,Y] = meshgrid(linspace(-3,2,100),linspace(-2,2,100));
                Th = x0(sys.params.n)*ones(size(X));
                gs_B = msubs(obj.B{ii},x,[X(:),Y(:),Th(:)]');
                [~,H] = contour(X,Y,reshape(double(gs_B),100,100),[0 0],'LineColor',color,'LineWidth',3);
                %set(H,'LineColor',color)
                hold on
            end
            drawnow
        end
        
        function u = execute(obj,x)
            %
            % look up the command input given the state
            
            gamma = 1;
            
            teval = 0.02;
            xs = msspoly('xs',sys.params.n);
            
            % Use a weighted Euclidean distance to resolve the current state of the TVLQR controller.
            weights = 1 * (~obj.sys.params.isCyclic) + 0.2 * obj.sys.params.isCyclic;
            
            t_trials = getTimeVec(obj.x0);
            
            minDelta = inf;
            
            for i = 1:1:length(obj.B)
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
            
            % select the best-fit parameters and compute the TVLQR feedback controller 
            K = double(obj.K,teval);
            [x0,t0] = double(obj.x0,teval);
            u0 = double(obj.u0,teval);
                        
            x0pp = PPTrajectory(foh(t0,x0)); 
            x0pp = setOutputFrame(x0pp,p.getStateFrame);
            
            K = K(end-length(u0)+1:end,:);
            
            u_test = [
                (K*(x+obj.sys.params.isCyclic*obj.sys.params.limsNonRegState(1) - x0))'
                (K*(x - x0))'
                (K*(x+obj.sys.params.isCyclic*obj.sys.params.limsNonRegState(2) - x0))'
                ]';
            % u_idx = abs(u_test(end,:)) == min(abs(u_test(end,:)));  
            u_ctrl = u_test(:,minIdx);
            u_pre = u0 + u_ctrl(:,1);
            
            % Compute the polynomial approximation of the drake plant
            [poly] = obj.sys.getSystemPoly(x0pp, c, V);
            
            fs = poly.getPolyDynamics(t_trials(minIdx));
            if (poly.getNumInputs > 0)   % zero all inputs
                fs = subs(fs,poly.getInputFrame.getPoly,zeros(poly.getNumInputs,1));
            end
            gs = [0;0;1];
            
            % Compute time derivatives and Lie derivatives of B
            Bs  = B{minIdx};
            Bsm = B{minIdx-1};
            if minIdx > 1
                dBdt = (Bs - Bsm)/(t_trials(minIdx) - t_trials(minIdx-1));  % take an approximate derivative - valid for small delta t's
            else
                dBdt = 0;
            end
            a = diff(Bs,xs)*fs + dBdt;
            b = diff(Bs,xs)*gs;

            % Compute the TV-CBF control law
            if b == 0
                u_cbf = 0;
            else
                u_cbf = -1/(b'*b)*(a + sqrt(a^2 + gamma^2*b'*b))*b;
            end
            
            u = u_cbf + u_pre;
            
        end
        
        function disp(obj)
            fprintf('Barrier Function Atomic Controller with parameters\n');
        end
        
    end
    
end