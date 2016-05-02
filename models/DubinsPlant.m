classdef DubinsPlant < DrakeSystem
    % Defines the dynamics for the Dubins car/unicycle model.
    
    properties
        sysparams;
    end
    
    methods
        function obj = DubinsPlant(sysparams)
            obj = obj@DrakeSystem(3,0,1,3,0,1);
            obj.sysparams = sysparams;
            
            obj = setOutputFrame(obj,getStateFrame(obj));  % allow full state feedback
            
            obj = setInputLimits(obj,-obj.sysparams.max_w,obj.sysparams.max_w);
        end
        
        function [xdot, df, d2f, d3f] = dynamics(obj,t,x,u)
            theta = x(3);
            xdot = [obj.sysparams.v*cos(theta);  obj.sysparams.v*sin(theta); u(1)];
            
            if (nargout>1)
                [df,d2f,d3f]= dynamicsGradients(obj,t,x,u,nargout-1);
            end
        end
        
        function [xdot] = dynamicsWaypointSteering(obj,t,x,xyPath,gotopt1)
            
            global gotopt
            
            [u,gotopt] = obj.steerToXYWaypoints(x,xyPath,gotopt);
            
            xdot = dynamics(obj,t,x,u);
            
        end
        
        function y = output(obj,t,x,u)
            y = x;
        end
        
        function y = state2SEconfig(obj,t,x,u)
            if size(x,1) == 1
                x = x';
            end
            y = x;
        end
        
        function x = SEconfig2state(obj,t,y,u)
            if size(y,1) == 1
                y = y';
            end
            x = y;
        end
        
        function [idxRegStates, idxNonRegStates, H] = getRegNonRegStates(obj,varargin)
            
            n = obj.sysparams.n;
            epsilon = 1e-6;
            
            if ~isempty(varargin)
                t = varargin{1};
                x = varargin{2};
                u = varargin{3};
            else
                t = [];
                x = zeros(n,1);
                u = [];
            end
            
            idxRegStates = [];
            if nargout == 3
                y = obj.state2SEconfig(t,x,u);
                H = zeros(length(y),n);
            else
                H = [];
            end
            
            for i = 1:n
                % Find the output at the baseline
                tryState = x;
                y = obj.state2SEconfig(t,tryState,u);
                
                % Perturb the states, one by one
                tryState(i) = tryState(i) + epsilon;
                yPert = obj.state2SEconfig(t,tryState,u);
                
                if (yPert(1) - y(1)) ~= 0 || (yPert(2) - y(2)) ~= 0
                    % flag any states that are coupled to x-y position (first two coordinates of an SE(2) system).
                    idxRegStates = [idxRegStates; i];
                end
                
                if nargout == 3
                    H(:,i) = (yPert - y)/epsilon;
                end
            end
            idxNonRegStates = setdiff(1:n,idxRegStates);
        end
        
        function x = getInitialState(obj)
            x = [0 0 0]';
        end
        
        function [u, waypointIndex] = steerToXYWaypoints(obj,x,xyWaypoints,waypointIndex)
            
            p = obj;
            
            e = obj.sysparams.e;
            closeEnough = obj.sysparams.closeEnough;
            
            xr = x(1); yr = x(2); th = x(3);
            
            if waypointIndex <= size(xyWaypoints,1)
                desXR = xyWaypoints(waypointIndex,1);  desYR = xyWaypoints(waypointIndex,2);
            else
                desXR = xyWaypoints(end,1);  desYR = xyWaypoints(end,2);
            end
            
            euclidDist2NextPt = sqrt((xr - desXR)^2 + (yr - desYR)^2);
            
            % If robot is within acceptance radius, index to next waypoint
            if abs(euclidDist2NextPt) < abs(closeEnough)
                waypointIndex = waypointIndex + 1;
                %     disp(['waypointIndex = ',num2str(waypointIndex)])
            end
            
            Vx = desXR - xr;
            Vy = desYR - yr;
            [~, uW] = obj.feedbackLin(Vx, Vy, th, e);
            
            % Limit on velocity
            %           uV = 0.08;%max(uV,0);  % Limits velocity on (1,0)
            %           uW = max(min(uW,0.13),-0.13);  % Limits angular velocity on (-1,1)
            
            % TODO: apply angular rate limits using the Drake approach here
            % p = setInputLimits(p,-obj.max_w,obj.max_w);
            
            uW = max(min(uW, obj.sysparams.max_w), -obj.sysparams.max_w);  % Limits angular velocity on (-1,1)
            
            u = uW;
            
        end
        
        function [c, V] = computeTVLQR(obj, utraj, xtraj)
            %
            % Compute the time-varying LQR controller
            
            p = obj;
            
            Q = obj.sysparams.Q;
            R = obj.sysparams.R;
            Qf = obj.sysparams.Qf;
            
            % Set input limits
            p = setInputLimits(p,-Inf,Inf);
            
            % Do tvlqr
            [c,V] = tvlqr(p,xtraj,utraj,Q,R,Qf);
        end
        
        function [cmdV, cmdW] = feedbackLin(obj, Vx, Vy, thetaR, e)
            % FEEDBACKLIN: Apply feedback linearization to convert desired inertial
            % velocities to local robot velocities.
            %
            %   [CMDV,CMDW] = FEEDBACKLIN(VX,VY,THETAR,E)
            %
            %   INPUTS
            %       Vx          desired x-velocity in global frame (m/s)
            %       Vy          desired y-velocity in global frame (m/s)
            %       thetaR      current robot angle relative to inertial frame (deg)
            %       e           feedback linearization distance (in meters)
            %
            %   OUTPUTS
            %       cmdV        forward velocity command (m/s)
            %       cmdW        angular velocity command (rad/s)
            %
            
            
            VxyI = [Vx; Vy];
            
            % Rotation from intertial to body frame
            RIB = [ cos(thetaR) -sin(thetaR);
                sin(thetaR)  cos(thetaR)
                ];
            
            % Rotation from body to inertial frame
            RBI = inv(RIB);
            
            % Perform feedback linearization
            VW = [1 0;0 1/e]*RBI*VxyI;
            
            % Output commanded V and w
            cmdV = VW(1);
            cmdW = VW(2);
        end
        
        function [poly, p] = getSystemPoly(obj, xtraj, c, V)
            %
            % Compute the polynomial approximation of the drake plant
            %  xtraj and V must be a PPTrajectory, c must be an AffineSystem (Drake)
            
            p = obj;
            
            % Set input limits
            p = setInputLimits(p,-Inf,Inf);
            
            % compute the polynomial approximation
            poly = taylorApprox(feedback(p,c),xtraj,[],3);
            
            num_xc = poly.getNumContStates();
            if (isa(V,'Trajectory'))
                V1 = V.inFrame(poly.getStateFrame);
                Q = eye(num_xc);
                V0 = tvlyap(poly,V1,Q,Q);
            else
                V0 = V;
            end
            
            poly = poly.inStateFrame(V0.getFrame); % convert system to Lyapunov function coordinates
            
        end
        
        function [utraj,xtraj]=runDircol(p,x0,xf,tf0)
            
            N = 21;
            prog = DircolTrajectoryOptimization(p,N,[.1 1]);
            prog = prog.addInputConstraint(BoundingBoxConstraint(.5*p.umin,.5*p.umax),1:N);
            prog = prog.addStateConstraint(ConstantConstraint(x0),1);
            prog = prog.addStateConstraint(ConstantConstraint(xf),N);
            prog = prog.addRunningCost(@cost);
            prog = prog.addFinalCost(@finalCost);
            
            function [g,dg] = cost(dt,x,u)
                R = 0;
                g = u'*R*u;
                %g = sum((R*u).*u,1);
                %dg = [zeros(1,1+size(x,1)),2*u'*R];
                dg = zeros(1, 1 + size(x,1) + size(u,1));
            end
            
            function [h,dh] = finalCost(t,x)
                h = t;
                dh = [1,zeros(1,size(x,1))];
            end
            
            traj_init.x = PPTrajectory(foh([0,tf0],[x0,xf]));
            info = 0;
            while (info~=1)
                tic
                [xtraj,utraj,z,F,info] = prog.solveTraj(tf0);
                toc
            end
            
        end
        
    end
    
end
