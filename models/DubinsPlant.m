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
                [df,d2f,d3f]= dynamicsGradientsDubins(obj,t,x,u,nargout-1);
            end
        end
        
        function y = output(obj,t,x,u)
            y = x;
        end
            
        function x = getInitialState(obj)
            x = [0 0 0]';
        end
        
        function y = state2SEconfig(obj,t,x,u)
            if size(x,1) == 1
                x = x';
            end
            y = x;
        end

        % NOTE: The following two methods are optional.
        function [xdot] = dynamicsWaypointSteering(obj,t,x,xyPath,gotopt1)
            
            global gotopt
            
            [u,gotopt] = obj.steerToXYWaypoints(x,xyPath,gotopt);
            
            xdot = dynamics(obj,t,x,u);
            
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
        
    end
    
    methods (Access = private)
        
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
        
    end
    
end
