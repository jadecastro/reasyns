
classdef SystemDynamics
    
    properties
        polyMdlFun;
        mdlFun;
        ctrlFun;
        drakeplant;
        params;
        
        H;
%         isCyclic;
%         limsNonRegState;
    end
    
    methods
        function obj = SystemDynamics(polyMdlFun,mdlFun,ctrlFun,drakeplant,H,params)
            % Constructor
            
            obj.polyMdlFun = polyMdlFun;
            obj.mdlFun = mdlFun;
            obj.ctrlFun = ctrlFun;
            obj.drakeplant = drakeplant;
            obj.H = H;
            obj.params = params;
        end
        
        function y = state2config(obj,x)
            %
            y(1:length(obj),obj.params.m) = 0;
            if ~isempty(obj.H)
                for i = 1:length(obj)
                    y(i,:) = obj.H*x(i,:);
                end
            else
                error('H not defined. Only linear state-output relationships are supported.')
            end
        end
        
        function [c, V] = computeTVLQR(obj, utraj, xtraj)
            %
            % Compute the time-varying LQR controller
            
            Q = obj.params.ctrloptions.Q;
            R = obj.params.ctrloptions.R;
            Qf = obj.params.ctrloptions.Qf;
            
            % Declare system model
            p = obj.drakeplant;
            
            % Set input limits
            p = setInputLimits(p,-Inf,Inf);
            
            % Do tvlqr
            [c,V] = tvlqr(p,xtraj,utraj,Q,R,Qf);
        end
        
        function [poly] = getSystemPoly(obj, xtraj, c, V)
            %
            % Compute the polynomial approximation of the drake plant
            %  xtraj and V must be a PPTrajectory, c must be an AffineSystem (Drake) 
            
            % Declare system model
            p = obj.drakeplant;
            
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
        
        function x = simulate(obj,x,u,deltat)
            %
            % simulate forward in time over a given time step in a sample-and-hold fashion with respect to the input.
            
            obj.mdlFun()
        end
        
    end
end