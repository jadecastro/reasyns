classdef TwoDOFManipulatorPlant < Manipulator
    
    properties
        % parameters from Spong95 (except inertias are now relative to the
        % joints)
        % axis)
        l1 = 1; l2 = 2;
        m1 = 1; m2 = 1;
        g = 9.81;
        b1=.1;  b2=.1;
        %    b1=0; b2=0;
        lc1 = .5; lc2 = 1;
        Ic1 = .083;  Ic2 = .33;
        
        xG
        uG
        
        sysparams;
    end
    
    methods
        function obj = TwoDOFManipulatorPlant()
            obj = obj@Manipulator(2,1);
            obj = setInputLimits(obj,-1,1);
            
            obj = setInputFrame(obj,CoordinateFrame('AcrobotInput',1,'u',{'tau'}));
            obj = setStateFrame(obj,CoordinateFrame('AcrobotState',4,'x',{'theta1','theta2','theta1dot','theta2dot'}));
            obj = setOutputFrame(obj,obj.getStateFrame);
            
            obj.xG = Point(obj.getStateFrame,[pi;0;0;0]);
            obj.uG = Point(obj.getInputFrame,0);
            
            obj = setParamFrame(obj,CoordinateFrame('AcrobotParams',6,'p',...
                { 'b1','b2','lc1','lc2','Ic1','Ic2' }));
            %      obj = setParamFrame(obj,CoordinateFrame('AcrobotParams',10,'p',...
            %        { 'l1','l2','m1','m2','b1','b2','lc1','lc2','I1','I2' }));
            %      obj = setParamFrame(obj,CoordinateFrame('AcrobotParams',1,'p',...
            %        { 'lc2' }));
            obj = setParamLimits(obj,zeros(obj.getParamFrame.dim,1));
            
            obj.sysparams.n = getNumStates(obj); % number of states
            obj.sysparams.m = getNumInputs(obj); % number of inputs
            obj.sysparams.limsNonRegState = [-100 -100; 100 100];
            obj.sysparams.stateLimits = [-pi -100 -pi -100; pi 100 pi 100];
            obj.sysparams.isCyclic = [1 0 1 0]';
            obj.sysparams.isOutputLinear = 0;
            
        end
                
        function [H,C,B] = manipulatorDynamics(obj,q,qd)
            % keep it readable:
            m1=obj.m1; m2=obj.m2; l1=obj.l1; g=obj.g; lc1=obj.lc1; lc2=obj.lc2; b1=obj.b1; b2=obj.b2;
            I1 = obj.Ic1 + obj.m1*obj.lc1^2; I2 = obj.Ic2 + obj.m2*obj.lc2^2;
            m2l1lc2 = m2*l1*lc2;  % occurs often!
            
            c = cos(q(1:2,:));  s = sin(q(1:2,:));  s12 = sin(q(1,:)+q(2,:));
            
            h12 = I2 + m2l1lc2*c(2);
            H = [ I1 + I2 + m2*l1^2 + 2*m2l1lc2*c(2), h12; h12, I2 ];
            
            C = [ -2*m2l1lc2*s(2)*qd(2), -m2l1lc2*s(2)*qd(2); m2l1lc2*s(2)*qd(1), 0 ];
            G = g*[ m1*lc1*s(1) + m2*(l1*s(1)+lc2*s12); m2*lc2*s12 ];
            
            % accumate total C and add a damping term:
            C = C*qd + G + [b1;b2].*qd;
            
            B = [0; 1];
            % B = eye(2);
        end
        
        % todo: also implement sodynamics here so that I can keep the
        % vectorized version?
        
        function [f,df,d2f,d3f] = dynamics(obj,t,x,u)
            f = dynamics@Manipulator(obj,t,x,u);
            if (nargout>1)
                [df,d2f,d3f]= dynamicsGradientsTwoDOFManipulator(obj,t,x,u,nargout-1);
            end
        end
        
        function prepareModelGradients(obj)
            taylorOrder = 3;
            generateGradients('dynamics',taylorOrder,'dynamicsGradientsTwoDOFManipulator',obj,0,randn(obj.sysparams.n,1),0);
        end

        
        function x = getInitialState(obj)
            x = 0.1*randn(4,1);
        end

        function y = output(obj)
            y = x;
        end
        
        function y = state2SEconfig(obj,t,x,u)
            y(1) = obj.lc1*cos(x(1)) + obj.lc2*cos(x(3));
            y(2) = obj.lc1*sin(x(1)) + obj.lc2*sin(x(3));
            y = y';
        end
        
        function n = getNumPositions(obj)
            n = 2;
        end
        
        function n = getNumVelocities(obj)
            n = 2;
        end
        
    end
    
end
