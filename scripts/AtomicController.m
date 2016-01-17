
classdef AtomicController
    
    properties
        % system model
        sys 
        
        % base point trajectories
        x0
        u0
        
        % controller 
        K;
    end
    
    properties(SetAccess = private)
        abs_tol = 1e-6;
    end
    
    methods
        function obj = AtomicController(x0,u0,K,sys)
            % Constructor
            
            if ~isa(x0,'Traject') || ~isa(u0,'Traject') || ~isa(K,'Traject')
                error('One or more of the arguments is not of type Traject!')
            end
            t = getTimeVec(x0);
            tlen = length(t);
            if tlen ~= length(getTimeVec(x0)) || tlen ~= length(getTimeVec(u0)) || tlen ~= length(getTimeVec(K))
                warning('One or more time vectors have dimensions that don''t match')
            end
            
            obj.x0 = x0;
            obj.u0 = u0;
            obj.K = K;
            
            obj.sys = sys;
        end
        
%         function display(obj)
%             % 
%             fprintf('Atomic Controller with parameters\n');
%             
%             [m, n] = size(obj);
%             if (m > 1) | (n > 1)
%                 fprintf('atomic controller consists of %dx%d array of ellipsoids.\n', m, n);
%                 return;
%             end
%             
%             fprintf('x0:\n'); disp(obj.x0);
%             fprintf('u0:\n'); disp(obj.u0);
%             fprintf('K:\n'); disp(obj.K);
%             fprintf('sys:\n'); disp(obj.sys);
%             
%             if isempty(obj)
%                 fprintf('Empty atomic controller.\n');
%                 return;
%             end
%             
%             return;
%         end
        
        function ac = normalizerho(ac)
            %
            error('normalizerho: Expecting an overloaded function.')
        end
        
        function Eproj = projection(ac)
            %
            error('projection: Expecting an overloaded function.')
        end
        
%         function ac = interp(ac,ts)
%             %
%             for i = 1:Nsteps;
%                 if any(ac.isCyclic)
%                     ith(i) = 1;
%                     if ac(i).x0(nn) > pi
%                         ith(i) = 2;
%                     elseif ac(i).x0(nn) < -pi
%                         ith(i) = 3;
%                     end
%                 else
%                     ith(i) = 1;
%                 end
%                 
%                 indx = find(min(abs(ac(i).t - ts{ith(i)})) == abs(ac(i).t - ts{ith(i)}),1,'first');
%                 ac(i).K = ac(i).K{ith(i)}(t(i));
%                 ac(i).x0 = ac(i).x0(i,:) + isCyclic'*thSgn(ith(i))*2*pi;
%             end
%         end
        
        function union
        end
        
        function intersect
        end
        
        function res = isinside(ac, X, s)
            %
            error('isinternal: Expecting an overloaded function.')
        end
        
%         function sample
%         end
        
        function acres = merge(acobj1,acobj2)
            %
            xn = timecat(acobj1.x0,acobj2.x0);
            un = timecat(acobj1.u0,acobj2.u0);
            Kn = timecat(acobj1.K,acobj2.K);
            acres = AtomicController(xn,un,Kn,acobj1.sys);
        end
        
        function plot(ac,fignum,color)
            %
            %TODO: treat both 2d and 3d cases
            Eproj = projection(ac,sys);
            figure(fignum)
            plot(Eproj)
        end
        
%         function subsasgn
%             %
%             error('wip')
%         end
        
    end
    
end