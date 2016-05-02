
classdef AtomicController
    
    properties
        % system model
        sys ;
        
        % FSM state pair
        pre;
        post;
        flagTransitionSet = false;
        
        % base point trajectories
        x0;
        u0;
        
        % controller 
        K;
    end
    
    properties(SetAccess = private)
        abs_tol = 1e-6;
    end
    
    methods
        function obj = AtomicController(x0,u0,K,sys,varargin)
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
            
            if ~isempty(varargin)
                obj.setTransition(varargin);
            end
        end
        
        function obj = setTransition(obj,varargin)
            obj.pre = varargin{1};
                if length(varargin) > 1
                    obj.post = varargin{2};
                else
                    obj.post = obj.pre;
                end
            obj.flagTransitionSet = true;
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
        
        function obj = normalizerho(obj)
            %
            error('normalizerho: Expecting an overloaded function.')
        end
        
        function Eproj = projection(obj)
            %
            error('projection: Expecting an overloaded function.')
        end
        
%         function obj = interp(obj,ts)
%             %
%             for i = 1:Nsteps;
%                 if any(obj.isCyclic)
%                     ith(i) = 1;
%                     if obj(i).x0(nn) > pi
%                         ith(i) = 2;
%                     elseif obj(i).x0(nn) < -pi
%                         ith(i) = 3;
%                     end
%                 else
%                     ith(i) = 1;
%                 end
%                 
%                 indx = find(min(abs(obj(i).t - ts{ith(i)})) == abs(obj(i).t - ts{ith(i)}),1,'first');
%                 obj(i).K = obj(i).K{ith(i)}(t(i));
%                 obj(i).x0 = obj(i).x0(i,:) + isCyclic'*thSgn(ith(i))*2*pi;
%             end
%         end
        
        function union
        end
        
        function intersect
        end
        
        function res = isinside(obj, X, s)
            %
            error('isinternal: Expecting an overloaded function.')
        end
        
%         function sample
%         end
        
        function obj = merge(obj1,obj2)
            %
            xn = timecat(obj1.x0,obj2.x0);
            un = timecat(obj1.u0,obj2.u0);
            Kn = timecat(obj1.K,obj2.K);
            obj = AtomicController(xn,un,Kn,obj1.sys);
        end
        
        function plot(obj,fignum,color)
            %
            %TODO: treat both 2d and 3d cases
            Eproj = projection(obj,sys);
            figure(fignum)
            plot(Eproj)
        end
        
%         function subsasgn
%             %
%             error('wip')
%         end
        
    end
    
end