
classdef PolytopeAC < AtomicController
    
    properties
        % polytope properties
        vert;
    end
    
    methods
        function obj = PolytopeAC(x0,u0,K,vert,sys)
            % Constructor
            
            obj = obj@AtomicController(x0,u0,K,sys);
            t = getTimeVec(obj.x0);
            tlen = length(t);
            if length(vert) ~= tlen
                error('Length of vert must match that of t.')
            end            
            if ~isa(vert,'Traject')
                error('One or more of the arguments is not of type Traject!')
            end
            
            obj.vert = vert;
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
        
        function Pproj = projection(obj,sys)
            %
            n = obj.x0.pp.dim;
            if ~isempty(sys.H)
                nproj = length(sys.H);
                if ~isequal(sys.H,eye(nproj))
                    error('PROJECTION: arbitrary basis vectors not yet supported...')
                end
                for i = 1:length(obj)
                    p = polytope(obj(i));
                    for j = 1:length(p)
                        % TODO: extend to arbitrary basis vectors rather than just the existing frame
                        Pproj(j,i) = projection(p(j),1:nproj);
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

        function p = polytope(obj)
            % Output an array of polytope objects based on parameters from
            % the atomic controller
            for i = 1:size(obj,1)
                for j = 1:size(obj,2)
                    t = getTimeVec(obj(i,j).x0);
                    vert = ppval(obj(i,j).vert.pp,t);
                    for k = 1:size(vert,3)
                        p(k) = polytope(vert(:,:,k));
                    end
                end
            end
        end
        
        function union
        end
        
        function intersect
        end
        
        function res = isinternal(obj, X, s, varargin)
            %
            [res1,res2] = isinside(polytope(obj), X);
            if s == 'u'
                res = isequal(res2,1:length(obj));
            elseif s == 'i'
                res = res1;
            else
                error('ISINTERNAL: unrecognized setting for s.')
            end
            
            if ~isempty(varargin)
                sys = varargin{1};
                %NB: not complete.
                if ~res
                    for i = find(sys.params.isCyclic)'
                        X(i) = X(i) + 2*pi;
                        [res1,res2] = isinside(polytope(obj), X);
                        if s == 'u'
                            res = isequal(res2,1:length(obj));
                        elseif s == 'i'
                            res = res1;
                        end
                        X(i) = X(i) - 2*pi;
                    end
                end
                if ~res
                    for i = find(sys.params.isCyclic)'
                        X(i) = X(i) - 2*pi;
                        [res1,res2] = isinside(polytope(obj), X);
                        if s == 'u'
                            res = isequal(res2,1:length(obj));
                        elseif s == 'i'
                            res = res1;
                        end
                        X(i) = X(i) + 2*pi;
                    end
                end
            end
                
        end
        
        function res = isinside(obj,regobj,sys)
            %
            reg = regobj.p;
            for i = 1:length(obj)
                res(i) = (reg == intersect(polytope(obj(i)),reg));
                
            end
            res = ~any(res);
        end
        
%         function sample
%         end
        
        function obj = append(obj)
            error('wip')
        end
        
        function obj = merge(obj,acobj2)
            %
            actmp = merge@AtomicController(obj,acobj2);
            vert_n = timecat(obj.vert,acobj2.vert);
            obj = QuadraticAC(actmp.x0,actmp.u0,actmp.K,vert_n,actmp.sys);
        end
        
        function plot(obj,sys,fignum,threedflag,color)
            %
            %TODO: treat both 2d and 3d cases
            
            if exist('threedflag')
                error('wip')
            else
                Pproj = projection(obj,sys);
                figure(fignum)
                plot(Pproj,'b')
            end
            drawnow
        end
        
        function disp(obj)
            fprintf('Polytope Atomic Controller with parameters\n');
        end
        
    end
    
end
