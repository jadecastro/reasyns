
classdef polytopeAC < atomiccontroller
    
    properties
        % polytope properties
        vert;
    end
    
    methods
        function obj = polytopeAC(x0,u0,K,vert,sys)
            % Constructor
            
            obj = obj@atomiccontroller(x0,u0,K,sys);
            t = getTimeVec(obj.x0);
            tlen = length(t);
            if length(vert) ~= tlen
                error('Length of vert must match that of t.')
            end            
            if ~isa(vert,'traject')
                error('One or more of the arguments is not of type traject!')
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
        
        function Pproj = projection(ac,sys)
            %
            n = ac.x0.pp.dim;
            if ~isempty(sys.H)
                nproj = length(sys.H);
                if ~isequal(sys.H,eye(nproj))
                    error('PROJECTION: arbitrary basis vectors not yet supported...')
                end
                for i = 1:length(ac)
                    p = polytope(ac(i));
                    for j = 1:length(p)
                        % TODO: extend to arbitrary basis vectors rather than just the existing frame
                        Pproj(j,i) = projection(p(j),1:nproj);
                    end
                end
            else
                error('H not defined. Only linear state-output relationships supported.')
            end
        end
        
%         function newobj = interp(ac,newobj,ts)
%             
%             newobj = ac.interp(ac,newobj,ts);
%             
%             for i = 1:Nsteps;
%                 if any(ac.isCyclic)
%                     ith(i) = 1;
%                     if ac(i).x(nn) > pi
%                         ith(i) = 2;
%                     elseif ac(i).x(nn) < -pi
%                         ith(i) = 3;
%                     end
%                 else
%                     ith(i) = 1;
%                 end
%                 
%                 indx = find(min(abs(ac(i).t - ts{ith(i)})) == abs(ac(i).t - ts{ith(i)}),1,'first');
%                 ac(i).P = ac(i).P{ith(i)};
%             end
%         end

        function p = polytope(ac)
            % Output an array of polytope objects based on parameters from
            % the atomic controller
            for i = 1:size(ac,1)
                for j = 1:size(ac,2)
                    t = getTimeVec(ac(i,j).x0);
                    vert = ppval(ac(i,j).vert.pp,t);
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
        
        function res = isinternal(acobj, X, s, varargin)
            %
            [res1,res2] = isinside(polytope(acobj), X);
            if s == 'u'
                res = isequal(res2,1:length(acobj));
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
                        [res1,res2] = isinside(polytope(acobj), X);
                        if s == 'u'
                            res = isequal(res2,1:length(acobj));
                        elseif s == 'i'
                            res = res1;
                        end
                        X(i) = X(i) - 2*pi;
                    end
                end
                if ~res
                    for i = find(sys.params.isCyclic)'
                        X(i) = X(i) - 2*pi;
                        [res1,res2] = isinside(polytope(acobj), X);
                        if s == 'u'
                            res = isequal(res2,1:length(acobj));
                        elseif s == 'i'
                            res = res1;
                        end
                        X(i) = X(i) + 2*pi;
                    end
                end
            end
                
        end
        
        function res = isinside(acobj,regobj,sys)
            %
            reg = regobj.p;
            for i = 1:length(acobj)
                res(i) = (reg == intersect(polytope(acobj(i)),reg));
                
            end
            res = ~any(res);
        end
        
%         function sample
%         end
        
        function ac = append(ac)
            error('wip')
        end
        
        function acres = merge(acobj1,acobj2)
            %
            actmp = merge@atomiccontroller(acobj1,acobj2);
            vert_n = timecat(acobj1.vert,acobj2.vert);
            acres = quadraticAC(actmp.x0,actmp.u0,actmp.K,vert_n,actmp.sys);
        end
        
        function plot(ac,sys,fignum,threedflag,color)
            %
            %TODO: treat both 2d and 3d cases
            
            if exist('threedflag')
                error('wip')
            else
                Pproj = projection(ac,sys);
                figure(fignum)
                plot(Pproj,'b')
            end
            drawnow
        end
        
        function disp(ac)
            fprintf('Polytope Atomic Controller with parameters\n');
        end
        
    end
    
end