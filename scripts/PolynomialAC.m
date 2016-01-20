
classdef PolynomialAC < AtomicController
    
    properties
        rho;
        V;
    end
    
    methods
        function obj = PolynomialAC(x0,u0,K,rho,V,sys)
            % Constructor
            
            obj = obj@AtomicController(x0,u0,K,sys);
            if ~isa(rho,'Traject')% || ~isa(V,'Traject')
                error('One or more of the arguments is not of type Traject!')
            end
            t = getTimeVec(obj.x0);
            tlen = length(t);
            if tlen ~= length(getTimeVec(rho)) % || ~all(t == getTimeVec(V))
                warning('One or more time vectors have dimensions that don''t match')
            end
%             if length(V) ~= tlen
%                 error('Length of V must match that of t.')
%             end
            
            obj.rho = rho;
            obj.V = V;
        end
        
        function obj = normalizerho(obj)
            %
            error('wip.')
        end
        
        % TODO: generalize to arbitrary polynomials
        function Eproj = projection(obj)
            %
            error('Work in progress.')
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
%                 obj(i).rho = obj(i).rho{ith(i)};
%                 obj(i).V = obj(i).V{ith(i)}(t(i));
%             end
%         end
        
        
        function obj = merge(obj1,obj2)
            %
            actmp = merge@AtomicController(obj1,obj2);
            rhon = timecat(obj1.rho,obj2.rho);
            t = getTimeVec(actmp.x0);  t1 = getTimeVec(obj1.x0);  t2 = getTimeVec(obj2.x0);
            [~,indx1] = intersect(t1,t);
            [~,indx2] = intersect(t2,setdiff(t,t1));
            Vn = [obj1.V(indx1) obj2.V(indx2)];
            obj = PolynomialAC(actmp.x0,actmp.u0,actmp.K,rhon,Vn,actmp.sys);
        end
        
        function plot(obj,fignum,color)
            %
            %TODO: treat both 2d and 3d cases
            Eproj = projection(obj);
            figure(fignum)
            plot(Eproj)
        end
        
    end
    
end