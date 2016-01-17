
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
            if length(V) ~= tlen
                error('Length of V must match that of t.')
            end
            
            obj.rho = rho;
            obj.V = V;
        end
        
        function ac = normalizerho(ac)
            %
            error('wip.')
        end
        
        % TODO: generalize to arbitrary polynomials
        function Eproj = projection(ac)
            %
            error('Work in progress.')
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
%                 ac(i).rho = ac(i).rho{ith(i)};
%                 ac(i).V = ac(i).V{ith(i)}(t(i));
%             end
%         end
        
        function union
            %
            error('Work in progress.')
        end
        
        function intersect
            %
            error('Work in progress.')
        end
        
        function res = isinside(ac, X, s)
            %
            error('Work in progress.')
        end
        
%         function sample
%             %
%             error('Work in progress.')
%         end
        
        function append
            %
            error('Work in progress.')
        end
        
        function acres = merge(acobj1,acobj2)
            %
            actmp = merge@AtomicController(acobj1,acobj2);
            rhon = timecat(acobj1.rho,acobj2.rho);
            t = getTimeVec(actmp.x0);  t1 = getTimeVec(acobj1.x0);  t2 = getTimeVec(acobj2.x0);
            [~,indx1] = intersect(t1,t);
            [~,indx2] = intersect(t2,setdiff(t,t1));
            Vn = [acobj1.V(indx1) acobj2.V(indx2)];
            acres = PolynomialAC(actmp.x0,actmp.u0,actmp.K,rhon,Vn,actmp.sys);
        end
        
        function plot(ac,fignum,color)
            %
            %TODO: treat both 2d and 3d cases
            Eproj = projection(ac);
            figure(fignum)
            plot(Eproj)
        end
        
    end
    
end