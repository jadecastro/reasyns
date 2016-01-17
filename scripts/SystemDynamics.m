
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
        
        function y = state2config(sysobj,x)
            %
            y(1:length(sysobj),sysobj.params.m) = 0;
            if ~isempty(sysobj.H)
                for i = 1:length(sysobj)
                    y(i,:) = sysobj.H*x(i,:);
                end
            else
                error('H not defined. Only linear state-output relationships supported.')
            end
        end
        
        function simulate(sysobj)
            %
            
        end
        
    end
end