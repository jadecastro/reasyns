
classdef systemdynamics
    
    properties
        n;
        m;
        H;
        isCyclic;
        f;
    end
    
    methods
        function obj = systemdynamics
            
        end
        
        function y = state2config(x)
            y(1:length(obj),obj.m) = 0;
            if ~isempty(obj.H)
                for i = 1:length(obj)
                    y(i,:) = obj.H*x(i,:);
                end
            else
                error('H not defined. Only linear state-output relationships supported.')
            end
        end
        

        
        
    end
end