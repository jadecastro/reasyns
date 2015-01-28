classdef threeRegAutData < MapData
    
    properties(setAccess = private, getAccess = public)
        aut
        isReactive
    end
    
    methods
        function obj = threeRegAutData
            % Automaton definition
            aut.q{1} = 1;
            aut.q{2} = 2;
            aut.q{3} = 3;
            aut.q{4} = 2;
            aut.trans{1} = [1 2];
            aut.trans{2} = [2 3];
            aut.trans{3} = [2 1];
            aut.trans{4} = [3 4];
            aut.trans{5} = [4 1];
            %aut.trans{6} = [1 1];
            
            isReactive(1) = false;
            isReactive(2) = true;
            isReactive(3) = true;
            isReactive(4) = false;
            isReactive(5) = false;
        end
    end
end