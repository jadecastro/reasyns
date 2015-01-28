function [V_opt, policy_opt] = solveRL(paramSet, model);

numAct = size(model.P,3);

V(1:model.stateCount) = 0;

while norm(delta) > 1;
    for i = 1:length(model.R)
        for j = 1:numAct
            tmp(j) = sum(model.P(i,:,j));
        end
        Vnew = model.R(i) + model.gamma*max(tmp*V(i));
        delta(i) = Vnew - V(i);
        V(i) = Vnew;
    end
end

V_opt = V;
policy_opt = [];
