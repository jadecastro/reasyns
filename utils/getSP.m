function [policy,states] = getSP(paramSet,model)

numAct = size(model.P,3);

s = randi(model.stateCount);

policy = [];
states = s;
while s ~= model.stateCount
    a = randi(numAct);
    s = getNextState(paramSet,stateCount,s,a);
    policy = [policy; a];
    states = [states; s];
end

