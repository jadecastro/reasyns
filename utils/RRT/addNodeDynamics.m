function [qNew,t,Xk,Uk,nearI] = addNodeDynamics(vBound,vObs1,vObs2,qGoal,node,gaussWeight,M,stepSize,Hout,n,isCyclic,modelType)
% ADDNODE: Adds a node according to a Gaussian sampling scheme which is
% stepSize away from the nearest point in the tree.

maxTrials = 100;
weightDist = ones(1,n) + 10*isCyclic';

isect = false;

[qNew,t,Xk,Uk,nearI] = deal([]);

trialMat = repmat([0 2*pi -2*pi],n,1).*repmat(isCyclic,1,3);

for i = 1:maxTrials %while isect
%     i
    % Bias the sampling toward qGoal
%     if rand < gaussWeight
        qRand = randn(1,n)*chol(M) + qGoal;
%     else
%         qRand = min(vBound) + (max(vBound) - min(vBound)).*rand(1,n);
%     end
    
    % Set qNear as the nearest coordinate in V (in configuration space) to qRand
    for j = 1:size(node,1) 
        if any(isCyclic)
            qRandMat = repmat(qRand',1,3)+trialMat;
            for k = 1:3
                distTmp(k) = norm(weightDist.*(qRandMat(:,k)' - node(j,:)));
            end
            dist(j) = min(distTmp);
        else
            dist(j) = norm(weightDist.*(qRand - node(j,:)));
        end
    end
    nearI = find(dist == min(dist),1);
    qNear = node(nearI,:);  
    
    % Set qNewStar as a stepSize away from qNear in the direction of qRand
    %qNewStar = qNear + stepSize/norm(qRand - qNear)*(qRand - qNear);
    
    % Construct motion plan between qNear and qNewStar; find qNew
%     [t,Xk,Uk] = genNominalTrajectory(qNear,[qNear(1:2); qNewStar(1:2)],Hout,n,isCyclic,modelType);
    [t,Xk,Uk] = genNominalTrajectory(qNear,[qNear(1:2); qRand(1:2)],Hout,n,isCyclic,modelType);
    
    if length(t) > 1
        
        clear xknrm
        xknrm(1) = 0;
        for k = 1:length(t)-1, xknrm(k+1) = norm(Xk(k+1,1:2) - Xk(k,1:2)); end
        indxRem = (cumsum(xknrm) > stepSize);
        t(indxRem) = [];
        Xk(indxRem,:) = [];
        Uk(indxRem,:) = [];
        
        qNew = Xk(end,:);
        
        % Test whether the path from qNear to qNew is in free space
        % isect = double(any(qNew < min(vBound)+radius) || any(qNew > max(vBound)-radius));
        isect = checkIntersection3(vBound,vObs1,vObs2,[],[],[],[],Xk,Hout,n,isCyclic,'finalPt');
        
        if ~isect, break, end
    end
end

if i == maxTrials, 
    disp('Cannot add node to tree')
end