function [qNew,nearI] = addNode(vObs,vBound,qGoal,node,gaussWeight,M,stepSize,radius)
% ADDNODE: Adds a node according to a Gaussian sampling scheme which is
% stepSize away from the nearest point in the tree.
%
%   [qNew,nearI] = addNode(vObs,vBound,qGoal,node,gaussWeight,M,stepSize,radius)
% 
%   INPUTS
%       vObs        cell array containing n-by-2 arrays of vertices defining each obstacle
%       vBound      4-by-2 array of vertices defining the boundary
%       qInit       initial point in Cartesian space (2-by-1)
%       qGoal       goal point in Cartesian space (2-by-1)
%       stepSize    step size
%       n           maximum number of nodes to add to tree
%
%   OUTPUTS
%       node        p-by-2 vector defining the node coordinates (node ID corresponds to the columns)
%       edge        q-by-2 vector defining the edges (pairs of node IDs)
% 

isect = false;
for i = 1:1000 %while isect
    % Bias the sampling toward qGoal
    if rand < gaussWeight
        qRand = randn(1,2)*chol(M) + qGoal;
    else
        qRand = (min(vBound)+radius) + (max(vBound) - min(vBound) - 2*radius).*rand(1,2);
    end
    
    % Set qNear as the nearest coordinate in V to qRand
    for j = 1:size(node,1), dist(j) = norm(qRand - node(j,:)); end
    nearI = find(dist == min(dist));
    qNear = node(nearI,:);
    
    % Set qNew as a stepSize away from qNear in the direction of qRand
    qNew = qNear + stepSize/norm(qRand - qNear)*(qRand - qNear);
    % Test whether the path from qNear to qNew is in free space
    isect = double(any(qNew < min(vBound)+radius) || any(qNew > max(vBound)-radius));
    for k = 1:100, qTest(k,:) = k/100*qNew + (1-k/100)*qNear; end
    for mm = 1:length(vObs), isect = isect + any(inpolygon(qTest(:,1),qTest(:,2),vObs{mm}(:,1),vObs{mm}(:,2))); end
    if radius == 0  % Checking the point criteria is enough
        if ~isect, break, end
    else  % Check if circle intersects any of the polygon line segments
        for alpha = 0:0.01:1
            qTest = alpha*qNew + (1-alpha)*qNear;
            for mm = 1:length(vObs)
                for j = 1:size(vObs{mm},1)-1
                    a = (vObs{mm}(j+1,1)-vObs{mm}(j,1))^2 + (vObs{mm}(j+1,2)-vObs{mm}(j,2))^2;
                    b = -2*((qTest(1)-vObs{mm}(j,1))*(vObs{mm}(j+1,1)-vObs{mm}(j,1)) + (qTest(2)-vObs{mm}(j,2))*(vObs{mm}(j+1,2)-vObs{mm}(j,2)));
                    c = qTest(1)^2 + qTest(2)^2 + vObs{mm}(j,1)^2 + vObs{mm}(j,2)^2 - 2*(qTest(1)*vObs{mm}(j,1) + qTest(2)*vObs{mm}(j,2)) - radius^2;
                    u = -b/2/a;
                    if (b^2-4*a*c >= 0 && u >= 0 && u <= 1) || (norm(qTest - vObs{mm}(j,:)) < radius) || (norm(qTest - vObs{mm}(j+1,:)) < radius), isect = isect+1; end
%                     for k = 1:100,
%                         qTestLine(k,:) = k/100*vObs{mm}(j,:) + (1-k/100)*vObs{mm}(j+1,:);
%                         if (norm(qTest - qTestLine(k,:)) < radius), isect = isect+1; end
%                     end
                end
            end
        end
        if ~isect, break, end
    end
end

if i == 1000, 
    i 
    warning('addNode: max number of iterations reached.')
    qNew = [];
    nearI = NaN;
end