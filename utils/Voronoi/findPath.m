function [path] = findPath(qInit,qGoal,node,edge,xyNode)%,xyInit,xyGoal)
% FINDPATH: Outputs an optimal path in a polygonal environment using
% Djikstra's algorithm using a graph defined by a collection of contiguous
% convex hulls defined in vCell.
%
%   [path] = findPath(vCell,node,edge,xyNode,xyInit,xyGoal)
% 
%   INPUTS
%       vCell       p-dimensional cell array containing m-by-2 arrays of vertices defining the grid cell
%       node        p-by-1 vector defining the nodes 
%       edge        q-by-2 vector defining the edges (pairs of nodes)
%       xyNode      p-by-2 vector defining the coordinates of nodes
%       xyInit      initial point in Cartesian space (2-by-1)
%       xyGoal      goal point in Cartesian space (2-by-1)
%
%   OUTPUTS
%       path        vector defining the waypoints from xyInit to xyGoal
% 
% 
%   Cornell University
%   MAE 4180/5180 CS 3758: Autonomous Mobile Robots
%   Homework #7
%   DeCastro, Jonathan


for i = 1:length(edge)
    dist(i) = norm(xyNode(edge(i,1),:) - xyNode(edge(i,2),:));
end

% for i = 1:length(vCell)
%     if inpolygon(xyInit(1,1),xyInit(1,2),vCell{i}(:,1),vCell{i}(:,2)), qInit = i; end
%     if inpolygon(xyGoal(1,1),xyGoal(1,2),vCell{i}(:,1),vCell{i}(:,2)), qGoal = i; end
% end

% initialize costs
cost(1:length(node)) = 1e9;
cost(qInit) = 0;

q = qInit;

aliveIndx = node;
tentIndx = [];
deadIndx = [];
path = [];

% for count = 1:4 
while ~isempty(aliveIndx) && ismember(qGoal,aliveIndx)
    
    % Find edge costs to each successor
    if q == qGoal
        edgeIndx = find(any(ismember(edge,q),2));
    else
        edgeIndx = find(any(ismember(edge,q),2) & ~any(ismember(edge,deadIndx),2));
    end
    if ~isempty(tentIndx)
        for k = 1:length(tentIndx)-1
            if ~any(~any(ismember(edge(edgeIndx,:),deadIndx),2))  % If no transitions to alive nodes, pick next highest successor
                aliveIndx = setdiff(aliveIndx,q);  % update alive index
                deadIndx = setdiff(node,aliveIndx);
                tentIndx(1) = [];
                q = tentIndx(1);
                if q == qGoal
                    edgeIndx = find(any(ismember(edge,q),2));
                else
                    edgeIndx = find(any(ismember(edge,q),2) & ~any(ismember(edge,deadIndx),2));
                end
            else
                break
            end
        end
    end
        
    % Update costs
    clear nodeIndx
    for m = 1:length(edgeIndx), 
        nodeIndx(m) = edge(edgeIndx(m),edge(edgeIndx(m),:)~=q); 
    end
    if isempty(edgeIndx)
        keyboard
    end
    for j = 1:length(nodeIndx)
        if cost(nodeIndx(j)) > cost(q) + dist(edgeIndx(j));
            cost(nodeIndx(j)) = cost(q) + dist(edgeIndx(j));
        end
    end
    
    path = [path; q];
    
    % Prune finalized nodes from the graph
    aliveIndx = setdiff(aliveIndx,q);  % update alive index
    deadIndx = setdiff(node,aliveIndx);
    [costSort,iSort] = sort(cost(nodeIndx));
    tentIndx = nodeIndx(iSort);
    q = tentIndx(1);  % Choose lowest cost successor
end

