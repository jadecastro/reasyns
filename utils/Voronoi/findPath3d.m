function [path, edgePath, pathDist] = findPath3d(qInit,qGoal,node,edge,xyzNode)
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

xyNode = xyzNode(1:2);
for i = 1:length(edge)
    dist(i) = norm(xyNode(edge(i,1),:) - xyNode(edge(i,2),:));
end

% % Hard-coded example
% qInit = 1;
% qGoal = 8;
% node = [1;2;3;4;5;6;7;8];
% edge = [1 2;1 3;1 4;2 3;2 5;2 7;3 4;3 5;3 6;4 6;5 6;5 8;6 8;7 8];
% dist = [2;5;4;2;7;12;1;4;3;4;1;5;7;3];

% initialize costs
cost(1:length(node)) = 1e9;
cost(qInit) = 0;

q = qInit;

aliveIndx = node;
deadIndx = [];
edgePath = [];
edgeBest = [];
% for count = 1:4 
while ~isempty(aliveIndx) && ismember(qGoal,setdiff(aliveIndx,q))
    
    % Find edge costs to each successor
    if q == qGoal
        edgeIndx = find(any(ismember(edge,q),2));
    else
        edgeIndx = find(any(ismember(edge,q),2) & ~any(ismember(edge,deadIndx),2));
    end
        
    % Update costs
    clear nodeIndx
    for m = 1:length(edgeIndx), 
        nodeIndx(m) = edge(edgeIndx(m),edge(edgeIndx(m),:)~=q); 
    end
    for j = 1:length(nodeIndx)
        if cost(nodeIndx(j)) > cost(q) + dist(edgeIndx(j));
            cost(nodeIndx(j)) = cost(q) + dist(edgeIndx(j));
            edgeBest(any(ismember(edgeBest,nodeIndx(j)),2),:) = [];
            edgeBest = [edgeBest; q nodeIndx(j)];
        end
    end
    
    % Prune finalized nodes from the graph
    aliveIndx = setdiff(aliveIndx,q);  % update alive index
    deadIndx = [q; deadIndx];
    [~,iSort] = sort(cost(aliveIndx));
    
    q1 = aliveIndx(iSort(1));  % Choose lowest cost successor
    while isempty(find(any(ismember(edge,q1),2) & ~any(ismember(edge,deadIndx),2), 1))
        deadIndx(1) = [];
    end
    q = q1;
    tmp = edgeBest(any(ismember(edgeBest,q),2),:);
    edgePath = [edgePath; tmp];
end
edgePath = sort(edgePath,2);

% Extract the path
pathDist = 0;
path = qGoal;
while path(1) ~= qInit
    tmp = find(any(ismember(edgePath(:,:),path(1)),2));
    path = [edgePath(tmp(1),edgePath(tmp(1),:)~=path(1)); path];
    edgePath(tmp(1),:) = [];
    pathDist = pathDist + dist(any(edge==path(1),2) & any(edge==path(2),2));
end


