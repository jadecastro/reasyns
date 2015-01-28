function [node,edge,xyNode,Vx,Vy,X,Y] = buildVoronoi(Xreg,dX)
% node: n-by-1 vector of nodes in the Voronoi decomposition
% edge: m-by-2 array defining connectivity of nodes
% xyNode: n-by-2 array defining the locations of the nodes
% Vx: cell array of raw vertex x-component data
% Vx: cell array of raw vertex y-component data
% X: cell array of obstacle raw x-component data
% Y: cell array of obstacle raw y-component data

if iscell(Xreg)
    xreg = Xreg;
else
    xreg{1} = Xreg;
end

for i = 1:length(xreg)
    % Do for external edges
    v1{i} = discretizePolytope(xreg{i}.extVert',dX);
    indx = 0;
    for j = 1:size(v1{i},1)
        indx = indx+1;
        X{i}(indx) = v1{i}(j,1);
        Y{i}(indx) = v1{i}(j,2);
    end
    % Do for internal edges
    v1{i} = discretizePolytope(xreg{i}.intVert',dX);
    for j = 1:size(v1{i},1)
        indx = indx+1;
        X{i}(indx) = v1{i}(j,1);
        Y{i}(indx) = v1{i}(j,2);
    end
    
    [Vx{i},Vy{i}] = voronoi(X{i},Y{i});
    
    [Vx{i},Vy{i}] = removeIntersects(xreg{i},Vx{i},Vy{i});
    
    xyNode = [];
    node = [];
    nodeIdx = 0;
    for j = 1:size(Vx{i},2)
        tmpNode1 = [Vx{i}(1,j) Vy{i}(1,j)];
        if isempty(xyNode)
            xyNode1 = [inf inf];
        else
            xyNode1 = xyNode;
        end
        if ~any(all([abs(tmpNode1(1)-xyNode1(:,1))<=1e-6 abs(tmpNode1(2)-xyNode1(:,2))<=1e-6],2))
            xyNode = [xyNode; tmpNode1];
            nodeIdx = nodeIdx+1;
            node = [node; nodeIdx];
        end
        tmpNode2 = [Vx{i}(2,j) Vy{i}(2,j)];
        if ~any(all([abs(tmpNode2(1)-xyNode(:,1))<=1e-6 abs(tmpNode2(2)-xyNode(:,2))<=1e-6],2))
            xyNode = [xyNode; tmpNode2];
            nodeIdx = nodeIdx+1;
            node = [node; nodeIdx];
        end
        node1 = find(all([abs(tmpNode1(1)-xyNode(:,1))<=1e-6 abs(tmpNode1(2)-xyNode(:,2))<=1e-6],2));
        node2 = find(all([abs(tmpNode2(1)-xyNode(:,1))<=1e-6 abs(tmpNode2(2)-xyNode(:,2))<=1e-6],2));
        try
        edge(j,:) = [node1 node2];    
        catch
            keyboard
        end
    end
    
    figure(10)
    hold on
    plot(Vx{i},Vy{i},'-',X{i},Y{i},'.')
end

edge = edge(edge(:,1) ~= edge(:,2),:);
