function [edge,node,edge1,xyNode1,goalNode] = getEdges(Xinit,Xgoal,dX,node,edgeTmp,xyNode,ctrlType)

global t x

edgeTmpD = sort(edgeTmp,2,'Descend');

if strmatch(ctrlType,'stay')
    % This control mode has the robot visit all nodes with more than two edges, essentially traveling the "highway" parts of the roadmap
    parent = [];
    goalNode = [];
    for i = 1:size(node,1)
        if sum(any(edgeTmp == node(i),2)) > 2
            goalNode = [goalNode; i];
        end
    end
    goalNode = sort(goalNode);
    
    edge1 = edgeTmp;
    edge2 = edgeTmpD;
    xyNode1 = xyNode;
    
elseif strmatch(ctrlType,'leave')
    % Find verices in the goal region and add them to the roadmap
    ve = discretizePolytope(Xgoal.extVert',dX);
    vi = discretizePolytope(Xgoal.intVert',dX);
    verte = [];
    verti = [];
    % Do any of the external/internal faces of the initial region intersect an external face point of the goal region?
    for i = 1:size(ve,1)
        extext = double(subs(Xinit.ext,x,[ve(i,1);ve(i,2);0;0]));
        intext = double(subs(Xinit.int,x,[ve(i,1);ve(i,2);0;0]));
        if any(intext <= 0) || all(extext <= 0);
            verte = [verte; ve(i,:)];
        end
    end
    % Do any of the external/internal faces of the initial region intersect an internal face point of the goal region?
    for i = 1:size(vi,1)
        extint = double(subs(Xinit.ext,x,[vi(i,1);vi(i,2);0;0]));
        intint = double(subs(Xinit.int,x,[vi(i,1);vi(i,2);0;0]));
        if all(extint <= 0) || any(intint <= 0);
            verti = [verti; vi(i,:)];
        end
    end
    
    %% Find the shortest path from any node on the roadmap to the goal region.
    % Note: the last step in the path is selected using a greedy strategy.. to find the true shortest path, need to modify 'findPath'
    
    parente = []; parenti = [];
    goale = []; goali = [];
    if ~isempty(verte)
        for j = 1:size(verte,1)
            for i = 1:length(node)
                diste(i) = sum((xyNode(i,:) - verte(j,:)).^2);
            end
            indx = find(diste == min(diste),1,'first');
            parentTmp(j) = node(indx);
            distGoal(j) = diste(indx);
        end
        parente = unique(parentTmp);
        for i = 1:length(parente)
            indxSet = find(parente(i) == parentTmp);
            indxTmp = find(distGoal(indxSet) == min(distGoal(indxSet)),1,'first');
            goale(i) = indxSet(indxTmp(1));
        end
    end
    if ~isempty(verti)
        for j = 1:size(verti,1)
            for i = 1:length(node)
                disti(i) = sum((xyNode(i,:) - verti(j,:)).^2);
            end
            indx = find(disti == min(disti),1,'first');
            parentTmp(j) = node(indx);
            distGoal(j) = disti(indx);
        end
        parenti = unique(parentTmp);
        for i = 1:length(parenti)
            indxSet = find(parenti(i) == parentTmp);
            goali(i) = indxSet(distGoal(indxSet) == min(distGoal(indxSet)));
        end
    end
    goal = [goale';goali'+length(goale)];
    vert = [verte;verti];
    
    % Now, populate the edge structure
    parent = [parente';parenti'];
    child = node(end)+(1:length(parent))';
    node = [node; child];
    if all(node' ~= 1:length(node))
        error('node vector must be the countable sequence from 1 to length(node)')
    end
    
    edge1 = [edgeTmp; parent child];
    edge2 = [edgeTmpD; child parent];
    xyNode1 = [xyNode;vert(goal,:)];
    goalNode = child;
end

xyzNode = [xyNode1 nan(size(xyNode1),1)];
for i = 1:size(edge1,1)
    edge(2*i-1).val = edge1(i,:);
    edge(2*i).val = edge2(i,:);
    edge(2*i-1).barrier = [];
    edge(2*i).barrier = [];
    edge(2*i-1).barrier = [];
    edge(2*i).barrier = [];
    edge(2*i-1).xyz = [xyzNode(edge1(i,1),:); xyzNode(edge1(i,2),:)];
    edge(2*i).xyz = [xyzNode(edge2(i,1),:); xyzNode(edge2(i,2),:)];
end
