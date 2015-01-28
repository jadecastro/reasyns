function [path,i] = buildRRT3(vBound,vObs1,vObs2,ellBound1,ellBound2,ellC1,ellC2,qInit,qGoal,stepSize,n,radius)
% BUILDRRT: Computes a path from qInit to qGoal using RRT assuming a
% circular robot.
%
%   [node,edge,xyNode] = buildRRT(vObs,vBound,qInit,qGoal,n,nNbrs,connectRad,type,radius)
% 
%   INPUTS
%       vObs        cell array containing n-by-2 arrays of vertices defining each obstacle
%       vBound      4-by-2 array of vertices defining the boundary
%       qInit       initial point in Cartesian space (2-by-1)
%       qGoal       goal point in Cartesian space (2-by-1)
%       stepSize    step size
%       n           maximum number of nodes to add to tree
%       radius      robot radius (m)
%
%   OUTPUTS
%       path        vector defining the waypoints from qInit to qGoal
% 


% Initialize the set V
node = qInit;
edge = [];            
path.node = [];
path.q = [];

% Choose random points
gaussWeight = 0.99;  % Weight [0-1] on Gaussian sampling biasing at qGoal
M = 0.0002*eye(2);    % Covariance matrix

for i = 1:n
    % Plotting
    plotNewNode(node,edge)
    
    % Test whether any elements in V are line-of-sight with qGoal
    qGoal(3) = atan2(qGoal(2)-node(end,2),qGoal(1)-node(end,1));  % Theta determined by angle of connecting line segment
    Npts = 10;
    for k = 1:Npts, qTest(k,:) = k/Npts*node(end,:) + (1-k/Npts)*qGoal; end
    isect = checkIntersection3(vBound,vObs1,vObs2,ellBound1,ellBound2,ellC1,ellC2,qTest);

%     for alpha = 0:0.01:1
%         qTest = alpha*node(end,:) + (1-alpha)*qGoal;
%         for mm = 1:length(vObs)
%             for j = 1:size(vObs{mm},1)-1
%                 a = (vObs{mm}(j+1,1)-vObs{mm}(j,1))^2 + (vObs{mm}(j+1,2)-vObs{mm}(j,2))^2;
%                 b = -2*((qTest(1)-vObs{mm}(j,1))*(vObs{mm}(j+1,1)-vObs{mm}(j,1)) + (qTest(2)-vObs{mm}(j,2))*(vObs{mm}(j+1,2)-vObs{mm}(j,2)));
%                 c = qTest(1)^2 + qTest(2)^2 + vObs{mm}(j,1)^2 + vObs{mm}(j,2)^2 - 2*(qTest(1)*vObs{mm}(j,1) + qTest(2)*vObs{mm}(j,2)) - radius^2;
%                 u = -b/2/a;
%                 if (b^2-4*a*c >= 0 && u >= 0 && u <= 1) || (norm(qTest - vObs{mm}(j,:)) < radius) || (norm(qTest - vObs{mm}(j+1,:)) < radius), isect = isect+1; end
%                 
%             end
%         end
%     end
    
    if ~isect,
        node = [node; qGoal];
        newI = size(node,1);
        edge = [edge; [i newI]];
        
        plotNewNode(node,edge)
        
        % Move backward from qGoal to qInit via edges to find the path
        nodeI = edge(end,2);
        path.node = nodeI;
        path.q = qGoal;
        while path.node(1) ~= 1
            path.node = [edge(nodeI-1,1); path.node];
            path.q = [node(edge(nodeI-1,1),:); path.q];
            nodeI = edge(nodeI-1,1);
        end
        return
    end
    
    % If not, go fish
    [qNew,nearI] = addNode3(vBound,vObs1,vObs2,ellBound1,ellBound2,ellC1,ellC2,qGoal(1:2),node,gaussWeight,M,stepSize,radius);
    
    % Update the set of nodes and edges
    node = [node; qNew];
    newI = size(node,1);
    edge = [edge; [nearI newI]];
    
end

end

function plotNewNode(node,edge)
% Plot the new node in the current figure

figure(3)
hold on
plot(node(end,1),node(end,2),'r.')
if ~isempty(edge), plot([node(edge(end,1),1) node(edge(end,2),1)],[node(edge(end,1),2) node(edge(end,2),2)],'r'), end
drawnow

end
