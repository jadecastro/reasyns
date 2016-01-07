function [path] = buildCtrlSpaceRRT(vBound,vObs1,vObs2,ellBndInv11,ellBndInv21,ellCInv11,ellCInv21,qInit,qGoal,stepSize,sys,reg,ac,options)
% BUILDRRT: Computes a path from qInit to qGoal using RRT assuming a
% circular robot.

maxNodes = 500;
sampleSkipColl = 2;

Hout = sys.params.H;
n = sys.params.n;
isCyclic = sys.params.isCyclic;

% Initialize the set V
node = qInit;
Xk = qInit;
t = 0;
nodeXU = [];
edge = [];            
path.node = [];
path.q = [];
path.t = [];
path.x = [];
path.u = [];

% Choose random points
gaussWeight = 0.9;  % Weight [0-1] on Gaussian sampling biasing at qGoal
% M = 0.1*eye(n);    % Covariance matrix
M = diag([0.1, 0.1, 1]);    % Covariance matrix

%isect1 = true;
isect2 = true;
inGoalSet = false;

for j = 1:length(ellCInv11)
    tmpEllCInv11{j,1} = ellCInv11{j};
end

figure(3), clf
hold on
axis equal

for i = 1:maxNodes
    i
    plotNewNode(node,edge) % uncomment to plot- will slow things down
    
    % Test whether any elements in V are in the goal region (TODO: for now-- must be at least two points in the vector)
    if i > 1 && edge(end,1) ~= 1  % want at least 1 segment!
        %isect1 = [];
        isect2 = [];
        inGoalSet = false;
        tsum = 0;
        for ii = 1:length(nodeXU.t), tsum = tsum + length(nodeXU.t{ii}); end
        tmp = tsum - length(nodeXU.t{i-1}) + 1;
        for k = 1:sampleSkipColl:length(t)
            %isect1(k) = checkIntersection3(vBound,vObs1,vObs2,ellBndInv11,ellBndInv21,ellCInv11,ellCInv21,Xk(k,:),Hout,n,isCyclic,'allPts',ac,reg,sys);
%             isect2(k) = checkIntersection4(ellBndInv11,Xk(k,:),Hout,n,isCyclic);
            isect2(k) = isinternal(ac,Xk(k,:)','u');
            if isect2(k)
                disp('intersection with goal set!')
%                 keyboard
            else
                disp('no intersection.')
            end
            %isectBnd = checkIntersection4(ellBndInv11,Xk,Hout,n,isCyclic)  % isect == true --> for all m,n there is some n for which the point is outside for all m, where ell{m,n}.
            %isectC = checkIntersection4(tmpEllCInv11,Xk,Hout,n,isCyclic)  % isect == true --> for all m,n there is some n for which the point is outside for all m, where ell{m,n}.
            %isect = isectBnd && isectC;
            % TODO: which is point in the curve which we consider the last one inside?
%             if ~isect1(k) && ~isect2(k) && k ~= length(t) && k > 10
            if isect2(k) % && k ~= length(t) && k > 10
                inGoalSet = true;
                nodeXU.t{i-1}(k+1:end) = [];
                nodeXU.x{i-1}(k+1:end,:) = [];
                nodeXU.u{i-1}(k+1:end,:) = [];
                disp('intersection with goal set found.')
                k
                break,
            end
        end
    end
    %     isect = all(isect1 & isect2);

    %     if ~isect,
    if inGoalSet,
        disp('reached goal set')
        plotNewNode(node,edge)
        
        % Move backward from qGoal to qInit via edges to find the path
        nodeI = edge(end,2);
        while nodeI ~= 1
            tlast = nodeXU.t{nodeI-1}(end);
            path.node = [edge(nodeI-1,1); path.node];
            path.q = [node(edge(nodeI-1,1),:); path.q];
            path.t = [nodeXU.t{nodeI-1}(2:end); path.t+tlast];
            path.x = [nodeXU.x{nodeI-1}(2:end,:); path.x];
            path.u = [nodeXU.u{nodeI-1}(2:end,:); path.u];
            tmpt = nodeXU.t{nodeI-1}(1);
            tmpx = nodeXU.x{nodeI-1}(1,:);
            tmpu = nodeXU.u{nodeI-1}(1,:);
            nodeI = edge(nodeI-1,1);
        end
        path.t = [tmpt; path.t];
        path.x = [tmpx; path.x];
        path.u = [tmpu; path.u];
        length(path.t)
        return
    end
    
    % If not, go fish
    [qNew,t,Xk,Uk,nearI] = addNodeDynamics(vBound,vObs1,vObs2,qGoal,node,gaussWeight,M,stepSize,ac,reg,sys,options);
    if isempty(t) || isempty(qNew)
        disp('add node failed.')
        break
    end
    
    % Update the set of nodes and edges
%     qNew
    node = [node; qNew];
    nodeXU.t{i} = t;
    nodeXU.x{i} = Xk;
    nodeXU.u{i} = Uk;
    newI = size(node,1);
    edge = [edge; [nearI newI]];
end
disp('cannot construct RRT.')

end

function plotNewNode(node,edge)
% Plot the new node in the current figure

figure(3)
hold on
plot(node(end,1),node(end,2),'r.')
if ~isempty(edge), plot([node(edge(end,1),1) node(edge(end,2),1)],[node(edge(end,1),2) node(edge(end,2),2)],'r'), end
drawnow

end
