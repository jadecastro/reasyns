function [path] = buildReachabilityRRT(vBound,vObs1,vObs2,ellBndInv11,ellBndInv21,ellCInv11,ellCInv21,qInit,qGoal,sys,reg,ac,options)
% BUILDREACHABILITYRRT: Computes a path from qInit to qGoal using RRT assuming a circular robot. 
% Implements the reachability-based rejection sampling technique of Shkolnik, Walter, and Tedrake, '09 

debug = true;

if isstruct(reg)
    regUnion = reg.union;
    regGoal = reg.goal;
else
    regUnion = reg;
    regGoal = reg;
end

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

% initialize the reachability tree. NB: no collision checking here 
% format: [(id of the parent node on the RRT)  (q of the child node in the reachability tree)]
if sys.sysparams.m > 1
    warning('m>1 is not yet supported; for now, the RRT will only build a reachability tree for the first input.')
end

qReach = computeNodesForAllMaximalActions(sys,qInit,options.pathLengthRRT,options);
nodeReach = [];
newNodeCount = size(qReach,1);
for i = 1:newNodeCount
    nodeReach = [nodeReach; [1 qReach(i,:)]];
end

%TODO: should we add additional goal containment checks with the reachability tree? note: this could seriously hurt the number of necessary computations! 

%isect1 = true;
isect2 = true;
inGoalSet = false;

for j = 1:length(ellCInv11)
    tmpEllCInv11{j,1} = ellCInv11{j};
end

figure(3)
plot(regUnion)
hold on
axis equal

QinvBallTest = inv(sys.sysparams.Qf);

for i = 1:options.maxNodes
    disp(['RRT iteration: ',num2str(i)]);
    
    %     plotNewNode(node,edge,'k') % uncomment to plot- will slow things down
    plotNewReachNode(sys,node,nodeReach,newNodeCount,'g') % uncomment to plot- will slow things down
    
    % Test whether any elements in V are in the goal region (TODO: for now-- must be at least two points in the vector)
    if i > 1 && edge(end,1) ~= 1  % want at least 1 segment!
        %isect1 = [];
        isect2 = [];
        inGoalSet = false;
        tsum = 0;
        for ii = 1:length(nodeXU.t), tsum = tsum + length(nodeXU.t{ii}); end
        tmp = tsum - length(nodeXU.t{i-1}) + 1;
        for k = 1:options.sampleSkipColl:length(t)
            %isect1(k) = checkIntersection(vBound,vObs1,vObs2,ellBndInv11,ellBndInv21,ellCInv11,ellCInv21,Xk(k,:),Hout,n,isCyclic,'allPts',ac,reg,sys);
            %             isect2(k) = checkIntersection4(ellBndInv11,Xk(k,:),Hout,n,isCyclic);
            isect2(k) = true;
         
            for acidx = 1:length(ac)
                isect2(k) = isect2(k) & isinternal(ac(acidx),Xk(k,:)','u');
                if debug && norm(Xk(k,1:2)' - qGoal(1:2)') < sys.sysparams.closeEnough
                    isect2(k) = true;
                end
            end
            %             if isect2(k)
            %                 disp('found an intersection with the goal set!')
            %                 keyboard
            %             else
            %                 disp('no intersection.')
            %             end
            
            %isectBnd = checkIntersection4(ellBndInv11,Xk,Hout,n,isCyclic)  % isect == true --> for all m,n there is some n for which the point is outside for all m, where ell{m,n}.
            %isectC = checkIntersection4(tmpEllCInv11,Xk,Hout,n,isCyclic)  % isect == true --> for all m,n there is some n for which the point is outside for all m, where ell{m,n}.
            %isect = isectBnd && isectC;
            % TODO: which is point in the curve which we consider the last one inside?
            
            %             if ~isect1(k) && ~isect2(k) && k ~= length(t) && k > 10
            if isect2(k) % && k > 10 && k ~= length(t) 
                
                ballTest = ellipsoid(Xk(k,:)',QinvBallTest);  % construct a ball representing the expected funnel level set at the final point along the trajectory
                
                if ~isempty(ac)
                    acceptCriterion = ac.funnelContainsEllipsoid(sys,ballTest,100);
                elseif ~isempty(regGoal)
                    % acceptCriterion = regGoal.isinside(sys,Xk(k,:)');
                    acceptCriterion = regGoal.regionContainsEllipsoid(sys,ballTest);
                else
                    error('Unhandled exception. acceptCriterion must be defined by either supplying a goal atomic controller or a goal region.')
                end
                
                if acceptCriterion
                    inGoalSet = true;
                    nodeXU.t{i-1}(k+1:end) = [];
                    nodeXU.x{i-1}(k+1:end,:) = [];
                    nodeXU.u{i-1}(k+1:end,:) = [];
                    break
                end
            end
        end
    end
    %     isect = all(isect1 & isect2);

    %     if ~isect,
    if inGoalSet,
        %         plotNewNode(node,edge,'r')
        
        % Move backward from qGoal to qInit via edges to find the path
        path.node = [];
        path.q = [];
        path.t = [];
        path.x = [];
        path.u = [];
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
        if length(path.node) > 3
            disp('terminating RRT: intersection with goal set.')
            return
        else
%             disp('nevermind... continuing.')
        end
    end
    
    % If not, go fish
    [qNew,t,Xk,Uk,nearI,nodeReach] = addNodeDynamicsReachability(vBound,vObs1,vObs2,qGoal,node,nodeReach, ...
        options.gaussWeight,options.M,options.pathLengthRRT,ac,regUnion,sys,options);
    if isempty(t) || isempty(qNew)
%         disp('add node failed.')
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
    
    % update the reachability tree with the id of the RRT parent and qNew
    qReach = computeNodesForAllMaximalActions(sys,qNew,options.pathLengthRRT,options);
    newNodeCount = 0;
    for j = 1:size(qReach,1)
        if ~checkIntersection(vBound,vObs1,vObs2,[],[],[],[],qReach,'finalPt',[],regUnion,sys);
            newNodeCount = newNodeCount + 1;
            nodeReach = [nodeReach; [newI qReach(j,:)]];
        end
    end
    if isempty(nodeReach)
%         disp('No nodes left in the reachability graph! Cannot expand the tree.')
        break
    end
end
disp('failed to construct RRT.')

end

function qNode = computeNodesForAllMaximalActions(sys,qNear,stepSize,options)

nInputs = length(sys.umin);
if nInputs ~= length(sys.umax), error('the properties umin and umax should be of the same size.'); end

[~,bVec] = buildBinaryCube([],[],nInputs);

qNode = [];
for i = 1:size(bVec,1)
    for j = 1:nInputs
        U0(j) = 0;
        if bVec(i,j) == -1
            U0(j) = sys.umin(j);
        elseif bVec(i,j) == 1
            U0(j) = sys.umax(j);
        end
    end
    
    xNew = computeOpenLoopTrajectory(sys,qNear,U0,stepSize,options);
    t = getTimeVec(xNew);
    qNode = [qNode; double(xNew,t(end))'];
end

end

function plotNewNode(node,edge,color)
% Plot the new node in the current figure

figure(3)
hold on
plot(node(end,1),node(end,2),[color,'.'])
if ~isempty(edge), plot([node(edge(end,1),1) node(edge(end,2),1)],[node(edge(end,1),2) node(edge(end,2),2)],color), end
drawnow

end

function plotNewReachNode(sys,node,nodeReach,newNodeCount,color)
% Plot the two most recent reachability nodes in the current figure (NB: does not update the display due to pruning of the existing tree!) 

figure(3)
hold on

for i = 1:newNodeCount
    q1 = sys.state2SEconfig([],node(nodeReach(end-i+1,1),:),[]);
    q2 = sys.state2SEconfig([],nodeReach(end-i+1,2:end),[]);
    plot(q2(1),q2(2),[color,'o'])
    plot([q1(1) q2(1)],[q1(2) q2(2)],color)
end
drawnow

end
