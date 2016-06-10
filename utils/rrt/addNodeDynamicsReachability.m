function [qNew,t,Xk,Uk,nearI,nodeReach] = addNodeDynamicsReachability(vBound,vObs1,vObs2,qGoal,node,nodeReach,gaussWeight,M,stepSize,ac,reg,sys,options)
% ADDNODE: Adds a node according to a Gaussian sampling scheme which is
% stepSize away from the nearest point in the tree.

n = sys.sysparams.n;
isCyclic = sys.sysparams.isCyclic;

maxTrials = 200;
weightDist = ones(1,n) + 0.01*isCyclic';

isect = false;

[qNew,t,Xk,Uk,nearI,nearReachI] = deal([]);

trialMat = repmat([0 2*pi -2*pi],n,1).*repmat(isCyclic,1,3);

for i = 1:maxTrials %while isect
    try
        if mod(i,10) == 0, disp(['  RRT sample iteration: ',num2str(i)]); end
        
        % Bias the sampling toward qGoal
        if rand < gaussWeight
            qRand = randn(1,n)*chol(M) + qGoal;
        else
            qRand = [min(vBound) -pi] + [(max(vBound) - min(vBound)) 2*pi].*rand(1,n);
        end
        
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
        for j1 = 1:size(nodeReach,1)
            if any(isCyclic)
                qRandMat = repmat(qRand',1,3)+trialMat;
                for k = 1:3
                    distTmp(k) = norm(weightDist.*(qRandMat(:,k)' - nodeReach(j1,2:end)));
                end
                distReach(j1) = min(distTmp);
            else
                distReach(j1) = norm(weightDist.*(qRand - nodeReach(j1,2:end)));
            end
        end
        [minDist,~] = min(dist);
        [minDistReach,nearReachI] = min(distReach);
        if minDistReach <= minDist  % if the random point is closer to the reachability tree than the RRT, accept it!
            
            % set qNear to the parent of the node on the reachability tree. NB: we defer pruning the reachability tree to the end, i.e. after we've verified collision-freeness 
            qNear = node(nodeReach(nearReachI,1),:);
            
            % Construct motion plan between qNear and qNewStar; find qNew
            
            % Below, we supply two options: 
            %   (1) if an appropriate steering control law is available,
            %   then we compute that controller via computeTrajectory.
            %   (2) if no steering control law, then form a cube over all
            %   control limits (kept constant over the specified time
            %   interval), then choose the combination that ultimately
            %   minimizes the distance to the desired sample.
            
            if (ismethod(sys,'dynamicsWaypointSteering') && ismethod(sys,'steerToXYWaypoints'))
                [u0,x0] = computeTrajectory(sys,qNear,[qNear(1:2); qRand(1:2)],options,[],'RRT',stepSize);
            else                
                [u0,x0] = computeTrajectoryFromAllMaximalActions(sys,qNear,qRand(1:2),stepSize,options);
            end
            
            [Xk,t] = downsampleUniformly(x0,options.sampSkipColl);
            t = t';
            Xk = Xk';
            Uk = double(u0,t);
            
            if length(t) > 1
                
                figure(3)
                plot(Xk(:,1),Xk(:,2),'k')
                drawnow
                
                qNew = Xk(end,:);
                
                % Test whether the path from qNear to qNew is in free space
                isect = checkIntersection(vBound,vObs1,vObs2,[],[],[],[],Xk,'finalPt',[],reg,sys);
                
                if ~isect, break, end
            end
        end
    catch ME
        disp(ME.message)
%         disp('trajectory too small; retrying')
%         keyboard
    end
    
end

if i == maxTrials, 
    qNew = [];
    disp(['Cannot add node ',num2str(node),' to tree'])
else
    % Since qNew is collision free, save the id of the parent and remove the child node from the reachability tree 
    nearI = nodeReach(nearReachI,1);
    nodeReach(nearReachI,:) = [];
end

end

function [u0,x0] = computeTrajectoryFromAllMaximalActions(sys,qNear,xyRand,stepSize,options)

nInputs = length(sys.umin);
if nInputs ~= length(sys.umax), error('umin and umax should be of the same size.'); end

[~,bVec] = buildBinaryCube([],[],nInputs);

for i = 1:size(bVec,1)
    for j = 1:nInputs
        U0(j) = 0;
        if bVec(i,j) == -1
            U0(j) = sys.umin(j);
        elseif bVec(i,j) == 1
            U0(j) = sys.umax(j);
        end
    end
    
    xTest{i} = computeOpenLoopTrajectory(sys,qNear,U0,stepSize,options);
    U0trial(i,:) = U0';
    
    distTrial(i) = norm(sys.state2SEconfig([],xTest{i}(end,:),[]) - xyRand);
end

[~,iMin] = min(distTrial);

x0 = xTest{iMin};
u0 = repmat(U0trial(iMin,:),size(x0,1),1);

end
