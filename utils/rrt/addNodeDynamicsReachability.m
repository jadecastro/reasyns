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
        i
        % Bias the sampling toward qGoal
        if rand < gaussWeight
            qRand = randn(1,n)*chol(M) + qGoal;
        else
            qRand = [min(vBound) -pi] + [(max(vBound) - min(vBound)) 2*pi].*rand(1,n);
        end
        figure(3)
        hold on
        plot(qRand(1),qRand(2),'m.')
        
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
            
            % Set qNewStar as a stepSize away from qNear in the direction of qRand
            %qNewStar = qNear + stepSize/norm(qRand - qNear)*(qRand - qNear);
            
            % Construct motion plan between qNear and qNewStar; find qNew
            %     [t,Xk,Uk] = genNominalTrajectory(qNear,[qNear(1:2); qNewStar(1:2)],Hout,n,isCyclic,modelType);
            %     [t,Xk,Uk] = genNominalTrajectory(qNear,[qNear(1:2); qRand(1:2)],Hout,n,isCyclic,modelType);
            [u0,x0] = computeTrajectory(sys,qNear,[qNear(1:2); qRand(1:2)],options,[],'RRT',stepSize);
            
            [Xk,t] = downsampleUniformly(x0,options.sampSkipColl);
            t = t';
            Xk = Xk';
            Uk = double(u0,t);
            
            if length(t) > 1
                
                %                 clear xknrm indxRem
                %                 xknrm(1) = 0;
                %                 for k = 1:length(t)-1, xknrm(k+1) = norm(Xk(k+1,1:2) - Xk(k,1:2)); end
                %                 indxRem = (cumsum(xknrm) > stepSize);
                %                 if any(indxRem)
                %                     t(indxRem) = [];
                %                     Xk(indxRem,:) = [];
                %                     Uk(indxRem,:) = [];
                %                 end
                
                plot(Xk(:,1),Xk(:,2),'k')
                drawnow
                
                qNew = Xk(end,:);
                
                % Test whether the path from qNear to qNew is in free space
                % isect = double(any(qNew < min(vBound)+radius) || any(qNew > max(vBound)-radius));
                isect = checkIntersection(vBound,vObs1,vObs2,[],[],[],[],Xk,'finalPt',[],reg,sys);
                
                if ~isect, break, end
            end
        end
    catch ME
        disp(ME.message)
        disp('trajectory too small; retrying')
%         keyboard
    end
    
end

if i == maxTrials, 
    qNew = [];
    disp('Cannot add node to tree')
else
    % Since qNew is collision free, save the id of the parent and remove the child node from the reachability tree 
    nearI = nodeReach(nearReachI,1);
    nodeReach(nearReachI,:) = [];
end