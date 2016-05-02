function [ac,c,isectIdx] = computeAtomicController(u0,x0,sys,reg,ellToCompose,options,varargin)
% reg must contain the field 'init'; if it also includes 'goal' and
% 'union', then additional appropriate funnel containment checks will be made.
% 
% If only two output parameters are supplied, then the function checks if
% invariant violations occur and returns failure if so.

% If isectIdx is specified as an additional output parameter, then do not
% return failure when the funnel leaves the invariant; instead return the
% funnel and the time indices at which the violations occur (isectIdx).

debug = true;

rho_if = varargin{1};

[ac, c] = recursiveControllerComputation(u0,x0,sys,reg,ellToCompose,options,rho_if);

if isfield(reg,'goal')
    regGoal = reg.goal;
elseif isfield(reg,'init')
    regGoal = reg.init;
end

disp('Checking the computed funnel...')

% Determine if the final ellipse is invariant to the goal region
% NB: this should be alreay covered by the trajectory checks, but do it anyway for modularity.
[H,K] = double(regGoal.p);
hpp = hyperplane(H',K');
ellProj = projection(ac,sys);
if intersect(ellProj(end),hpp,'u')
    error('Final computed ellipse of the computed funnel is not contained inside the final region.  Adjust rho_if or Qf.')
end

% Also determine if there is enough volume within the initial region for
% the final ellipse to fit 
% NB: assumes we are using the same controller/dynamics throughout.  Also
% assumes the final ellipse is a ball.
[H,K] = double(reg.init.p);
hpp = hyperplane(H',K');
t = ac.x0.getTimeVec;

isectIdx = [];
if ~debug
    Ntrials = 100;
    isContained = false;
    for i = 1:Ntrials
        % Try a certain number of randomly-computed funnels, concentrated about
        % the first point in the next funnel.  If none are found to be in the
        % interior, declare failure.
        
        if i == 1
            centerTest = double(ac.x0,t(1));
        else
            centerTest = double(ac.x0,t(1)) + 10*sys.sysparams.Qrand*randn(sys.sysparams.n,1);
        end
        
        [~,~,H] = sys.getRegNonRegStates([],centerTest,[]);
        ballTestProj = projection(ballTest,H(1:2,:)');
        
        if ~intersect(ballTestProj,hpp,'u')  % if contained within the region, check that it is also contained within the funnel.
            
            ballTest = ellipsoid(centerTest,inv(sys.sysparams.Qf));
            
            isContained = ac.funnelContainsEllipsoid(sys,ballTest,100);
            
            if isContained
                break;
            end
        end
    end
    if ~isContained
        error('Computed funnel does not contain enough volume to compose with an incoming funnel.  Either adjust rho_if or Qf or the feedback controller parameters.')
    end
    
    % Lastly, check for invariance to the region pair
    if isfield(reg,'goal')
        [res1, isectIdx1] = isinside(ac,reg.init,sys);
        [res2, isectIdx2] = isinside(ac,reg.goal,sys);
        
        % Do a simple intersection test first; do a more expensive
        % pair-containment test for any ellipses that fail the test
        isectIdxTest = intersect(isectIdx1,isectIdx2);
        isectIdx = [];
        for i = isectIdxTest;
            % TODO: could this go into any of the existing classes, i.e. Region?
            [center,Qell] = double(ellProj(i));
            xyPoint = ellipsoidrand(center,Qell,20);  % we know this is 2-D, so can get away with a handful of boundary points
            for j = 1:size(xyPoint,1)
                if ~isinside(vertcat(reg.init.p),xyPoint) || ~isinside(vertcat(reg.goal.p),xyPoint)
                    isectIdx = [isectIdx; i];
                    break
                end
            end
        end
        
        % since 'atomiccontroller/isinside' is not really doing containment
        % checking, need to manually add in any that intersect with either one
        % of the regions.
        isectIdx = sort([isectIdx; setdiff(isectIdx1,isectIdxTest); setdiff(isectIdx2,isectIdxTest)]);
        
    elseif isfield(reg,'init')
        [res, isectIdx] = isinside(ac,reg.init,sys);
    end
    
    if ~isempty(isectIdx)
        if nargout < 3
            error('Computed funnel is not contained within the region pair.  Adjust rho_if or any of the feedback controller parameters.')
        else
            warning('Computed funnel is not contained within the region pair.  Adjust rho_if or any of the feedback controller parameters.')
        end
    end
end

end

function [ac,c] = recursiveControllerComputation(u0,x0,sys,reg,ellToCompose,options,varargin)
%

ac = [];
c = [];

if isfield(reg,'union')
    regInvariant = reg.union;
else
    regInvariant = reg.init;
end 
    
try
    rho_if = varargin{1};
    if length(varargin) == 1, endSegmentFlag = true; else endSegmentFlag = varargin{2}; end
    [ac,c] = computeAtomicControllerSegment(u0,x0,sys,regInvariant,ellToCompose,rho_if,options,endSegmentFlag);

catch ME
    if strcmp(ME.identifier, 'Drake:PolynomialTrajectorySystem:InfeasibleRho')
        if length(x0) > 10
            % bisect the (sub) trajectory
            [x01, x02] = bisect(x0);
            [u01, u02] = bisect(u0);
            
            figidx = 20;
            
            disp('two')
            %TODO: ensure containment of sequenced funnels - reverse ordering 
            [ac2, c2] = recursiveControllerComputation(u02,x02,sys,reg,ellToCompose,options,rho_if,endSegmentFlag);
            
            t0 = ac2.P.getTimeVec();
            
            sys.sysparams.Qf = double(ac2.P,t0(1));
            
            tmp = projection(ac2,sys);
            figure(figidx)
            plot(tmp(1),'r')
            hold on
            %             Qf = getMaximalQ(funnelI2,Xk2(1+Noverlap,:));
            %keyboard
            
            disp('one')
            [ac1, c1] = recursiveControllerComputation(u01,x01,sys,reg,ellToCompose,options,rho_if,false);
            
            t0 = ac1.P.getTimeVec();
            
            sys.sysparams.Qf = double(ac1.P,t0(1));
            
            tmp = projection(ac1,sys);
            figure(figidx)
            plot(tmp(end),'b')
            hold on
            
            %collect into one atomiccontroller
            ac = merge(ac1,ac2); 
            c = [c1; c2];
        else
            error('Rho infeasible after segmenting the trajecory to its maximal extent. No smaller segments allowed.')
        end
    else
        rethrow(ME)
    end
end

end

function Q = getMaximalQ(funnel,X)

N = 100;
n = length(X);
H = eye(2);
isCyclic = [zeros(n-1,1); 1];

Qi = 1e3*eye(length(X));
Ei = ellipsoid(X',inv(Qi));
E = Ei;
Ep = Ei;

for j = 1:length(funnel.t)
    ellBndInv{1}(j,1).x = funnel.x(j,:)';
    ellBndInv{1}(j,1).P = funnel.P(:,:,j)/funnel.rho(j); 
end

% starting with a tiny rho, iteratively check containment while
% increating rho
for k = 1:100
    [xbar,invQ] = double(E);
    qTest = ellipsoidrand(xbar,invQ,N);
    clear badIndx isect
    [isect,badIndx] = checkIntersection5([],[],ellBndInv,[],qTest,H,n,isCyclic);
    if isect
        break
    end
    Ep = E;
    E = shape(E,1.1);  % increase by 10% and try again..
end
[xbar,invQ] = double(Ep);
Q = inv(invQ);
% if debugFlg
%     plot(projection(Ep,[H; zeros(n-length(H),length(H))]),'g');
%     drawnow
%     %             keyboard
% end

end


