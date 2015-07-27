function [ac] = computePolytopeAtomicControllerDubins(u0,x0,sys,ac,reg,varargin)

flowstarDir = '/home/jon/Software/flowstar-1.2.0';
verifTimeStep = 0.004;

% Q = ctrloptions.Q;
% R = ctrloptions.R;
% Qf = ctrloptions.Qf;

% [~,tk] = double(x0);
% if nargin > 4
%     [xRed,~] = downsampleUniformly(x0,sampSkip);
%     [uRed,tRed] = downsampleUniformly(u0,sampSkip);
% else
%     [xRed,~] = double(x0);
%     [uRed,tRed] = double(u0);    
% end

% NB: For polytope-based reach set computations, we will forego including a feedback controller

% xtraj = PPTrajectory(foh(tRed,xRed)); % should we be using xRed, uRed here?
% utraj = PPTrajectory(foh(tRed,uRed(2,:)));
% 
% % Declare Dubins car model
% p = DubinsPlant();
% 
% % Set input limits
% p = setInputLimits(p,-Inf,Inf);
% 
% % Do tvlqr
% utraj = setOutputFrame(utraj,p.getInputFrame);
% xtraj = setOutputFrame(xtraj,p.getStateFrame);
% [c,V] = tvlqr(p,xtraj,utraj,Q,R,Qf);
% poly = taylorApprox(feedback(p,c),xtraj,[],3);
% ts = xtraj.getBreaks();

% Options for subsequent flowstar computations
options = struct();
options.max_iterations = 5; % Maximum number of iterations to run for

% Form each trajectory into an array; extract endpoint values
t = x0.getTimeVec();
xi = ppval(x0.pp,0);
xf = ppval(x0.pp,t(end));

% Find which hyperplane for which our reach set will be aligned - assume every region is convex for now
polyIdx = 1;    % index of the polytope
[Hp,Kp] = double(reg(1).p(polyIdx));
d(1:length(Kp)) = Inf;
t_ac = ac.x0.getTimeVec();
for idx = 1:length(Kp)
    H(idx) = hyperplane(Hp(idx,:)',Kp(idx));
end
idxfound = false;
for i = 2:length(t_ac)
    x1 = ppval(ac.x0.pp,t_ac(i-1));
    x2 = ppval(ac.x0.pp,t_ac(i));
    for vidx = 1:size(reg.v,1)
        x3 = reg.v(vidx,:);
        x4 = reg.v(mod(vidx,size(reg.v,1))+1,:);
        if intersectPoint(x1(1),x1(2),x2(1),x2(2),x3(1),x3(2),x4(1),x4(2))
            idxfound = true;
            break
        end
    end
    if idxfound
        break
    end
end
i
% find the hyperplane index from the found vertex index
idxfound = false;
for idx = 1:length(H)
    if contains(H(idx),x3') && contains(H(idx),x4')
        idxfound = true;
        break
    end
end
if ~idxfound
    error('did not find a hyperplane from the given vertices!')
end

% Decide the orientation of the initial set (a box) based on collision of the reach set with the estimated nearest hyperplane.
% theta = xi(3);  % placeholder
boxAlignmentAngle = atan2(-Hp(idx,1)/Hp(idx,2),1);
boxAlignmentAngle = boxAlignmentAngle + pi/2;  % convert from tangent to hyperplane normal
boxAlignmentAngle = boxAlignmentAngle + pi;  % TODO: choose whether or not to add pi based on containment of initial box

rot = [cos(boxAlignmentAngle), -sin(boxAlignmentAngle); sin(boxAlignmentAngle), cos(boxAlignmentAngle)]'; 

thetaInitRelative = xi(3) - boxAlignmentAngle;  % the angle specified to flowstar which determines the initial angle of the Dubin's car model, relative to the specified orientation of the initial box.
thetaFinalRelative = xf(3) - boxAlignmentAngle;  

initStateBoxRelative = [-0.01 0.01; -0.01 0.01; -0.02+thetaInitRelative 0.02+thetaInitRelative];
finalStateRelative = [rot*(xf(1:2) - xi(1:2)); thetaFinalRelative];

% plots
figure(2), clf, hold on, axis equal
plot(reg(1),'y')
initVert = [
    (rot'*[initStateBoxRelative(1,1) initStateBoxRelative(2,1)]' + xi(1:2))'; 
    (rot'*[initStateBoxRelative(1,1) initStateBoxRelative(2,2)]' + xi(1:2))';
    (rot'*[initStateBoxRelative(1,2) initStateBoxRelative(2,2)]' + xi(1:2))';
    (rot'*[initStateBoxRelative(1,2) initStateBoxRelative(2,1)]' + xi(1:2))';
    ];
plot(region(initVert),'r')
plot(x0,'g',2)
ax = axis;
plot(hyperplane(Hp(idx,:)',Kp(idx)))
axis(ax); drawnow


% Configure input file for flowstar
infname = '/home/jon/Software/flowstar-1.2.0/dubins_car_template.model';
outfname = '/home/jon/Software/flowstar-1.2.0/dubins_car_new.model';
writeToFlowstarModelFile(infname, outfname, initStateBoxRelative, finalStateRelative)


% External call to flowstar - expects output to be written to an m-file (dubins.m)
disp(' Computing reach set...');
[callSuccess, logDump] = system('source ~/.bashrc && cd /home/jon/Software/flowstar-1.2.0 && ./flowstar < dubins_car_new.model');
if callSuccess
    error(['External call to flowstar failed.  ',logDump])
else
    disp(' Reach set successfully generated.');
end

% Next, retrieve the data through a single call to the written m-file
run('/home/jon/Software/flowstar-1.2.0/outputs/dubins.m')

figure(2)
for i = 1:length(X)
    tmp(:,1) = X{i};
    tmp(:,2) = Y{i};
    for j = 1:size(tmp,1)
        tmp(j,:) = (rot'*tmp(j,:)' + xi(1:2))';
    end
    X{i} = tmp(:,1);
    Y{i} = tmp(:,2);
    
    if length(X) < 1000 || ~rem(i,10)
        plot([X{i}],[Y{i}],'Color',[0.5 0.5 0.5])
    end
    vert(:,:,i) = [[X{i}],[Y{i}]];
    
    % construct new trajectories that we can use (NB: for the current
    % setup, only need the initial and final values
    theta_tmp = i/length(X)*xf(3) + (length(X) - i)/length(X)*xi(3);  % placeholder until we figure out how to bring in theta from the verifier.
    x_new(i,:) = [mean(X{i}) mean(Y{i}) theta_tmp];  % extremely dumb first pass to get a "trajectory" of some kind brought in, since flowstar doesn't work with trajectories.
    u_new(i,:) = [0.05 1/0.1*(-(xf(1) - x_new(i,:))*sin(theta_tmp) + (xf(2) - x_new(i,:))*cos(theta_tmp))]; % this is also just a placeholder. Better to do an feval with the function in the 'sys' object.
end

tVerif = 0:verifTimeStep:verifTimeStep*(length(X)-1);
vert0 = traject(tVerif,vert);

% assemble the atomic controller, re-sampling the original vectors
if verLessThan('matlab','7.15')
    zeromatrix = repmat([0 0 0],[1,1,length(tVerif)]);
else
    zeromatrix = repmat([0 0 0],1,1,length(tVerif));
end
K0 = traject(tVerif,zeromatrix);   % For polytope-verified atomic controllers, let's deactivate the feedback controller
x0 = traject(tVerif,x_new');
u0 = traject(tVerif,u_new');

ac = polytopeAC(x0,u0,K0,vert0,sys);

% ac = resample(ac); % resample back to the original resolution

% plot the projection
% figure(5)
% hold on
% options.inclusion = 'projection';
% options.plotdims = [1 2];
% plotFunnel(Vxframe,options);
% fnplt(xtraj,[1 2]);
% axis equal

figure(5)
hold on
plot(ac,sys,5)

% Create a 1-D hyperplane that separates the reach set from the plane for the next region
d = [];
clear jmin dmintmp
for i = 1:size(vert,3)
    for j = 1:size(vert,1)-1
        dmintmp(j) = dist2hyperplane(vert(j,:,i),H(idx));
    end
    [dtmp,jmin] = sort(dmintmp);
    d = [d; dtmp(1:2)];
end
[dmin, imin] = min(d);
xmin = vert(jmin(1:2),:,imin);
mmin = (xmin(1,2) - xmin(2,2))/(xmin(1,1) - xmin(2,1));
hmin = [-mmin 1];
kmin = xmin(1,2) - mmin*xmin(1,1);
Hmin = hyperplane(hmin',kmin);

% Now, let's use the hyperplane to split the region into two new regions
% first, find the points...
point = [];
Hidx = [];
for idx = 1:length(H)
    m1 = -Hp(idx,1)/Hp(idx,2);
    k1 = Kp(idx)/Hp(idx,2);
    pointtmp = [(k1 - kmin)/(mmin - m1); (mmin*k1 - m1*kmin)/(mmin - m1)];
    if isinside(reg.p,pointtmp)
        point = [point; pointtmp'];
        Hidx = [Hidx; idx];
    end
end
if size(point,1) < 2
    error('did not find any in-region intersections with the hyperplane!')
end
if size(point,1) > 2
    warning('interections should have no more than two entries!!')
end

% contruct the new regions...
vidxsav = [];
for idx = Hidx';  % these must already be ordered as per the region polynomial (see double() call above)
    for vidx = 1:size(reg.v,1)
        if contains(H(idx),reg.v(vidx,:)') && contains(H(idx),reg.v(mod(vidx,size(reg.v,1))+1,:)')
            break
        end
    end
    vidxsav = [vidxsav; vidx];
end
reg1 = region([reg.v(1:vidxsav(1),:); point; reg.v(vidxsav(2):end,:)]);
reg2 = region([point(1,:); reg.v(vidxsav(1):vidxsav(2),:); point(2,:)]);

figure(10), hold on
plot(reg1,'g')
plot(reg2,'m')
ax = axis;
plot(Hmin)
axis(ax); drawnow

% 
% % Tests to make sure simulated trajectories stay inside computed funnel
% doTest = 0;
% if doTest
%     
%     % Create closed loop system with optimized controller
%     sysCl = feedback(p,c);
%     V0 = Vtraj.getPoly(0); % Get inlet of funnel
%     xinit = getLevelSet(decomp(V0),V0,[0;0;0]);
%     
%     Vsall = [];
%     figure(30)
%     hold on
%     grid on
%     for j = 1:5
%         Vs = [];
%         x0 = 0.95*xinit(:,j) + xtraj.eval(0); % Simulate from 0.95*boundary of funnel
%         xsim = sysCl.simulate([0 ts(end)],x0);
%         for k = 1:length(ts)
%             Vs = [Vs, Vxframe.eval(ts(k),xsim.eval(ts(k)))];
%         end
%         Vsall = [Vsall, Vs];
%         
%         plot(ts,Vs)
%         plot(ts,ones(1,length(ts)),'ro')
%         drawnow;
%     end
%     
%     if ~(all(Vsall < 1))
%         success = false;
%         disp('A Trajectory left the funnel!!')
%     else
%         success = true;
%         disp('All simulated trajectories stayed inside funnel.')
%     end
%     
% end

