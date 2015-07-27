function [v,w,errorMsg,acLastData] = executeControllersSingleStep(aut,sys,ac_trans,ac_inward,x,t,currReg,nextReg,acLastData)

persistent t_offset t_base tend ell t_trials

global u_sav x_sav x0_sav

delayTime = -0.1;

timeBaseFlag = false;
useLastAC = true;

x(3) = x(3) + 0.2;

prevCtrl = false;

if isempty(t_base)
    t_base = clock;
end

v = 0;  w = 0.001;
errorMsg = 'nothing to see here...';

try
    if size(x,2) > size(x,1), x = x'; end
    
    trans = vertcat(aut.trans{:});
    currStateVec = find([aut.q{:}]==currReg);
    nextStateVec = find([aut.q{:}]==nextReg);
    currTrans = find(ismember(trans(:,1),currStateVec) & ismember(trans(:,2),nextStateVec), 1);
    if isempty(currTrans)
        error('you have specified an invalid transition!')
    end
    currState = trans(currTrans,1);  nextState = trans(currTrans,2);
    
    % identify the funnel we're in and get the ellipse index
    iTrans = false;
    iIn = false;
    if isinternal(ac_trans{currTrans}, x,'u',sys)
%     if isinternal(ac_trans{currTrans}, x,'u')
        iTrans = currTrans;
        ac = ac_trans{currTrans};
    end
    if isempty(iTrans) && ~isempty(ac_inward)
        if isinternal(ac_inward{currState}, x,'u',sys)
%         if isinternal(ac_inward{currState}, x,'u')
            iIn = currState;
            ac = ac_inward{currState};
        end
    end
    
    % if no funnels in the new region/transition OR have temporarily left a funnel, then keep activating the previous one.
    currStateOld = [];  nextStateOld = [];
    if ~iTrans && ~iIn
        if useLastAC
            % propagate the appropriate old states on to the next iteration
            if ~isempty(acLastData{5}) && ~isempty(acLastData{6})
                currStateOld = acLastData{5};  nextStateOld = acLastData{6};
            else
                currStateOld = acLastData{3};  nextStateOld = acLastData{4};
            end
            currTransOld = trans(:,1) == currStateOld & trans(:,2) == nextStateOld;
            if ~isempty(acLastData{1})
                ac = ac_trans{currTransOld};
            else
                ac = ac_inward{currStateOld};
            end
            iTrans = acLastData{1};  iIn = acLastData{2};
            prevCtrl = true;
            disp('WARNING: no funnels found! Using the previous controller.')
        else
            iTrans = currTrans;
            ac = ac_trans{currTrans};
            disp('WARNING: no funnels found! Forcing an (unverified) controller.')
        end
    end
    
    
    % determine the new time offset if something has changed
    if timeBaseFlag
        acData = {iTrans, iIn, currState, nextState, currStateOld, nextStateOld};
        for i = 1:length(acData)-2,
            cmp(i) = acData{i}==acLastData{i};
        end
        cmp
        if (~all(cmp) || ~(all(ismember([acLastData{1:4}],[acData{1:4}])) && all([acLastData{3:4}] == [acData{3:4}]))) && ~prevCtrl
            t_trials = getTimeVec(ac.x0);
            tend = t_trials(end);
            t_base = clock
            ell = ellipsoid(ac);
            
            for i = 1:length(ell)
                if isinternal(ell(i),x,'u')
%                 if isinternal(ell(i),x,'u')
                    t_offset = t_trials(i);
                    break
                end
            end
        end
        
        acLastData = acData;
        
        if isempty(t)  % if time is not given, we compute it here
            t = etime(clock,t_base) + delayTime
        end
        t_offset
        teval = min(tend,t + t_offset);
    else
        acData = {iTrans, iIn, currState, nextState, currStateOld, nextStateOld};
        for i = 1:length(acData)-2,
            cmp(i) = acData{i}==acLastData{i};
        end
        cmp
        if ~all(cmp) || ~(all(ismember([acLastData{1:4}],[acData{1:4}])) && all([acLastData{3:4}] == [acData{3:4}]))
            t_base = clock
            ell = ellipsoid(ac);
            t_trials = getTimeVec(ac.x0);
        end
        acLastData = acData;
        
        minDelta = inf;
        weights = [1;1;0.2];
        for i = 1:length(ell)
            xtmp = ppval(ac.x0.pp,t_trials(i));
            testMinDist = min([norm(weights.*(xtmp(1:3) - (x(1:3)+[0 0 2*pi]'))) norm(weights.*(xtmp(1:3) - x(1:3))) norm(weights.*(xtmp(1:3) - (x(1:3)-[0 0 2*pi]')))]);
            if testMinDist < minDelta
                teval = t_trials(i);
                minDelta = testMinDist;
            end
        end
        if isempty(t)  % if time is not given, we compute it here
            t = etime(clock,t_base)
        end
    end
    
    % compute a command
    K = ppval(ac.K.pp,teval);
    x0 = ppval(ac.x0.pp,teval)
    u0 = ppval(ac.u0.pp,teval);
    
    u_test = [(K*(x+[0 0 2*pi]' - x0)) (K*(x - x0)) (K*(x-[0 0 2*pi]' - x0))];
    u_ctrl = u_test(:,abs(u_test(2,:)) == min(abs(u_test(2,:))));
    u = u0 + 1*u_ctrl;
%     u(2) = max(min(u(2),16),-16);
    
    v = u(1);  w = u(2);
    
    u_sav = [u_sav; t u'];
    x_sav = [x_sav; t x'];
    x0_sav = [x0_sav; t x0'];
    
catch ME
    rethrow(ME)
    errorMsg = ME.message;
end

