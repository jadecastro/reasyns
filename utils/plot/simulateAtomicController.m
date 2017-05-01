function simulateAtomicController(ac_trans, sys, varargin)

if isempty(ac_trans) || isempty(sys)
    error('simulate: ''ac_trans'' and ''sys'' are required but are empty.')
end

figure(4)
clf
if nargin > 2
    % Plot the workspace.
    reg = varargin{1};
    plot(reg);
end
hold on
axis equal

iTrans = 13;  % The transition of the automaton corresponding to the atomic controller we wish to simulate.
doRandom = true;  % If true, then the simulation is initialized over random 
                  % initial conditions within the initial ellipse.
                  % Otherwise, simulate with a deterministic initial
                  % condition.
T = 0.01;  % Time step.

if doRandom
    % Populate the initial state with random samples inside the initial ellipse.
    ellInit = ac_trans{iTrans}.ell(1);
    [c, Q] = double(ellInit);
    xRand = sampleEllipsoidBoundary(zeros(length(c),1), Q, 1);
    randNormDistFromCenter = rand(1,length(c));
    xInitialOffset = randNormDistFromCenter.*xRand;
else
    % Use a determinisitic intitial condition.
    xInitialOffset = [0.01; 0.05; 0; 0]';
end
xk = ac_trans{iTrans}.x0.double(0)' + xInitialOffset;  % Initial state.
if ~isinternal(ellInit, xk')
    error('simulateAtomicController: The supplied initial condition is not within the initial ellipse.')
end

t = getTimeVec(ac_trans{iTrans}.x0);
Tend = floor(t(end)/T);  % Time horizon over which to simulate.

% Plot the funnel.
ac_trans{iTrans}.plot(sys,4);

% Simulate it!
tVect = [];
uVect = [];
xVect = [];
for k = 1:Tend
    tk = T*(k - 1);
    tVect = [tVect; tk];
    
    % Execute one step of the controller.
    uk = ac_trans{iTrans}.evaluate(xk');
    uVect = [uVect; uk];
    
    % Evaluate the dynamics.
    [~, xint] = ode45(@(t, x) sys.dynamics(t, x, uk), [tk, tk + T], xk);
    xk = xint(end,:);
    xVect = [xVect; xk];
end

figure(4)
plot(xVect(:,1), xVect(:,2), 'k', 'LineWidth', 2)
