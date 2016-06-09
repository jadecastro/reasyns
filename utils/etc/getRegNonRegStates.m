function [idxRegStates, idxNonRegStates, H] = getRegNonRegStates(obj,varargin)

n = obj.sysparams.n;
epsilon = 1e-6;

if ~isempty(varargin)
    t = varargin{1};
    x = varargin{2};
    u = varargin{3};
else
    t = [];
    x = zeros(n,1);
    u = [];
end

idxRegStates = [];
if nargout == 3
    y = obj.state2SEconfig(t,x,u);
    H = zeros(length(y),n);
else
    H = [];
end

for i = 1:n
    % Find the output at the baseline
    tryState = x;
    y = obj.state2SEconfig(t,tryState,u);
    
    % Perturb the states, one by one, in the positive direction
    tryState(i) = tryState(i) + epsilon;
    yPert = obj.state2SEconfig(t,tryState,u);
    
    if (yPert(1) - y(1)) ~= 0 || (yPert(2) - y(2)) ~= 0
        % flag any states that are coupled to x-y position (first two coordinates of an SE(2) system).
        idxRegStates = [idxRegStates; i];
    end
    
    if nargout == 3
        H(:,i) = (yPert - y)/epsilon;
    end
end
idxNonRegStates = setdiff(1:n,idxRegStates);