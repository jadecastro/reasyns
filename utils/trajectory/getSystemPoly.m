function [poly, p] = getSystemPoly(p, xtraj, c, V)
%
% Compute the polynomial approximation of the drake plant
%  xtraj and V must be a PPTrajectory, c must be an AffineSystem (Drake)

% Set input limits
p = setInputLimits(p,-Inf,Inf);

% compute the polynomial approximation
poly = taylorApprox(feedback(p,c),xtraj,[],3);

num_xc = poly.getNumContStates();
if (isa(V,'Trajectory'))
    V1 = V.inFrame(poly.getStateFrame);
    Q = eye(num_xc);
    V0 = tvlyap(poly,V1,Q,Q);
else
    V0 = V;
end

poly = poly.inStateFrame(V0.getFrame); % convert system to Lyapunov function coordinates
