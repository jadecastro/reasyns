function [c, V] = computeTVLQR(p, utraj, xtraj)
%
% Compute the time-varying LQR controller

Q = p.sysparams.Q;
R = p.sysparams.R;
Qf = p.sysparams.Qf;

% Set input limits
p = setInputLimits(p,-Inf,Inf);

% Do tvlqr
[c,V] = tvlqr(p,xtraj,utraj,Q,R,Qf);