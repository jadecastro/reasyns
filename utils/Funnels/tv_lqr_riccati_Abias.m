function [ts,Ss] = tv_lqr_riccati_Abias(tspan,A,B,Q,R,Qf,T)
    % tspan can either be a pair [tstart tend] or a vector of values

    treatAsTimeInvariant = false;
    
    t0 = tspan(1);
        
    n = size(A(t0),1); m = size(B(t0),2);
    if size(A(t0)) ~= [n n], error('A must be n-by-n'); end
    if size(Qf) ~= [n n], error('Qf must be n-by-n'); end
    if size(Q(t0)) ~= [n n], error('Q must be n-by-n'); end
    if size(B(t0)) ~= [n m], error('B must be n-by-m'); end
    if size(R(t0)) ~= [m m], error('R must be m-by-m'); end
    
    if ~any(any(isnan(Qf)))
        [ts,Ss] = matrixODE(@ode45,...
            @(t,S) -(inv(T)*A(t)*T)'*S-S*(inv(T)*A(t)*T)-Q(t)+S*inv(T)*B(t)*inv(R(t))*B(t)'*inv(T)'*S,...
            fliplr(tspan),Qf);
    elseif treatAsTimeInvariant
        %         Ss = reshape(repmat(Q(0),1,length(tspan)),n,n,length(tspan));  % this was the (incorrect) implementation used in the VDP model
        Ss = []; % solve time-invariant LQR problem
    end
%     [ts,Ss] = matrixODE(@ode45,...
%                         @(t,S) -(inv(T)*A(t)*T + ...
%                         0e13*eps*[1 0 0;0 1 0;0 0 0])'*S-S*(inv(T)*A(t)*T + ...
%                         0e13*eps*[1 0 0;0 1 0;0 0 0])-Q(t)+S*inv(T)*B(t)*inv(R(t))*B(t)'*inv(T)'*S,...
%                         fliplr(tspan),Qf);
end