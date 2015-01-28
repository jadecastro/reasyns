% SOSDEMO2 --- Lyapunov Function Search 
% Section 3.2 of SOSTOOLS User's Manual
% 

pvar x0 x1 x2;
vars = [x0; x1; x2];

Va = 434.6-632.98*x0+353.24*x0^2-101.16*x1+84.867*x1^2-88.782*x1*x0-35.901*x2+30.105*x2^2-46.73*x2*x0+9.1231*x2*x1;
Vdot_a = [-870.01+1003.7*x0-90.915*x0^2+558.9*x1+97.887*x1^2-732.63*x1*x0-141*x2+371.18*x2^2+19.742*x2^3-0.52465*x2^4-104.56*x2*x0-1.8182*x2*x0^2-400.34*x2^2*x0-7.234*x2^2*x0^2+8.4104*x2^3*x0+352.17*x2*x1-175.15*x2*x1^2-288.41*x2^2*x1-45.807*x2^2*x1^2-25.107*x2^3*x1+95.089*x2*x1*x0+365.41*x2^2*x1*x0];

Lxmonom = monomials([x0 x1 x2],0:2);

% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vars);


[prog,L] = sospolyvar(prog,Lxmonom,'wscoeff');

% =============================================
% Next, define SOSP constraints

% Constraint
prog = sosineq(prog,-(Vdot_a+L*(double(rho)-Va)));


% =============================================
% And call solver
prog = sossolve(prog);

% =============================================
% Finally, get solution
SOLV = sosgetsol(prog,L)


