

syms t ddth ddphi ddpsi ddu dth dphi dpsi du th phi psi u Ja Jt md mc R l g

dx = dpsi*R*cos(phi);
dy = dpsi*R*sin(phi);

T = 0.5*Ja*(dth^2 + dphi^2*cos(th)^2) + 0.5*Jt*(dphi*sin(th) + dpsi)^2 + ...
    0.5*md*(R^2*dth^2 + 2*R*(dy*cos(phi) - dx*sin(phi))*dth*cos(th) + ...
    (dx - R*dphi*sin(th)*cos(phi))^2 + (dy - R*dphi*sin(th)*sin(phi))^2) + ...
    0.5*mc*((R+l)^2*(dth + du/(R+l))^2 + 2*(R+l)*(dy*cos(phi) - dx*sin(phi))*(dth+du/(R+l))*cos(th) + ...
    (dx - ((R+l)*sin(th)+u*cos(th))*dphi*cos(phi))^2 + (dy - ((R+l)*sin(th)+u*cos(th))*dphi*sin(phi))^2);

V = md*g*R*cos(th) + mc*g*((R+l)*cos(th)+u*sin(th));
L = T-V;

dLddTh = diff(L,dth);
dLddPhi = diff(L,dphi);
dLddPsi = diff(L,dpsi);
dLddU = diff(L,du);

dLdTh = diff(L,th);
dLdPhi = diff(L,phi);
dLdPsi = diff(L,psi);
dLdU = diff(L,u);

% Form EOM from Lagrange
lagEqn = Lagrange(L,[th dth ddth phi dphi ddphi psi dpsi ddpsi u du ddu]);

% Now, parameterize the model and form the 
mc = 3;
md = 1;
R = 1;
l = 4;
Ja = md*R^2/2;
Jt = md*R^2/4;
g = 1;

homEqn = subs(lagEqn);
thEqn = matlabFunction(homEqn(1));
phiEqn = matlabFunction(homEqn(2));
psiEqn = matlabFunction(homEqn(3));
uEqn = matlabFunction(homEqn(4));

