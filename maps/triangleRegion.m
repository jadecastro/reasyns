
% Set up map
pix2m = 1/100;

vBnd = [-3 3; 3 3; 3 -3; -3 -3];
vReg{1} = [vBnd; 1 0; 2 1; 2 -1];
vReg{2} = [1 0; 2 1; 2 -1];

epsilon = 0.05;

% Note on usage: adjacent faces must appear consecutively
Xbnd(1) = -(x(1) + 10) + epsilon*x(2) + eps*sum(x(3:n));
Xbnd(2) = -(x(2) + 10) + epsilon*x(1) + eps*sum(x(3:n));
Xbnd(3) = -(-x(1) + 10) + epsilon*x(2) + eps*sum(x(3:n));
Xbnd(4) = -(-x(2) + 10) + epsilon*x(1) + eps*sum(x(3:n));

X1(1) = -(0.7*x(1) + x(2) + 1.0 + 0.0) + eps*sum(x(3:n));
X1(2) = -(0.7*x(1) - x(2) + 1.0 - 0.0) + eps*sum(x(3:n));
X1(3) = -(-x(1) + 1.0) + eps*sum(x(2:n));

% Define region 1
Xreg{1}.ext = -Xbnd;
Xreg{1}.int = X1;

% Define region 2
Xreg{2}.ext = -X1;
Xreg{2}.int = [];

% Find vertices
for i = 1:length(Xreg)
    Xreg{i}.extVert = [];
    Xreg{i}.intVert = [];
    polytmp = subs(Xreg{i}.ext,[x(3);x(4)],[0;0]);
    for k = 1:length(Xreg{i}.ext)
        kp1 = mod(k,length(Xreg{i}.ext))+1;
        Xreg{i}.extVert(:,k) = newton(polytmp([k kp1])',[x(1);x(2)]);
    end
    polytmp = subs(Xreg{i}.int,[x(3);x(4)],[0;0]);
    for k = 1:length(Xreg{i}.int)
        kp1 = mod(k,length(Xreg{i}.int))+1;
        Xreg{i}.intVert(:,k) = newton(polytmp([k kp1])',[x(1);x(2)]);
    end
end
