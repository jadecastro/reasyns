
% Set up map
pix2m = 1/100;

vBnd = [-3 3; 3 3; 3 -3; -3 -3];
vReg{1} = [vBnd; 1 0; 2 1; 2 -1];
vReg{2} = [1 0; 2 1; 2 -1];

epsilon = 0.0001;

% Note: adjacent faces must appear consecutively
Xbnd(1) = -(x(1) + 0.5) + epsilon*sum(x(2:n));
Xbnd(2) = -(x(2) + 1) + epsilon*sum([x(1)]);
Xbnd(3) = -(-x(1) + 1) + epsilon*sum(x(2:n));
% Xbnd(3) = -(-x(1) + 30)^2 + epsilon*sum(x(2:n)^2);
Xbnd(4) = -(-x(2) + 1) + epsilon*sum([x(1)]);

X1(1) = (0.7*x(1) + x(2) + 1.0 + 0.0) ;
X1(2) = (0.7*x(1) - x(2) + 1.0 - 0.0) ;
X1(3) = (-x(1) + 1.0) + epsilon*sum(x(2:n));

% Define region 1
Xreg{1}.ext = -Xbnd;
Xreg{1}.int = X1;

% Define region 2
Xreg{2}.ext = -X1;
Xreg{2}.int = [];

% Find vertices
if n == 4
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
else
    for i = 1:length(Xreg)
        Xreg{i}.extVert = [];
        Xreg{i}.intVert = [];
        polytmp = Xreg{i}.ext;
        for k = 1:length(Xreg{i}.ext)
            kp1 = mod(k,length(Xreg{i}.ext))+1;
            Xreg{i}.extVert(:,k) = newton(polytmp([k kp1])',[x(1);x(2)]);
        end
        polytmp = Xreg{i}.int;
        for k = 1:length(Xreg{i}.int)
            kp1 = mod(k,length(Xreg{i}.int))+1;
            Xreg{i}.intVert(:,k) = newton(polytmp([k kp1])',[x(1);x(2)]);
        end
    end
end

