function [varargout] = Ellipse_plot(A, C, varargin)
%
%  Ellipse_Plot(A,C,N) plots a 2D ellipse or a 3D ellipsoid 
%  represented in the "center" form:  
%               
%                   (x-C)' A (x-C) <= 1
%
%  A and C could be the outputs of the function: "MinVolEllipse.m",
%  which computes the minimum volume enclosing ellipsoid containing a 
%  set of points in space. 
% 
%  Inputs: 
%  A: a 2x2 or 3x3 matrix.
%  C: a 2D or a 3D vector which represents the center of the ellipsoid.
%  N: the number of grid points for plotting the ellipse; Default: N = 20. 
%
%  Example:
%  
%       P = rand(3,100);
%       t = 0.001;
%       [A , C] = MinVolEllipse(P, t)
%       figure
%       plot3(P(1,:),P(2,:),P(3,:),'*')
%       hold on
%       Ellipse_plot(A,C)
%  
%
%  Nima Moshtagh
%  nima@seas.upenn.edu
%  University of Pennsylvania
%  Feb 1, 2007
%  Updated: Feb 3, 2007

%%%%%%%%%%%  start  %%%%%%%%%%%%%%%%%%%%%%%%%%%

% See if the user wants a different value for N.
%----------------------------------------------

color = [0.95 0.9 0.95];

if nargin >= 3,
    noplot = varargin{1};
else
    noplot = 0;
end

if nargin >= 4,
    N = varargin{2}; % Default value for grid
else
    N = 1;%20; % Default value for grid
end 

% check the dimension of the inputs: 2D or 3D
%--------------------------------------------
if length(C) == 3,
    Type = '3D';
elseif length(C) == 2,
    Type = '2D';
else
    display('Cannot plot an ellipse with more than 3 dimensions!!');
    return
end

% "singular value decomposition" to extract the orientation and the
% axes of the ellipsoid
[U D V] = svd(A);

if strcmp(Type, '2D'),
    % get the major and minor axes
    %------------------------------------
    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));

    theta = [0:1/N:2*pi+1/N];

    % Parametric equation of the ellipse
    %----------------------------------------
    state(1,:) = a*cos(theta); 
    state(2,:) = b*sin(theta);

    % Coordinate transform 
    %----------------------------------------
    X = V * state;
    X(1,:) = X(1,:) + C(1);
    X(2,:) = X(2,:) + C(2);
    
    varargout{1} = X(1,:);
    varargout{2} = X(2,:);
    
elseif strcmp(Type,'3D'),
    % generate the ellipsoid at (0,0,0)
    %----------------------------------
    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));
    c = 1/sqrt(D(3,3));
    warning off
    rmpath('C:\Users\Jon\Dropbox\Research\CreatePath\ParametricVerification\ellipsoids')
    rmpath('C:\Users\jad455\Dropbox\Research\CreatePath\ParametricVerification\ellipsoids')
    rmpath('C:\Cornell\Research\CreatePath\ParametricVerification\ellipsoids')
    [X,Y,Z] = ellipsoid(0,0,0,a,b,c,N);
    addpath('C:\Cornell\Research\CreatePath\ParametricVerification\ellipsoids')
    addpath('C:\Users\jad455\Dropbox\Research\CreatePath\ParametricVerification\ellipsoids')
    addpath('C:\Users\Jon\Dropbox\Research\CreatePath\ParametricVerification\ellipsoids')
    warning on
    
    %  rotate and center the ellipsoid to the actual center point
    %------------------------------------------------------------
    XX = zeros(N+1,N+1);
    YY = zeros(N+1,N+1);
    ZZ = zeros(N+1,N+1);
    for k = 1:length(X),
        for j = 1:length(X),
            point = [X(k,j) Y(k,j) Z(k,j)]';
            P = V * point;
            XX(k,j) = P(1)+C(1);
            YY(k,j) = P(2)+C(2);
            ZZ(k,j) = P(3)+C(3);
        end
    end
    
    varargout{1} = XX;
    varargout{2} = YY;
    varargout{3} = ZZ;
end


if ~noplot
    % Plot the ellipse
    %----------------------------------------
    if strcmp(Type,'2D'),
        fill(X(1,:),X(2,:),color);
        hold on;
        %     plot(C(1),C(2),'r*');
        %    axis equal, grid
        
    else
        H = surf(XX,YY,ZZ);
        set(H,'LineStyle','none','FaceColor',color)
        %    axis equal
        %    hidden off
    end
end
