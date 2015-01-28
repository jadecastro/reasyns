function[cmdV, cmdW] = feedbackLinCar(Vx, Vy, theta, phi, l, e)
% FEEDBACKLIN: Apply feedback linearization to convert desired inertial
% velocities to local robot velocities.
% 
%   [CMDV,CMDW] = FEEDBACKLIN(VX,VY,THETAR,E) 
% 
%   INPUTS
%       Vx          desired x-velocity in global frame (m/s)
%       Vy          desired y-velocity in global frame (m/s)
%       thetaR      current robot angle relative to inertial frame (deg)
%       e           feedback linearization distance (in meters)
% 
%   OUTPUTS
%       cmdV        forward velocity command (m/s)
%       cmdW        angular velocity command (rad/s)
% 


VxyI = [Vx; Vy];

A = [ cos(theta)-tan(phi)*(sin(theta)+e*sin(theta+phi)/l) -e*sin(theta+phi);
      sin(theta)+tan(phi)*(cos(theta)+e*cos(theta+phi)/l) e*cos(theta+phi)
      ];
    
VW = inv(A)*[Vx ; Vy];

% Output commanded V and w
cmdV = VW(1);
cmdW = VW(2);

