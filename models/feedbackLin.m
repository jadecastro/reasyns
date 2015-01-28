function[cmdV, cmdW] = feedbackLin(Vx, Vy, thetaR, e)
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

% Rotation from intertial to body frame
RIB = [ cos(thetaR) -sin(thetaR); 
        sin(thetaR)  cos(thetaR)
        ];

% Rotation from body to inertial frame
RBI = inv(RIB);

% Perform feedback linearization
VW = [1 0;0 1/e]*RBI*VxyI;

% Output commanded V and w
cmdV = VW(1);
cmdW = VW(2);

