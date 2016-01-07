function H = drawCar(x)

scl = 0.16;  % NB: axle displacement = 0.5 with this geometry

bodyCoords = scl*[[-0.8 -0.8 0.8 0.8];[-1.2 1.2 1.2 -1.2]];
wheel1Coords = scl*[[-0.8 -1.2 -1.2 -0.8];[0.4 0.4 -0.4 -0.4]+0.7];
wheel2Coords = scl*[[0.8 1.2 1.2 0.8];[0.4 0.4 -0.4 -0.4]+0.7];
wheel3Coords = scl*[[-0.8 -1.2 -1.2 -0.8];[0.4 0.4 -0.4 -0.4]-0.7];
wheel4Coords = scl*[[0.8 1.2 1.2 0.8];[0.4 0.4 -0.4 -0.4]-0.7];

angleBC = atan2(bodyCoords(1,:),bodyCoords(2,:)) - x(3) + pi/2;
magBC = sqrt(bodyCoords(1,:).^2+bodyCoords(2,:).^2);
rotBodyCoords = [magBC.*sin(angleBC); magBC.*cos(angleBC)];
angleW1 = atan2(wheel1Coords(1,:),wheel1Coords(2,:)) - x(3) + pi/2;
magW1 = sqrt(wheel1Coords(1,:).^2+wheel1Coords(2,:).^2);
rotWheel1Coords = [magW1.*sin(angleW1); magW1.*cos(angleW1)];
angleW2 = atan2(wheel2Coords(1,:),wheel2Coords(2,:)) - x(3) + pi/2;
magW2 = sqrt(wheel2Coords(1,:).^2+wheel2Coords(2,:).^2);
rotWheel2Coords = [magW2.*sin(angleW2); magW2.*cos(angleW2)];
angleW3 = atan2(wheel3Coords(1,:),wheel3Coords(2,:)) - x(3) + pi/2;
magW3 = sqrt(wheel3Coords(1,:).^2+wheel3Coords(2,:).^2);
rotWheel3Coords = [magW3.*sin(angleW3); magW3.*cos(angleW3)];
angleW4 = atan2(wheel4Coords(1,:),wheel4Coords(2,:)) - x(3) + pi/2;
magW4 = sqrt(wheel4Coords(1,:).^2+wheel4Coords(2,:).^2);
rotWheel4Coords = [magW4.*sin(angleW4); magW4.*cos(angleW4)];

H1 = fill(rotBodyCoords(1,:)+x(1),rotBodyCoords(2,:)+x(2),'y');
H2 = fill(rotWheel1Coords(1,:)+x(1),rotWheel1Coords(2,:)+x(2),'k');
H3 = fill(rotWheel2Coords(1,:)+x(1),rotWheel2Coords(2,:)+x(2),'k');
H4 = fill(rotWheel3Coords(1,:)+x(1),rotWheel3Coords(2,:)+x(2),'k');
H5 = fill(rotWheel4Coords(1,:)+x(1),rotWheel4Coords(2,:)+x(2),'k');
H = [H1;H2;H3;H4;H5];

