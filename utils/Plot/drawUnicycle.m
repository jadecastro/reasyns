function H = drawUnicycle(x,z)

scl = 0.20;

bodyCoords = scl*[[-1 -1 -0.5 0.5 1 1];[-0.8 0.75 1.3 1.3 0.75 -0.8]];
wheel1Coords = scl*[[-1 -1.4 -1.4 -1];[0.5 0.5 -0.5 -0.5]];
wheel2Coords = scl*[[1 1.4 1.4 1];[0.5 0.5 -0.5 -0.5]];

angleBC = atan2(bodyCoords(1,:),bodyCoords(2,:)) - x(3) + pi/2;
magBC = sqrt(bodyCoords(1,:).^2+bodyCoords(2,:).^2);
rotBodyCoords = [magBC.*sin(angleBC); magBC.*cos(angleBC)];
angleW1 = atan2(wheel1Coords(1,:),wheel1Coords(2,:)) - x(3) + pi/2;
magW1 = sqrt(wheel1Coords(1,:).^2+wheel1Coords(2,:).^2);
rotWheel1Coords = [magW1.*sin(angleW1); magW1.*cos(angleW1)];
angleW2 = atan2(wheel2Coords(1,:),wheel2Coords(2,:)) - x(3) + pi/2;
magW2 = sqrt(wheel2Coords(1,:).^2+wheel2Coords(2,:).^2);
rotWheel2Coords = [magW2.*sin(angleW2); magW2.*cos(angleW2)];

if isempty(z)
    H1=fill(rotBodyCoords(1,:)+x(1),rotBodyCoords(2,:)+x(2),'c');
    H2=fill(rotWheel1Coords(1,:)+x(1),rotWheel1Coords(2,:)+x(2),'k');
    H3=fill(rotWheel2Coords(1,:)+x(1),rotWheel2Coords(2,:)+x(2),'k');
else
    H1=fill3(rotBodyCoords(1,:)+x(1),rotBodyCoords(2,:)+x(2),z*ones(size(rotBodyCoords(1,:))),'c');
    H2=fill3(rotWheel1Coords(1,:)+x(1),rotWheel1Coords(2,:)+x(2),z*ones(size(rotWheel1Coords(1,:))),'k');
    H3=fill3(rotWheel2Coords(1,:)+x(1),rotWheel2Coords(2,:)+x(2),z*ones(size(rotWheel2Coords(1,:))),'k');
end
H = [H1;H2;H3];
