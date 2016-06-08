% Aeolis Sat
clear
close all
clc

c = cube([-1 1], [-1 1], [-1 1]);

b1 = cylin(.2, .2, 1);
rotate(b1, [0 1 0], 90, [0 0 0]);
translate(b1, 1, 0, 0)
b2 = cylin(.2, .2, 1);
rotate(b2, [0 1 0], -90, [0 0 0]);
translate(b2, -1, 0, 0)

c1 = cone(.2, .2, 1)
rotate(c1, [0 1 0], -90, [0 0 0]);
translate(c1, 3, 0, 0)

c2 = cone(.2, .2, 1)
rotate(c2, [0 1 0], 90, [0 0 0]);
translate(c2, -3, 0, 0)

p1 = cube([2 7], [-1 1], [-.1 .1]);
%translate(p1, 2, 0, 0)
p2 = cube([-7 -2], [-1 1], [-.1 .1]);

panelR = merge(b1, p1, c1)
panelL = merge(b2, p2, c2)

xlabel 'x'
ylabel 'y'
zlabel 'z'
light
lighting gouraud
cameratoolbar('ResetCameraAndSceneLight')
axis equal