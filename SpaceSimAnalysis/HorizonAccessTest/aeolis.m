function aeolis_h = aeolis()
% Aeolis Sat
% clear
% close all
% clc

s = 10;
c = cube([-1 1]*s, [-1 1]*s, [-1 1]*s);

b1 = cylin(.2*s, .2*s, 1*s);
rotate(b1, [0 1 0], 90, [0 0 0]);
translate(b1, 1*s, 0, 0)
b2 = cylin(.2*s, .2*s, 1*s);
rotate(b2, [0 1 0], -90, [0 0 0]);
translate(b2, -1*s, 0, 0)

c1 = cone(.2*s, .2*s, 1*s)
rotate(c1, [0 1*s 0], -90, [0 0 0]);
translate(c1, 3*s, 0, 0)

c2 = cone(.2*s, .2*s, 1*s)
rotate(c2, [0 1 0], 90, [0 0 0]);
translate(c2, -3*s, 0, 0)

p1 = cube([2*s 7*s], [-1*s 1*s], [-.1*s .1*s]);
%translate(p1, 2, 0, 0)
p2 = cube([-7 -2]*s, [-1 1]*s, [-.1 .1]*s);

aeolis_h = merge(b1, p1, c1, b2, p2, c2, c);
%panelL = merge(b2, p2, c2);

% xlabel 'x'
% ylabel 'y'
% zlabel 'z'
% light
% lighting gouraud
% cameratoolbar('ResetCameraAndSceneLight')
% axis equal