
clear
clc
close all

figure
set(gcf, 'renderer', 'opengl');

hold on
grid on

% Load Scedule Data
schedData = xlsread('schedData');
% Plot position of Target
plot3(schedData(1,6), schedData(1,7),schedData(1,8), 'r*')
plot3(schedData(1,18), schedData(1,19), schedData(1,20), 'g*')

% Load s/c position and Velocity Data
propData = xlsread('propData');
% Plot s/c trajectory
plot3(propData(:,2), propData(:,3), propData(:,4))

% Define Constants
Re = 6378;
mu = 398601;

% Create Earth
[sx sy sz] = sphere(360);
sx = sx*Re;
sy = sy*Re;
sz = sz*Re;
I = imread('earthbare.jpg');
earth = surf(sx,sy,sz)
% Rotate Earth so aligned with GMT
rotate(earth, [0 0 1], schedData(2,1));
set(earth, 'edgeColor', 'none', 'CData', I, 'Facecolor', 'texturemap');

% Set lighting
lig = camlight('infinite');
set(lig, 'position', [0 1 0]);
lighting gouraud

xlabel 'x'
ylabel 'y'
zlabel 'z'
cameratoolbar('show')
cameratoolbar('ResetCameraAndSceneLight')

r = propData(1,2:4)';
v = propData(1,5:7)';

% Create Horizon View Cone
r_n = norm(r);
betaH = asin(Re/r_n);
Dh = r_n*cos(betaH);
y = Dh * cos(betaH);
x = Dh * sin(betaH);
[cx cy cz] = cylinder(0:x);
cz = cz * y;
viewCone = surf(cx,cy,cz, 'edgecolor', 'none', 'facecolor', [.6 .1 .1], 'facealpha', .4)
% rotate cone and move into proper location
rotVec = cross(r/r_n, [0 0 1]);
rotAngle = 180 - 180/pi*acos(r'*[0; 0; 1]/r_n);
rotate(viewCone, rotVec, rotAngle)
translate(viewCone, r(1), r(2), r(3));

% create s/c
run aeolis
translate(aeolis_h, r(1), r(2), r(3));
rotate(aeolis_h, v, 90, r)

% calculate angular momentum vector
h = cross(r, v);
h_hat = h/norm(h);

% calculate the orbital period
T = 2*pi/sqrt(mu)*(Re+300)^(3/2);
dir = 1;
simStep = 20;
for simTime = (propData(1,1)+simStep):simStep:propData(end,1)%2:length(propData)
% reply = '2';
% i = 2;
% while (reply ~= 'q' || reply ~= 'Q')
%     reply = input('Direction?','s');
%     if (reply == '.')
%         i = i + 1;
%         dir = 1;
%         if (i > length(propData))
%             i = length(propData);
%         end
%     elseif (reply == ',')
%         i = i - 1;
%         dir = -1;
%         if (i < 2)
%             i = 2;
%         end
%     end
    % determine how to translate s/c
    r_old = interp1(propData(:,1), propData(:,2:4), simTime-simStep)';
    %propData(i-dir,2), propData(i-dir,3), propData(i-dir,4)]';
    r = interp1(propData(:,1), propData(:,2:4), simTime)';
    %[propData(i,2), propData(i,3), propData(i,4)]';
    d_r_sat = (r - r_old);
    % determine s/c and view cone rotation angle
%    t_old = propData(i - dir,1);
%    t = propData(i,1);
    dt = simStep;%(t - t_old);
    ang = 360*dt/T;
    
    % Move s/c and view cone
    translate(aeolis_h, d_r_sat(1), d_r_sat(2), d_r_sat(3));
%     translate(panelL, d_r_sat(1), d_r_sat(2), d_r_sat(3));
%     translate(panelR, d_r_sat(1), d_r_sat(2), d_r_sat(3));
    translate(viewCone, d_r_sat(1), d_r_sat(2), d_r_sat(3));
    % determine horizon angle for viewing
    alpha = asin(Re/norm(r)) + 2*pi/180;
    norm_l = sqrt(norm(r)^2 - Re^2);
    lr = norm_l*cos(alpha);
    lv = norm_l*sin(alpha);
    v = interp1(propData(:,1), propData(:,5:7), simTime)';
    %propData(i,5:7)';
    l = lr*(-r/norm(r)) + lv*(v/norm(v));
    % camera target location
    targ = r + l;
    

    
    % rotate s/c and view cone
    rotate(aeolis_h, h_hat, ang, r);
    rotate(viewCone, h_hat, ang, r);
    
    % camera position vector
    p = r*1.01 + 20*h_hat;
    
    campos(p);
    camup(r);
    camva(30);
    camtarget(targ);
    drawnow
%    pause;
end