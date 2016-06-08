% script used to test the Horizon propagator
% comparing to Matlab ode45 solver
% Eric Mehiel
% 12/5/06

clear
clc
close all

horizonPropData = xlsread('propData.xls');
ht = horizonPropData(:,1);
hp = horizonPropData(:,2:4);
hv = horizonPropData(:,5:7);

initPos = [6678.0,0.0,0.0]';
initVel = [0.0,7.2599,2.6424]';
options = odeset('reltol', 1e-6, 'abstol', 1e-6);
[t, y] = ode45('eoms', [0 5400], [initPos; initVel], options);

plot(t,y(:,1:3), 'b');
hold on
plot(ht, hp, 'r');
figure
plot(t,y(:,4:6),'b');
hold on
plot(ht, hv, 'r');

[t, y] = ode45('eoms', [0:10:5400], [initPos; initVel], options);

hpi = interp1(ht, hp, 0:10:5400, 'spline');
hvi = interp1(ht, hv, 0:10:5400, 'spline');
% Plot the percent error of both position and velocity given a spline fit
% to Horizon Propagator Data
figure
plot(0:10:5400,(hpi - y(:,1:3))/norm(initPos), 'b');
hold on
plot(0:10:5400, (hvi - y(:,4:6))/norm(initVel), 'r');