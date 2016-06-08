clear
clc
close all

% Name = "t185"
% Type = "gt"
% latitude = "-0.203448476"
% longitude = "-127.0248349"
% altitude = "0"
% Priority = "5"
% Value = "8"
% MinQuality = "3"
% DesiredCapTime = "3600"
% NonzeroValCapTime = "7200">

startJD = 2454368.649079210;
        
Re = 6378;
posData = dlmread('Asset1_PositionData.txt', '\t', 0, 1)';
velData = dlmread('Asset1_VelocityData.txt', '\t', 0, 1)';
schedData = dlmread('Asset1_TargetData.txt', '\t', 2, 1)';
load info

taskStart = schedData(:,2);
targLoc = schedData(:,5:6);

for i = 1:length(posData)
   
    v = velData(i,2:4);
    r = posData(i,2:4);
    rn = norm(r);
    t = posData(i,1);
    
    [lat, lon, alt] = eci2lla(startJD + t/86400, r);
    sysPLat(i) = lat;
    sysPLon(i) = lon;
    sysPAlt(i) = alt;

%     b = asin(Re/norm(r));
%     h_hat = cross(r, v)/norm(cross(r, v));
%     ph = h_hat*sin(b);
%     pr = -r/rn*cos(b);
%     p_hat = ph + pr;
%     pn = norm(r)*cos(b);
%     p = pn*p_hat;
%     tm = r + p;
%     if (t >24600)
%         y = 7;
%     end
%         
%     [lat, lon, alt] = eci2lla(startJD + t/86400,tm);
%     tlat1(i) = lat;
%     tlon1(i) = lon;
%     talt1(i) = alt;
%     
%     ph = -h_hat*sin(b);
%     pr = -r/rn*cos(b);
%     p_hat = ph + pr;
%     pn = norm(r)*cos(b);
%     p = pn*p_hat;
%     tm = r + p;
%     [lat, lon, alt] = eci2lla(startJD + t/86400,tm);
%     tlat2(i) = lat;
%     tlon2(i) = lon;
%     talt2(i) = alt;
    
end

t = 0;
ct = 0;
for i = 1:length(taskStart)
    ind = find(posData(:,1) == taskStart(i));
    if(Sheet1{i+1, 2} == 'FacilityTarget')
        ct = ct+1;
        sysCLat(ct) = sysPLat(ind);
        sysCLon(ct) = sysPLon(ind);
        sysCAlt(ct) = sysPAlt(ind);
        targetCLat(ct) = targLoc(i,1);
        targetCLon(ct) = targLoc(i,2);
        targetCAlt(ct) = 0;
    else
        t = t+1;
        sysLat(t) = sysPLat(ind);
        sysLon(t) = sysPLon(ind);
        sysAlt(t) = sysPAlt(ind);
        targetLat(t) = targLoc(i,1);
        targetLon(t) = targLoc(i,2);
        targetAlt(t) = 0;
    end
end

targetLat = -1024/180*[targetLat, targetLat] + 1024/2;
targetLon = 2048/360*[targetLon, targetLon+360] + 1023;
targetAlt = [targetAlt, targetAlt];

targetCLat = -1024/180*[targetCLat, targetCLat] + 1024/2;
targetCLon = 2048/360*[targetCLon, targetCLon+360] + 1023;
targetCAlt = [targetCAlt, targetCAlt];

sysPLat = -1024/180*[sysPLat, sysPLat] + 1024/2;
sysPLon = 2048/360*[sysPLon, sysPLon+360] + 1023;
sysPAlt = [sysPAlt, sysPAlt];

sysCLat = -1024/180*[sysCLat, sysCLat] + 1024/2;
sysCLon = 2048/360*[sysCLon, sysCLon+360] + 1023;
sysCAlt = [sysCAlt, sysCAlt];

sysLat = -1024/180*[sysLat, sysLat] + 1024/2;
sysLon = 2048/360*[sysLon, sysLon+360] + 1023;
sysAlt = [sysAlt, sysAlt];

% tlat1 = 768/180*[tlat1, tlat1] + 768/2;
% tlon1 = 1024/360*[tlon1, tlon1+360] - 500;
% talt1 = [talt1, talt1];
% 
% tlat2 = 768/180*[tlat2, tlat2] + 768/2;
% tlon2 = 1024/360*[tlon2, tlon2+360] - 500;
% talt2 = [talt2, talt2];

earth = imread('land_ocean_ice_2048.jpg');
image(earth);
%figure
hold on
grid on

ind = find(abs(diff(sysPLon)) > 100);
ind = [0, ind];

plot3(sysPLon(ind(end)+1:end), sysPLat(ind(end)+1:end), sysPAlt(ind(end)+1:end))
plot3(sysLon, sysLat, sysAlt, '.')

for i = 1:length(ind)-1
    plot3(sysPLon(ind(i)+1:ind(i+1)), sysPLat(ind(i)+1:ind(i+1)), sysPAlt(ind(i)+1:ind(i+1)))
end

for i = 1:length(targetLon)
    line([targetLon(i), sysLon(i)], [targetLat(i), sysLat(i)], [targetAlt(i) sysAlt(i)], 'color', 'g')
end

for i = 1:length(targetCLon)
    line([targetCLon(i), sysCLon(i)], [targetCLat(i), sysCLat(i)], [targetCAlt(i) sysCAlt(i)], 'color', 'y')
end

plot(targetLon, targetLat, 'g.', 'MarkerSize', 12)
plot(targetCLon, targetCLat, 'y.', 'MarkerSize', 12)

% plot3(tlon1, tlat1, talt1, 'y.', 'markerSize', 6)
% plot3(tlon2, tlat2, talt2, 'y.', 'markerSize', 6)
% xlim([180 540])
% ylim([-90, 90])
% zlim([-10, 700])
%set(gca, 'yDir', 'normal')
set(gca, 'xTick', [0:2048/12:2048], 'xTickLabel', [-180:30:180])
set(gca, 'yTick', [0:1024/12:1024], 'yTickLabel', [90:-15:-90])
