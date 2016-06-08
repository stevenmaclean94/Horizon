clear
clc
close

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
        
Re = 6378;
propData = xlsread('propData.xls');
schedData = xlsread('schedData1.xls');

taskStart = schedData(:,2);
targ = schedData(:,6:8);
scpos = schedData(:,18:20);

for i = 1:length(scpos)
    
    tPos = schedData(i,6:8);
    r = schedData(i,18:20);
    rn = norm(r);
    start = schedData(i,2);
    
    [lat, lon, alt] = eci2lla2(start,tPos);
    targetLat(i) = lat;
    if lon < 0;
        lon = lon+360;
    end
    targetLon(i) = lon;
    targetAlt(i) = alt;
    
    [lat, lon, alt] = eci2lla2(start,r);
    sysLat(i) = lat;
    if lon < 0;
        lon = lon+360;
    end
    sysLon(i) = lon;
    sysAlt(i) = alt;
    
    v = interp1(propData(:,1), propData(:,5:7), (start-2454680)*86400);
    b = asin(Re/norm(r));
    h_hat = cross(r, v)/norm(cross(r, v));
    ph = h_hat*sin(b);
    pr = -r/rn*cos(b);
    p_hat = ph + pr;
    pn = norm(r)*cos(b);
    p = pn*p_hat;
    tm = r + p;
    [lat, lon, alt] = eci2lla2(start,tm);
    tlat1(i) = lat;
    if lon < 0;
        lon = lon+360;
    end
    tlon1(i) = lon;
    talt1(i) = alt;
    
    ph = -h_hat*sin(b);
    pr = -r/rn*cos(b);
    p_hat = ph + pr;
    pn = norm(r)*cos(b);
    p = pn*p_hat;
    tm = r + p;
    [lat, lon, alt] = eci2lla2(start,tm);
    tlat2(i) = lat;
    if lon < 0;
        lon = lon+360;
    end
    tlon2(i) = lon;
    talt2(i) = alt;
end

for i = 1:length(propData)
    [lat, lon, alt] = eci2lla2(2454680 + propData(i,1)/86400, propData(i,2:4));
    sysPLat(i) = lat;
    if lon < 0;
        lon = lon+360;
    end
    sysPLon(i) = lon;
    sysPAlt(i) = alt;
end

targetLat = 768/180*[targetLat, targetLat] + 768/2;
targetLon = 1024/360*[targetLon, targetLon+360] - 500;
targetAlt = [targetAlt, targetAlt];

sysPLat = 768/180*[sysPLat, sysPLat] + 768/2;
sysPLon = 1024/360*[sysPLon, sysPLon+360] - 500;
sysPAlt = [sysPAlt, sysPAlt];

sysLat = 768/180*[sysLat, sysLat] + 768/2;
sysLon = 1024/360*[sysLon, sysLon+360] - 500;
sysAlt = [sysAlt, sysAlt];

tlat1 = 768/180*[tlat1, tlat1] + 768/2;
tlon1 = 1024/360*[tlon1, tlon1+360] - 500;
talt1 = [talt1, talt1];

tlat2 = 768/180*[tlat2, tlat2] + 768/2;
tlon2 = 1024/360*[tlon2, tlon2+360] - 500;
talt2 = [talt2, talt2];

earth = imread('earthbareDisp.jpg');
image(earth);%figure
hold on
grid on

ind = find(abs(diff(sysPLon)) > 100);
ind = [0, ind];

for i = 1:length(ind)-1
plot3(sysPLon(ind(i)+1:ind(i+1)), sysPLat(ind(i)+1:ind(i+1)), sysPAlt(ind(i)+1:ind(i+1)))
end

for i = 1:length(targetLon)
    line([targetLon(i), sysLon(i)], [targetLat(i), sysLat(i)], [0 300], 'color', 'g')
end
plot(targetLon, targetLat, 'r.')

plot3(sysPLon(ind(end)+1:end), sysPLat(ind(end)+1:end), sysPAlt(ind(end)+1:end))
plot3(sysLon, sysLat, 300*ones(size(sysAlt)), '.')

plot3(tlon1, tlat1, talt1, 'k.', 'markerSize', 2)
plot3(tlon2, tlat2, talt2, 'k.', 'markerSize', 2)
% xlim([180 540])
% ylim([-90, 90])
% zlim([-10, 700])
set(gca, 'yDir', 'normal')
set(gca, 'xTick', [0:1024/6:1024], 'xTickLabel', [-180:60:180])
set(gca, 'yTick', [0:768/6:768], 'yTickLabel', [-90:30:90])